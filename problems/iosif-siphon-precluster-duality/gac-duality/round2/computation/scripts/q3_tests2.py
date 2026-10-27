"""Second batch of tests of agent Q3's candidates (round2/ec/candidates.md, update of 19:28).

Environment: Python 3.12 (WSL Fedora 44), standard library + msolve 0.10.1 (certified real root isolation) for the
explicit examples; exact rational arithmetic elsewhere.  Networks: agent P's generator (enumerate_ec.gen_networks,
imported unchanged) for the WF-EC families given on the command line; first-round records for Q3-CR.
Usage (from round2/computation/): python scripts/q3_tests2.py [NPROC]

  Ustar  for WF-EC N (families (3,3,4), (4,2,4), (2,4,4), weighted) and 2 random rational kappa and 3 random
         positive points k of N* (not steady states): det(tI - J*) = t^(r-s) q(t) with q Hurwitz (exact);
  GS     for the same N, 3 random rational x > 0 and rates k: G = D(phi)(x) A_k is Metzler and G u(k) = 0, with u(k)
         the positive kernel vector of A_k (exact);
  CR     (conjecture Q3-CR) every first-round network that is connected with (H1), (H2) and s = #clusters - 1 has
         no critical siphon in N nor in N*, and is absolutely consistent when decided by k12_abscons.verdict;
  B      Q3's examples: NM (k = (1,2,1), class x1+x2+x3 = 5/2) and NB (k = 1, class 2x1+x2+x3 = 6) have two
         positive steady states in that class (msolve count); NT* is not certified absolutely consistent (exact verdict NOT1, or OPEN when d >= 2) and
         has no positive steady state for kappa = (1, 2, 1) in a sampled class (msolve);
  K8s    Q3's pair M3, M4: both (H1)-(H3), weakly reversible, pure, rows and columns of Gamma pairwise linearly
         independent; for every relabelling of M4's reactions that matches the doubling classes and the minimal
         T-traps (minimal siphons of the dual) of M3, P differs (recomputed with crn.py), and there is one.
"""
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import permutations
from multiprocessing import Pool

import e2lib
from e2lib import CNet, crn, ecore
import q3_tests
import k12_abscons
import ms_dual

sys.path.insert(0, e2lib.PSCRIPTS)
import enumerate_ec as EE  # noqa: E402


def ustar_gs(n, rx, seed):
    rng = random.Random(seed)
    N = CNet(n, rx)
    G = N.Gamma
    r = len(rx)
    s = N.rank()
    Ye = [[rx[j][0][i] for j in range(r)] for i in range(n)]
    ok_u = True
    for _ in range(2):
        kap = [Fraction(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(n)]
        for _ in range(3):
            k = [Fraction(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(r)]
            mono = [kap[i] * ms_dual.mono(k, Ye[i]) for i in range(n)]
            J = [[sum(Fraction(G[i][j]) * mono[i] * Ye[i][l] / k[l] for i in range(n)) for l in range(r)]
                 for j in range(r)]
            cp = q3_tests.charpoly(J)
            if any(c != 0 for c in cp[s + 1:]) or not q3_tests.hurwitz(cp[:s + 1]):
                ok_u = False
    # GS
    srcs = sorted({a for a, b in rx})
    m = len(srcs)
    ok_g = True
    for _ in range(3):
        x = [Fraction(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(n)]
        kk = [Fraction(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(r)]
        A = [[sum(kk[t] * (rx[t][1][p] - rx[t][0][p]) for t in range(r) if rx[t][0] == srcs[c]) for c in range(m)]
             for p in range(n)]
        K = crn.nullspace(A, m)
        if len(K) != 1:
            ok_g = False
            continue
        u = K[0]
        if u[0] < 0:
            u = [-v for v in u]
        if not all(v > 0 for v in u):
            ok_g = False
            continue
        Dphi = [[srcs[c][p] * ms_dual.mono(x, srcs[c]) / x[p] for p in range(n)] for c in range(m)]
        Gm = [[sum(Dphi[c][p] * A[p][c2] for p in range(n)) for c2 in range(m)] for c in range(m)]
        if any(Gm[c][c2] < 0 for c in range(m) for c2 in range(m) if c != c2):
            ok_g = False
        if any(sum(Gm[c][c2] * u[c2] for c2 in range(m)) != 0 for c in range(m)):
            ok_g = False
    return ok_u, ok_g


def worker(args):
    n, m, rmax, pidx = args
    perms = list(permutations(range(n)))
    seen = set()
    out = []
    for rx in EE.gen_networks(n, m, rmax, 0, pidx):
        rx = tuple(sorted(rx))
        can = EE.canonical(rx, perms)
        if can != rx or can in seen:
            continue
        seen.add(can)
        if not e2lib.petri_connected(n, rx):
            continue
        N = CNet(n, rx)
        if not N.is_consistent() or not N.is_conservative() or N.rank() != e2lib.n_sources(rx) - 1:
            continue
        bal = all(len({x for x in a if x > 0}) <= 1 for a, b in rx)
        u, g = ustar_gs(n, rx, hash(rx) & 0xffff)
        out.append((n, rx, bal, u, g))
    return out


def cr_check(rc):
    n, rx = rc["n"], rc["rx"]
    N = rc["N"]
    if not (N["cons"] and N["consist"] and e2lib.petri_connected(n, rx)):
        return None
    ncl = len(rc["N"]["CLr"])
    if rc["s"] != ncl - 1:
        return None
    v = k12_abscons.verdict(n, rx)[0]
    return dict(n=n, rx=rx, ec=e2lib.is_EC(rx), csN=bool(rc["N"]["CS"]), csD=bool(rc["D"]["CS"]), v=v)


def examples():
    ok = True
    F = Fraction

    def count(Xn, Xrx, k, x0):
        return ms_dual.msolve_count(len(Xrx), e2lib.dual_rx(Xn, Xrx), k, x0)
    NM = (((0, 0, 2), (1, 1, 0)), ((0, 0, 3), (0, 1, 2)), ((1, 2, 0), (0, 0, 3)))
    NB = (((0, 0, 1), (0, 1, 0)), ((0, 0, 2), (1, 0, 0)), ((1, 1, 0), (0, 0, 3)))
    NT = (((1, 1, 0, 0), (0, 2, 0, 0)), ((1, 1, 0, 0), (0, 0, 2, 1)), ((0, 0, 2, 1), (2, 0, 0, 0)))
    cm = count(3, NM, [F(1), F(2), F(1)], [F(5, 6)] * 3)
    cb = count(3, NB, [F(1)] * 3, [F(1), F(2), F(2)])
    vt = k12_abscons.verdict(3, e2lib.dual_rx(4, NT))
    ct = ms_dual.msolve_count(4, NT, [F(1), F(2), F(1), F(1)], [F(1), F(2), F(3)])
    okb = cm[0] == 2 and cb[0] == 2 and vt[0] in ("NOT1", "OPEN") and ct[0] == 0
    ok &= okb
    print(f"{'PASS' if okb else 'FAIL'} B: NM two positive steady states {cm}; NB {cb}; NT* verdict {vt}; NT* "
          f"positive steady states for kappa = (1,2,1,1) in the class of (1,2,3): {ct}")
    # K8s
    M3 = (((0, 0, 0, 1), (0, 0, 2, 0)), ((0, 0, 0, 1), (0, 1, 1, 0)), ((0, 0, 2, 0), (0, 2, 0, 0)),
          ((0, 1, 1, 0), (2, 0, 0, 0)), ((0, 2, 0, 0), (0, 0, 0, 1)), ((2, 0, 0, 0), (0, 0, 2, 0)))
    M4 = (((0, 0, 0, 1), (0, 1, 1, 0)), ((0, 0, 2, 0), (0, 0, 0, 1)), ((0, 0, 2, 0), (0, 2, 0, 0)),
          ((0, 1, 1, 0), (2, 0, 0, 0)), ((0, 2, 0, 0), (0, 0, 0, 1)), ((2, 0, 0, 0), (0, 0, 2, 0)))

    def props(rx):
        N = CNet(4, rx)
        D = N.dual()
        pure = not any(a[i] and b[i] for a, b in rx for i in range(4))
        dbl = {}
        for j, (a, b) in enumerate(rx):
            dbl.setdefault(a, set()).add(j)
        mt = [ecore.to_set(S) for S in ecore.minimal_masks(ecore.siphon_masks(D))]
        P = [frozenset(b) for b in N.preclusters("onestep")]
        return dict(H=N.is_conservative() and N.is_consistent() and N.rows_distinct(), WR=N.weakly_reversible(),
                    pure=pure, strong=N.rows_pairwise_independent() and N.columns_pairwise_independent(),
                    dbl=frozenset(frozenset(v) for v in dbl.values()), mt=frozenset(mt), P=frozenset(P))
    p3, p4 = props(M3), props(M4)
    match = 0
    differ = 0
    for pi in permutations(range(6)):
        rel = lambda fam: frozenset(frozenset(pi[j] for j in B) for B in fam)  # noqa: E731
        if rel(p4["dbl"]) == p3["dbl"] and rel(p4["mt"]) == p3["mt"]:
            match += 1
            differ += rel(p4["P"]) != p3["P"]
    okk = all(p[k] for p in (p3, p4) for k in ("H", "WR", "pure", "strong")) and match > 0 and differ == match
    ok &= okk
    print(f"{'PASS' if okk else 'FAIL'} K8s: M3, M4 (H1)-(H3)/WR/pure/strong: {[p3[k] for k in ('H', 'WR', 'pure', 'strong')]}"
          f" {[p4[k] for k in ('H', 'WR', 'pure', 'strong')]}; relabellings matching doubling classes and minimal "
          f"T-traps: {match}, with different P: {differ}; P(M3) = {sorted(map(sorted, p3['P']))}, "
          f"P(M4) = {sorted(map(sorted, p4['P']))}")
    return ok


def main():
    t0 = time.time()
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    print(f"# q3_tests2.py  Python {platform.python_version()} on {platform.system()}; exact; msolve 0.10.1; "
          f"nproc={nproc}")
    ok = True
    # Ustar and GS
    for fam in ((3, 3, 4), (4, 2, 4), (2, 4, 4)):
        n, m, rmax = fam
        nparts = sum(1 for _ in EE.set_partitions(list(range(n))))
        with Pool(nproc) as pool:
            res = [q for ch in pool.map(worker, [(n, m, rmax, p) for p in range(nparts)]) for q in ch]
        bu = [q for q in res if not q[3]]
        bg = [q for q in res if not q[4]]
        ok &= not bu and not bg
        print(f"{'PASS' if not bu else 'FAIL'} Ustar family {fam}: {len(res)} WF-EC networks "
              f"({sum(not q[2] for q in res)} unbalanced), J* Hurwitz on im Gamma^T at 6 random points each "
              f"({len(bu)} failures)")
        for q in bu[:2]:
            print(f"     failure: {e2lib.describe(q[0], q[1])}")
        print(f"{'PASS' if not bg else 'FAIL'} GS family {fam}: G = Dphi A_k Metzler with G u(k) = 0 at 3 random "
              f"(x, k) each ({len(bg)} failures)")
        for q in bg[:2]:
            print(f"     failure: {e2lib.describe(q[0], q[1])}")
    # CR
    recs = e2lib.load(e2lib.E_FILES)
    with Pool(nproc) as pool:
        cr = [q for q in pool.map(cr_check, recs, chunksize=500) if q is not None]
    bad = [q for q in cr if q["csN"] or q["csD"] or q["v"] == "NOT1"]
    nonec = [q for q in cr if not q["ec"]]
    vc = {}
    for q in cr:
        vc[q["v"]] = vc.get(q["v"], 0) + 1
    print(f"{'HOLDS' if not bad else 'REFUTED'} CR (conjecture): {len(cr)} connected networks with (H1), (H2), "
          f"s = #clusters - 1 ({len(nonec)} not EC); critical siphon in N or N*, or not absolutely consistent: "
          f"{len(bad)}; verdicts {vc}")
    bad.sort(key=lambda q: (q["n"] + len(q["rx"]), q["rx"]))
    for q in bad[:3]:
        print(f"     counterexample: {e2lib.describe(q['n'], q['rx'])} (EC={q['ec']}, crit N={q['csN']}, "
              f"crit N*={q['csD']}, verdict {q['v']})")
    ok &= examples()
    print(f"{'PASS' if ok else 'FAIL'}: all claims (CR is a conjecture test)")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
