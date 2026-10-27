"""Exact decision of persistence of the dual N* when s = 1 (task 2 of BRIEF-E2; Q1 task 1).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (e2lib -> ecore -> crn).
Usage (from round2/computation/):  python scripts/s1_persistence.py [--nsample K]
Input: every first-round record file (~/gacE-data/*.pkl.gz); networks with s = 1 and (H1), (H2) for N, so that N*
is conservative and consistent and its compatibility classes are compact segments.

Setting.  s = 1: Gamma = a b^T (a in Z^n primitive, b in Q^r).  N* has species k_1..k_r and reactions i = 1..n with
source Ye[i,:] and reaction vector a_i b.  On the class k0 + R b the dynamics is t' = g(t),
g = sum_i kappa_i a_i k^{Ye[i,:]}.  The end t+ (species with b_j < 0 vanish) and the end t- (b_j > 0) are reached
with zero set any nonempty J of {b_j < 0}, resp. {b_j > 0}, and any positive values x on the other coordinates.
At the end of side sigma (+1 for t+, -1 for t-), a reaction i is 'bad' if sign(a_i) = sigma (it pushes towards the
end) and 'good' if sign(a_i) = -sigma.  Near the end, with the J-coordinates of order eps, the term i is of order
eps^{d_i(J)}, d_i(J) = sum_{j in J} Ye[i,j], with coefficient kappa_i |a_i| beta_i x^{e_i}
(beta_i = prod_{j in J} |b_j|^{Ye[i,j]}, e_i = Ye[i, J^c]).  The class is non-persistent iff g points to the end on
a whole neighbourhood of it, i.e. iff the lowest-order coefficient has the bad sign (g is a polynomial).

Decision (proved in notes.md, Prop. E2.1):
  ALL  (persistent for every kappa): for every side and every J, min_good d(J) < min_bad d(J).
  NONE (non-persistent for every kappa): some side and J with min_bad d(J) < min_good d(J), or a tie in which some
       tied bad exponent e_i lies outside the convex hull of the tied good exponents.
  SOME (persistent for some kappa, not for others): otherwise, when the rate LP below is feasible.  Certificates:
       kappa_1 satisfying, at every tie, the AM-GM domination inequalities (each bad monomial is bounded by a convex
       combination of good monomials with the same weights, coefficientwise strictly smaller), and kappa_2 = bad
       rates large at one tie (lowest-order coefficient > 0 at x = 1).
  HYPERPLANE: all nontrivial dual reactions share one source y; then g = (sum_i kappa_i a_i) k^y and N* is
       persistent iff sum_i kappa_i a_i = 0 (a hyperplane of rate constants), non-persistent otherwise.
  UNDECIDED: ties only, rate LP infeasible (not decided by this script).
Claims checked (PASS/FAIL): F factorisation Gamma = a b^T exact; L structural lemma (bad-at-t+ reactions contain all
of {b_j<0} in their source, and dually); A ADLS(N*) => ALL; B drain species in N* => NONE; W WF-EC N => ALL
(Cor. P of agent P); K both certificates of every SOME network re-verified with the exact kappas; T sampled exact
evaluation of g near the ends agrees with the leading-order prediction.
"""
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import combinations

import e2lib
from e2lib import CNet, endo


def sgn(x):
    return (x > 0) - (x < 0)


def factor(n, rx):
    G = CNet(n, rx).Gamma
    r = len(rx)
    j0 = next(j for j in range(r) if any(G[i][j] for i in range(n)))
    a = [G[i][j0] for i in range(n)]
    from math import gcd
    g = 0
    for x in a:
        g = gcd(g, abs(x))
    a = [x // g for x in a]
    i0 = next(i for i in range(n) if a[i])
    if a[i0] < 0:
        a = [-x for x in a]
    b = [Fraction(G[i0][j], a[i0]) for j in range(r)]
    ok = all(G[i][j] == a[i] * b[j] for i in range(n) for j in range(r))
    return a, b, ok


def analyse(n, rx):
    """Returns dict(cat, cert, ties, lemma_ok, factor_ok, a, b)."""
    r = len(rx)
    a, b, fok = factor(n, rx)
    Ye = [[rx[j][0][i] for j in range(r)] for i in range(n)]
    lemma = True
    for i in range(n):
        for j in range(r):
            if a[i] * b[j] < 0 and Ye[i][j] == 0:  # reaction i consumes k_j net, so k_j must be in its source
                lemma = False
    none_cert = None
    ties = []
    for side in (1, -1):
        pool = [j for j in range(r) if sgn(b[j]) == -side]
        bad = [i for i in range(n) if sgn(a[i]) == side]
        good = [i for i in range(n) if sgn(a[i]) == -side]
        for size in range(1, len(pool) + 1):
            for J in combinations(pool, size):
                d = {i: sum(Ye[i][j] for j in J) for i in bad + good}
                dbad = min((d[i] for i in bad), default=None)
                dgood = min((d[i] for i in good), default=None)
                if dbad is None:
                    continue
                if dgood is None or dbad < dgood:
                    if none_cert is None:
                        none_cert = ("order", side, J, dbad, dgood)
                    continue
                if dgood < dbad:
                    continue
                Jc = [j for j in range(r) if j not in J]

                def coef(i):
                    c = Fraction(abs(a[i]))
                    for j in J:
                        c *= abs(b[j]) ** Ye[i][j]
                    return c
                tb = [i for i in bad if d[i] == dbad]
                tg = [i for i in good if d[i] == dgood]
                gexp = sorted({tuple(Ye[i][j] for j in Jc) for i in tg})
                rows = []
                for i in tb:
                    e = tuple(Ye[i][j] for j in Jc)
                    if e in gexp:
                        lam = [Fraction(int(g == e)) for g in gexp]
                    else:
                        lam = e2lib.in_convex_hull(list(e), [list(g) for g in gexp]) if Jc else None
                    if lam is None:
                        if none_cert is None:
                            none_cert = ("vertex", side, J, e)
                        break
                    rows.append((i, coef(i), lam))
                else:
                    ties.append(dict(side=side, J=J, tb=tb, tg=tg, gexp=gexp,
                                     gcoef={i: coef(i) for i in tg},
                                     gidx={i: gexp.index(tuple(Ye[i][j] for j in Jc)) for i in tg},
                                     rows=rows))
    out = dict(a=a, b=b, factor_ok=fok, lemma_ok=lemma, ties=ties)
    srcs = {tuple(Ye[i]) for i in range(n) if a[i] != 0}
    if len(srcs) == 1:
        # every nontrivial dual reaction has the same source y: g = (sum_i kappa_i a_i) k^y, so N* is persistent
        # iff sum_i kappa_i a_i = 0 (g identically 0); otherwise the orbit drifts monotonically to an end.
        out["cat"], out["cert"] = "HYPERPLANE", dict(normal=a)
        return out
    if none_cert is not None:
        out["cat"], out["cert"] = "NONE", none_cert
        return out
    if not ties:
        out["cat"], out["cert"] = "ALL", None
        return out
    # rate LP: kappa_i >= 1, and for each tie and each good exponent k:
    #   sum_rows kappa_i c_i lam_k - sum_{good i at k} kappa_i c_i <= -1
    act = [i for i in range(n) if a[i] != 0]
    pos = {i: t for t, i in enumerate(act)}
    ineqs = [([-int(t == u) for u in range(len(act))], -1) for t in range(len(act))]
    for T in ties:
        for k in range(len(T["gexp"])):
            v = [Fraction(0)] * len(act)
            for i, c, lam in T["rows"]:
                v[pos[i]] += c * lam[k]
            for i in T["tg"]:
                if T["gidx"][i] == k:
                    v[pos[i]] -= T["gcoef"][i]
            ineqs.append((v, -1))
    kap = e2lib.fm_point(ineqs, len(act))
    T0 = ties[0]
    sg = sum(T0["gcoef"].values())
    mb = min(c for i, c, lam in T0["rows"])
    M = sg / mb + 1
    kappa2 = [Fraction(1) if sgn(a[i]) != T0["side"] else M for i in range(n)]
    if kap is None:
        out["cat"], out["cert"] = "UNDECIDED", dict(kappa2=kappa2, tie=T0)
        return out
    kappa1 = [Fraction(0)] * n
    for i in act:
        kappa1[i] = kap[pos[i]]
    for i in range(n):
        if a[i] == 0:
            kappa1[i] = Fraction(1)
    out["cat"], out["cert"] = "SOME", dict(kappa1=kappa1, kappa2=kappa2, tie=T0)
    return out


def leading(n, rx, a, b, kappa, side, J, x):
    """Exact lowest-order coefficient (and its order) of sigma-signed g at the end of side `side`, zero set J,
    other coordinates x (dict j -> value).  Returns (order, coefficient of g)."""
    r = len(rx)
    Ye = [[rx[j][0][i] for j in range(r)] for i in range(n)]
    terms = {}
    for i in range(n):
        if a[i] == 0:
            continue
        dd = sum(Ye[i][j] for j in J)
        c = kappa[i] * a[i]
        for j in range(r):
            if j in J:
                c *= abs(b[j]) ** Ye[i][j]
            else:
                c *= x[j] ** Ye[i][j]
        terms[dd] = terms.get(dd, 0) + c
    for dd in sorted(terms):
        if terms[dd] != 0:
            return dd, terms[dd]
    return None, Fraction(0)


def g_at(n, rx, a, b, kappa, side, J, x, eps):
    r = len(rx)
    k = [Fraction(0)] * r
    for j in range(r):
        base = Fraction(0) if j in J else x[j]
        k[j] = base - side * eps * b[j]
    tot = Fraction(0)
    for i in range(n):
        m = kappa[i] * a[i]
        for j in range(r):
            m *= k[j] ** rx[j][0][i]
        tot += m
    return tot


def verify_certs(n, rx, res):
    """Re-verify a SOME certificate exactly: kappa1 satisfies the domination inequalities at every tie and every
    non-tie end is strictly fine (structure, independent of kappa); kappa2 gives a positive bad coefficient."""
    a, b = res["a"], res["b"]
    k1, k2 = res["cert"]["kappa1"], res["cert"]["kappa2"]
    if any(k <= 0 for k in k1) or any(k <= 0 for k in k2):
        return False
    for T in res["ties"]:
        for kk in range(len(T["gexp"])):
            v = sum(k1[i] * c * lam[kk] for i, c, lam in T["rows"])
            v -= sum(k1[i] * T["gcoef"][i] for i in T["tg"] if T["gidx"][i] == kk)
            if not v < 0:
                return False
    T = res["cert"]["tie"]
    r = len(rx)
    x = {j: Fraction(1) for j in range(r) if j not in T["J"]}
    dd, c = leading(n, rx, a, b, k2, T["side"], T["J"], x)
    return c != 0 and sgn(c) == T["side"]


def main():
    t0 = time.time()
    nsample = 3
    if "--nsample" in sys.argv:
        nsample = int(sys.argv[sys.argv.index("--nsample") + 1])
    print(f"# s1_persistence.py  Python {platform.python_version()} on {platform.system()}; exact arithmetic")
    recs = e2lib.load(e2lib.E_FILES)
    recs = [rc for rc in recs if rc["s"] == 1 and rc["N"]["cons"] and rc["N"]["consist"]]
    print(f"# {len(recs)} distinct networks with s = 1 and (H1), (H2) for N (all first-round record files)")
    CL = ["all", "H", "WR", "WR0"]
    cats = ["ALL", "SOME", "HYPERPLANE", "NONE", "UNDECIDED"]
    count = {c: {k: 0 for k in cats} for c in CL}
    hd = {k: 0 for k in cats}
    smallest = {}
    fails = {k: [] for k in "FLABWKT"}
    napp = {k: 0 for k in "FLABWKT"}
    rng = random.Random(11)
    none_kinds = {"order": 0, "vertex": 0}
    for rc in recs:
        n, rx = rc["n"], rc["rx"]
        res = analyse(n, rx)
        cat = res["cat"]
        for c in e2lib.classes(rc):
            count[c][cat] += 1
        if e2lib.H(rc) and e2lib.Hd(rc):
            hd[cat] += 1
        if cat == "NONE":
            none_kinds[res["cert"][0]] += 1
        key = (n + len(rx), len(rx), rx)
        for c in e2lib.classes(rc):
            kk = (c, cat)
            if kk not in smallest or key < smallest[kk][0]:
                smallest[kk] = (key, rc, res)
        napp["F"] += 1
        if not res["factor_ok"]:
            fails["F"].append(rx)
        napp["L"] += 1
        if not res["lemma_ok"]:
            fails["L"].append(rx)
        if not rc["D"]["CS"]:  # no critical siphon in N*
            napp["A"] += 1
            if cat != "ALL":
                fails["A"].append(rx)
        if endo.drain_species(len(rx), e2lib.dual_rx(n, rx)):
            napp["B"] += 1
            if cat != "NONE":
                fails["B"].append(rx)
        if e2lib.wfec(rc):
            napp["W"] += 1
            if cat != "ALL":
                fails["W"].append(rx)
        if cat == "SOME":
            napp["K"] += 1
            if not verify_certs(n, rx, res):
                fails["K"].append(rx)
        # T: sampled exact evaluation near the ends
        a, b = res["a"], res["b"]
        r = len(rx)
        for _ in range(nsample):
            kappa = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(n)]
            side = rng.choice((1, -1))
            pool = [j for j in range(r) if sgn(b[j]) == -side]
            if not pool:
                continue
            J = tuple(sorted(rng.sample(pool, rng.randint(1, len(pool)))))
            x = {j: Fraction(rng.randint(1, 30), rng.randint(1, 30)) for j in range(r) if j not in J}
            dd, c = leading(n, rx, a, b, kappa, side, J, x)
            if c == 0:
                continue
            napp["T"] += 1
            val = g_at(n, rx, a, b, kappa, side, J, x, Fraction(1, 10 ** 40))
            if sgn(val) != sgn(c):
                fails["T"].append(rx)
            # consistency with the category: ALL => never bad sign; NONE says nothing about this sample
            if cat == "ALL" and sgn(c) == side:
                fails["T"].append(("ALL but bad sign", rx))
    print()
    print(f"{'class':6s} " + " ".join(f"{k:>10s}" for k in cats) + "      total")
    for c in CL:
        tot = sum(count[c].values())
        print(f"{c:6s} " + " ".join(f"{count[c][k]:>10d}" for k in cats) + f" {tot:>10d}")
    print(f"H & H* " + " ".join(f"{hd[k]:>10d}" for k in cats) + f" {sum(hd.values()):>10d}")
    print(f"# NONE certificates: {none_kinds['order']} by a strictly lower bad order, {none_kinds['vertex']} by a "
          f"Newton-polytope vertex at a tie")
    print()
    names = {"F": "factorisation Gamma = a b^T exact", "L": "structural lemma (net consumers contain the species)",
             "A": "ADLS(N*) (no critical siphon) => ALL", "B": "drain species in N* => NONE",
             "W": "WF-EC N => ALL (Cor. P)", "K": "SOME certificates re-verified exactly",
             "T": "sampled exact g near the ends agrees with leading order (and ALL never has a bad sign)"}
    ok = True
    for k in "FLABWKT":
        st = "PASS" if not fails[k] else "FAIL"
        ok &= not fails[k]
        print(f"{st} {k}: {names[k]}  ({napp[k]} applicable, {len(fails[k])} failures)")
        for fx in fails[k][:3]:
            print(f"     failure: {fx}")
    print()
    print("# smallest example per class and category (by n + r, then r)")
    for c in CL:
        for cat in cats:
            if (c, cat) not in smallest:
                continue
            key, rc, res = smallest[(c, cat)]
            n, rx = rc["n"], rc["rx"]
            print(f"[{c} / {cat}] n={n} r={len(rx)} WR={rc['N']['WR']} def={rc['N']['deficiency']} "
                  f"H={e2lib.H(rc)} H*={e2lib.Hd(rc)} a={res['a']} b={[str(x) for x in res['b']]}")
            print(f"   N : {e2lib.describe(n, rx)}")
            print(f"   N*: {e2lib.describe(len(rx), e2lib.dual_rx(n, rx), [f'k{j+1}' for j in range(len(rx))])}")
            if cat == "NONE":
                print(f"   certificate: {res['cert']}  (side +1: species k_j with b_j<0 vanish; J 0-based)")
            if cat in ("SOME", "UNDECIDED"):
                ce = res["cert"]
                if cat == "SOME":
                    print(f"   kappa1 (persistent) = {[str(x) for x in ce['kappa1']]}")
                print(f"   kappa2 (non-persistent) = {[str(x) for x in ce['kappa2']]}; tie side {ce['tie']['side']}, "
                      f"J = {ce['tie']['J']}")
    # explicit generic SOME example (constructed by hand, outside the families: molecularity 4-5)
    ex = (((1, 0, 1, 1), (2, 1, 0, 0)), ((1, 2, 1, 1), (0, 1, 2, 2)), ((1, 1, 0, 2), (0, 0, 1, 3)))
    resx = analyse(4, ex)
    okx = resx["cat"] == "SOME" and verify_certs(4, ex, resx)
    Nx = CNet(4, ex)
    print(f"{'PASS' if okx else 'FAIL'} S: constructed example is SOME (generic: persistent for an open set of rates, "
          f"not for another), certificates verified; N (H1)={Nx.is_conservative()} (H2)={Nx.is_consistent()} "
          f"s={Nx.rank()} WR={Nx.weakly_reversible()}")
    print(f"   N : {Nx.describe()}")
    print(f"   N*: {e2lib.describe(3, e2lib.dual_rx(4, ex), ['k1', 'k2', 'k3'])}")
    print(f"   category {resx['cat']}, ties {[(T['side'], T['J']) for T in resx['ties']]}")
    if resx["cat"] == "SOME":
        print(f"   kappa1 (persistent) = {[str(x) for x in resx['cert']['kappa1']]}, kappa2 (not) = "
              f"{[str(x) for x in resx['cert']['kappa2']]}")
    ok &= okx
    print()
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
