"""Tests of agent Q2's candidates (round2/preclusters/candidates.md) on first-round families Q2 did not run, with an
independent implementation of kappa_V.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic.  Imports Q2's library
round2/preclusters/scripts/pclust.py (unchanged) only for its lower bound pi_plus (pi+ <= CF, Q2's Lemma 3.4).
Usage (from round2/computation/): python scripts/q2_tests.py [--limit K] [--nproc P] FILE...
  FILE = record name in ~/gacE-data (first round) or a path (e.g. ~/gacE2-data/sample_5_3_5_8_s3.pkl.gz).

Independent kappa_V (definition of candidates.md): A_V = Gamma diag(V) Q (n x m, Q = reaction-to-source indicator);
classes c ~ c' iff e_c - e_c' lies in the row space of A_V, computed here as: u_c = u_c' for every u in ker A_V
(row space = (ker)^perp).  V = positive integer combinations of the flux rays (generic points of C0).
Claims, on every side X (X = N, and X = N* when N has (H1)) with (H2):
  T1 (Q2-1, soundness of pi+)  pi+(X) refines kappa_V(X) for 2 random V;
  T2 (Q2-5)  X EC, connected, (H1), (H2)  =>  kappa_V(X) is one block (2 random V);
  T3 (Q2-6)  X weakly reversible with deficiency 0  =>  kappa_V(X) = linkage classes (2 random V);
  T3b (report, Q2's K5 check)  same hypothesis  =>  pi+(X) = linkage classes (CF determined by Q2's lower bound);
  T4 (Q2-8)  dim ker A_V <= number of blocks of kappa_V;
  T5 (Q2-3)  M* (full coupling, Q2's full_coupling) refines pi (Q2's flat_closure), and pi = M* iff kappa_V = M*
             (generic V) iff every doubling pair is fully coupled;
  T6 (Q2-10) N with (H1), CF(N*) determined (pi+ = meet of 3 kappa_V): every block of CF(N*) lies in a minimal
             siphon of N iff N has no twin pair and every ray-MIPS of N lies in a minimal siphon.
Reports:
  R12 (Q2-12) N weakly reversible: CF(N) merges two linkage classes (certified by pi+ not refining LC) / CF refines
      LC (certified by the meet of 3 kappa_V refining LC) / undetermined; smallest certified merge, its deficiency;
  R13 (Q2-13) N with (H1): every minimal siphon of N lies in one block of CF(N*): holds (inside a pi+(N*) block),
      fails (inside no block of the meet of 3 kappa_V(N*)), undetermined; counts for H and WR, smallest failures.
"""
import os
import sys
import time
import random
import platform
from fractions import Fraction
from multiprocessing import Pool

import e2lib
from e2lib import CNet, crn, ecore

QDIR = os.path.normpath(os.path.join(e2lib.HERE, "..", "..", "preclusters", "scripts"))
sys.path.insert(0, QDIR)
import pclust as pc  # noqa: E402


def norm(P):
    return tuple(sorted(tuple(sorted(b)) for b in P))


def refines(P, Q):
    return all(any(set(b) <= set(c) for c in Q) for b in P)


def meet(P, Q):
    out = []
    for b in P:
        for c in Q:
            x = set(b) & set(c)
            if x:
                out.append(x)
    return norm(out)


def sources(X):
    srcs = []
    cl = []
    for t, (a, b) in enumerate(X.reactions):
        if a not in srcs:
            srcs.append(a)
        cl.append(srcs.index(a))
    return srcs, cl


def my_kappa(X, V):
    srcs, cl = sources(X)
    m = len(srcs)
    G = X.Gamma
    A = [[sum(Fraction(V[t]) * G[p][t] for t in range(X.r) if cl[t] == c) for c in range(m)] for p in range(X.n)]
    K = crn.nullspace(A, m)
    key = [tuple(u[c] for u in K) for c in range(m)]
    groups = {}
    for t in range(X.r):
        groups.setdefault(key[cl[t]], []).append(t)
    return norm(groups.values()), len(K)


def rand_flux(X, rng):
    rays = X.flux_rays()
    V = [Fraction(0)] * X.r
    for v in rays:
        w = rng.randint(1, 10 ** 6)
        for t in range(X.r):
            V[t] += w * v[t]
    return V


def lc_part(X):
    cs = X.complexes()
    idx = {c: i for i, c in enumerate(cs)}
    parent = list(range(len(cs)))

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for a, b in X.reactions:
        parent[find(idx[a])] = find(idx[b])
    g = {}
    for t, (a, b) in enumerate(X.reactions):
        g.setdefault(find(idx[a]), []).append(t)
    return norm(g.values())


def side(n, rx, rng, want_pi=True):
    X = CNet(n, rx)
    o = {"fails": []}
    if not X.is_consistent():
        return None
    Ks = []
    for _ in range(3):
        K, dk = my_kappa(X, rand_flux(X, rng))
        Ks.append(K)
        if dk > len(K):
            o["fails"].append("T4")
    o["meet"] = Ks[0]
    for K in Ks[1:]:
        o["meet"] = meet(o["meet"], K)
    pip = norm(pc.pi_plus(X)) if want_pi else None
    o["pi"] = pip
    if pip is not None and not all(refines(pip, K) for K in Ks):
        o["fails"].append("T1")
    ec = e2lib.is_EC(rx) and e2lib.petri_connected(n, rx) and X.is_conservative()
    o["ec"] = ec
    if ec and not all(len(K) == 1 for K in Ks):
        o["fails"].append("T2")
    # T5 (Q2-3): M* refines pi; pi = M* iff kappa_V = M* (generic V) iff every doubling pair is fully coupled
    Ms = norm(pc.full_coupling(X))
    fl = norm(pc.flat_closure(X))
    if not refines(Ms, fl):
        o["fails"].append("T5")
    dbl_ok = all(any(t in b and u in b for b in Ms) for t in range(X.r) for u in range(t + 1, X.r)
                 if rx[t][0] == rx[u][0])
    if not ((fl == Ms) == dbl_ok == (Ks[0] == Ms)):
        o["fails"].append("T5")
    o["WR"] = X.weakly_reversible()
    o["def"] = X.deficiency()
    o["LC"] = lc_part(X)
    if o["WR"] and o["def"] == 0:
        if not all(K == o["LC"] for K in Ks):
            o["fails"].append("T3")
        if pip is not None and pip != o["LC"]:
            o["fails"].append("T3b")
    return o


def work(args):
    n, rx, seed, h, wr = args
    rng = random.Random(seed)
    res = {"n": n, "rx": rx, "H": h, "WRN": wr}
    N = CNet(n, rx)
    res["N"] = side(n, rx, rng)
    if N.is_conservative():
        res["D"] = side(len(rx), e2lib.dual_rx(n, rx), rng)
        # R13: minimal siphons of N vs blocks of CF(N*)
        D = res["D"]
        if D is not None:
            mins = [ecore.to_set(S) for S in ecore.minimal_masks(ecore.siphon_masks(N))]
            st = "holds"
            for Sg in mins:
                if any(set(Sg) <= set(b) for b in D["pi"]):
                    continue
                if not any(set(Sg) <= set(b) for b in D["meet"]):
                    st = "fails"
                    break
                st = "undet"
            res["R13"] = st
            # T6 (Q2-10), when CF(N*) is determined by the bounds (pi+ = meet of kappa_V)
            if D["pi"] == D["meet"]:
                ms = [set(ecore.to_set(S)) for S in ecore.minimal_masks(ecore.siphon_masks(N))]
                lhs = all(any(set(b) <= S for S in ms) for b in D["pi"])
                rows = [tuple(a[i] for a, b in rx) for i in range(n)]
                twins = any(rows[i] == rows[j] and any(rows[i]) for i in range(n) for j in range(i + 1, n))
                mips = [set(b) for b in N.mips_ray()]
                rhs = (not twins) and all(any(B <= S for S in ms) for B in mips)
                res["T6"] = (lhs == rhs)
    return res


def main():
    t0 = time.time()
    argv = sys.argv[1:]
    limit = int(argv[argv.index("--limit") + 1]) if "--limit" in argv else 3000
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 4
    files = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or not argv[i - 1].startswith("--"))]
    print(f"# q2_tests.py  Python {platform.python_version()} on {platform.system()}; exact; limit {limit} per file; "
          f"files {[os.path.basename(f) for f in files]}")
    rng = random.Random(2)
    jobs = []
    for fn in files:
        path = os.path.expanduser(fn) if "/" in fn else fn
        recs = e2lib.load([path])
        rng.shuffle(recs)
        for rc in recs[:limit]:
            jobs.append((rc["n"], rc["rx"], rng.randint(0, 10 ** 9), e2lib.H(rc), rc["N"]["WR"]))
    print(f"# {len(jobs)} networks")
    with Pool(nproc) as pool:
        res = pool.map(work, jobs, chunksize=10)
    fails = {k: [] for k in ("T1", "T2", "T3", "T4", "T3b", "T5", "T6")}
    napp = {k: 0 for k in fails}
    r12 = {"merge": [], "refines": 0, "undet": 0}
    r13 = {}
    r13_fail = {"H": [], "WR": []}
    for q in res:
        for sd in ("N", "D"):
            o = q.get(sd)
            if o is None:
                continue
            napp["T1"] += 1
            napp["T4"] += 1
            napp["T5"] += 1
            napp["T2"] += bool(o["ec"])
            napp["T3"] += bool(o["WR"] and o["def"] == 0)
            napp["T3b"] += bool(o["WR"] and o["def"] == 0)
            for f in o["fails"]:
                fails[f].append((sd, q["n"], q["rx"]))
        if "T6" in q:
            napp["T6"] += 1
            if not q["T6"]:
                fails["T6"].append(("N", q["n"], q["rx"]))
        o = q["N"]
        if o is not None and o["WR"]:
            if not refines(o["pi"], o["LC"]):
                r12["merge"].append(q)
            elif refines(o["meet"], o["LC"]):
                r12["refines"] += 1
            else:
                r12["undet"] += 1
        if "R13" in q:
            for cl in ("all",) + (("H",) if q["H"] else ()) + (("WR",) if q["H"] and q["WRN"] else ()):
                r13.setdefault(cl, {}).setdefault(q["R13"], 0)
                r13[cl][q["R13"]] += 1
            if q["R13"] == "fails" and q["H"]:
                r13_fail["H"].append(q)
                if q["WRN"]:
                    r13_fail["WR"].append(q)
    ok = True
    names = {"T1": "pi+ refines kappa_V (Q2-1 soundness)", "T2": "EC+connected+(H1)+(H2) => kappa_V one block (Q2-5)",
             "T3": "WR, deficiency 0 => kappa_V = linkage classes (Q2-6)", "T4": "dim ker A_V <= #kappa_V (Q2-8)",
             "T3b": "(report) WR, deficiency 0 => pi+ = linkage classes (Q2's K5 form; not part of Q2-6)",
             "T5": "M* refines pi, and pi = M* iff kappa_V = M* iff doubling pairs fully coupled (Q2-3)",
             "T6": "blocks of CF(N*) in minimal siphons iff no twins and every ray-MIPS in a minimal siphon (Q2-10; "
                   "CF determined cases)"}
    for k in ("T1", "T2", "T3", "T4", "T5", "T6", "T3b"):
        if k != "T3b":
            ok &= not fails[k]
        lab = ("PASS" if not fails[k] else "FAIL") if k != "T3b" else ("HOLDS" if not fails[k] else "REFUTED")
        print(f"{lab} {k}: {names[k]}  ({len(fails[k])} failures over {napp[k]} sides)")
        for f in fails[k][:2]:
            print(f"     failure side {f[0]}: {e2lib.describe(f[1], f[2])}")
    key = lambda q: (q["n"] + len(q["rx"]), len(q["rx"]), q["rx"])  # noqa: E731
    r12["merge"].sort(key=key)
    defs = {}
    for q in r12["merge"]:
        defs[q["N"]["def"]] = defs.get(q["N"]["def"], 0) + 1
    print(f"# R12 (WR networks): CF merges two linkage classes (certified) {len(r12['merge'])} (by deficiency "
          f"{dict(sorted(defs.items()))}), CF refines LC (certified) {r12['refines']}, undetermined {r12['undet']}")
    for q in r12["merge"][:3]:
        print(f"     merge: {e2lib.describe(q['n'], q['rx'])}  (deficiency {q['N']['def']}, H={q['H']}, "
              f"pi+ = {q['N']['pi']}, LC = {q['N']['LC']})")
    print(f"# R13 (every minimal siphon of N inside one block of CF(N*)): {r13}")
    for cl in ("H", "WR"):
        r13_fail[cl].sort(key=key)
        for q in r13_fail[cl][:2]:
            print(f"     fails ({cl}): N = {e2lib.describe(q['n'], q['rx'])}")
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
