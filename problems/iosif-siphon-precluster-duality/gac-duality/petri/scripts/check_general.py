"""Theorems A, B and Proposition H of petri/notes.md on general (not necessarily EC) small networks (agent P).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rationals (crn.py, petri.py).
Usage: python check_general.py N MAXMOL RMAX [NPROC].  Runtime printed at the end.

Networks: the manuscript's enumeration families (enumerate_small.py: distinct non-loop reactions between nonzero
complexes of molecularity <= MAXMOL, every species occurring, up to species permutation), restricted to (H1), (H2).
Claims (PASS/FAIL, exit 1 on failure):
  A   every minimal siphon contains at most one species of each same-consumer class, every minimal T-trap at most
      one reaction of each same-reactant-set class;
  B   Q3 species side  <=> no twin pair and every species in a minimal siphon; reaction side likewise;
  H   every precluster of N* with >= 2 elements is a union of ray-MIPS blocks of N; every precluster of N with
      >= 2 elements is a union of fully coupled classes (MIPS of N*);
  G   minimal siphons = supports of extreme conservation laws  ==>  every minimal siphon is a union of MIPS blocks
      and every MIPS block lies in a minimal siphon.
Reported: K8 (is P(N*) determined by the twin classes and the minimal siphons? groups with several values),
counts with (H3), counts of K9(iv) "every support of an extreme conservation law is a union of minimal
siphons", and of K6 "(no critical siphon in N) == (no critical siphon in N*)".
"""
import os
import sys
import time
import platform
from itertools import combinations, permutations
from multiprocessing import Pool

from petri import Network, critical_siphons, refines, twin_classes, ec_sets

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                                                 "scripts")))
from enumerate_small import complexes, canonical  # noqa: E402


def union_of_blocks(B, part):
    return all((b <= B) or not (b & B) for b in part)


def evaluate(N):
    D = N.dual()
    cl, rep = {}, {}
    S, Sd = N.minimal_siphons(), D.minimal_siphons()
    cons = {}
    for i in range(N.n):
        cons.setdefault(frozenset(j for j in range(N.r) if N.reactions[j][0][i] > 0), set()).add(i)
    reac = {}
    for j in range(N.r):
        reac.setdefault(frozenset(i for i in range(N.n) if N.reactions[j][0][i] > 0), set()).add(j)
    cl["A"] = all(len(R & c) <= 1 for R in S for c in cons.values()) and \
        all(len(U & c) <= 1 for U in Sd for c in reac.values())
    P, Pd = N.preclusters(), D.preclusters()
    twins = [b for b in twin_classes(N) if len(b) >= 2]
    dbl = [b for b in ec_sets(N) if len(b) >= 2]
    qa = refines(Pd, S)
    qb = refines(P, Sd)
    ea = (not twins) and bool(S) and set().union(*S) == set(range(N.n))
    eb = (not dbl) and bool(Sd) and set().union(*Sd) == set(range(N.r))
    cl["B"] = (qa == ea) and (qb == eb)
    M, Md = N.mips_ray(), D.mips_ray()
    cl["H"] = all(union_of_blocks(b, M) for b in Pd if len(b) >= 2) and \
        all(union_of_blocks(b, Md) for b in P if len(b) >= 2)
    sup = sorted({frozenset(i for i in range(N.n) if z[i] > 0) for z in N.conservation_rays()},
                 key=lambda x: (len(x), sorted(x)))
    if sup == S:
        cl["G"] = all(union_of_blocks(R, M) for R in S) and all(any(b <= R for R in S) for b in M)
    rep["K9iv"] = all(set().union(*([R for R in S if R <= X] or [set()])) == set(X) for X in sup)
    rep["K6"] = (not critical_siphons(N)) == (not critical_siphons(D))
    rep["H3"] = N.rows_distinct()
    key = lambda part: tuple(sorted(tuple(sorted(b)) for b in part))
    rep["_K8A"] = ((N.n, key(twins), key(S)), key(Pd))
    rep["_K8B"] = ((N.r, key(dbl), key(Sd)), key(P))
    return cl, rep


def worker(args):
    n, m, rmax, first = args
    C = complexes(n, m)
    R = sorted((a, b) for a in C for b in C if a != b)
    perms = list(permutations(range(n)))
    out = []
    rest = R[first + 1:]
    for k in range(0, rmax):
        for tail in combinations(rest, k):
            rx = (R[first],) + tail
            if not all(any(a[i] or b[i] for a, b in rx) for i in range(n)):
                continue
            if canonical(n, rx, perms) != rx:
                continue
            N = Network(n, rx)
            if not N.is_consistent() or not N.is_conservative():
                continue
            out.append((rx,) + evaluate(N))
    return out


def main():
    n, m, rmax = (int(a) for a in sys.argv[1:4])
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    t0 = time.time()
    print(f"# check_general.py n={n} m={m} r<={rmax}; Python {platform.python_version()} on {platform.system()}")
    C = complexes(n, m)
    nR = len(C) * (len(C) - 1)
    tot, h3 = 0, 0
    claims, reps, first = {}, {}, {}
    k8 = {}
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, [(n, m, rmax, f) for f in range(nR)]):
            for rx, cl, rep in chunk:
                tot += 1
                h3 += rep["H3"]
                for k, v in cl.items():
                    a, b = claims.get(k, (0, 0))
                    claims[k] = (a + v, b + (not v))
                    if not v and k not in first:
                        first[k] = rx
                if rep["H3"]:
                    for kk in ("_K8A", "_K8B"):
                        g, val = rep[kk]
                        k8.setdefault(kk, {}).setdefault(g, {}).setdefault(val, rx)
                for k in ("K9iv", "K6"):
                    a, b = reps.get(k, (0, 0))
                    reps[k] = (a + rep[k], b + (not rep[k]))
                    if not rep[k] and (k not in first or len(rx) < len(first[k])):
                        first[k] = rx
    print(f"# networks with (H1), (H2): {tot}; of which (H3): {h3}")
    bad = 0
    for k in sorted(claims):
        a, b = claims[k]
        bad += b > 0
        print(f"{'PASS' if b == 0 else 'FAIL'}  {k}: {a}/{a + b}" + (f"  first failure {first[k]}" if b else ""))
    for k in sorted(reps):
        a, b = reps[k]
        print(f"REPORT  {k}: holds {a}/{a + b}" + (f"  smallest failure {first[k]}" if b else ""))
    for kk, side in (("_K8A", "species side: (n, twin classes, minimal siphons) -> P(N*)"),
                     ("_K8B", "reaction side: (r, doubling classes, minimal T-traps) -> P(N)")):
        groups = k8.get(kk, {})
        amb = [(g, v) for g, v in groups.items() if len(v) > 1]
        print(f"REPORT  K8 {side}: {len(groups)} groups, {len(amb)} with two different preclusterings (H3 class)")
        for g, v in amb[:2]:
            print(f"REPORT  K8 example {kk}: key={g}; networks {list(v.values())[:2]}")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
