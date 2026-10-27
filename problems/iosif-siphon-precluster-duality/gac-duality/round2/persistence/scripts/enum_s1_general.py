"""Enumerate ALL networks N (not necessarily weakly reversible) with rank s = 1, (H1), (H2) and optionally (H3), and
classify the persistence of N* by the face analysis of s1lib (agent Q1, round 2).

Environment: Python 3.12+ (run with Windows py -3, 3.14), standard library only, exact arithmetic.  Imports ecore.py (canonical form) and s1lib.py.
Usage: python enum_s1_general.py N_SPECIES MAX_MOLECULARITY MAX_REACTIONS [--noH3]

Reactions: y -> y + b a with a primitive (both signs: (H1)), complexes of molecularity <= m, b != 0; reaction sets
of size 2..R containing both signs of b ((H2) for rank one: some v > 0 with sum b_j v_j = 0).  (H3): entries of a
pairwise distinct.  One representative per orbit of species permutations.
Checks (PASS/FAIL): G1 the dual has rank 1 and is conservative; G2 Theorem S1-WR: no DOM face occurs in the dual
of a weakly reversible network; G3 Corollary S1-pure: if no reaction has a net product in its source, the dual is
persistent for all rates (ALL); G4 Corollary S1-two-species (n = 2): ALL iff every source contains strictly more of
the species its reaction consumes than of the species it produces; G5 no DOM face under the weaker hypothesis (E1)
(in every linkage class the two complexes extreme along a are sources); reports verdict counts (ALL / SOME / NONE / TWO-SIDED) split by WR / not WR, and the
smallest network with a DOM face (dual persistent for some but not all rate constants).
"""
import os
import sys
import time
import platform
from itertools import combinations, product
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "computation", "scripts")))
import ecore  # noqa: E402
from s1lib import classify, dual, show  # noqa: E402
from enum_s1 import primitive_dirs  # noqa: E402


def dirs_any(n, m, H3):
    out = []
    for a in product(range(-m, m + 1), repeat=n):
        if not (any(x > 0 for x in a) and any(x < 0 for x in a)):
            continue
        if H3 and len(set(a)) != n:
            continue
        G = 0
        for x in a:
            G = gcd(G, abs(x))
        if G != 1:
            continue
        if next(x for x in a if x != 0) < 0:
            continue
        out.append(a)
    return out


def main():
    t0 = time.time()
    args = sys.argv[1:]
    n, m, R = int(args[0]), int(args[1]), int(args[2])
    H3 = "--noH3" not in args
    print(f"# enum_s1_general.py n={n} m={m} R={R} H3={H3}")
    print(f"# Python {platform.python_version()} on {platform.system()} ({platform.machine()})")
    C = ecore.complexes(n, m, False)
    Cset = set(C)
    seen = set()
    count = {}
    smallest = {}
    fails = {"G1": 0, "G2": 0, "G3": 0, "G4": 0, "G5": 0}
    ne1 = 0
    for a in dirs_any(n, m, H3):
        rxs = []
        for y in C:
            for b in range(-m, m + 1):
                if b == 0:
                    continue
                yp = tuple(yi + b * ai for yi, ai in zip(y, a))
                if yp in Cset:
                    rxs.append((y, yp, b))
        for k in range(2, R + 1):
            for S in combinations(rxs, k):
                if not (any(t[2] > 0 for t in S) and any(t[2] < 0 for t in S)):
                    continue
                rx = [(t[0], t[1]) for t in S]
                if not all(any(yy[i] or zz[i] for yy, zz in rx) for i in range(n)):
                    continue
                can = ecore.canonical(n, rx)
                if can in seen:
                    continue
                seen.add(can)
                rx = list(can)
                N = ecore.CNet(n, rx)
                wr = N.weakly_reversible()
                p, drx = dual(n, rx)
                verdict, data = classify(p, drx)
                if verdict in ("not-rank-1", "not-conservative"):
                    fails["G1"] += 1
                    continue
                g, c, fa = data
                hasdom = any(f["type"] == "DOM" for f in fa)
                # G3 (Corollary S1-pure): no reaction has a net product in its source  =>  ALL
                noauto = all(src[i] == 0 for src, prd in rx for i in range(n) if prd[i] > src[i])
                if noauto and verdict != "ALL":
                    fails["G3"] += 1
                    print("FAIL G3", show(n, rx))
                # G4 (Corollary S1-two-species): for n = 2, ALL <=> every source has strictly more of its consumed
                # species than of its produced species
                if n == 2:
                    crit = all(src[[i for i in range(2) if prd[i] < src[i]][0]] > src[[i for i in range(2) if prd[i] > src[i]][0]] for src, prd in rx)
                    if crit != (verdict == "ALL"):
                        fails["G4"] += 1
                        print("FAIL G4", show(n, rx), verdict)
                # G5 (Theorem S1-WR under the weaker hypothesis (E1)): in every linkage class the two complexes
                # that are extreme along a are sources  =>  no DOM face
                srcs = {u for u, v in rx}
                comp = {}
                for u, v in rx:
                    for w_ in (u, v):
                        comp.setdefault(w_, w_)
                def find(x):
                    while comp[x] != x:
                        comp[x] = comp[comp[x]]
                        x = comp[x]
                    return x
                for u, v in rx:
                    comp[find(u)] = find(v)
                classes = {}
                for w_ in list(comp):
                    classes.setdefault(find(w_), []).append(w_)
                k0 = next(i for i in range(n) if a[i] != 0)
                e1 = all(min(L, key=lambda y: y[k0] * a[k0]) in srcs and max(L, key=lambda y: y[k0] * a[k0]) in srcs
                         for L in classes.values())
                if e1 and hasdom:
                    fails["G5"] += 1
                    print("FAIL G5", show(n, rx))
                if e1:
                    ne1 += 1
                if wr and hasdom:
                    fails["G2"] += 1
                    print("FAIL G2", show(n, rx))
                key = (verdict, wr)
                count[key] = count.get(key, 0) + 1
                size = (len(rx), sum(sum(u) + sum(v) for u, v in rx))
                tag = verdict + (" WR" if wr else " nonWR") + (" withDOM" if hasdom else "")
                if tag not in smallest or size < smallest[tag][0]:
                    smallest[tag] = (size, rx, [(f["Z"], f["sigma"], f["type"]) for f in fa if f["type"] != "IN"])
    print(f"# networks: {len(seen)}; satisfying (E1) (extreme complexes of every linkage class are sources): {ne1}")
    for (v, wr), x in sorted(count.items()):
        print(f"{v:10s} {'WR' if wr else 'nonWR':6s} {x}")
    for tag, (size, rx, fl) in sorted(smallest.items()):
        print(f"# smallest {tag}: N = {show(n, rx)} | N* = {show(*dual(n, rx), names=[f'k{j+1}' for j in range(len(rx))])} | non-IN faces {fl}")
    for k, v in fails.items():
        print(f"{k}: {'PASS' if v == 0 else 'FAIL'} ({v} failures)")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if any(fails.values()) else 0)


if __name__ == "__main__":
    main()
