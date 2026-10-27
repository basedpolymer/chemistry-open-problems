"""Enumerate weakly reversible networks N with stoichiometric rank s = 1 and (H1)-(H3), and classify the
persistence of the dual N* for all / some / no rate constants (agent Q1, round 2).

Environment: Python 3.12+ (run with Windows py -3, 3.14), standard library only, exact arithmetic.  Imports ../../../computation/scripts/ecore.py
(canonical form, strongly connected patterns) without editing it.
Usage: python enum_s1.py N_SPECIES MAX_MOLECULARITY MAX_REACTIONS [--zero] [--out FILE]

Generation.  For s = 1 all reaction vectors are multiples of one primitive vector a; the complexes of every linkage
class lie on one line y0 + Z a.  For every primitive a with entries of both signs (H1) and pairwise distinct
entries (H3, since the rows of Gamma are a_i b^T), every set of >= 2 complexes of molecularity <= m on one such
line, and every labelled strongly connected digraph on it, we form linkage classes; networks are disjoint unions of
linkage classes (distinct complexes) with <= R reactions in which every species occurs.  (H2) holds because the
network is weakly reversible.  One representative per orbit of the species permutations (ecore.canonical).

Claims checked (PASS/FAIL, exit status 1 on FAIL):
  E1  every generated network is WR, has rank 1 and satisfies (H1)-(H3);
  E2  the dual N* has rank 1 and is conservative (it is consistent too);
  E3  face formula: for every face of N* and 3 seeded random (rates, boundary point) pairs, the lowest-order
      coefficient of the exact endpoint expansion of F equals eps^mu * L with L given by the face formula, and the
      endpoint behaviour agrees with the face type (IN repelling, OUT attracting);
  E5  Theorem S1-WR: no face of the dual of a WR network of rank one is of type DOM;
  E4  verdict consistency: ALL <=> every face IN (by construction) and no drain-species of E (endo.drain_species)
      occurs in an ALL dual (a drain species is an OUT singleton face).
Reports: counts of verdicts ALL / SOME / NONE / TWO-SIDED by (n, r, molecularity), face-type statistics, the
smallest example of each verdict.
"""
import os
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import combinations, product
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "computation", "scripts")))
import ecore  # noqa: E402
import endo  # noqa: E402
import s1lib  # noqa: E402
from s1lib import classify, dual, show, faces, face_analysis, endpoint_expansion, lowest, leading_L  # noqa: E402


def primitive_dirs(n, m):
    """Primitive a in Z^n with both signs, pairwise distinct entries, |a_i| <= m, first nonzero entry > 0
    (a and -a give the same lines)."""
    out = []
    for a in product(range(-m, m + 1), repeat=n):
        if not (any(x > 0 for x in a) and any(x < 0 for x in a)):
            continue
        if len(set(a)) != n:
            continue
        G = 0
        for x in a:
            G = gcd(G, abs(x))
        if G != 1:
            continue
        first = next(x for x in a if x != 0)
        if first < 0:
            continue
        out.append(a)
    return out


def line_key(y, a):
    n = len(y)
    return tuple(y[i] * a[k] - y[k] * a[i] for i in range(n) for k in range(i + 1, n))


def lc_instances_s1(n, m, R, zero):
    """Per direction a: list of linkage-class instances (mask over the global complex list, reactions)."""
    C = ecore.complexes(n, m, zero)
    pats = {k: ecore.sc_patterns(k, R) for k in range(2, R + 1)}
    per_dir = []
    for a in primitive_dirs(n, m):
        inst = []
        lines = {}
        for idx, y in enumerate(C):
            lines.setdefault(line_key(y, a), []).append(idx)
        for L in lines.values():
            if len(L) < 2:
                continue
            for k in range(2, min(R, len(L)) + 1):
                for sub in combinations(L, k):
                    mask = 0
                    for t in sub:
                        mask |= 1 << t
                    for pat in pats[k]:
                        inst.append((mask, tuple((C[sub[i]], C[sub[j]]) for i, j in pat)))
        inst.sort(key=lambda t: len(t[1]))
        if inst:
            per_dir.append((a, inst))
    return per_dir


def networks(n, m, R, zero):
    seen = set()
    for a, inst in lc_instances_s1(n, m, R, zero):
        for first in range(len(inst)):
            for rx in ecore.wr_networks_from(first, inst, R, n):
                can = ecore.canonical(n, rx)
                if can in seen:
                    continue
                seen.add(can)
                yield can


def main():
    t0 = time.time()
    args = sys.argv[1:]
    n, m, R = int(args[0]), int(args[1]), int(args[2])
    zero = "--zero" in args
    print(f"# enum_s1.py n={n} m={m} R={R} zero={zero}")
    print(f"# Python {platform.python_version()} on {platform.system()} ({platform.machine()})")
    rng = random.Random(20261027)
    fails = {"E1": 0, "E2": 0, "E3": 0, "E4": 0, "E5": 0}
    singgap = []
    nout = 0
    ndrain = 0
    count = {}
    ftype = {}
    smallest = {}
    examples = []
    total = 0
    for rx in networks(n, m, R, zero):
        rx = list(rx)
        N = ecore.CNet(n, rx)
        # rank 1 within one linkage class is guaranteed, but two classes on different directions cannot occur
        # (we only combine instances; check the global rank)
        if N.rank() != 1:
            continue
        if not (N.weakly_reversible() and N.is_conservative() and N.is_consistent() and N.rows_distinct()):
            fails["E1"] += 1
            continue
        total += 1
        p, drx = dual(n, rx)
        verdict, data = classify(p, drx)
        if verdict in ("not-rank-1", "not-conservative"):
            fails["E2"] += 1
            print("FAIL E2", show(n, rx))
            continue
        g, c, fa = data
        for f in fa:
            ftype[f["type"]] = ftype.get(f["type"], 0) + 1
            for _ in range(3):
                kappa = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(len(drx))]
                xplus = [Fraction(0) if i in f["Z"] else Fraction(rng.randint(1, 9), rng.randint(1, 9)) for i in range(p)]
                poly = endpoint_expansion(p, drx, g, c, kappa, xplus, f["sigma"])
                k, lead = lowest(poly)
                L = leading_L(p, drx, f, g, c, kappa, xplus)
                ok = True
                if L != 0:
                    ok = (k == f["mu"] and lead == L)
                else:
                    ok = (k is None or k > f["mu"])
                if f["type"] == "IN":
                    ok = ok and f["sigma"] * L < 0
                if f["type"] == "OUT":
                    ok = ok and f["sigma"] * L > 0
                if not ok:
                    fails["E3"] += 1
                    if fails["E3"] < 5:
                        print("FAIL E3", show(n, rx), f["Z"], f["sigma"], k, lead, L)
        sing_ok = all(f["type"] == "IN" for f in fa if len(f["Z"]) == 1)
        if sing_ok and verdict != "ALL":
            singgap.append(rx)
        hasdom = any(f["type"] == "DOM" for f in fa)
        if hasdom:
            fails["E5"] += 1
            print("FAIL E5 (DOM face in the dual of a WR network)", show(n, rx))
        if any(f["type"] == "OUT" for f in fa):
            nout += 1
        drains = endo.drain_species(p, drx)
        if drains:
            ndrain += 1
        if verdict == "ALL" and drains:
            fails["E4"] += 1
            print("FAIL E4", show(n, rx))
        r = len(rx)
        key = (verdict, r)
        count[key] = count.get(key, 0) + 1
        size = (r, sum(sum(a) + sum(b) for a, b in rx))
        if verdict not in smallest or size < smallest[verdict][0]:
            smallest[verdict] = (size, rx)
        examples.append((verdict, rx))
    print(f"# networks (WR, s = 1, (H1)-(H3)): {total}")
    verdicts = sorted({v for v, _ in count})
    rs = sorted({r for _, r in count})
    print("# verdict counts by number of reactions r:")
    print("verdict      " + " ".join(f"r={r:<6}" for r in rs) + " total")
    for v in verdicts:
        row = [count.get((v, r), 0) for r in rs]
        print(f"{v:12s} " + " ".join(f"{x:<8}" for x in row) + f" {sum(row)}")
    print("# face types over all duals:", ftype)
    print(f"# duals with an OUT face (rate-independent extinction in every class reaching it): {nout}; with a drain species (E.4): {ndrain}")
    print(f"# networks whose singleton faces are all IN but which are not ALL (a face with |Z| >= 2 decides): {len(singgap)}")
    if singgap:
        sg = min(singgap, key=lambda rx: (len(rx), sum(sum(a) + sum(b) for a, b in rx)))
        print(f"#   smallest: N = {show(n, sg)}")
    for v, (size, rx) in sorted(smallest.items()):
        print(f"# smallest {v}: N = {show(n, rx)}   N* = {show(*dual(n, rx), names=[f'k{j+1}' for j in range(len(rx))])}")
    for k, v in fails.items():
        print(f"{k}: {'PASS' if v == 0 else 'FAIL'} ({v} failures)")
    out = None
    if "--out" in args:
        out = args[args.index("--out") + 1]
        with open(out, "w") as fh:
            for v, rx in examples:
                fh.write(f"{v}\t{rx}\n")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if any(fails.values()) else 0)


if __name__ == "__main__":
    main()
