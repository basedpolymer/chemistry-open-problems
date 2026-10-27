"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: independent re-count of Q1's rank-one table
(round2/persistence/notes.md, section 2) with an independent generator and an independent face classifier
(v3lib.face_types).  Does not import any author code.

Family (n, m, R): weakly reversible networks with s = 1 and (H1)-(H3), complexes of molecularity <= m (zero complex
excluded: under (H1) it cannot share a linkage class with a nonzero complex), <= R reactions, disjoint union of
strongly connected linkage classes on pairwise distinct complexes, no self-loops or parallel reactions, every species
occurring; counted up to permutation of the species (complete canonical form: minimum over all permutations).

Claims: PASS if the counts (networks, ALL, NONE, NONE-with-an-OUT-face, DOM-only) equal Q1's table; R5 = no DOM face
in any dual (Theorem S1-WR).  Usage: python v3_rank1_enum.py n m R  (exit 1 on FAIL).
Environment: Python 3.12, standard library only, exact arithmetic.
"""
import sys
import time
import platform
from itertools import combinations, product, permutations
from math import gcd

import v3lib as L

Q1_TABLE = {  # (n, m, R): (networks, ALL, NONE, of which OUT face, some-not-all)
    (2, 2, 4): (12, 4, 8, 0, 0), (2, 3, 6): (460, 115, 345, 267, 0), (2, 4, 6): (3544, 381, 3163, 1833, 0),
    (2, 5, 5): (2742, 662, 2080, 1601, 0), (3, 2, 6): (12, 4, 8, 0, 0), (3, 3, 6): (426, 101, 325, 130, 0),
    (3, 4, 5): (829, 245, 584, 291, 0), (4, 3, 4): (1, 1, 0, 0, 0)}


def complexes(n, m):
    return [c for c in product(range(m + 1), repeat=n) if 0 < sum(c) <= m]


def directions(n, m):
    out = []
    for a in product(range(-m, m + 1), repeat=n):
        if not (any(x > 0 for x in a) and any(x < 0 for x in a)):
            continue
        if len(set(a)) != n:
            continue
        g = 0
        for x in a:
            g = gcd(g, abs(x))
        if g != 1:
            continue
        if next(x for x in a if x != 0) < 0:
            continue
        out.append(a)
    return out


def sc_digraphs(k, R):
    """All labelled strongly connected simple digraphs on vertices 0..k-1 with at most R edges (k >= 2)."""
    arcs = [(u, v) for u in range(k) for v in range(k) if u != v]
    out = []
    for e in range(k, min(R, len(arcs)) + 1):
        for E in combinations(arcs, e):
            # strong connectivity: reach from 0 forward and backward
            fw, bw = {0}, {0}
            ch = True
            while ch:
                ch = False
                for u, v in E:
                    if u in fw and v not in fw:
                        fw.add(v); ch = True
                    if v in bw and u not in bw:
                        bw.add(u); ch = True
            if len(fw) == k and len(bw) == k:
                out.append(E)
    return out


def main():
    n, m, R = map(int, sys.argv[1:4])
    t0 = time.time()
    print("# v3_rank1_enum.py n=%d m=%d R=%d  Python %s on %s; exact" % (n, m, R, platform.python_version(),
                                                                       platform.system()))
    C = complexes(n, m)
    scd = {k: sc_digraphs(k, R) for k in range(2, R + 1)}
    nets = set()
    for a in directions(n, m):
        lines = {}
        for y in C:
            key = tuple(y[i] * a[k] - y[k] * a[i] for i in range(n) for k in range(i + 1, n))
            lines.setdefault(key, []).append(y)
        inst = []  # (frozenset of complexes, tuple of reactions, #edges)
        for pts in lines.values():
            for k in range(2, min(len(pts), R) + 1):
                for S in combinations(pts, k):
                    for E in scd[k]:
                        inst.append((frozenset(S), tuple((S[u], S[v]) for u, v in E), len(E)))
        inst.sort(key=lambda x: x[2])

        def rec(start, used, rx, ne):
            if rx:
                spec = [any(y[i] for s, p in rx for y in (s, p)) for i in range(n)]
                if all(spec):
                    nets.add(L.canonical(n, list(rx)))
            for idx in range(start, len(inst)):
                S, E, e = inst[idx]
                if ne + e > R:
                    break
                if S & used:
                    continue
                rec(idx + 1, used | S, rx + E, ne + e)
        rec(0, frozenset(), (), 0)
    nets = sorted(nets)
    cnt = {'ALL': 0, 'NONE': 0, 'DOM-only': 0}
    out_face = 0
    dom_any = 0
    bad_h = 0
    smallest = {}
    for key in nets:
        Rx = [tuple(r) for r in key]
        if not (L.is_wr(Rx) and L.rank(L.Gamma(n, Rx)) == 1 and L.conservative(n, Rx) and L.H3(n, Rx)):
            bad_h += 1
            continue
        nd, Rd = L.dual(n, Rx)
        v, ft = L.verdict_rank_one(nd, Rd)
        cnt[v] += 1
        if v == 'NONE' and any(t == 'OUT' for _, _, t in ft):
            out_face += 1
        if any(t == 'DOM' for _, _, t in ft):
            dom_any += 1
        if v not in smallest or (len(Rx), sum(map(sum, (s for s, p in Rx)))) < smallest[v][0]:
            smallest[v] = ((len(Rx), sum(map(sum, (s for s, p in Rx)))), L.show(n, Rx))
    total = len(nets) - bad_h
    print("networks %d (generator sanity failures %d): ALL %d, NONE %d (with an OUT face %d), DOM-only %d; "
          "duals with a DOM face: %d" % (total, bad_h, cnt['ALL'], cnt['NONE'], out_face, cnt['DOM-only'], dom_any))
    for v, (k_, s_) in sorted(smallest.items()):
        print("  smallest %s: %s" % (v, s_))
    fails = 0
    exp = Q1_TABLE.get((n, m, R))
    if exp:
        got = (total, cnt['ALL'], cnt['NONE'], out_face, cnt['DOM-only'])
        ok = got == exp
        print(("PASS" if ok else "FAIL") + ": counts equal Q1's table %s (got %s)" % (exp, got))
        fails += 0 if ok else 1
    ok = dom_any == 0 and bad_h == 0
    print(("PASS" if ok else "FAIL") + ": R5 no DOM face in any dual (Theorem S1-WR) and generator sanity")
    fails += 0 if ok else 1
    print("# runtime %.1f s" % (time.time() - t0))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
