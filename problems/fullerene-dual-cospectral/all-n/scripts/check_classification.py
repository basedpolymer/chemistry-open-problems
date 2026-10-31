#!/usr/bin/env python3
"""check_classification.py -- every certified cospectral pair with n <= 138 is explained (PASS/FAIL, exact).

Input: ../../certificates/cospectral_pairs_T.json (890 classes of the exhaustive scan) and
../../certificates/tube_family_members.json (cap rotations of Theorem 2 of ../../result.md).

For each class {G1, G2}:
  * clique partitions of G1, G2: perfect matchings M of the pentagons by edges of T such that the edges of T are
    partitioned by M and by one colour class of the triangles of T - M;
  * if a partition exists on both sides, H = (intersection graph of the lines) is computed for each, and the class is a
    COLOUR TWIN when some H1 is isomorphic to some H2 (nauty) and the two colour classes of that H give back G1, G2;
  * otherwise the class must be one of the cap rotations of Theorem 2.
Claims:
  K1  882 classes are colour twins, each with exact identity N_b^T N_b = N_w^T N_w;
  K2  the 8 others are exactly the cap-rotation pairs with n = 44, 50, 68, 74, 92, 98, 116, 122;
  K3  no pair of a class has cospectral T^6 (exact charpoly).
Runtime about 4 min.  Environment as in check_family.py.
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from twin import *

CERT = os.path.join(ROOT, 'certificates')
FAILED = []


def report(name, ok, detail=''):
    print('%s  %s%s' % ('PASS' if ok else 'FAIL', name, ('  ' + detail) if detail else ''), flush=True)
    if not ok:
        FAILED.append(name)


def incid(rows, V):
    return flint.fmpz_mat(len(rows), V, [1 if v in f else 0 for f in rows for v in range(V)])


def main():
    pairs = json.load(open(os.path.join(CERT, 'cospectral_pairs_T.json')))
    tube = json.load(open(os.path.join(CERT, 'tube_family_members.json')))
    tubepairs = {}
    for r in tube:
        if r['nonisomorphic']:
            tubepairs.setdefault(r['n'], []).append((r['graph6_delta0'], r['graph6_delta1']))
    ntwin = 0
    others = []
    idok = True
    t6ok = True
    for c in pairs:
        a1 = from_graph6(c['graph6_dual'][0]); a2 = from_graph6(c['graph6_dual'][1])
        r1 = rotation_system(a1); r2 = rotation_system(a2)
        m = len(a1)
        C1 = clique_partitions(r1); C2 = clique_partitions(r2)
        hex1 = [v for v in range(m) if len(a1[v]) == 6]
        hex2 = [v for v in range(m) if len(a2[v]) == 6]
        t6ok &= charpoly(a1, hex1) != charpoly(a2, hex2)
        found = False
        if C1 and C2:
            H1 = [line_graph(L) for M, L in C1]
            H2 = [line_graph(L) for M, L in C2]
            can = labelg_canon(H1 + H2)
            s1, s2 = can[:len(H1)], can[len(H1):]
            for i, x in enumerate(s1):
                if x in s2:
                    rh = rotation_system(H1[i])
                    F = faces_of(rh)
                    col = two_colour_faces(F)
                    t0, t1 = colour_class_graph(F, col, 0), colour_class_graph(F, col, 1)
                    cc = labelg_canon([t0, t1, a1, a2])
                    if {cc[0], cc[1]} == {cc[2], cc[3]} and cc[2] != cc[3]:
                        V = len(H1[i])
                        Nb = incid([f for f, k in zip(F, col) if k == 0], V)
                        Nw = incid([f for f, k in zip(F, col) if k == 1], V)
                        idok &= (Nb.transpose() * Nb) == (Nw.transpose() * Nw)
                        found = True
                        break
        if found:
            ntwin += 1
        else:
            others.append(c)
    report('K1 colour twins: %d of %d classes, with N_b^T N_b = N_w^T N_w' % (ntwin, len(pairs)), ntwin == 882 and idok)
    ok = sorted(c['n'] for c in others) == [44, 50, 68, 74, 92, 98, 116, 122]
    for c in others:
        cc = labelg_canon([from_graph6(c['graph6_dual'][0]), from_graph6(c['graph6_dual'][1])])
        match = False
        for g0, g1 in tubepairs.get(c['n'], []):
            c2 = labelg_canon([from_graph6(g0), from_graph6(g1)])
            match |= set(c2) == set(cc)
        ok &= match
    report('K2 the remaining classes are the cap rotations of Theorem 2, n = %s' % sorted(c['n'] for c in others), ok)
    report('K3 in no class are the two T^6 graphs cospectral', t6ok)
    if FAILED:
        sys.exit(1)


if __name__ == '__main__':
    main()
