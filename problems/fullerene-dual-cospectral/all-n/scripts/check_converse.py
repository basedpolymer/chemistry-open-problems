#!/usr/bin/env python3
"""check_converse.py -- cross-check of Theorem A against the scan, n <= NMAX (default 100).

Goes through ALL isomers (buckygen) and finds those with a clique partition (perfect matching M of the pentagons by edges
of T, all merged quadrilaterals of T - M in one colour class).  For each such T it builds H, the two colour classes and
the twin T'.  Claims:
  X1  every T with a clique partition has a twin T' which is again one of the isomers, and T' is cospectral with T;
  X2  the set of unordered pairs {T, T'} with T != T' equals the set of A(T)-cospectral classes of
      ../../certificates/cospectral_pairs_T.json for that n, minus the cap rotations of n = 44, 50, 68, 74, 92, 98.
Prints per n: isomers, isomers with a clique partition, self-twins, twin pairs.  Runtime about 20 min for NMAX = 100.
"""
import collections
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from twin import *

CERT = os.path.join(ROOT, 'certificates')


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    pairs = json.load(open(os.path.join(CERT, 'cospectral_pairs_T.json')))
    bad = False
    for n in [20] + list(range(24, nmax + 1, 2)):
        p = subprocess.run(['buckygen', '-q', '%dd' % n], capture_output=True, check=True)
        hits = []
        total = 0
        for rot in read_planar_code(p.stdout):
            total += 1
            m = len(rot)
            pent = [v for v in range(m) if len(rot[v]) == 5]
            ps = set(pent)
            # every pentagon needs a pentagon neighbour
            if not all(any(w in ps for w in rot[v]) for v in pent):
                continue
            C = clique_partitions(rot)
            if C:
                hits.append((rot, C))
        twins_found = set()
        selft = 0
        allg = []
        recs = []
        for rot, C in hits:
            for M, L in C:
                h = line_graph(L)
                rh = rotation_system(h)
                F = faces_of(rh)
                col = two_colour_faces(F)
                t0, t1 = colour_class_graph(F, col, 0), colour_class_graph(F, col, 1)
                recs.append((rot, t0, t1))
        if recs:
            can = labelg_canon([g for r in recs for g in r])
            for i in range(len(recs)):
                c, c0, c1 = can[3 * i], can[3 * i + 1], can[3 * i + 2]
                if c not in (c0, c1):
                    bad = True
                    print('FAIL: T is not a colour class of its own H, n =', n)
                if c0 == c1:
                    selft += 1
                else:
                    twins_found.add(frozenset((c0, c1)))
            for rot, t0, t1 in recs:
                if charpoly(t0) != charpoly(t1):
                    bad = True
                    print('FAIL: colour classes not cospectral, n =', n)
        scan = [c for c in pairs if c['n'] == n]
        scan_set = set()
        for c in scan:
            cc = labelg_canon([from_graph6(c['graph6_dual'][0]), from_graph6(c['graph6_dual'][1])])
            scan_set.add(frozenset(cc))
        missing = twins_found - scan_set
        extra = scan_set - twins_found
        ok = not missing and len(extra) == (1 if n in (44, 50, 68, 74, 92, 98) else 0)
        bad |= not ok
        print('n=%d isomers=%d with_partition=%d self_twin_records=%d twin_pairs=%d scan_classes=%d %s'
              % (n, total, len(hits), selft, len(twins_found), len(scan), 'PASS' if ok else 'FAIL'), flush=True)
    print('PASS  X1, X2 up to n = %d' % nmax if not bad else 'FAIL  X1/X2')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
