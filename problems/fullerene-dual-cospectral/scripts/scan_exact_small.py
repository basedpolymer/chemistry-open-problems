#!/usr/bin/env python3
"""scan_exact_small.py -- exhaustive EXACT scan of small n (pure python-flint, no hashing).

For every C_n isomer produced by buckygen (triangulation = dual T_n), compute the integer characteristic
polynomials of the adjacency matrices of T_n, T_n^6 (induced on degree-6 vertices = hexagons) and T_n^5
(induced on degree-5 vertices = pentagons), and report, per n, the number of isomers, the number of
spectra shared by >= 2 isomers ("non-unique spectra", as in Bille-Buchstaber-Spodarev Table 1) and the
number of isomers involved.

Usage: python3 scan_exact_small.py BUCKYGEN_BINARY nmin nmax [outfile.json]
Environment: python3.13, python-flint 0.9. Runtime: ~1 min for n <= 80.
"""
import json
import subprocess
import sys
import time
from collections import defaultdict
import flint
sys.path.insert(0, __import__('os').path.dirname(__file__))
from planarcode import read_planar_code, to_graph6


def charpoly_of(adj, verts):
    idx = {v: i for i, v in enumerate(verts)}
    k = len(verts)
    M = [[0] * k for _ in range(k)]
    for v in verts:
        for w in adj[v]:
            if w in idx:
                M[idx[v]][idx[w]] = 1
    if k == 0:
        return (1,)
    return tuple(int(c) for c in flint.fmpz_mat(M).charpoly().coeffs())


def main():
    bg, nmin, nmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    out = sys.argv[4] if len(sys.argv) > 4 else None
    results = {}
    for n in range(nmin, nmax + 1, 2):
        if n == 22:
            continue
        t0 = time.time()
        data = subprocess.run([bg, '-q', '%dd' % n], capture_output=True, check=True).stdout
        cls = {'T': defaultdict(list), 'T6': defaultdict(list), 'T5': defaultdict(list)}
        count = 0
        for i, adj in enumerate(read_planar_code(data)):
            count += 1
            m = len(adj)
            assert m == n // 2 + 2
            V5 = [v for v in range(m) if len(adj[v]) == 5]
            V6 = [v for v in range(m) if len(adj[v]) == 6]
            cls['T'][charpoly_of(adj, list(range(m)))].append(i)
            cls['T6'][charpoly_of(adj, V6)].append(i)
            cls['T5'][charpoly_of(adj, V5)].append(i)
        r = {'n': n, 'isomers': count}
        for key in ('T', 'T6', 'T5'):
            groups = [g for g in cls[key].values() if len(g) > 1]
            r[key + '_nonunique_spectra'] = len(groups)
            r[key + '_isomers_in_nonunique'] = sum(len(g) for g in groups)
            if key in ('T', 'T6'):
                r[key + '_groups'] = groups
        r['seconds'] = round(time.time() - t0, 2)
        results[n] = r
        print(n, count, 'T:', r['T_nonunique_spectra'], 'T6:', r['T6_nonunique_spectra'],
              'T5:', r['T5_nonunique_spectra'], '(isomers', r['T5_isomers_in_nonunique'], ')',
              'groups T:', r['T_groups'][:10], 'T6:', r['T6_groups'][:10], r['seconds'], 's', flush=True)
    if out:
        with open(out, 'w') as f:
            json.dump(results, f, indent=1)


if __name__ == '__main__':
    main()
