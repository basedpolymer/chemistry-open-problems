#!/usr/bin/env python3
"""t6_ds_all_graphs.py -- Conjecture 1(a) of the published version, read literally: "the graph T_n^6 is uniquely
determined by its spectrum" in the usual sense of spectral graph theory (no non-isomorphic graph, among ALL
graphs, is cospectral with it).

For every C_n isomer with n in [NMIN, NMAX] (T_n^6 has k = n/2 - 10 vertices), we enumerate with nauty's geng
all graphs on k vertices with the same number of edges e (cospectral graphs have equal k and e) and count the
graphs with the same integer characteristic polynomial as T_n^6 (python-flint).  If at least two
non-isomorphic graphs (geng outputs one graph per isomorphism class) share the characteristic polynomial,
T_n^6 has a cospectral mate and is not determined by its spectrum.  The mate is printed in graph6 together
with the certificate: equal integer characteristic polynomials and different nauty canonical forms.

Usage: python3 t6_ds_all_graphs.py NMIN NMAX [out.json]
Environment: python3.13, python-flint 0.9, nauty 2.8.8 (geng, labelg), buckygen 1.1. Runtime: ~1-3 min for n <= 40.
"""
import json
import os
import subprocess
import sys

import flint

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from planarcode import read_planar_code, to_graph6, from_graph6  # noqa: E402
from scan_driver import BUCKYGEN, LABELG, SRC, induced_adj  # noqa: E402
from spiral import canonical_spiral  # noqa: E402

GENG = os.path.join(SRC, 'tools', 'nauty2_8_8', 'geng')


def cp_of(adj):
    k = len(adj)
    if k == 0:
        return (1,)
    M = [[0] * k for _ in range(k)]
    for v in range(k):
        for w in adj[v]:
            M[v][w] = 1
    return tuple(int(c) for c in flint.fmpz_mat(M).charpoly().coeffs())


def main():
    nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    out = sys.argv[3] if len(sys.argv) > 3 else None
    report = []
    cache = {}
    for n in range(nmin, nmax + 1, 2):
        if n == 22:
            continue
        data = subprocess.run([BUCKYGEN, '-q', '%dd' % n], capture_output=True, check=True).stdout
        gs = list(read_planar_code(data))
        sp = [canonical_spiral(g) for g in gs]
        order = sorted(range(len(gs)), key=lambda i: sp[i])
        rank = {i: r + 1 for r, i in enumerate(order)}
        nonds = 0
        for i in order:
            g = gs[i]
            V6 = [v for v in range(len(g)) if len(g[v]) == 6]
            h = induced_adj(g, V6)
            k = len(h)
            e = sum(len(a) for a in h) // 2
            cp = cp_of(h)
            key = (k, e)
            if key not in cache:
                if k == 0:
                    cache[key] = {}
                else:
                    r = subprocess.run([GENG, '-q', str(k), '%d:%d' % (e, e)], capture_output=True, check=True)
                    d = {}
                    for line in r.stdout.decode().split():
                        d.setdefault(cp_of(from_graph6(line)), []).append(line)
                    cache[key] = d
            same = cache[key].get(cp, []) if k else []
            mates = []
            if len(same) >= 2:
                own = subprocess.run([LABELG, '-q'], input=(to_graph6(h) + '\n').encode(), capture_output=True,
                                     check=True).stdout.decode().split()[0]
                canon = subprocess.run([LABELG, '-q'], input=('\n'.join(same) + '\n').encode(), capture_output=True,
                                       check=True).stdout.decode().split()
                mates = [s for s, c in zip(same, canon) if c != own]
                assert len(mates) == len(same) - 1, 'T6 must be one of the geng graphs'
                nonds += 1
            report.append({'n': n, 'atlas': rank[i], 'spiral': list(sp[i]), 'T6_vertices': k, 'T6_edges': e,
                           'T6_graph6': to_graph6(h) if k else '', 'charpoly_low_to_high': [str(c) for c in cp],
                           'cospectral_mates_graph6': mates})
            if mates:
                print('n=%d C%d:%d T6=%s (k=%d,e=%d) has %d cospectral mate(s): %s' %
                      (n, n, rank[i], to_graph6(h), k, e, len(mates), ' '.join(mates[:3])), flush=True)
        print('n=%d: %d isomers, %d with T6 not determined by its spectrum (among all graphs)' % (n, len(gs), nonds),
              flush=True)
    if out:
        with open(out, 'w') as f:
            json.dump(report, f, indent=1)


if __name__ == '__main__':
    main()
