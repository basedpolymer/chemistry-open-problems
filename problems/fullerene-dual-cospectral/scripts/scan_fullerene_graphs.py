#!/usr/bin/env python3
"""scan_fullerene_graphs.py -- side question (Fowler 2005: "As yet, no cospectral fullerenes are known"): are there
two non-isomorphic fullerene graphs C_n (the carbon skeletons, cubic, n vertices; Hueckel spectrum) with the
same adjacency spectrum?  Exhaustive and complete, same method as scan_driver.py: buckygen -d writes the
fullerene graphs, dualspecG16 (GENERIC mode, traces mod 2^16 up to K) hashes A(C_n), every key collision is
settled with the exact integer characteristic polynomial (python-flint).

Usage: python3 scan_fullerene_graphs.py NMIN NMAX [K]   -> certificates/fullerene_graph_scan.json (appended)
Environment: as scan_driver.py.  Runtime: seconds for n <= 80, ~20 min for n = 110.
"""
import json
import os
import subprocess
import sys
import time

import flint

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from planarcode import read_planar_code  # noqa: E402
from scan_driver import BUCKYGEN, PCSELECT, SRC, oeis_counts  # noqa: E402

BIN = os.path.join(SRC, 'build', 'dualspecG16')
OUT = os.path.join(ROOT, 'certificates', 'fullerene_graph_scan.json')


def cp(adj):
    m = len(adj)
    return tuple(int(c) for c in flint.fmpz_mat(m, m, [1 if w in adj[v] else 0 for v in range(m)
                                                        for w in range(m)]).charpoly().coeffs())


def main():
    nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    K = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    oe = oeis_counts()
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for n in range(nmin, nmax + 1, 2):
        if n == 22 or str(n) in res:
            continue
        t0 = time.time()
        p1 = subprocess.Popen(['nice', '-n', '10', BUCKYGEN, '-q', '-d', '%dd' % n], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
        p2 = subprocess.run(['nice', '-n', '10', BIN, '64', str(K)], stdin=p1.stdout, capture_output=True)
        p1.stdout.close()
        p1.wait()
        assert p2.returncode == 0
        lines = p2.stdout.decode().splitlines()
        coll = [tuple(map(int, l.split()[2:4])) for l in lines if l.startswith('C ')]
        summ = [l for l in lines if 'SUMMARY' in l][-1]
        graphs = int(summ.split('graphs=')[1].split()[0])
        classes = []
        if coll:
            idx = sorted({i for c in coll for i in c})
            fn = os.path.join(SRC, 'scanraw', 'fg%03d_idx.txt' % n)
            with open(fn, 'w') as f:
                f.write('\n'.join(map(str, idx)) + '\n')
            q1 = subprocess.Popen([BUCKYGEN, '-q', '-d', '%dd' % n], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            q2 = subprocess.run([PCSELECT, fn], stdin=q1.stdout, capture_output=True)
            q1.stdout.close()
            q1.wait()
            sel = dict(zip(idx, read_planar_code(q2.stdout)))
            groups = {}
            for i, j in coll:
                groups.setdefault(i, {i}).add(j)
            for first, mem in groups.items():
                parts = {}
                for x in sorted(mem):
                    parts.setdefault(cp(sel[x]), []).append(x)
                classes += [v for v in parts.values() if len(v) > 1]
        res[str(n)] = {'n': n, 'isomers': graphs, 'oeis': oe.get(n), 'K': K, 'key_collisions': len(coll),
                       'cospectral_classes_exact': classes, 'seconds': round(time.time() - t0, 1)}
        print(n, graphs, graphs == oe.get(n), 'collisions', len(coll), 'cospectral classes', classes, flush=True)
        with open(OUT, 'w') as f:
            json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
