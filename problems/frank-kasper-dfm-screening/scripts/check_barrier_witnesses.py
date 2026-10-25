#!/usr/bin/env python3
"""Certify eight relaxed witnesses, including auxiliary versus centroid cones.

Python 3.12+, numpy, python-flint. About one minute, one process.
Default uses committed integer meshes, catalogue data and the stage CSV. --build rounds
exploratory meshes saved by run_barrier.py to a dyadic grid, then verifies them.
Upper bounds for explicit partitions, never claims of local/global optimality.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import csv
import hashlib
import json
import math

import numpy as np
from audit_certificates import (require, topology_checks, exact_energy,
                                dot, cross, sub, matvec)

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['cat:N55', 'cat:N59', 'cat:N60', 'cat2:Y02', 'cat2:Y03',
         'cat2:Y04', 'cat2:Y05', 'cat2:Y06']


def centroid_check(c):
    """Exact signed cone test at the actual volume centroids (not mesh vertices)."""
    L, X, C = (c[k].astype(object).tolist() for k in ('L', 'X', 'C'))
    n = int(c['ncell'])
    volumes, firsts, cones = [0]*n, [[0]*3 for _ in range(n)], []
    for ids, shifts, f in zip(c['tri'], c['tsh'], c['tface']):
        f = int(f)
        ps = [tuple(a+b for a, b in zip(X[int(i)], matvec(L, tuple(map(int, s)))))
              for i, s in zip(ids, shifts)]
        off = matvec(L, tuple(map(int, c['fsj'][f])))
        for i, sign, shift in ((int(c['fi'][f]), 1, (0, 0, 0)),
                               (int(c['fj'][f]), -1, off)):
            origin = tuple(a+b for a, b in zip(C[i], shift))
            a, b, d = [sub(p, origin) for p in ps]
            determinant = sign*dot(a, cross(b, d))
            require(determinant > 0, 'nonpositive auxiliary cone')
            volumes[i] += determinant
            for k in range(3):
                firsts[i][k] += determinant*(a[k]+b[k]+d[k])
            normal = tuple(sign*(x+y+z) for x, y, z in
                           zip(cross(b, d), cross(d, a), cross(a, b)))
            cones.append((i, determinant, normal))
    minimum, negative = None, 0
    for i, determinant, normal in cones:
        value = Q(determinant)-Q(dot(firsts[i], normal), 4*volumes[i])
        minimum = value if minimum is None else min(minimum, value)
        negative += value < 0
    scaled = minimum*Q(n, 6*abs(dot(L[0], cross(L[1], L[2]))))
    return dict(negative_centroid_cones=negative,
                minimum_centroid_cone_over_mean_cell_volume=str(scaled))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--build', action='store_true')
    args = ap.parse_args()
    data_path = ROOT/'certificates'/'catalogue_topology.npz'
    manifest = json.loads((ROOT/'certificates'/'catalogue_topology_manifest.json').read_text())
    require(hashlib.sha256(data_path.read_bytes()).hexdigest() == manifest['npz_sha256'],
            'catalogue hash mismatch')
    dest = ROOT/'certificates'/'meshes'/'barrier'
    dest.mkdir(exist_ok=True)
    rows = []
    with (ROOT/'data'/'barrier_continuation.csv').open(newline='') as f:
        stages = list(csv.DictReader(f))
    require(len(stages) == 48 and {r['name'] for r in stages} == set(NAMES),
            'expected 48 continuation stages for eight names')
    for name in NAMES:
        selected = [r for r in stages if r['name'] == name]
        require(sorted(Q(r['mu']) for r in selected) == sorted(Q(1, 10**k) for k in range(1, 7)),
                'barrier weight schedule mismatch: '+name)
        require(all(int(r['level']) == 0 and float(r['min_cone']) > 0
                    and math.isfinite(float(r['F'])) for r in selected),
                'invalid saved continuation diagnostic: '+name)
    with np.load(data_path, allow_pickle=False) as topo:
        for name in NAMES:
            stem = name.replace(':', '_')+'_L0.npz'
            path = dest/stem
            if args.build:
                with np.load(ROOT/'work'/'barrier'/stem, allow_pickle=False) as raw:
                    c = {k: (np.rint(raw[k]*2**44).astype(np.int64)
                             if k in ('X', 'C', 'L') else raw[k].astype(np.int64))
                         for k in raw.files}
                c['K'] = np.array(44, dtype=np.int64)
                np.savez_compressed(path, **c)
            with np.load(path, allow_pickle=False) as data:
                c = dict(data)
            require(all(a.dtype.kind in 'iu' for a in c.values()), 'noninteger data')
            require(int(c['K']) == 44, 'wrong grid exponent')
            ident = ('II_' if name.startswith('cat2:') else 'I_')+name.split(':')[1]
            tc = topology_checks(c, name, expected_tetrahedra=topo[ident+'_tets'])
            low, high, geom = exact_energy(c)
            final = next(r for r in stages if r['name'] == name and Q(r['mu']) == Q(1, 10**6))
            require(abs(float(final['F'])-float(high)) < 1e-9,
                    'saved final stage disagrees with rounded witness: '+name)
            cc = centroid_check(c)
            if name == 'cat2:Y04':
                require(cc['negative_centroid_cones'] > 0,
                        'Y04 no longer supplies the claimed centroid counterexample')
            upper_int = -((-high.numerator*10**9)//high.denominator)
            upper = f'{upper_int//10**9}.{upper_int%10**9:09d}'
            row = dict(name=name, status='PASS', file='barrier/'+stem,
                       sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                       mesh_upper_9dp=upper, topology=tc, geometry=geom, **cc)
            rows.append(row)
            print('PASS', name, 'F <=', upper, '; negative centroid cones:',
                  cc['negative_centroid_cones'], flush=True)
    (ROOT/'certificates'/'barrier_witness_report.json').write_text(
        json.dumps(dict(catalogue_sha256=manifest['npz_sha256'],
                        continuation_stages_checked=len(stages), claims=rows), indent=2)+'\n')
    print('PASS all eight witnesses: exact topology, positive auxiliary cones, volume and energy')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('FAIL barrier witness verification:', exc, flush=True)
        raise SystemExit(1)
