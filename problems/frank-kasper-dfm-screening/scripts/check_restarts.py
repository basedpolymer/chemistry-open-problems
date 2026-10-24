#!/usr/bin/env python3
"""Numerical symmetry-breaking restarts; NOT an optimality certificate.

Python 3.12+, numpy/scipy/numba/gemmi. One process; several minutes.
No crystal symmetry is imposed on perturbed foam vertices or box strains.
The file is flushed after each completed level. Invalid final meshes are marked.
"""
from pathlib import Path
import argparse
import csv
import time

import dfm
import numpy as np
from dfm import registry, mesh, relax, energy
from run_dfm import initial_mesh

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--names', nargs='+', default=['sigma', 'cat:E11', 'H', 'mu'])
    ap.add_argument('--seeds', nargs='+', type=int, default=[1701, 2903])
    ap.add_argument('--maxlevel', type=int, default=2)
    args = ap.parse_args()
    target = ROOT/'data'/'symmetry_breaking_restarts.csv'
    with target.open('w', newline='') as f:
        out = csv.DictWriter(f, fieldnames=['name', 'seed', 'level', 'F', 'min_cone', 'valid_numerically',
                                           'gmax', 'stop_reason', 'seconds'])
        out.writeheader()
        for name in args.names:
            for seed in args.seeds:
                M = initial_mesh(registry.get(name))
                rng = np.random.default_rng(seed)
                X = M.X.copy()
                scale = energy.dfm_values(M)['R0']
                noise = rng.standard_normal(X.shape)*scale*.015
                for _ in range(20):
                    M.X = X+noise
                    if relax.cone_check(M) > 1e-12:
                        break
                    noise *= .5
                else:
                    raise RuntimeError('no valid perturbed start')
                strain = rng.standard_normal((3, 3))*.01
                T = np.eye(3)+(strain+strain.T)/2
                M.L = T@M.L
                M.X = M.X@T.T
                M.ref = M.ref@T.T
                for level in range(args.maxlevel+1):
                    started = time.time()
                    if level:
                        M = mesh.refine(M)
                    infos = relax.relax_cycles(M, cycles=2, maxiter=6000, box=True)
                    value = energy.dfm_values(M)
                    cone = float(relax.cone_check(M))
                    valid = bool(np.isfinite(value['F']) and cone > 0 and np.linalg.det(M.L) > 0)
                    row = {'name': name, 'seed': seed, 'level': level, 'F': f"{value['F']:.12f}",
                           'min_cone': f'{cone:.5e}', 'valid_numerically': valid, 'gmax': infos[-1]['gmax'],
                           'stop_reason': infos[-1]['msg'], 'seconds': round(time.time()-started, 2)}
                    out.writerow(row)
                    f.flush()
                    print('NUMERICAL', row, flush=True)
                    if not valid:
                        break


if __name__ == '__main__':
    main()
