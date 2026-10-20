"""Random-parameter counterexample search: for each sampled (alpha, z, q), minimise the total reboiler vapor of
every configuration (multi-start ipopt) and compare with the Petlyuk value.

Environment: WSL Fedora, Python 3, pyomo 6.10, ipopt 3.14; multiprocessing with W workers (default 5).
Runtime: about 1-2 min per parameter point with 4 workers (N = 4, 152 configurations).
Floating point exploration only (not a certificate).

Usage: python search_random.py N dist qmode npoints seed outfile [workers]
  dist  in {uniform, hard, tiny, mixed}
  qmode in {1, 01, sub, sup}  (q = 1; q ~ U[0,1]; q ~ U[1,1.6] subcooled; q ~ U[-0.6,0] superheated)
Optional environment variables: PERMISSIVE=1 (TG literal model), IDS_FILE=<file with configuration ids>.
Output: one JSON line per parameter point: params, Petlyuk reboiler vapor, per-config best ratio V_c/V_P.
"""
from __future__ import annotations

import json
import math
import multiprocessing as mpc
import random
import sys
import time

from hsmodel import all_configurations, petlyuk_value, petlyuk_peaks


def sample(N, dist, qmode, rng):
    def ratios():
        if dist == 'uniform':
            return [math.exp(rng.uniform(math.log(1.05), math.log(10.0))) for _ in range(N - 1)]
        if dist == 'hard':  # mixture of very close and very spread volatilities
            return [rng.choice([math.exp(rng.uniform(math.log(1.01), math.log(1.3))),
                                math.exp(rng.uniform(math.log(3), math.log(30)))]) for _ in range(N - 1)]
        return [math.exp(rng.uniform(math.log(1.02), math.log(20.0))) for _ in range(N - 1)]
    rat = ratios()
    al = [1.0]
    for r in reversed(rat):
        al.insert(0, al[0] * r)
    if dist == 'tiny' or (dist == 'mixed' and rng.random() < 0.5):
        w = [rng.choice([math.exp(rng.uniform(math.log(1e-3), math.log(3e-2))), rng.random() + 0.05])
             for _ in range(N)]
    else:
        w = [rng.expovariate(1.0) + 1e-3 for _ in range(N)]
    s = sum(w)
    z = [x / s for x in w]
    if qmode == '1':
        q = 1.0
    elif qmode == '01':
        q = rng.choice([0.0, 1.0, rng.random(), rng.random()])
    elif qmode == 'sub':
        q = rng.uniform(1.0, 1.6)
    elif qmode == 'sup':
        q = rng.uniform(-0.6, 0.0)
    else:
        raise ValueError(qmode)
    return al, z, q


_CFGS = None


def _init(N):
    global _CFGS
    import os
    _CFGS = all_configurations(N)
    if os.environ.get('IDS_FILE'):  # optional restriction to a list of configuration ids
        ids = set(int(x) for x in open(os.environ['IDS_FILE']).read().split())
        _CFGS = [c for c in _CFGS if c['id'] in ids]


def evaluate_point(args):
    import os
    from nlp import solve_config
    perm = os.environ.get('PERMISSIVE') == '1'  # optional: TG literal (permissive) model
    al, z, q, seed = args
    t0 = time.time()
    vt, vp = petlyuk_value(al, z, q)
    vp = float(vp)
    P = [float(x) for x in petlyuk_peaks(al, z, q)]
    res = {}
    for c in _CFGS:
        best, _ = solve_config(c, al, z, q, starts=2, seed=seed + c['id'], permissive=perm)
        r = best[0] / vp if best else None
        if r is not None and r < 1 + 2e-3:
            best2, _ = solve_config(c, al, z, q, starts=6, seed=seed + 7919 + c['id'], vscale=3.0,
                                    permissive=perm)
            if best2 and best2[0] / vp < r:
                r = best2[0] / vp
        res[c['id']] = r
    return {'alpha': al, 'z': z, 'q': q, 'peaks': P, 'VP': vp, 'ratios': res, 'time': time.time() - t0}


def main():
    N = int(sys.argv[1])
    dist, qmode = sys.argv[2], sys.argv[3]
    npts, seed = int(sys.argv[4]), int(sys.argv[5])
    out = sys.argv[6]
    W = int(sys.argv[7]) if len(sys.argv) > 7 else 5
    rng = random.Random(seed)
    tasks = []
    for t in range(npts):
        al, z, q = sample(N, dist, qmode, rng)
        tasks.append((al, z, q, seed * 100000 + t * 1000))
    with mpc.Pool(W, initializer=_init, initargs=(N,)) as pool, open(out, 'a') as fh:
        for r in pool.imap_unordered(evaluate_point, tasks):
            vals = [v for v in r['ratios'].values() if v is not None]
            r['min_ratio'] = min(vals)
            r['argmin'] = [k for k, v in r['ratios'].items() if v is not None and v <= min(vals) + 1e-7]
            fh.write(json.dumps(r) + '\n')
            fh.flush()
            flag = 'CANDIDATE' if r['min_ratio'] < 1 - 1e-6 else ''
            print(f"min_ratio={r['min_ratio']:.9f} n_argmin={len(r['argmin'])} t={r['time']:.1f} {flag}",
                  flush=True)


if __name__ == '__main__':
    main()
