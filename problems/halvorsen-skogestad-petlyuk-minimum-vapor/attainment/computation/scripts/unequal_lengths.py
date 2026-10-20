#!/usr/bin/env python3
"""unequal_lengths.py -- Agent C.  Analysis Theorem A.5 claims that every section may have its own number of stages
n_sigma, only min n_sigma >= n_0 being required.  Test: for N = 3 (alpha = (4,2,1), z = 1/3, q = 1, pattern BT) and
N = 4 (alpha = (8,4,2,1), z = 1/4, q = 1, pattern BBT) solve the exact finite-stage system with
 (a) random lengths n_sigma in [n, 4n] (5 seeds), (b) one section at n and all others at 8n (every section in turn),
from the limit configuration; check convergence, every stage balance, positivity, S > 0, the distance to the limit
(should be governed by min n_sigma) and sigma_min(DG_n).
Prints PASS/FAIL; writes runs/unequal_lengths.json.  Numerical evidence only.
Usage: python unequal_lengths.py
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env).  Runtime: several minutes.
"""
import os
import sys
import json
import time
import random
import mpmath as mp
from petnet import (Problem, Design, limit_config, System, reconstruct, stage_check, impurities, svals,
                    solve_audited)

HERE = os.path.dirname(os.path.abspath(__file__))
fails = []
out = []
t0 = time.time()


def setup(N):
    if N == 3:
        pb = Problem([4, 2, 1], [1, 1, 1], 1)
        des = Design(pb, 0.05, eta_nb_default=0.025, eta_bin={(0, 1): ('B', 0.003), (1, 2): ('T', 0.003)})
    else:
        pb = Problem([8, 4, 2, 1], [1, 1, 1, 1], 1)
        des = Design(pb, 0.05, eta_nb_default=0.01,
                     eta_bin={(0, 1): ('B', 0.003), (1, 2): ('B', 0.003), (2, 3): ('T', 0.003)})
    return pb, des, limit_config(des)


def run(N, nst, label):
    nmax = max(nst.values())
    nmin = min(nst.values())
    res = solve_audited(lambda: setup(N), nst, int(60 + 0.9 * nmax))
    sn, v, nF, chk = res['sn'], res['v'], res['nF'], res['chk']
    vl = sn.start()
    dmaj = max(abs(v[j] - vl[j]) for j in range(len(v)) if vl[j] != 0)
    sm = svals(sn.jacobian(v))[-1]
    imp = impurities(sn, v)
    ok = res['phys']
    print(f"  [{'PASS' if ok else 'FAIL'}] N={N} {label}: min n {nmin}, max n {nmax}, |G_n| {mp.nstr(nF, 2)}, "
          f"stage rel {mp.nstr(chk['worst_rel'], 2)}, min x {mp.nstr(chk['minx'], 3)}, min S {mp.nstr(chk['minS'], 4)}, "
          f"|v_n - v_inf| {mp.nstr(dmaj, 3)}, sigma_min {mp.nstr(sm, 4)}, impurities {[mp.nstr(x, 3) for x in imp]}",
          flush=True)
    out.append(dict(N=N, label=label, nmin=nmin, nmax=nmax, ok=ok, dist=float(dmaj), sigma_min=float(sm),
                    impurities=[float(x) for x in imp]))
    if not ok:
        fails.append(f"N={N} {label}")


for N, nbase in [(3, 24), (4, 24)]:
    mp.mp.dps = 60
    pb, des, lim = setup(N)
    secs = list(pb.secs)
    run(N, {s: nbase for s in secs}, f"uniform n = {nbase}")
    for seed in range(5):
        rnd = random.Random(seed)
        run(N, {s: rnd.randint(nbase, 4 * nbase) for s in secs}, f"random lengths in [{nbase}, {4 * nbase}], seed {seed}")
    for s0 in secs:
        run(N, {s: (nbase if s == s0 else 8 * nbase) for s in secs}, f"{s0} at {nbase}, others at {8 * nbase}")
with open(os.path.join(HERE, '..', 'runs', 'unequal_lengths.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
sys.exit(1 if fails else 0)
