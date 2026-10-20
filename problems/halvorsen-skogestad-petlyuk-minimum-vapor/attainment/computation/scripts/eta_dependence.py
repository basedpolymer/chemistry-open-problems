#!/usr/bin/env python3
"""eta_dependence.py -- Agent C.  How the finite-stage threshold and constants depend on the forcing eta
(N = 3, alpha = (4,2,1), z = 1/3, q = 1, delta = 0.05, prefractionator slack 0.025, pattern BT, eta_1 = eta_2 = eta).
For eta in 1e-2 .. 1e-6 and n in 16 .. 512: is the zero of G_n found by Newton from v_inf (or by continuation from
the next larger n) a physical steady state; the observed threshold n0(eta); C_n = |v_n - v_inf| / rho^n.
Expected from the analysis: n0(eta) ~ log(1/eta)/log(1/rho) and C bounded for fixed eta, growing like 1/eta.
Prints a table and PASS/FAIL (PASS = for every eta the largest n is physical); writes runs/eta_dependence.json.
Numerical evidence only.  Usage: python eta_dependence.py
Environment: Python 3.12, mpmath 1.3, python-flint (WSL Fedora sci env).  Runtime: a few minutes."""
import os
import sys
import json
import time
import math
import mpmath as mp
from petnet import Problem, Design, limit_config, solve_audited

HERE = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
fails = []
out = []
NS = [16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]
print("eta      | n0 observed | C_n = |v_n - v_inf| / rho^n at n = 64, 128, 256, 512 | product impurities at n = 512")
for eta in ['1e-2', '3e-3', '1e-3', '1e-4', '1e-5', '1e-6']:
    def build(eta=eta):
        pb = Problem([4, 2, 1], [1, 1, 1], 1)
        des = Design(pb, mp.mpf('0.05'), eta_nb_default=mp.mpf('0.025'),
                     eta_bin={(0, 1): ('B', mp.mpf(eta)), (1, 2): ('T', mp.mpf(eta))})
        return pb, des, limit_config(des)
    mp.mp.dps = 60
    pb, des, lim = build()
    rho = max(max([float(lim['roots'][s][lim['kstar'][s]] / lim['roots'][s][j]) for j in range(lim['kstar'][s])] +
                  [float(lim['roots'][s][j] / lim['roots'][s][lim['kstar'][s]]) for j in range(lim['kstar'][s] + 1, 3)])
              for s in pb.secs)
    rows = {}
    vprev = None
    for n in sorted(NS, reverse=True):
        dps = int(40 + 0.8 * (n + 1))
        res = solve_audited(build, n, dps)
        how = 'v_inf'
        if not res['phys'] and vprev is not None:
            res2 = solve_audited(build, n, dps, start=vprev)
            if res2['phys']:
                res, how = res2, 'continued'
        v = res['v']
        vl = res['sn'].start()
        d = max(abs(v[j] - vl[j]) for j in range(len(v)) if vl[j] != 0)
        rows[n] = dict(physical=res['phys'], how=how, dist=float(d), C=float(d) / rho ** n,
                       imp=[float(x) for x in __import__('petnet').impurities(res['sn'], v)])
        if res['phys']:
            vprev = [str(x) for x in v]
    n0 = None
    for n in sorted(NS, reverse=True):
        if rows[n]['physical']:
            n0 = n
        else:
            break
    if n0 is None or not rows[max(NS)]['physical']:
        fails.append(eta)
    out.append(dict(eta=eta, rho=rho, n0=n0, rows={str(k): v for k, v in rows.items()}))
    print(f"{eta:8} | {str(n0):>11} | " + ", ".join(f"{rows[n]['C']:.3g}" for n in [64, 128, 256, 512]) +
          f" | {', '.join(f'{x:.2e}' for x in rows[512]['imp'])}   (rho = {rho:.4f})", flush=True)
with open(os.path.join(HERE, '..', 'runs', 'eta_dependence.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
sys.exit(1 if fails else 0)
