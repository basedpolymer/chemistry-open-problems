#!/usr/bin/env python3
"""forcing_study.py -- Agent C.  Role of the forced impurities (N = 3, alpha = (4,2,1), z = 1/3, q = 1, delta = 0.05,
prefractionator slack 0.025).

Part A (limit system vs eta): for eta_A = eta_C = eta in a list, the smallest singular value of DG_inf(v_inf) (my
variables), of its major and trace blocks, and the pinned junction ratios xi = A/(A+B) at J[1,2], zeta = B/(B+C) at
J[2,3].  Expected: sigma_min > 0 for every eta > 0 and -> 0 as eta -> 0 (the unforced limit system is singular).
Part B (finite n, no forcing): at fixed n the forced solution is continued to eta = 0 (same chart); sigma_min(DG_n)
at eta = 0 is compared with eta = 0.003 as n grows.  Expected (network.md 10.4): geometric decay without forcing,
bounded below with forcing.
Prints PASS/FAIL; writes runs/forcing_study.json.  Numerical evidence only.
Usage: python forcing_study.py
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env).  Runtime: about 2 minutes.
"""
import os
import sys
import json
import time
import mpmath as mp
from petnet import Problem, Design, limit_config, System, svals, classify, submatrix, impurities, ratio_form_rows

HERE = os.path.dirname(os.path.abspath(__file__))
fails = []
out = {}
t0 = time.time()
mp.mp.dps = 60
pb = Problem([4, 2, 1], [1, 1, 1], 1)


def design(eta):
    return Design(pb, 0.05, eta_nb_default=0.025, eta_bin={(0, 1): ('B', eta), (1, 2): ('T', eta)})


print("Part A: limit Jacobian versus the forcing eta (eta_A = eta_C = eta)")
rowsA = []
for eta in ['0.1', '0.03', '0.01', '0.003', '0.001', '1e-4', '1e-5', '1e-6', '1e-8']:
    des = design(mp.mpf(eta))
    lim = limit_config(des)
    sl = System(des, lim, None)
    v = sl.start()
    cl = classify(sl, v)
    sv = svals(cl['J'])
    sm = svals(submatrix(cl['J'], cl['mr'], cl['mc']))[-1]
    st = svals(submatrix(cl['J'], cl['tr'], cl['tc']))[-1]
    J1, J3 = lim['X'][(0, 1)], lim['X'][(1, 2)]
    xi, ze = J1[0] / (J1[0] + J1[1]), J3[1] / (J3[1] + J3[2])
    sr = svals(ratio_form_rows(sl, v, cl['J']))[-1]
    rowsA.append(dict(eta=eta, sigma_min=float(sv[-1]), sigma_min_ratio_form=float(sr), major=float(sm), trace=float(st),
                      xi=float(xi), zeta=float(ze)))
    print(f"  eta {eta:>6}: sigma_min(DG_inf) product form {mp.nstr(sv[-1], 5):>11}, ratio form (L) {mp.nstr(sr, 5):>9}; "
          f"major block {mp.nstr(sm, 5):>11}  trace block {mp.nstr(st, 5):>11}  xi {mp.nstr(xi, 8)}  zeta {mp.nstr(ze, 8)}")
ok = all(r['sigma_min'] > 0 for r in rowsA) and rowsA[-1]['sigma_min'] < rowsA[0]['sigma_min']
ratios = [rowsA[k]['sigma_min'] / float(rowsA[k]['eta']) for k in range(len(rowsA))]
print(f"  sigma_min/eta: {[f'{x:.4g}' for x in ratios]}")
print(f"[{'PASS' if ok else 'FAIL'}] DG_inf invertible for every eta > 0 tested, sigma_min -> 0 as eta -> 0")
if not ok:
    fails.append('A')
out['A'] = rowsA

print("Part B: finite n, sigma_min(DG_n) with forcing (eta = 0.003) and without (eta = 0)")
print("  (unforced solution: Agent N's stage solver stagenet.py (double precision, uniform start), then polished by")
print("   Newton in the junction variables; both solutions re-checked stage by stage; a continuation eta -> 0 in the")
print("   chart variables is NOT used: it lands on non-physical zeros of G_n, see scripts/spurious_zeros.py)")
import numpy as np
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'total', 'network', 'scripts'))
from ftc_fixed_flows import build_ftc
from stagenet import solve_lse
from petnet import reconstruct, stage_check
rowsB = []
for n in [8, 16, 24, 32, 40, 48]:   # n = 64: the unforced system is singular to ~1e-16 and a double-precision start is useless
    mp.mp.dps = max(60, int(40 + 0.8 * n))
    pbn = Problem([4, 2, 1], [1, 1, 1], 1)
    globals()['pb'] = pbn
    des = design(mp.mpf('0.003'))
    lim = limit_config(des)
    sn = System(des, lim, n)
    v, nF, _ = sn.newton(sn.start())
    chk = stage_check(sn, v, reconstruct(sn, v))
    phys_f = (not chk['negs']) and chk['minS'] > 0
    sf = svals(sn.jacobian(v))[-1]
    # unforced: stagenet
    net, P, H, C = build_ftc([4.0, 2.0, 1.0], [1 / 3, 1 / 3, 1 / 3], 1.0, n, 0.05)
    Xa, ra = net.solve(maxit=100)
    xi_, rel = solve_lse(net, np.log(np.maximum(Xa, 1e-300)))
    Xs = np.exp(xi_)
    nm = lambda s_: net.names.index(s_)
    al = np.array([4.0, 2.0, 1.0])
    def netflow(sec, bu):
        V = [f for d, f in net.vap[bu] if d == nm(f"{sec}.0")][0]
        L = [f for d, f in net.liq[nm(f"{sec}.0")] if d == bu][0]
        y = al * Xs[bu] / (al @ Xs[bu])
        return V * y - L * Xs[nm(f"{sec}.0")]
    umap = {(0, 1): 'J01', (0, 2): 'J02', (1, 2): 'J12', (0, 0): 'P0', (1, 1): 'P1', (2, 2): 'P2'}
    Xd = {u: [mp.mpf(float(x)) for x in Xs[nm(umap[u])]] for u in pbn.units}
    Wd = {}
    for s_ in pbn.secs:
        bu, tu = pbn.sec_ends(s_)
        Wd[s_] = [mp.mpf(float(x)) for x in netflow(str(s_), nm(umap[bu]))]
    s0 = System(design(mp.mpf(0)), lim, n)
    v0, nF0, _ = s0.newton(s0.pack(Xd, Wd))
    chk0 = stage_check(s0, v0, reconstruct(s0, v0))
    phys_u = (not chk0['negs']) and chk0['minS'] > 0 and nF0 < mp.mpf(10) ** (-(mp.mp.dps - 15))
    su = svals(s0.jacobian(v0))[-1]
    X0, W0 = s0.unpack(v0)
    J1 = X0[(0, 1)]
    imp = impurities(s0, v0)
    rowsB.append(dict(n=n, forced=float(sf), forced_physical=phys_f, unforced=float(su), unforced_physical=phys_u,
                      stagenet_residual=float(ra), xi_unforced=float(J1[0] / (J1[0] + J1[1])),
                      imp_unforced=[float(x) for x in imp]))
    print(f"  n {n:3d}: sigma_min forced {mp.nstr(sf, 5):>10} ({'physical' if phys_f else 'NOT physical'})   "
          f"unforced {mp.nstr(su, 5):>10} ({'physical' if phys_u else 'NOT physical'}; xi at J[1,2] = "
          f"{mp.nstr(J1[0] / (J1[0] + J1[1]), 6)}, impurities {[mp.nstr(x, 3) for x in imp]})", flush=True)
okB = all(r['forced_physical'] and r['unforced_physical'] for r in rowsB)
f_ok = min(r['forced'] for r in rowsB[2:]) > 0.5 * rowsB[-1]['forced']
u_dec = rowsB[-1]['unforced'] < 1e-4 * rowsB[-1]['forced'] and all(
    rowsB[k + 1]['unforced'] < rowsB[k]['unforced'] for k in range(1, len(rowsB) - 1))
ok = okB and f_ok and u_dec
print(f"[{'PASS' if ok else 'FAIL'}] forced: sigma_min bounded below; unforced: sigma_min decreasing geometrically (physical solutions)")
if not ok:
    fails.append('B')
out['B'] = rowsB
with open(os.path.join(HERE, '..', 'runs', 'forcing_study.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
sys.exit(1 if fails else 0)
