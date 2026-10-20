#!/usr/bin/env python3
"""Agent D: cross-check of the limit configuration against independent finite-stage numerics (Agent N,
total/network/network.md §11.3, script total/network/scripts/ftc3_forced.py, floating point).

N's run: alpha = (4, 2, 1), z = 1/3, q = 1, R = max H + 0.05, prefractionator slack = circulation delta/2 = 0.025
added to V and L of both prefractionator sections (= Astra's deformation with eta_[1,3] = 0.025, h = 0), binary
columns forced with eta_A = eta_C = 0.003, pattern (B, T) (A down the bottom of [A,B], C up the top of [B,C]).
In the vapor network of design.md the binary vapors are determined by R and the non-binary vapors, so N's flows are
exactly F(delta = 0.05, eta_[1,3] = 0.025, eta_1 = eta_2 = 0.003, pi = BT).
Reported finite-stage limits (n = 64..96): xi = A/(A+B) at J[1,2] -> 0.696154, zeta = B/(B+C) at J[2,3] -> 0.553931,
side-product impurity -> 1.768e-2.  This script computes the same quantities from the LIMIT configuration and
compares (agreement to the 6 printed digits = PASS).

Usage: python crosscheck_N3.py        Environment: Python 3 + mpmath.  Runtime: < 1 s.
"""
import sys
import mpmath as mp
from design_lib import build_spec, limit_config

mp.mp.dps = 50
S = build_spec([4, 2, 1], [mp.mpf(1) / 3] * 3, 1, mp.mpf('0.05'), {(0, 2): mp.mpf('0.025')},
               {0: mp.mpf('0.003'), 1: mp.mpf('0.003')}, 'BT')
S = limit_config(S)
X1, X2 = S['X'][('J', 0, 1)], S['X'][('J', 1, 2)]
xi = X1[0] / (X1[0] + X1[1])
zeta = X2[1] / (X2[1] + X2[2])
imp = 1 - S['prods'][1]['x'][1]
print(f"limit xi   = {mp.nstr(xi, 12)}   (N, n = 64..96: 0.696154)")
print(f"limit zeta = {mp.nstr(zeta, 12)}   (N, n = 64..96: 0.553931)")
print(f"side-product impurity = {mp.nstr(imp, 8)}   (N: 1.768e-2)")
print(f"R = {mp.nstr(S['R'], 10)}, H = {[mp.nstr(h, 8) for h in S['H']]}")
for key in [('T', 0, 1), ('B', 0, 1), ('T', 1, 2), ('B', 1, 2)]:
    s = S['secs'][key]
    print(key, 'V =', mp.nstr(s['V'], 8), 'L =', mp.nstr(s['L'], 8), 'w =', [mp.nstr(v, 6) for v in s['w']],
          'modes =', [mp.nstr(m['val'], 6) for m in s['modes']], 'pinch =', mp.nstr(s['pinch']['val'], 6))
ok = abs(xi - mp.mpf('0.696154')) < 6e-7 and abs(zeta - mp.mpf('0.553931')) < 6e-7 and abs(imp - mp.mpf('0.01768')) < 6e-6

# identity of the total flows: rebuild N's flows (HS point of total/network/scripts/ftc_fixed_flows.py, imported
# read-only, + N's circulation weights on the two monotone paths reboiler -> condenser, + N's forcing) and compare
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'total', 'network', 'scripts'))
try:
    from ftc_fixed_flows import attainment_point
    mp.mp.dps = 50
    d0, b0, VT0, VB0, H0, C0 = attainment_point([4.0, 2.0, 1.0], [1 / 3, 1 / 3, 1 / 3], 1.0)
    dl, eA, eC = mp.mpf('0.05'), mp.mpf('0.003'), mp.mpf('0.003')
    w = {('T', 0, 2): dl / 2, ('B', 0, 2): dl / 2, ('T', 0, 1): dl, ('B', 0, 1): dl / 2, ('T', 1, 2): dl / 2,
         ('B', 1, 2): dl}
    dev = mp.mpf(0)
    for (kind, i, m), c in w.items():
        if kind == 'T':
            V = VT0[(i, m)] + c
            L = VT0[(i, m)] - mp.fsum(d0[(i, m)].values()) + c
        else:
            V = VB0[(i, m)] + c
            L = VB0[(i, m)] + mp.fsum(b0[(i, m)].values()) + c
        if (i, m) == (0, 1):
            L += eA
        if (i, m) == (1, 2):
            L -= eC
        s = S['secs'][(kind, i, m)]
        dev = max(dev, abs(V - s['V']), abs(L - s['L']))
    # N's attainment_point solves its feed roots by float bisection, so agreement is to ~1e-15
    print(f"max |V, L (N's flows) - V, L (F(0.05, (0.025; 0.003, 0.003), BT))| = {mp.nstr(dev, 3)}")
    ok = ok and dev < 1e-12
except ImportError as ex:
    print("N's script not available:", ex)
print("PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
