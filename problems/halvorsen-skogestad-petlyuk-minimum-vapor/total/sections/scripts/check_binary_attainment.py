#!/usr/bin/env python3
"""
check_binary_attainment.py (Agent S) -- checks for Proposition 10 of sections.md (two components: a column with
n stages per section and V_B = s > H_1 has steady states with purities -> 1).
  (1) sympy: with 1-q eliminated by the feed equation, f(s) = alpha s (s+1-q-z) - (s+1-q)(s+1-z) has the roots
      H_1 = (1-z)/(theta-1) and s_- = 1/(alpha-1) - z theta/(alpha-theta), and
      H_1 - s_- - z alpha (theta-1)/((alpha-theta)(alpha-1)) = (1-z)/(theta-1) - (1-z)/(alpha-1) > 0  (exact identity).
  (2) numerics (illustration, not part of the proof): the shooting construction of the proof (bisection on
      eta = x_B), for a few data and n; the resulting profile is re-checked as a steady state by the independent
      network solver stagenet.py (residual), and R, purities are printed.
Prints PASS/FAIL; exits non-zero on failure. Environment: Python 3, sympy, numpy, scipy. Runtime: seconds.
"""
import sys
import numpy as np
import sympy as sp
from stagenet import Net

ok = True
a, z, th, s = sp.symbols('alpha z theta s', positive=True)
omq = a * z / (a - th) + (1 - z) / (1 - th)          # 1 - q from the feed equation
qq = 1 - omq
H1 = (1 - z) / (th - 1)
sm = 1 / (a - 1) - z * th / (a - th)
f = a * s * (s + 1 - qq - z) - (s + 1 - qq) * (s + 1 - z)
c1 = sp.simplify(f.subs(s, H1)) == 0
c2 = sp.simplify(f.subs(s, sm)) == 0
c3 = sp.simplify(H1 - sm - z * a * (th - 1) / ((a - th) * (a - 1)) - ((1 - z) / (th - 1) - (1 - z) / (a - 1))) == 0
print('f(H1) = 0:', c1, '| f(s_-) = 0:', c2, '| identity for H1 - s_-:', c3)
ok &= c1 and c2 and c3


import mpmath as mp
mp.mp.dps = 120


def shoot(alpha, zz, q, sv, nT, nB, eta):
    """shooting profiles of the proof, in 120-digit arithmetic"""
    alpha, zz, q, sv, eta = (mp.mpf(v) for v in (alpha, zz, q, sv, eta))
    e = lambda x: alpha * x / (1 + (alpha - 1) * x)
    einv = lambda y: y / (alpha - (alpha - 1) * y)
    D, B = zz, 1 - zz
    VT = sv + 1 - q; LT = VT - D; LB = sv + 1 - zz
    xD = 1 - (1 - zz) * eta / zz
    xb = [eta]
    for k in range(nB + 1):
        xb.append((sv * e(xb[-1]) + B * eta) / LB)
    xt = [einv(xD)]
    for k in range(nT):
        xt.append(einv((LT * xt[-1] + D * xD) / VT))
    return xb, xt, xD


def solve_column(alpha, zz, q, sv, n):
    eta1 = mp.mpf('1e-100')
    while True:
        xb, xt, _ = shoot(alpha, zz, q, sv, n, n, eta1)
        if xb[-1] > xt[-1] or eta1 > 0.5:
            break
        eta1 *= 2
    lo, hi = mp.mpf(0), eta1
    for it in range(400):
        mid = (lo + hi) / 2
        xb, xt, _ = shoot(alpha, zz, q, sv, n, n, mid)
        if xb[-1] > xt[-1]:
            hi = mid
        else:
            lo = mid
    return hi


for (alpha, zz, q) in [(2.0, 0.5, 1.0), (1.5, 0.3, 0.0), (3.0, 0.7, 0.5), (1.2, 0.5, 1.3)]:
    from scipy.optimize import brentq
    gfe = lambda t: alpha * zz / (alpha - t) + (1 - zz) / (1 - t) - (1 - q)
    tht = brentq(gfe, 1 + 1e-12, alpha - 1e-12, xtol=1e-15)
    h1 = float((1 - zz) / (tht - 1))
    sv = 1.02 * h1
    for n in [20, 40, 80, 160]:
        eta = solve_column(alpha, zz, q, sv, n)
        xb, xt, xD = shoot(alpha, zz, q, sv, n, n, eta)
        # rebuild the column in stagenet and check that the shooting profile is a steady state
        net = Net([alpha, 1.0], [zz, 1 - zz], q)
        VT = sv + 1 - q; LT = VT - zz; LB = sv + 1 - zz
        cond = net.add_node('mix'); tops = [net.add_node('stage') for _ in range(n)]
        fst = net.add_node('stage'); bots = [net.add_node('stage') for _ in range(n)]; reb = net.add_node('stage')
        chain = [cond] + tops + [fst] + bots + [reb]
        net.add_stream('F', fst, 'F', 1.0)
        for k in range(len(chain) - 1):
            below = k >= n + 1
            net.add_stream(chain[k], chain[k + 1], 'L', LB if below else LT)
            net.add_stream(chain[k + 1], chain[k], 'V', sv if below else VT)
        net.add_stream(cond, ('P', 'D'), 'L', zz); net.add_stream(reb, ('P', 'B'), 'L', 1 - zz)
        C = np.zeros((len(chain), 2))
        xb = [float(v) for v in xb]; xt = [float(v) for v in xt]; xD = float(xD); eta = float(eta)
        C[0] = [xD, 1 - xD]
        for k in range(n):
            C[1 + k] = [xt[k], 1 - xt[k]]
        C[1 + n] = [xt[n], 1 - xt[n]]
        for k in range(n):
            xk = xb[n - k]            # bottom stages listed top-down: stage n ... 1
            C[2 + n + k] = [xk, 1 - xk]
        C[-1] = [xb[0], 1 - xb[0]]
        net._b = net.balances(); net._prep()
        res = np.max(np.abs(net.residual(C.ravel())))
        match = abs(xb[-1] - xt[-1])
        ok &= match < 1e-12
        print(f'alpha={alpha} z={zz} q={q} H1={h1:.6f} s={sv:.6f} n={n}: x_B={eta:.3e} 1-x_D={1 - xD:.3e} '
              f'feed-stage mismatch={match:.1e} steady-state residual={res:.1e}')
        ok &= res < 1e-9
print('PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
