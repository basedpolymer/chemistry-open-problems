#!/usr/bin/env python3
"""check_triangle_free_drg.py -- Theorem 13'' of result.md.

For every distance-regular graph with a_1 = 0 (triangle-free), valency k >= 2, n vertices, diameter >= 2,
Proposition 10 simplifies (b_1 = k - 1; the Gamma^2 terms cancel) to

   Kf(G) - Kf(G - v) = [ (k-1)(n-1)^2 - 2k(n-k-1) Kf(G)/n ] / [ k((k-2)n + 2) ]      at every vertex v,

so G is Kirchhoff Soltes iff Kf(G) = n(k-1)(n-1)^2 / (2k(n-k-1)).
Checks: (1) symbolic derivation from the radial Green functions (sympy) for D = 2 and for general D
(only Gamma_0, Gamma_1, Gamma_2 enter; Gamma_2 - Gamma_1 depends on k_1 b_1 = k(k-1) only);
(2) agreement with the exact rational DRG formula (drg_formula.py) on all a_1 = 0 arrays of drg_families
(cycles, hypercubes, odd graphs, generalized polygons, triangle-free SRGs, cubic DRGs);
(3) hypercubes: Kf(Q_d) = 2^d S_1/2 with S_1 = sum_i C(d,i)/i, and Delta > 0 for every d >= 2:
    d <= 4 exactly; d >= 5 from S_1 < 4n/(d+1) and (d^2 - 1)(n-1)^2 > 4dn(n-d-1);
(4) cycles: Delta(C_n) = -n(n-1)(n-5)/12.
PASS/FAIL; exit 1 on failure.  Runtime ~1 min.
"""
import os
import sys
from fractions import Fraction as Fr
from math import comb

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drg_formula import drg_kf_pair  # noqa: E402

ok = True


def check(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name)
    ok &= bool(cond)


# (1) symbolic: a_1 = 0 so b_1 = k - 1; Gamma_1 = Gamma_0 - r_1/2, Gamma_2 = Gamma_1 - (r_2 - r_1)/2 with
# r_1 = 2(1 - 1/n)/k, r_2 - r_1 = 2(1 - (1 + k)/n)/(k b_1); Gamma^2: psi_0 - psi_1 = Gamma_0/k,
# psi_1 - psi_2 = (Gamma_0 + k Gamma_1)/(k b_1)  (layer currents; valid for any diameter >= 2)
n, k, g, h = sp.symbols('n k g h', positive=True)
b1 = k - 1
r1 = 2 * (1 - 1 / n) / k
r21 = 2 * (1 - (1 + k) / n) / (k * b1)
G0 = g
G1 = G0 - r1 / 2
G2 = G1 - r21 / 2
H0 = h
H1 = H0 - G0 / k
H2 = H1 - (G0 + k * G1) / (k * b1)
m0, m1, m2 = 2 * (G0 - G1), G0 - G1, G2 - 2 * G1 + G0
cross = n * (G1 - G0) ** 2
w0, w2 = 2 * H0 - 2 * H1 + cross, H2 - 2 * H1 + H0 + cross
s = n * (G0 - G1)
t = -n * (H1 - H0) - n * n * G0 * (G1 - G0)
z22 = (n - 1) + n * n * G0
wsum = w0 + w2 * (k - 1)
tr = (k - 1) * (w0 - w2) / (m0 - 1 - m2) + (2 * s * k * t - z22 * wsum) / (s * s * k)
Delta = sp.simplify(n * n * G0 - (n - 1) * (2 * n * G0 - 1 - tr))
closed = ((k - 1) * (n - 1) ** 2 - 2 * k * (n - k - 1) * g * n) / (k * ((k - 2) * n + 2))   # Kf = n^2 g
check('symbolic: Delta = [(k-1)(n-1)^2 - 2k(n-k-1)Kf/n] / [k((k-2)n+2)]  (h = Gamma^2_vv cancels)',
      sp.simplify(Delta - closed) == 0)


def closed_value(nn, kk, kf):
    return (Fr((kk - 1) * (nn - 1) ** 2) - Fr(2 * kk * (nn - kk - 1)) * kf / nn) / (kk * ((kk - 2) * nn + 2))


# (2) agreement with the exact formula on many a_1 = 0 arrays
arrays = []
for nn in range(4, 60):
    D = nn // 2
    arrays.append(([2] + [1] * (D - 1), [1] * (D - 1) + ([2] if nn % 2 == 0 else [1])))
for d in range(2, 40):
    arrays.append(([d - i for i in range(d)], [i + 1 for i in range(d)]))
for kk in range(2, 30):
    b, c = [], []
    for i in range(kk):
        j = i // 2
        b.append(kk + 1 - j if i % 2 == 0 else kk - j)
    for i in range(1, kk + 1):
        c.append((i + 1) // 2)
    arrays.append((b, c))
for q in range(2, 30):
    for m in (3, 4, 6):
        arrays.append(([q + 1] + [q] * (m - 1), [1] * (m - 1) + [q + 1]))
for kk in range(2, 40):
    for mu in range(1, kk + 1):
        if (kk * (kk - 1)) % mu == 0:
            arrays.append(([kk, kk - 1], [1, mu]))
for b, c in [([3, 2, 2, 1], [1, 1, 2, 3]), ([3, 2, 2, 1], [1, 1, 1, 2]), ([3, 2, 1, 1, 1], [1, 1, 1, 2, 3]),
             ([3, 2, 2, 1, 1], [1, 1, 2, 2, 3]), ([3, 2, 2, 2, 1, 1, 1], [1, 1, 1, 1, 1, 1, 3]),
             ([3, 2, 2, 2, 2, 1, 1, 1], [1, 1, 1, 1, 2, 2, 2, 3])]:
    arrays.append((b, c))
bad = 0
for b, c in arrays:
    nn, kk, kf, kfh = drg_kf_pair(b, c)
    if kf - kfh != closed_value(nn, kk, kf):
        bad += 1
check(f'closed form = exact DRG formula on {len(arrays)} triangle-free arrays', bad == 0)

# (3) hypercubes
allpos = True
for d in range(2, 121):
    nn = 2 ** d
    S1 = sum(Fr(comb(d, i), i) for i in range(1, d + 1))
    val = closed_value(nn, d, nn * S1 / 2)
    allpos &= val > 0
    if d >= 5:
        allpos &= S1 < Fr(4 * nn, d + 1)
check('hypercubes: Delta(Q_d) > 0 and S_1 < 4n/(d+1) (d >= 5) for 2 <= d <= 120 (exact)', allpos)
x = sp.symbols('x', nonnegative=True)
check('d >= 5: d^2 - 1 >= 4d (so (d^2-1)(n-1)^2 > 4d n(n-d-1) since n(n-d-1) < (n-1)^2)',
      all(cc >= 0 for cc in sp.Poly(sp.expand((x + 5) ** 2 - 1 - 4 * (x + 5)), x).all_coeffs()))
# (4) cycles
cyc_ok = all(closed_value(nn, 2, Fr(nn ** 3 - nn, 12)) == Fr(-nn * (nn - 1) * (nn - 5), 12) for nn in range(4, 300))
check('cycles: Delta(C_n) = -n(n-1)(n-5)/12 (zero only for n = 5)', cyc_ok)
print('PASS triangle-free DRG theorem' if ok else 'FAIL triangle-free DRG theorem')
sys.exit(0 if ok else 1)
