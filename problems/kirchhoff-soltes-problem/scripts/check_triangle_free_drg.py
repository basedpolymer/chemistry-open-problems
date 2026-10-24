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
(4) cycles: Delta(C_n) = -n(n-1)(n-5)/12;
(5) valency >= 5: Markowsky-Koolen bound r_D <= (1 + 94/101) r_1  =>  Delta > 0;
(7) Hamming graphs H(d,q): closed form linear in Kf and Delta > 0 for all d >= 2, q >= 2;
(6) valency 3 and 4 with diameter > 2: the classification lists quoted in Markowsky-Koolen (2010), Cases 2-3
    (including the four valency-4 graphs with a_1 = 1, local graph 2K_2): Delta != 0 (exact).
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
# (5) valency >= 5: Markowsky-Koolen (EJC 17 (2010) R78, Thm 1): r_D <= K r_1, K = 1 + 94/101, for every DRG of
# valency > 2.  With r_1 = 2(n-1)/(nk) and Kf <= C(n,2) r_D: Kf <= K (n-1)^2/k < n(k-1)(n-1)^2/(2k(n-k-1))
# iff 2K(n-k-1) < n(k-1), true for k >= 5 because 2K < 4 <= k-1.
K = sp.Rational(1) + sp.Rational(94, 101)
check('2K < 4 (so 2K(n-k-1) < 2K n < (k-1) n for k >= 5)', 2 * K < 4)
# (6) valency 3 and 4: classification lists (BCN Thm 7.5.1 for k = 3; Brouwer-Koolen 1999 for k = 4, D > 2),
# as quoted in Markowsky-Koolen 2010, Cases 2 and 3; D <= 2 graphs are complete or strongly regular (Thm 12).
cubic = [([3, 2, 1], [1, 2, 3]), ([3, 2, 2], [1, 1, 3]), ([3, 2, 2, 1], [1, 1, 2, 3]), ([3, 2, 2, 1], [1, 1, 1, 2]),
         ([3, 2, 2, 2], [1, 1, 1, 3]), ([3, 2, 1, 1, 1], [1, 1, 1, 2, 3]), ([3, 2, 2, 1, 1], [1, 1, 2, 2, 3]),
         ([3, 2, 2, 2, 2, 2], [1, 1, 1, 1, 1, 3]), ([3, 2, 2, 2, 1, 1, 1], [1, 1, 1, 1, 1, 1, 3]),
         ([3, 2, 2, 2, 2, 1, 1, 1], [1, 1, 1, 1, 2, 2, 2, 3])]
quartic_tf = [([4, 3, 1], [1, 3, 4]), ([4, 3, 2], [1, 2, 4]), ([4, 3, 2, 1], [1, 2, 3, 4]), ([4, 3, 3], [1, 1, 4]),
              ([4, 3, 3, 1], [1, 1, 3, 4]), ([4, 3, 3], [1, 1, 2]), ([4, 3, 3, 2, 2, 1, 1], [1, 1, 2, 2, 3, 3, 4]),
              ([4, 3, 3, 3], [1, 1, 1, 4]), ([4, 3, 3, 3, 3, 3], [1, 1, 1, 1, 1, 4])]
quartic_a1 = [([4, 2, 1], [1, 1, 4]), ([4, 2, 2], [1, 1, 2]), ([4, 2, 2, 2], [1, 1, 1, 2]), ([4, 2, 2, 2, 2, 2], [1, 1, 1, 1, 1, 2])]
res = []
for b, c in cubic + quartic_tf:
    nn, kk, kf, kfh = drg_kf_pair(b, c)
    res.append(kf - kfh)
for b, c in quartic_a1:        # a_1 = 1: local graph 2K_2, non-principal eigenvalues 1, -1, -1
    nn, kk, kf, kfh = drg_kf_pair(b, c, [Fr(1), Fr(-1), Fr(-1)])
    res.append(kf - kfh)
print('   Delta for the cubic (D>2), quartic triangle-free (D>2) and quartic a_1=1 (D>2) lists:',
      [round(float(x), 4) for x in res])
check('valency 3 and 4 DRGs with D > 2 (classification lists): Delta != 0', all(x != 0 for x in res))
# direct check of the a_1 = 1 formula on the line graph of the Petersen graph (array (4,2,1;1,1,4))
from kflib import kf_inv, delete  # noqa: E402
pe = [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
LP = {i: {j for j in range(15) if j != i and set(pe[i]) & set(pe[j])} for i in range(15)}
nn, kk, kf, kfh = drg_kf_pair([4, 2, 1], [1, 1, 4], [Fr(1), Fr(-1), Fr(-1)])
ex, exh = kf_inv(LP), kf_inv(delete(LP, 0))
check('L(Petersen): DRG formula with local graph 2K_2 equals direct exact values',
      kf == Fr(int(ex.p), int(ex.q)) and kfh == Fr(int(exh.p), int(exh.q)))
# (7) Hamming graphs H(d,q), d >= 2, q >= 2 (local graph d K_{q-1}): Delta is linear in Gamma_0 (Gamma^2 cancels)
dd, qq = sp.symbols('d q', positive=True)
kH = dd * (qq - 1)
b1H = (dd - 1) * (qq - 1)
a1H = qq - 2
r1H = 2 * (1 - 1 / n) / kH
r21H = 2 * (1 - (1 + kH) / n) / (kH * b1H)
G0 = g
G1 = G0 - r1H / 2
G2 = G1 - r21H / 2
H0 = h
H1 = H0 - G0 / kH
H2 = H1 - (G0 + kH * G1) / (kH * b1H)
m0, m1, m2 = 2 * (G0 - G1), G0 - G1, G2 - 2 * G1 + G0
cross = n * (G1 - G0) ** 2
w0, w1, w2 = 2 * H0 - 2 * H1 + cross, H0 - H1 + cross, H2 - 2 * H1 + H0 + cross
s_ = n * (G0 - G1)
t_ = -n * (H1 - H0) - n * n * G0 * (G1 - G0)
z22 = (n - 1) + n * n * G0
wsum = w0 + w1 * a1H + w2 * (kH - 1 - a1H)
trH = sum(mult * ((w0 - w2) + (w1 - w2) * th) / ((m0 - 1 - m2) + (m1 - m2) * th)
          for th, mult in ((qq - 2, dd - 1), (-1, dd * (qq - 2))))
trH += (2 * s_ * kH * t_ - z22 * wsum) / (s_ * s_ * kH)
DeltaH = n * n * G0 - (n - 1) * (2 * n * G0 - 1 - trH)
closedH = (((dd - 1) * (n - 1) ** 2 * (n * (kH - 1) + 1) + g * dd * n * (kH + 1 - n) * (n * (2 * kH - qq) + qq))
           / (dd * (n * (kH - 1) + 1) * (n * (kH - qq) + qq)))
check('Hamming: Delta = [(d-1)(n-1)^2(n(k-1)+1) - g d n (n-k-1)(n(2k-q)+q)] / [d(n(k-1)+1)(n(k-q)+q)], g = Kf/n^2',
      sp.simplify(DeltaH - closedH) == 0)
hamok = True
for d_ in range(2, 31):
    for q_ in range(2, 31):
        nn = q_ ** d_
        k_ = d_ * (q_ - 1)
        Sq = sum(Fr(comb(d_, i) * (q_ - 1) ** i, i) for i in range(1, d_ + 1))
        kf = nn * Sq / q_            # Kf = n sum_i mult_i/(q i) = n S_q / q
        gv = kf / nn ** 2
        val = (Fr((d_ - 1) * (nn - 1) ** 2 * (nn * (k_ - 1) + 1)) + gv * d_ * nn * (k_ + 1 - nn) * (nn * (2 * k_ - q_) + q_)) \
            / (d_ * (nn * (k_ - 1) + 1) * (nn * (k_ - q_) + q_))
        b, c = [(d_ - i) * (q_ - 1) for i in range(d_)], [i + 1 for i in range(d_)]
        th = None if q_ == 2 else [Fr(q_ - 2)] * (d_ - 1) + [Fr(-1)] * (d_ * (q_ - 2))
        if d_ <= 8 and q_ <= 8:
            n2, k2, a_, b_ = drg_kf_pair(b, c, th)
            hamok &= (a_ - b_ == val) and (a_ == kf)
        hamok &= val > 0
        if (d_ - 1) * (d_ + 1) * (q_ - 1) >= 4 * d_:
            hamok &= Sq < Fr(2 * q_ * nn, (d_ + 1) * (q_ - 1))
check('Hamming: Delta > 0 for 2 <= d, q <= 30 (exact), closed form = DRG formula (d, q <= 8), S_q bound', hamok)
check('Hamming: (d-1)(d+1)(q-1) >= 4d except (q=2, d<=4) and (q=3, d=2)',
      all(((d_ - 1) * (d_ + 1) * (q_ - 1) >= 4 * d_) or (q_ == 2 and d_ <= 4) or (q_ == 3 and d_ == 2)
          for d_ in range(2, 200) for q_ in range(2, 200)))
print('PASS triangle-free DRG theorem' if ok else 'FAIL triangle-free DRG theorem')
sys.exit(0 if ok else 1)
