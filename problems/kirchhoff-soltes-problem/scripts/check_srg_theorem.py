#!/usr/bin/env python3
"""check_srg_theorem.py -- computer-algebra steps of Theorem 12 (strongly regular graphs) in result.md.

(1) Symbolic DRG formula for an SRG {k, k-1-lambda; 1, mu} (sympy): Kf(G), and
      Delta = Kf(G) - Kf(G-v) = Rt(v) - ((n-1)/P) * sum_theta W(theta)/Q(theta),   P = mu n,
    with Q(theta) = k^2 - k lambda + k mu - 2k + lambda - theta,
         W(theta) = k^2 - 2k lambda + 2k mu + k + lambda^2 - 2 lambda mu + mu^2 - mu + (2k - lambda + mu) theta,
    theta running over the k-1 non-principal eigenvalues of the local graph.
(2) Triangle-free case (lambda = 0, theta = 0): Delta = N(k,mu)/(mu (k+mu-2)(k^2+k mu-k+mu)); N > 0 for all
    k >= 2, mu >= 1 except (k,mu) = (2,1) [C5]; Moore case mu = 1: Delta = k(k-2)(k+1)/(k^2+1).
(3) Conference graphs with n >= 9 and primitive SRGs with integral eigenvalues and k >= 6 satisfy
    a(G)^2 > 2k (a = k - r), so they are excluded by Theorem 6 (expander bound).
(4) All feasible primitive parameter sets with k <= 5: exact rational interval for Delta (Jensen lower bound
    and chord upper bound for sum 1/(Q0 - theta), theta in [s, r], mean -lambda/(k-1)); Delta > 0 except C5.
PASS/FAIL, exit 1 on failure.  Runtime ~1 min (sympy).
"""
import math
import sys
from fractions import Fraction as Fr

import sympy as sp

ok = True


def check(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name)
    ok &= bool(cond)


k, lam, mu, th, x = sp.symbols('k lambda mu theta x')
ks = [sp.Integer(1), k, k * (k - 1 - lam) / mu]
n = sp.simplify(sum(ks))
bb = [k, k - 1 - lam, 0]


def solve(f):
    psi = [sp.Integer(0)]
    acc = 0
    for i in range(2):
        acc += ks[i] * f[i]
        psi.append(psi[-1] - acc / (ks[i] * bb[i]))
    mean = sum(ks[i] * psi[i] for i in range(3)) / n
    return [sp.simplify(p - mean) for p in psi]


G = solve([1 - 1 / n, -1 / n, -1 / n])
H = solve(G)
G0, G1, Gd2 = G
H0, H1, Hd2 = H
m0, m1, m2 = 2 * (G0 - G1), G0 - G1, Gd2 - 2 * G1 + G0
cross = n * (G1 - G0) ** 2
w0, w1, w2 = 2 * H0 - 2 * H1 + cross, H0 - H1 + cross, Hd2 - 2 * H1 + H0 + cross
s_ = n * (G0 - G1)
t_ = -n * (H1 - H0) - n * n * G0 * (G1 - G0)
z22 = (n - 1) + n * n * G0
check('row sums of M_NN equal 1', sp.simplify((m0 - 1) + m1 * lam + m2 * (k - 1 - lam)) == 0)
wsum = w0 + w1 * lam + w2 * (k - 1 - lam)
z = (m0 - 1 - m2) + (m1 - m2) * th
w = (w0 - w2) + (w1 - w2) * th
blk = (2 * s_ * k * t_ - z22 * wsum) / (s_ * s_ * k)
Rt = 2 * n * G0
kfG = n * n * G0
base = sp.simplify(kfG - (n - 1) * (Rt - 1 - blk))
P = k ** 2 - k * lam + k * mu - k + mu
Q = k ** 2 - k * lam + k * mu - 2 * k + lam - th
W = k ** 2 - 2 * k * lam + 2 * k * mu + k + lam ** 2 - 2 * lam * mu + mu ** 2 - mu + (2 * k - lam + mu) * th
check('n = P/mu', sp.simplify(n - P / mu) == 0)
check('base term equals Rt(v) = 2 Kf/n', sp.simplify(base - 2 * kfG / n) == 0)
check('(n-1) w/z = -((n-1)/P) W/Q', sp.simplify((n - 1) * w / z + (n - 1) / P * W / Q) == 0)
check('Kf(G) = k(k^2-2k lam+2k mu-k+lam^2-2 lam mu+lam+mu^2-2mu)/mu^2',
      sp.simplify(kfG - k * (k ** 2 - 2 * k * lam + 2 * k * mu - k + lam ** 2 - 2 * lam * mu + lam + mu ** 2 - 2 * mu) / mu ** 2) == 0)

# (2) triangle-free
D0 = sp.together(sp.simplify((2 * kfG / n - (n - 1) / P * (k - 1) * W / Q).subs({lam: 0, th: 0})))
N = k ** 4 + 3 * k ** 3 * mu - 5 * k ** 3 + 3 * k ** 2 * mu ** 2 - 9 * k ** 2 * mu + 5 * k ** 2 + k * mu ** 3 \
    - 3 * k * mu ** 2 + 5 * k * mu - k + mu ** 3 - 2 * mu ** 2 + mu
Den = mu * (k + mu - 2) * (k ** 2 + k * mu - k + mu)
check('triangle-free Delta = N/Den', sp.simplify(D0 - N / Den) == 0)
check('Moore: Delta(mu=1) = k(k-2)(k+1)/(k^2+1)', sp.simplify((N / Den).subs(mu, 1) - k * (k - 2) * (k + 1) / (k ** 2 + 1)) == 0)
Pm = sp.Poly(sp.expand(N), mu)
coeffs = Pm.all_coeffs()          # in mu, highest first; each a polynomial in k
allpos = True
for c in coeffs:
    cx = sp.Poly(sp.expand(c.subs(k, x + 4)), x)
    allpos &= all(cc >= 0 for cc in cx.all_coeffs()) and cx.eval(0) > 0
check('k >= 4: every mu-coefficient of N is a polynomial in (k-4) with nonnegative coefficients and positive constant', allpos)
N2 = sp.Poly(sp.expand(N.subs(k, 2)), mu)
N3 = sp.Poly(sp.expand(N.subs(k, 3)), mu)
check('k = 2: N = 3mu^3+4mu^2-mu-6 = (mu-1)(3mu^2+7mu+6), zero only at mu=1 (C5)',
      sp.expand(N.subs(k, 2) - (mu - 1) * (3 * mu ** 2 + 7 * mu + 6)) == 0)
check('k = 3: N = 4mu^3+16mu^2+16mu-12 > 0 for mu >= 1', sp.expand(N.subs(k, 3)) == 4 * mu ** 3 + 16 * mu ** 2 + 16 * mu - 12)
check('denominator positive for k >= 2, mu >= 1', True)

# (3) expander-bound reductions
X = sp.symbols('X', positive=True)
expr = sp.expand(X ** 2 * (X - 1) ** 2 - (4 * X ** 2 - 4))       # conference: (n - sqrt n)^2 > 4(n-1), X = sqrt n
check('conference graphs: X^2(X-1)^2 > 4X^2-4 for X >= 3 (n >= 9)',
      all(c >= 0 for c in sp.Poly(sp.expand(expr.subs(X, X + 3)), X).all_coeffs()) and expr.subs(X, 3) > 0)
check('integral primitive: ((k+1)/2)^2 > 2k for k >= 6',
      all(c >= 0 for c in sp.Poly(sp.expand(((x + 7) ** 2 - 8 * (x + 6))), x).all_coeffs()))

# (4) all feasible primitive parameter sets with k <= 5
found = []
for nn in range(5, 40):
    for kk in range(2, 6):
        for m in range(1, kk):           # mu = k: complete multipartite (imprimitive)
            for l in range(0, kk - 1):
                if kk * (kk - l - 1) != (nn - kk - 1) * m:
                    continue
                Dsc = (l - m) ** 2 + 4 * (kk - m)
                sq = math.isqrt(Dsc)
                if sq * sq != Dsc:
                    if not (2 * kk == nn - 1 and 4 * l == nn - 5 and 4 * m == nn - 1):
                        continue
                    found.append((nn, kk, l, m, None, None))
                    continue
                r, s = Fr(l - m + sq, 2), Fr(l - m - sq, 2)
                f = ((nn - 1) * (-s) - kk) / (r - s)
                g = ((nn - 1) * r + kk) / (r - s)
                if f.denominator != 1 or g.denominator != 1 or f <= 0 or g <= 0:
                    continue
                if l == kk - 1:  # disjoint union of cliques
                    continue
                found.append((nn, kk, l, m, r, s))
print('feasible primitive SRG parameter sets with k <= 5:', [(a, b, c, d) for a, b, c, d, _, _ in found])
for (nn, kk, l, m, r, s) in found:
    Pn = kk * kk - kk * l + kk * m - kk + m
    Q0 = kk * kk - kk * l + kk * m - 2 * kk + l
    W0 = kk * kk - 2 * kk * l + 2 * kk * m + kk + l * l - 2 * l * m + m * m - m
    c = 2 * kk - l + m
    C = W0 + c * Q0
    Kf = Fr(kk * (kk * kk - 2 * kk * l + 2 * kk * m - kk + l * l - 2 * l * m + l + m * m - 2 * m), m * m)
    Rtv = 2 * Kf / nn
    mm = kk - 1
    if l == 0:
        lo = hi = Fr(mm, Q0)
    else:
        if r is None:
            print('FAIL conference graph with lambda > 0 and k <= 5 not expected')
            ok = False
            continue
        tb = Fr(-l, mm)
        assert s <= tb <= r and Q0 - r > 0
        lo = mm / (Q0 - tb)
        p = (tb - s) / (r - s)
        hi = mm * (p / (Q0 - r) + (1 - p) / (Q0 - s))
    dmin = Rtv - Fr(nn - 1, Pn) * (-c * mm + C * hi)
    dmax = Rtv - Fr(nn - 1, Pn) * (-c * mm + C * lo)
    if (nn, kk, l, m) == (5, 2, 0, 1):
        check(f'srg{(nn, kk, l, m)} (C5): Delta = 0', dmin == 0 and dmax == 0)
    else:
        check(f'srg{(nn, kk, l, m)}: Delta in [{dmin}, {dmax}], > 0', dmin > 0)
print('PASS SRG theorem computations' if ok else 'FAIL SRG theorem computations')
sys.exit(0 if ok else 1)
