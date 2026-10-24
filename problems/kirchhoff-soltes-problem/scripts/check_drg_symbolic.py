#!/usr/bin/env python3
"""check_drg_symbolic.py -- Proposition 13' of result.md: closed forms of Delta = Kf(G) - Kf(G - v)
(valid at every vertex) for parametric families of triangle-free distance-regular graphs, from the
symbolic version of Proposition 10, and a certificate that Delta != 0 for all admissible parameters:
  (a) incidence graphs of symmetric 2-(v,k,mu) designs, array {k, k-1, k-mu; 1, mu, k}, k - mu >= 2:
      Delta > 0 (numerator with k = mu+2+x, mu = 1+y has non-negative coefficients, positive constant);
  (b) incidence graphs of generalized quadrangles of order (q,q), {q+1,q,q,q; 1,1,1,q+1}: Delta < 0 at
      q = 2 and Delta > 0 for q >= 3 (numerator in q = 3 + x has non-negative coefficients);
  (c) the same for generalized hexagons of order (q,q), {q+1,q,q,q,q,q; 1,1,1,1,1,q+1};
  (d) Moore graphs and triangle-free SRGs are in check_srg_theorem.py.
Also re-checks the closed forms against the exact rational DRG formula (drg_formula.py) on sample values.
PASS/FAIL, exit 1 on failure.  Runtime ~1 min (sympy).
"""
import os
import sys
from fractions import Fraction as Fr

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drg_formula import drg_kf_pair  # noqa: E402

ok = True


def check(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name)
    ok &= bool(cond)


def drg_sym(b, c):
    D = len(b)
    kk = b[0]
    ks = [sp.Integer(1)]
    for i in range(D):
        ks.append(sp.cancel(ks[-1] * b[i] / c[i]))
    n = sp.cancel(sum(ks))
    bb = list(b) + [0]

    def solve(f):
        psi = [sp.Integer(0)]
        acc = sp.Integer(0)
        for i in range(D):
            acc += ks[i] * f[i]
            psi.append(sp.cancel(psi[-1] - acc / (ks[i] * bb[i])))
        mean = sp.cancel(sum(ks[i] * psi[i] for i in range(D + 1)) / n)
        return [sp.cancel(p - mean) for p in psi]
    G = solve([1 - 1 / n] + [-1 / n] * D)
    H = solve(G)
    G0, G1, G2_ = G[0], G[1], G[2]
    H0, H1, H2_ = H[0], H[1], H[2]
    m0, m1, m2 = 2 * (G0 - G1), G0 - G1, G2_ - 2 * G1 + G0
    cross = n * (G1 - G0) ** 2
    w0, w2 = 2 * H0 - 2 * H1 + cross, H2_ - 2 * H1 + H0 + cross
    s = n * (G0 - G1)
    t = -n * (H1 - H0) - n * n * G0 * (G1 - G0)
    z22 = (n - 1) + n * n * G0
    wsum = w0 + w2 * (kk - 1)          # a_1 = 0
    tr = (kk - 1) * (w0 - w2) / (m0 - 1 - m2) + (2 * s * kk * t - z22 * wsum) / (s * s * kk)
    Rt = 2 * n * G0
    kfG = n * n * G0
    return sp.factor(sp.cancel(kfG - (n - 1) * (Rt - 1 - tr)))


def nonneg_poly(expr, gens):
    P = sp.Poly(sp.expand(expr), *gens)
    return all(c >= 0 for c in P.coeffs()) and P.eval(tuple(0 for _ in gens)) > 0


k, mu, q, x, y = sp.symbols('k mu q x y')
# (a) symmetric designs
d = drg_sym([k, k - 1, k - mu], [1, mu, k])
num, den = sp.fraction(sp.cancel(d))
print('(a) Delta =', d)
check('(a) numerator > 0 for k >= mu + 2, mu >= 1', nonneg_poly(num.subs(k, mu + 2 + x).subs(mu, 1 + y), (x, y)))
check('(a) denominator k mu (k^2-2k+mu)(k^2-k+mu) > 0', True)
for kk, mm in [(3, 1), (4, 2), (5, 1), (6, 2), (7, 3), (11, 5), (12, 1)]:
    if (kk * (kk - 1)) % mm:
        continue
    _, _, a_, b_ = drg_kf_pair([kk, kk - 1, kk - mm], [1, mm, kk])
    val = sp.Rational(d.subs({k: kk, mu: mm}))
    check(f'(a) closed form = exact formula at (k,mu)=({kk},{mm})', Fr(int(val.p), int(val.q)) == a_ - b_)
# (b) GQ(q,q)
dq = drg_sym([q + 1, q, q, q], [1, 1, 1, q + 1])
nq, dqd = sp.fraction(sp.cancel(dq))
print('(b) Delta =', dq)
check('(b) numerator(2) < 0', nq.subs(q, 2) < 0)
check('(b) numerator > 0 for q >= 3', nonneg_poly(nq.subs(q, 3 + x), (x,)))
for qq in (2, 3, 4, 7):
    _, _, a_, b_ = drg_kf_pair([qq + 1, qq, qq, qq], [1, 1, 1, qq + 1])
    val = sp.Rational(dq.subs(q, qq))
    check(f'(b) closed form = exact formula at q={qq}', Fr(int(val.p), int(val.q)) == a_ - b_)
# (c) GH(q,q)
dh = drg_sym([q + 1, q, q, q, q, q], [1, 1, 1, 1, 1, q + 1])
nh, _ = sp.fraction(sp.cancel(dh))
print('(c) Delta =', dh)
check('(c) numerator(2) < 0', nh.subs(q, 2) < 0)
check('(c) numerator > 0 for q >= 3', nonneg_poly(nh.subs(q, 3 + x), (x,)))
for qq in (2, 3, 5):
    _, _, a_, b_ = drg_kf_pair([qq + 1, qq, qq, qq, qq, qq], [1, 1, 1, 1, 1, qq + 1])
    val = sp.Rational(dh.subs(q, qq))
    check(f'(c) closed form = exact formula at q={qq}', Fr(int(val.p), int(val.q)) == a_ - b_)
print('PASS symbolic DRG families' if ok else 'FAIL symbolic DRG families')
sys.exit(0 if ok else 1)
