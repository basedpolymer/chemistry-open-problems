#!/usr/bin/env python3
"""check_multipartite.py -- Theorem 11 of result.md (complete multipartite graphs).

For G = K_{n_1,...,n_k} (n = sum n_i >= 3, k >= 2, not a star) and v in a part of size n_j:
   Kf(G) - Kf(G - v) = phi(n_j) - S,
   phi(x) = (n + x - 2)/(n - x) + x(x-1)/((n-x)(n-1-x)),   S = sum_i n_i(n_i-1)/((n-n_i)(n-1-n_i)).
phi is strictly increasing, so the good vertices are the union of the parts of one size; if all vertices
are good all parts have the same size s and then Kf(G)-Kf(G-v) = [k(k-1)s^2-2ks+2]/[(k-1)s((k-1)s-1)] > 0.
Checks: (a) the formula against exact Kf (python-flint) for all partitions of n = 3..12 into >= 2 parts
(stars excluded); (b) the closed form for K_{k x s}, k,s <= 12; (c) monotonicity of phi (exact, symbolic).
PASS/FAIL, exit 1 on failure.  Runtime ~1 min.
"""
import os
import sys
from fractions import Fraction as Fr

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kflib import kf_inv, delete  # noqa: E402

ok = True


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield []
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest


def multipartite(parts):
    adj = {}
    lab = []
    for i, p in enumerate(parts):
        lab += [i] * p
    n = len(lab)
    for u in range(n):
        adj[u] = {w for w in range(n) if lab[w] != lab[u]}
    return adj, lab


def tofr(q):
    return Fr(int(q.p), int(q.q))


cnt = 0
for n in range(3, 13):
    for parts in partitions(n):
        if len(parts) < 2 or (len(parts) == 2 and min(parts) == 1):
            continue
        adj, lab = multipartite(parts)
        S = sum(Fr(p * (p - 1), (n - p) * (n - 1 - p)) for p in parts)
        phi = lambda x: Fr(n + x - 2, n - x) + Fr(x * (x - 1), (n - x) * (n - 1 - x))
        kG = tofr(kf_inv(adj))
        kfclosed = (len(parts) - 1) + n * sum(Fr(p - 1, n - p) for p in parts)
        if kG != kfclosed:
            ok = False
            print('FAIL Kf closed form', parts)
        seen = set()
        for v in range(n):
            j = lab[v]
            if j in seen:
                continue
            seen.add(j)
            d = kG - tofr(kf_inv(delete(adj, v)))
            cnt += 1
            if d != phi(parts[j]) - S:
                ok = False
                print('FAIL formula', parts, j)
print(f'(a) {cnt} (partition, part) cases checked exactly')
for kk in range(2, 13):
    for s in range(1, 13):
        if kk * s < 3:
            continue
        n = kk * s
        if kk == 2 and s == 1:
            continue
        S = kk * Fr(s * (s - 1), (n - s) * (n - 1 - s))
        phi = Fr(n + s - 2, n - s) + Fr(s * (s - 1), (n - s) * (n - 1 - s))
        closed = Fr(kk * (kk - 1) * s * s - 2 * kk * s + 2, (kk - 1) * s * ((kk - 1) * s - 1))
        if phi - S != closed or closed <= 0:
            ok = False
            print('FAIL K_{k x s}', kk, s)
print('(b) closed form for K_{k x s} checked, all positive')
x, n_ = sp.symbols('x n', positive=True)
phis = (n_ + x - 2) / (n_ - x) + x * (x - 1) / ((n_ - x) * (n_ - 1 - x))
dphi = sp.together(sp.diff(phis, x))
num, den = sp.fraction(dphi)
# substitute n = x + 2 + y (y >= 0: n - 1 - x >= 1 for a non-star part) and x = 1 + z (z >= 0)
y, zz = sp.symbols('y z', nonnegative=True)
numsub = sp.Poly(sp.expand(num.subs(n_, x + 2 + y).subs(x, 1 + zz)), y, zz)
allpos = all(c >= 0 for c in numsub.coeffs()) and numsub.eval((0, 0)) > 0
print('(c) phi\'(x) numerator after n = x+2+y, x = 1+z:', sp.factor(num.subs(n_, x + 2 + y).subs(x, 1 + zz)))
print(('PASS ' if allpos else 'FAIL ') + 'phi strictly increasing for 1 <= x <= n-2')
ok &= allpos
print('PASS complete multipartite theorem' if ok else 'FAIL complete multipartite theorem')
sys.exit(0 if ok else 1)
