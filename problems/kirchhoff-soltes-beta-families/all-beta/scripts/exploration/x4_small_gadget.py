#!/usr/bin/env python3
"""Exploration: K_N + uniform hubs; W = cliques K_{2q} with j leaves (residue j/q) + tree, glued at hub h_1."""
import sys, math
from fractions import Fraction as Fr
from sympy import factorint
def dr(a, b):
    N = a*b; B = b*(b-1)
    delta = Fr(-2*a, B) + Fr((b-1)*(b+1), b*b*(N+1)) + Fr((b+1)*(2*b-1), b**3*(N-1)) + Fr(b**3+b-1, b**3) - Fr(1, N)
    rho = Fr(-(b-2), B) + Fr(b-1, b*(N+1)) - Fr((2*b-1)*(b*b-2), b*b*(b-1)*(N-1)) + Fr(b+2, a*b*b)
    return delta, rho
b = int(sys.argv[1]); a0 = int(sys.argv[2]); a1 = int(sys.argv[3]); target = float(sys.argv[4])
best = None
for a in range(a0, a1+1):
    d, r = dr(a, b)
    D = d.denominator*r.denominator//math.gcd(d.denominator, r.denominator)
    fac = factorint(D)
    qs = [p**e for p, e in fac.items()]
    # max m allowed for the target proportion
    mmax = int((a-1)*b/target - a*(b+1))
    if mmax < 0: continue
    if sum(1 for q in qs if q > mmax/2): continue
    Dq = {q: D//q for q in qs}
    for m in range(0, mmax+1):
        tau = m*r - d
        if tau <= 0: continue
        A = tau.numerator*(D//tau.denominator)
        # A/D = sum c_q/q mod 1, c_q = A * (D/q)^{-1} mod q
        T_g = Fr(0); m_g = 0
        for q in qs:
            c = A*pow(Dq[q], -1, q) % q
            if c == 0: continue
            j = (c+1) % q
            T_g += 2 - Fr(1, q) + j*(1+Fr(1, q)); m_g += 2*q-1+j
        if m_g > m: continue
        L = m-m_g; Tt = tau-T_g
        assert Tt.denominator == 1
        Tt = int(Tt)
        if L <= Tt <= L*(L+1)//2 and (L > 0 or Tt == 0):
            n = a*(b+1)+m; prop = (a-1)*b/n
            print('b', b, 'a', a, 'm', m, 'm_g', m_g, 'L', L, 'Tt', Tt, 'n', n, 'prop', round(prop, 5), 'qs', qs, flush=True)
            if best is None or n < best[0]: best = (n, a, m)
            break
print('best', best)
