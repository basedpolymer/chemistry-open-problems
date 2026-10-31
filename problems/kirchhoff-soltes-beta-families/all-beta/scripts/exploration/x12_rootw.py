#!/usr/bin/env python3
"""Exploration: root w (hubless base vertex). Smallest instances with proportion above a target."""
import sys, math
from fractions import Fraction as Fr
from sympy import factorint
def dr(a, b):
    B = b*(b-1)
    delta = Fr(-2*a, B) + Fr(b**4-b**3-2*b+1, b**3*(b-1)) + Fr(2*(b-1)*(b+2), b*b*(a*b+2)) - Fr((b+1)*(2*b*b-2*b-1), b**3*(a*b+1)) + Fr(1, a*b*b*(b-1))
    rho = Fr(1, B) + Fr(2*(b-1), b*(a*b+2)) - Fr(b*b-3*b-1, b*b*(a*b+1)) - Fr(b*b-b-1, a*b*b*(b-1))
    return delta, rho
def gadget(q, c):
    k = q if q % 2 else 2*q
    w = 2 if q % 2 else 1
    Tbase = (k-1)*Fr(w, q)
    t = (c-(k-1)*w) % q
    j = (t*pow(w, -1, q)) % q if q > 1 else 0
    return k, j, Tbase+j*(1+Fr(w, q)), k-1+j
def solve(a, b, mmax):
    d, r = dr(a, b)
    D = d.denominator*r.denominator//math.gcd(d.denominator, r.denominator)
    qs = [p**e for p, e in factorint(D).items()]
    Dq = {q: D//q for q in qs}
    for m in range(0, mmax+1):
        tau = m*r-d
        if tau <= 0: continue
        A = tau.numerator*(D//tau.denominator)
        T_g = Fr(0); m_g = 0; gad = []
        for q in qs:
            c = A*pow(Dq[q], -1, q) % q
            if c == 0: continue
            k, j, Tq, mq = gadget(q, c)
            T_g += Tq; m_g += mq; gad.append((q, k, j))
        if m_g > m: continue
        L = m-m_g; Tt = tau-T_g
        Tt = int(Tt)
        if L <= Tt <= L*(L+1)//2 and (L > 0 or Tt == 0):
            return m, gad, L, Tt
    return None
if __name__ == '__main__':
    b = int(sys.argv[1]); a0 = int(sys.argv[2]); a1 = int(sys.argv[3]); target = float(sys.argv[4])
    for a in range(a0, a1+1):
        mmax = int(a*b/target-(a*b+1+a))
        if mmax < 0: continue
        d, r = dr(a, b)
        if d >= 0: continue
        D = d.denominator*r.denominator//math.gcd(d.denominator, r.denominator)
        qs = [p**e for p, e in factorint(D).items()]
        if sum(q-1 for q in qs if q % 2) > mmax: continue
        s = solve(a, b, mmax)
        if s:
            n = a*b+1+a+s[0]
            print(n, 'b', b, 'a', a, 'm', s[0], 'prop', round(a*b/n, 5), s[1:], flush=True)
