#!/usr/bin/env python3
"""Exploration: list (a,b) with cheap total clique cost: sum over prime powers q of D of 2q; compare with |delta|."""
import sys, math
from fractions import Fraction as Fr
from sympy import factorint
def dr(a, b):
    N = a*b; B = b*(b-1)
    delta = Fr(-2*a, B) + Fr((b-1)*(b+1), b*b*(N+1)) + Fr((b+1)*(2*b-1), b**3*(N-1)) + Fr(b**3+b-1, b**3) - Fr(1, N)
    rho = Fr(-(b-2), B) + Fr(b-1, b*(N+1)) - Fr((2*b-1)*(b*b-2), b*b*(b-1)*(N-1)) + Fr(b+2, a*b*b)
    return delta, rho
out = []
for b in range(2, 7):
    for a in range(3, 1200):
        d, r = dr(a, b)
        if d >= 0: continue
        D = d.denominator*r.denominator//math.gcd(d.denominator, r.denominator)
        qs = [p**e for p, e in factorint(D).items()]
        c = sum(2*q-1 for q in qs)
        n = a*(b+1)+c
        out.append((c/float(-d), b, a, c, float(d), n, (a-1)*b/n, qs))
out.sort()
for o in out[:15]:
    print(o)
