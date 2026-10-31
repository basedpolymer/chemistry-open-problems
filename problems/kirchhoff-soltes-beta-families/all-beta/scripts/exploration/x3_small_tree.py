#!/usr/bin/env python3
"""Exploration: K_N + uniform hubs (a blocks of size b), glue rooted tree W at hub h_1.
Solve m*C = A mod D for the smallest m; check tree transmission range m <= T <= m(m+1)/2."""
import sys, math
from fractions import Fraction as Fr
def dr(a, b):
    N = a*b; B = b*(b-1)
    # delta, rho partial fractions
    delta = Fr(-2*a, B) + Fr((b-1)*(b+1), b*b*(N+1)) + Fr((b+1)*(2*b-1), b**3*(N-1)) + Fr(b**3+b-1, b**3) - Fr(1, N)
    rho = Fr(-(b-2), B) + Fr(b-1, b*(N+1)) - Fr((2*b-1)*(b*b-2), b*b*(b-1)*(N-1)) + Fr(b+2, a*b*b)
    return delta, rho
bmax = int(sys.argv[1]); amax = int(sys.argv[2])
found = []
for b in range(3, bmax+1):
    for a in range(2, amax+1):
        d, r = dr(a, b)
        D = d.denominator*r.denominator//math.gcd(d.denominator, r.denominator)
        A = int(d*D); C = int(r*D)
        g = math.gcd(C, D)
        if A % g: continue
        Dg = D//g
        if Dg == 1:
            m0 = 0
        else:
            m0 = (A//g)*pow(C//g, -1, Dg) % Dg
        # m C = A mod D -> tau = m rho - delta integer
        for k in range(0, 3):
            m = m0 + k*Dg
            if m > a*b: break
            tau = m*r - d
            assert tau.denominator == 1
            T = int(tau)
            if m <= T <= m*(m+1)//2:
                n = a*(b+1)+m
                prop = (a-1)*b/n
                found.append((prop, b, a, m, T, n))
                print('b',b,'a',a,'m',m,'T',T,'n',n,'prop>=',round(prop,4), flush=True)
print(len(found))
