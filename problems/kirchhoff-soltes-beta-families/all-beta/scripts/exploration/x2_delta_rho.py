#!/usr/bin/env python3
"""Exploration: hub transmission, delta, rho symbolic; denominators; verification vs brute force."""
import sys, random
sys.path.insert(0, '/mnt/c/Users/demor/chemistry-open-problems/problems/kirchhoff-soltes-beta-families/scripts')
from x1_formula import build, Kf_formula
from fractions import Fraction as Fr
import kfexact as kx, sympy as sp

def Rt_hub(N, a, ap, Y, b1):
    x1 = Fr(1, b1)
    # derived by hand; check
    return x1*Fr((N+1)*(N+a-2), N) + Y*Fr(N+1, N) + Fr(N-ap, N+1) + Fr(ap-2, N)

random.seed(2); bad = 0
for t in range(30):
    r = random.randint(2, 5)
    blocks = [random.randint(1, 4) for _ in range(r)]
    hl = random.randint(0, 3)
    N = sum(blocks)+hl
    G = build(blocks, hl)
    a = len(blocks); ap = a+hl
    Y = sum(Fr(1, b) for b in blocks)
    f = Rt_hub(N, a, ap, Y, blocks[0])
    k = kx.rt(G, ('h', 0)); k = Fr(int(k.p), int(k.q))
    if f != k:
        bad += 1; print('MISMATCH', blocks, hl, f, k)
print('Rt mismatches', bad)

a, b, N = sp.symbols('a b N', positive=True)
def KfS(N, a, ap, Y):
    return Y*(N+1)*(N+a-1)/N + (N+a-2) + (1+ap)*(1-a)/(N+1) + a*(ap-2)/N
def RtS(N, a, ap, Y, b1):
    return (N+1)*(N+a-2)/(b1*N) + Y*(N+1)/N + (N-ap)/(N+1) + (ap-2)/N
BB = b*(b-1)
Y = a/b
Nv = a*b
delta = KfS(N, a, a, Y) - KfS(N-1, a, a, Y + 1/BB)
rho = RtS(N-1, a, a, Y+1/BB, b) - RtS(N, a, a, Y, b)
delta_ab = sp.factor(sp.together(delta.subs(N, Nv)))
rho_ab = sp.factor(sp.together(rho.subs(N, Nv)))
print('delta(a,b)=', delta_ab)
print('rho(a,b)=', rho_ab)
print('partial frac delta in a:', sp.apart(delta_ab, a))
print('partial frac rho in a:', sp.apart(rho_ab, a))
