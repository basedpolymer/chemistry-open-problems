#!/usr/bin/env python3
"""Exploration: base K_N, N=ab+1, a hubs on blocks of size b, one hubless vertex w = root.  delta, rho at w; also hub deletion."""
import sys, random
sys.path.insert(0, '/mnt/c/Users/demor/chemistry-open-problems/problems/kirchhoff-soltes-beta-families/scripts')
from fractions import Fraction as Fr
import networkx as nx, kfexact as kx, sympy as sp

def build(blocks, hubless=0):
    G = nx.Graph(); base = []
    for i, b in enumerate(blocks):
        for j in range(b):
            base.append(('u', i, j))
    for j in range(hubless):
        base.append(('w', j))
    for p in range(len(base)):
        for q in range(p+1, len(base)):
            G.add_edge(base[p], base[q])
    for i, b in enumerate(blocks):
        for j in range(b):
            G.add_edge(('h', i), ('u', i, j))
    return G
def KF(N, a, ap, Y): return Y*Fr((N+1)*(N+a-1), N) + (N+a-2) + Fr((1+ap)*(1-a), N+1) + Fr(a*(ap-2), N)
def RTH(N, a, ap, Y, b1): return Fr(1, b1)*Fr((N+1)*(N+a-2), N) + Y*Fr(N+1, N) + Fr(N-ap, N+1) + Fr(ap-2, N)
def RTW(N, a, ap, Y): return 2 + Y + (Y+a+ap-2)/N - Fr(1+ap, N+1)
def RTU(N, a, ap, Y, b1):
    """transmission of base vertex in a block of size b1 WITH hub"""
    return None
random.seed(3); bad = 0
for t in range(40):
    r = random.randint(2, 4); blocks = [random.randint(1, 4) for _ in range(r)]; hl = random.randint(1, 3)
    N = sum(blocks)+hl; G = build(blocks, hl); a = r; ap = r+hl; Y = sum(Fr(1, b) for b in blocks)
    k = kx.rt(G, ('w', 0)); k = Fr(int(k.p), int(k.q))
    if RTW(N, a, ap, Y) != k: bad += 1; print('bad', blocks, hl)
print('Rt(w) mismatches', bad)

a, b = sp.symbols('a b', positive=True)
def KfS(N, a, ap, Y): return Y*(N+1)*(N+a-1)/N + (N+a-2) + (1+ap)*(1-a)/(N+1) + a*(ap-2)/N
def RtW(N, a, ap, Y): return 2 + Y + (Y+a+ap-2)/N - (1+ap)/(N+1)
N = a*b+1; Y = a/b; ap = a+1; BB = b*(b-1)
delta = KfS(N, a, ap, Y) - KfS(N-1, a, ap, Y+1/BB)
rho = RtW(N-1, a, ap, Y+1/BB) - RtW(N, a, ap, Y)
print('delta', sp.apart(sp.together(delta), a))
print('rho', sp.apart(sp.together(rho), a))
