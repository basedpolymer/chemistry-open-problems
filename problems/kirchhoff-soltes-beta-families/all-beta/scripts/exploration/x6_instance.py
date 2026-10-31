#!/usr/bin/env python3
"""Find small solved instances (any proportion) and verify on the actual graph, exact rationals."""
import sys, math, time
from fractions import Fraction as Fr
from sympy import factorint
sys.path.insert(0, '/mnt/c/Users/demor/chemistry-open-problems/problems/kirchhoff-soltes-beta-families/scripts')
import kfexact as kx
import networkx as nx
from x1_formula import build

def dr(a, b):
    N = a*b; B = b*(b-1)
    delta = Fr(-2*a, B) + Fr((b-1)*(b+1), b*b*(N+1)) + Fr((b+1)*(2*b-1), b**3*(N-1)) + Fr(b**3+b-1, b**3) - Fr(1, N)
    rho = Fr(-(b-2), B) + Fr(b-1, b*(N+1)) - Fr((2*b-1)*(b*b-2), b*b*(b-1)*(N-1)) + Fr(b+2, a*b*b)
    return delta, rho

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
            j = (c+1) % q
            T_g += 2-Fr(1, q)+j*(1+Fr(1, q)); m_g += 2*q-1+j; gad.append((q, j))
        if m_g > m: continue
        L = m-m_g; Tt = tau-T_g
        Tt = int(Tt)
        if L <= Tt <= L*(L+1)//2 and (L > 0 or Tt == 0):
            return m, gad, L, Tt
    return None

def tree_depths(L, t):
    # depth multiset: path 1..k + leaves at depth 1 + one at depth 1+e
    if L == 0: return []
    k = 1
    while k+1 <= L and (k+1)*(k+2)//2+L-k-1 <= t: k += 1
    e = t-(k*(k+1)//2+L-k)
    return k, e

def add_W(G, root, gad, L, t):
    c = 0
    for q, j in gad:
        vs = [root]+[('c', c, i) for i in range(2*q-1)]
        for p in range(len(vs)):
            for r in range(p+1, len(vs)):
                G.add_edge(vs[p], vs[r])
        for i in range(j):
            G.add_edge(vs[1], ('l', c, i))
        c += 1
    if L:
        k, e = tree_depths(L, t)
        path = [root]+[('p', i) for i in range(k)]
        for i in range(k):
            G.add_edge(path[i], path[i+1])
        rest = L-k
        for i in range(rest):
            if i == 0 and e > 0:
                G.add_edge(path[e], ('q', i))   # depth e+1
            else:
                G.add_edge(root, ('q', i))
    return G

if __name__ == '__main__':
    b = int(sys.argv[1]); a0 = int(sys.argv[2]); a1 = int(sys.argv[3]); mfac = float(sys.argv[4]); nmax = int(sys.argv[5])
    res = []
    for a in range(a0, a1+1):
        s = solve(a, b, int(mfac*a))
        if s:
            n = a*(b+1)+s[0]
            res.append((n, a, s))
    res.sort(key=lambda z: z[0])
    for r in res[:8]:
        print(r)
    for n, a, s in res[:1]:
        if n > nmax:
            print('too big'); break
        m, gad, L, t = s
        G = build([b]*a)
        add_W(G, ('h', 0), gad, L, t)
        print('graph n', G.number_of_nodes(), 'expected', n, 'edges', G.number_of_edges())
        t0 = time.time()
        K = kx.kf(G); print('Kf(G)', K, time.time()-t0, flush=True)
        v = ('u', 1, 0)
        Gv = kx.minus(G, v)
        Kv = kx.kf(Gv); print('Kf(G-v)==Kf(G)?', Kv == K, time.time()-t0)
