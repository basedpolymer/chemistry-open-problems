#!/usr/bin/env python3
"""check_formulas.py -- certificates for Lemmas 3-9 of paper.tex (exact arithmetic, PASS/FAIL, exit code 1 on failure).

 F1  resistances of Lemma 4 (labelled R1-R8 below) against brute force on random graphs H(N; b_1..b_a), hubless vertices included
 F2  Kf and Rt of a hubless vertex (Lemma 5) against brute force, random graphs
 F3  delta, rho (Lemma 6): symbolic identity of rational functions in (a, b), and brute force on (a, b) grid
 F4  delta_h, rho_h (Lemma 6): the same
 F5  splice lemma (Lemma 3): random rooted graphs, brute force
 F6  gadget and clique data (Lemma 8): m, T, rho, delta of K_c with j leaves, brute force
 F7  tree lemma (Lemma 8): every (l, t) with l <= 45, by building the tree
 F8  factorisation identities and pairwise coprimality of Q_0..Q_4 (Lemma 9), symbolic and on a grid

Runtime: about 60 s, one core.  Environment: Fedora 44 under WSL, Python 3.12, sympy 1.14, python-flint 0.6, networkx 3.
"""
import os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts'))
from fractions import Fraction as Fr
from math import gcd
import networkx as nx
import sympy as sp
import kfexact as kx
from hubgraph import *

C = Checker()


def fr(z):
    return Fr(int(z.p), int(z.q))


random.seed(20261031)
# ---- F1: resistances -------------------------------------------------------------------------------------------
bad = tested = 0
for trial in range(60):
    r = random.randint(1, 5)
    blocks = [random.randint(1, 5) for _ in range(r)]
    hl = random.randint(0, 3)
    N = sum(blocks)+hl
    if N < 2:
        continue
    G = hub_graph(blocks, hl)
    x = [Fr(1, b) for b in blocks]
    def res(u, w):
        tu, tw = u[0], w[0]
        if tu == 'u' and tw == 'u':
            i, j = u[1], w[1]
            return Fr(2, N+1) if i == j else (2+(x[i]+x[j])/N)/(N+1)      # R1, R2
        if tu == 'h' and tw == 'u':
            i, j = u[1], w[1]
            return x[i]+(1-x[i])/(N+1) if i == j else x[i]*Fr(N+1, N)+(1+x[j]/N)/(N+1)      # R3, R4
        if tu == 'h' and tw == 'h':
            return (x[u[1]]+x[w[1]])*Fr(N+1, N)     # R5
        if tu == 'w' and tw == 'w':
            return Fr(2, N)     # R6
        if tu == 'w' and tw == 'u':
            return Fr(1, N)+(1+x[w[1]]/N)/(N+1)     # R7
        if tu == 'w' and tw == 'h':
            return Fr(1, N)+x[w[1]]*Fr(N+1, N)      # R8
        return res(w, u)
    nodes = list(G.nodes())
    for _ in range(30):
        u, w = random.sample(nodes, 2)
        tested += 1
        if fr(kx.resistance(G, u, w)) != res(u, w):
            bad += 1
C.check(bad == 0, 'F1 resistances R1-R8 on %d random pairs' % tested)

# ---- F2: Kf, Rt of hubless vertex --------------------------------------------------------------------------------
bad = tested = 0
for trial in range(120):
    r = random.randint(1, 6)
    blocks = [random.randint(1, 5) for _ in range(r)]
    hl = random.randint(0, 4)
    N = sum(blocks)+hl
    if N < 2:
        continue
    G = hub_graph(blocks, hl)
    a, ap = r, r+hl
    Y = sum(Fr(1, b) for b in blocks)
    tested += 1
    if fr(kx.kf(G)) != kf_closed(N, a, ap, Y):
        bad += 1
    if hl and fr(kx.rt(G, ('w', 0))) != rt_hubless_closed(N, a, ap, Y):
        bad += 1
C.check(bad == 0, 'F2 Kf and Rt(hubless) closed forms on %d random graphs' % tested)

# ---- F3, F4: delta, rho, delta_h, rho_h ------------------------------------------------------------------------------
a, b = sp.symbols('a b', positive=True)
B = b*(b-1)


def KfS(N, a, ap, Y):
    return Y*(N+1)*(N+a-1)/N + (N+a-2) + (1+ap)*(1-a)/(N+1) + a*(ap-2)/N


def RtS(N, a, ap, Y):
    return 2 + Y + (Y+a+ap-2)/N - (1+ap)/(N+1)


N = a*b+1
Y = a/b
d_def = KfS(N, a, a+1, Y) - KfS(N-1, a, a+1, Y+1/B)
r_def = RtS(N-1, a, a+1, Y+1/B) - RtS(N, a, a+1, Y)
d_cl = (-2*a/B + (b**4-b**3-2*b+1)/(b**3*(b-1)) + 2*(b-1)*(b+2)/(b**2*(a*b+2))
        - (b+1)*(2*b**2-2*b-1)/(b**3*(a*b+1)) + 1/(a*b**2*(b-1)))
r_cl = 1/B + 2*(b-1)/(b*(a*b+2)) - (b**2-3*b-1)/(b**2*(a*b+1)) - (b**2-b-1)/(a*b**2*(b-1))
C.check(sp.simplify(sp.together(d_def-d_cl)) == 0, 'F3 delta identity (rational functions in a, b)')
C.check(sp.simplify(sp.together(r_def-r_cl)) == 0, 'F3 rho identity (rational functions in a, b)')
dh_def = KfS(N, a, a+1, Y) - KfS(N, a-1, a+b, Y-1/b)
rh_def = RtS(N, a-1, a+b, Y-1/b) - RtS(N, a, a+1, Y)
dh_cl = (N+1)*(N+2*a-2)/(b*N) + 1 + (a*b-2*a-2*b)/(N+1) + (a-1)*(2-b)/N
rh_cl = -1/b + (b-2-1/b)/N - (b-1)/(N+1)
C.check(sp.simplify(sp.together(dh_def-dh_cl)) == 0, 'F4 delta_h identity')
C.check(sp.simplify(sp.together(rh_def-rh_cl)) == 0, 'F4 rho_h identity')
bad = tested = 0
for bb in range(2, 7):
    for aa in range(2, 7):
        G = hub_graph([bb]*aa, 1)
        K = fr(kx.kf(G))
        Rw = fr(kx.rt(G, ('w', 0)))
        Gv = kx.minus(G, ('u', 1, 0))
        Gh = kx.minus(G, ('h', 1))
        tested += 1
        ok = (K-fr(kx.kf(Gv)) == delta_closed(aa, bb)
              and fr(kx.rt(Gv, ('w', 0)))-Rw == rho_closed(aa, bb)
              and K-fr(kx.kf(Gh)) == delta_hub_closed(aa, bb)
              and fr(kx.rt(Gh, ('w', 0)))-Rw == rho_hub_closed(aa, bb))
        if not ok:
            bad += 1
C.check(bad == 0, 'F3-F4 delta, rho, delta_h, rho_h against brute force on %d pairs (a, b)' % tested)

# ---- F5: splice lemma ----------------------------------------------------------------------------------------------
bad = tested = 0
for trial in range(200):
    n1 = random.randint(3, 8)
    n2 = random.randint(2, 7)
    while True:
        Bg = nx.gnp_random_graph(n1, 0.6, seed=random.randint(0, 10**9))
        Wg = nx.gnp_random_graph(n2, 0.5, seed=random.randint(0, 10**9))
        if nx.is_connected(Bg) and nx.is_connected(Wg):
            break
    Gs = kx.splice([(Bg, 0), (Wg, 0)])
    for v in list(Bg.nodes())[1:]:
        Bv = kx.minus(Bg, v)
        if not nx.is_connected(Bv):
            continue
        dl = fr(kx.kf(Bg)-kx.kf(Bv))
        rh = fr(kx.rt(Bv, 0)-kx.rt(Bg, 0))
        lhs = fr(kx.kf(Gs)-kx.kf(kx.minus(Gs, (0, v))))
        rhs = dl-(n2-1)*rh+fr(kx.rt(Wg, 0))
        tested += 1
        if lhs != rhs:
            bad += 1
C.check(bad == 0, 'F5 splice lemma on %d random (B, W, v)' % tested)

# ---- F6: cliques and gadgets -----------------------------------------------------------------------------------------
bad = tested = 0
for q in range(1, 8):
    for jj in range(0, q+2):
        c = 2*q
        Z = nx.Graph()
        Z.add_node('x')
        vs = add_gadget(Z, 'x', c, jj, 'g')
        y = vs[0]
        tested += 1
        T = fr(kx.rt(Z, 'x'))
        ok = T == 2-Fr(1, q)+jj*(1+Fr(1, q)) and Z.number_of_nodes()-1 == 2*q-1+jj
        if c >= 3:                  # clique vertex different from x, y
            zz = vs[1] if len(vs) > 1 else None
            if zz is not None:
                Zz = kx.minus(Z, zz)
                dl = fr(kx.kf(Z)-kx.kf(Zz))
                rh = fr(kx.rt(Zz, 'x'))-T
                ok = ok and rh == Fr(2*(jj-1), c*(c-1)) and dl == 1+jj+Fr(2*jj, c*(c-1))
        if not ok:
            bad += 1
C.check(bad == 0, 'F6 gadget K_{2q} + j leaves: m, T, rho, delta on %d cases' % tested)
bad = 0
for c in range(3, 12):
    Z = nx.complete_graph(c)
    Zz = kx.minus(Z, 1)
    if fr(kx.kf(Z)-kx.kf(Zz)) != 1 or fr(kx.rt(Zz, 0)-kx.rt(Z, 0)) != Fr(-2, c*(c-1)):
        bad += 1
C.check(bad == 0, 'F6 K_c: delta = 1, rho = -2/(c(c-1)), 3 <= c <= 11')

# ---- F7: tree lemma ---------------------------------------------------------------------------------------------------
bad = tested = 0
for l in range(0, 46):
    for t in (range(l, l*(l+1)//2+1) if l else [0]):
        T = nx.Graph()
        T.add_node('r')
        add_tree(T, 'r', l, t)
        dist = nx.single_source_shortest_path_length(T, 'r')
        tested += 1
        if not nx.is_tree(T) or T.number_of_nodes() != l+1 or sum(dist.values()) != t:
            bad += 1
C.check(bad == 0, 'F7 rooted trees with l non-root vertices and depth sum t: all %d pairs, l <= 45' % tested)

# ---- F8: factorisations and coprimality ------------------------------------------------------------------------------
x = sp.symbols('x')
C.check(sp.expand(x*(x**2-3*x+1)+1-(x-1)*(x**2-2*x-1)) == 0 and sp.expand(x*(x**2-3*x+1)+2-(x-2)*(x**2-x-1)) == 0,
        'F8 factorisations of x(x^2-3x+1), +1, +2')
bad = tested = 0
for bb in range(2, 26):
    for ss in range(1, 41):
        xx, Q = q_list(bb, ss)
        if xx < 4:
            continue
        tested += 1
        if xx % 2 or min(Q) < 1:
            bad += 1
        for i in range(5):
            for j in range(i+1, 5):
                if gcd(Q[i], Q[j]) != 1:
                    bad += 1
        aa = (bb-1)*ss*Q[1]
        D = 1
        for q in Q:
            D *= q
        if (delta_closed(aa, bb)*D).denominator != 1 or (rho_closed(aa, bb)*D).denominator != 1:
            bad += 1
        if aa*bb != xx*Q[1] or aa*bb+1 != Q[4]*Q[2] or aa*bb+2 != (xx-2)*Q[3]:
            bad += 1
C.check(bad == 0, 'F8 Q_0..Q_4 pairwise coprime, denominators of delta and rho divide D, on %d pairs (b, s)' % tested)
C.finish()
