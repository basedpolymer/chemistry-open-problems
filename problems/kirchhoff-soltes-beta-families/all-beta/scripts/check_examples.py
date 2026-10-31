#!/usr/bin/env python3
"""check_examples.py -- complete good set of the explicit graphs E1-E4 of paper.tex, exact rationals, via Lemma 1.

check_direct.py tests these graphs on the actual Laplacian, one vertex per class.  Here every vertex type is covered:
for G = H_{a,b} + W the script evaluates, with exact fractions,
   block vertex:  Delta = delta - m rho + T                    (must be 0)
   hub:           Delta = delta_h - m rho_h + T                (must be nonzero)
   root w and every vertex y carrying leaves, inner tree vertices: cut vertices, never good
   vertex z of a branch Z of W, not a cut vertex: Delta = delta_Z(z) - m_R rho_Z(z) + Rt_R(w), with m_R = n - 1 - m_Z,
   Rt_R(w) = Rt_H(w) + T - T_Z; delta_Z and rho_Z are computed by brute force on Z (exact), each type once.
The tree of Lemma 5 is handled leaf by leaf.

Output: K(G) = the ab block vertices for E1-E4 and the exact proportion.
Runtime: 1 min, one core.  Environment: Fedora 44 under WSL, Python 3.12, python-flint 0.6, networkx 3.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts'))
from fractions import Fraction as Fr
import networkx as nx
import kfexact as kx
from hubgraph import *
from check_direct import EXAMPLES

C = Checker()


def fr(z):
    return Fr(int(z.p), int(z.q))


def gadget_T(c, j):
    return (c-1)*Fr(2, c)+j*(1+Fr(2, c))


for name, (a, b, gadgets, fillers, tree) in EXAMPLES.items():
    B = b*(b-1)
    N = a*b+1
    m = sum(c-1+j for c, j in gadgets)+sum((c-1)*y for c, y in fillers)+tree[0]
    n = N+a+m
    T = sum((gadget_T(c, j) for c, j in gadgets), Fr(0))+sum((y*(c-1)*Fr(2, c) for c, y in fillers), Fr(0))+tree[1]
    delta, rho = delta_closed(a, b), rho_closed(a, b)
    C.check(delta-m*rho+T == 0, '%s: Delta(block vertex) = 0 (a=%d, b=%d, m=%d, n=%d, T=%s)' % (name, a, b, m, n, T))
    dh = delta_hub_closed(a, b)-m*rho_hub_closed(a, b)+T
    C.check(dh != 0, '%s: Delta(hub) = %s > 0' % (name, float(dh)))
    RtH = rt_hubless_closed(N, a, a+1, Fr(a, b))
    ok = True
    for c, j in gadgets:
        Z = nx.Graph()
        Z.add_node('x')
        vs = add_gadget(Z, 'x', c, j, 'g')
        TZ = fr(kx.rt(Z, 'x'))
        assert TZ == gadget_T(c, j)
        mZ = Z.number_of_nodes()-1
        m_R = n-1-mZ
        RtR = RtH+T-TZ
        if len(vs) > 1:                     # clique vertex other than x, y
            z = vs[1]
            Zz = kx.minus(Z, z)
            D = fr(kx.kf(Z)-kx.kf(Zz))-m_R*(fr(kx.rt(Zz, 'x'))-TZ)+RtR
            ok = ok and D != 0
        if j:                               # a leaf
            z = ('g', 'leaf', 0)
            Zz = kx.minus(Z, z)
            D = fr(kx.kf(Z)-kx.kf(Zz))-m_R*(fr(kx.rt(Zz, 'x'))-TZ)+RtR
            ok = ok and D != 0
        else:                               # y is a clique vertex like the others
            z = vs[0]
            Zz = kx.minus(Z, z)
            D = fr(kx.kf(Z)-kx.kf(Zz))-m_R*(fr(kx.rt(Zz, 'x'))-TZ)+RtR
            ok = ok and D != 0
    C.check(ok, '%s: Delta != 0 on all non-cut gadget vertices' % name)
    # tree: leaves of every depth
    Tr = nx.Graph()
    Tr.add_node('x')
    add_tree(Tr, 'x', tree[0], tree[1])
    TZ = fr(kx.rt(Tr, 'x'))
    assert TZ == tree[1]
    m_R = n-1-tree[0]
    RtR = RtH+T-TZ
    ok = True
    seen = set()
    for z in list(Tr.nodes()):
        if z == 'x' or Tr.degree(z) != 1:
            continue
        d = nx.shortest_path_length(Tr, 'x', z)
        if d in seen:
            continue
        seen.add(d)
        Zz = kx.minus(Tr, z)
        D = fr(kx.kf(Tr)-kx.kf(Zz))-m_R*(fr(kx.rt(Zz, 'x'))-TZ)+RtR
        ok = ok and D != 0
    C.check(ok, '%s: Delta != 0 on tree leaves (depths %s)' % (name, sorted(seen)))
    print('%s: K(G) = the %d block vertices of %d, proportion %s = %.5f' % (name, a*b, n, Fr(a*b, n), a*b/n), flush=True)
C.finish()
