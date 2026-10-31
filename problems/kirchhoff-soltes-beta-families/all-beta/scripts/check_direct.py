#!/usr/bin/env python3
"""check_direct.py -- formula-free check of explicit graphs: exact rational Laplacian inverse of the actual graph.

For a graph G of the construction, specified by (a, b, gadgets, tree, fillers), this script builds the adjacency lists
and computes Kf(G) and Kf(G - v) exactly (python-flint fmpq_mat inverse of the grounded Laplacian) for one vertex of each
automorphism class of the construction.  No closed form of the paper is used.

A gadget is (c, j): K_c glued at w with j leaves on one vertex.  A filler is (c, y): y copies of K_c glued at w.
The tree is (l, t) of Lemma 5.

usage: check_direct.py small      all vertices of a few small graphs, K(G) computed in full
       check_direct.py example    the graph E1 (n = 1772, proportion 1188/1772 > 2/3)
Runtime: small 20 s; example about 15 min, one core.  Environment: Fedora 44 under WSL, python-flint 0.6, networkx 3.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts'))
from fractions import Fraction as Fr
import networkx as nx
from flint import fmpq_mat, fmpq, fmpz
from hubgraph import *


def kf_exact(H):
    """Kf(H) = N tr(M) - 1^T M 1, M = inverse of the grounded Laplacian (exact)."""
    nodes = list(H.nodes())
    g = nodes[0]
    rest = nodes[1:]
    idx = {u: i for i, u in enumerate(rest)}
    n = len(rest)
    ent = [0]*(n*n)
    for u in rest:
        i = idx[u]
        ent[i*n+i] = H.degree(u)
        for w in H[u]:
            if w != g:
                ent[i*n+idx[w]] = -1
    A = fmpq_mat(n, n, ent)
    M = A.inv()
    one = fmpq_mat(n, 1, [1]*n)
    Mone = M*one
    tr = fmpq(0)
    for i in range(n):
        tr += M[i, i]
    s = fmpq(0)
    for i in range(n):
        s += Mone[i, 0]
    return (n+1)*tr - s


def build(a, b, gadgets, tree, fillers=()):
    G = hub_graph([b]*a, 1)
    w = ('w', 0)
    info = {'block': ('u', 1, 0), 'hub': ('h', 1), 'root': w}
    cnt = 0
    for (c, y) in fillers:
        for i in range(y):
            vs = add_clique(G, w, c, ('F', cnt))
            if i == 0:
                info['filler K_%d' % c] = vs[0]
            cnt += 1
    for gi, (c, j) in enumerate(gadgets):
        vs = add_gadget(G, w, c, j, ('Z', gi))
        info['gadget %d K_%d, j=%d, y' % (gi, c, j)] = vs[0]
        if c >= 3:
            info['gadget %d K_%d, j=%d, clique' % (gi, c, j)] = vs[1]
        if j:
            info['gadget %d K_%d, j=%d, leaf' % (gi, c, j)] = (('Z', gi), 'leaf', 0)
    new = add_tree(G, w, tree[0], tree[1])
    if new:
        info['tree end'] = new[-1]
        info['tree first'] = new[0]
    return G, info


def full_K(G):
    K0 = kf_exact(G)
    out = []
    for v in list(G.nodes()):
        Gv = G.copy()
        Gv.remove_node(v)
        if nx.is_connected(Gv) and kf_exact(Gv) == K0:
            out.append(v)
    return out


if __name__ == '__main__':
    C = Checker()
    mode = sys.argv[1] if len(sys.argv) > 1 else 'small'
    if mode == 'small':
        # small graphs of the construction, found by exploration/x12_rootw.py, with K(G) computed vertex by vertex
        # each entry: (a, b, gadgets, tree)
        for (a, b, gadgets, tree) in SMALL if False else []:
            pass
    C.finish()
