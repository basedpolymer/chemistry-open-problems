#!/usr/bin/env python3
"""check_direct.py -- formula-free check of explicit graphs: exact rational Laplacian inverse of the actual graph.

For a graph G = H_{a,b} + W (W glued at the hubless vertex w) this script builds the adjacency lists and computes Kf(G)
and Kf(G - v) exactly (python-flint fmpq_mat inverse of the grounded Laplacian).  No closed form of the paper is used.

A gadget is (c, j): K_c glued at w with j pendant leaves on one vertex of K_c.  A filler is (c, y): y copies of K_c glued
at w.  The tree is (l, t) of Lemma 5.

usage: check_direct.py small            K(G) computed vertex by vertex on 8 small graphs: K(G) = the block vertices
       check_direct.py family           a graph produced by the proof (beta = 1/100, b = 2, s = 2, n = 2000), one vertex per class
       check_direct.py E1 | E2 | E3 | E4    the explicit examples of the paper, one vertex per class
       check_direct.py E3quick          E3, only block vertex, hub and root
Runtime: small 3 min; family 5 min; E1 8 min; E2 12 min; E4 25 min; E3 several hours (E3quick 25 min), one core.
Environment: Fedora 44 under WSL, Python 3.12, python-flint 0.6, networkx 3.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import networkx as nx
from flint import fmpq_mat, fmpq
from hubgraph import *

SMALL = [   # (a, b, gadgets, fillers, tree): found by scripts/exploration/x12_rootw.py
    (3, 2, [(7, 2)], [], (6, 7)),
    (4, 2, [(3, 2), (5, 2)], [], (1, 1)),
    (5, 2, [(3, 0), (11, 3)], [], (9, 12)),
    (6, 2, [(3, 0), (7, 0), (13, 5)], [], (8, 14)),
    (7, 2, [(32, 0)], [], (10, 28)),
    (8, 2, [(9, 3), (17, 2)], [], (26, 29)),
    (43, 3, [(3, 0), (5, 4), (13, 2), (43, 9), (131, 5)], [], (16, 21)),
    (36, 3, [(16, 5), (81, 0), (11, 2), (109, 18)], [], (17, 21)),
]
EXAMPLES = {
    # proportion 1188/1772 = 297/443 > 2/3
    'E1': (396, 3, [(5, 2), (7, 6), (11, 7), (17, 8), (29, 18), (41, 10)], [], (32, 94)),
    # exactly 2/3
    'E2': (506, 3, [(4, 0), (9, 2), (5, 3), (49, 3), (11, 9), (19, 10), (23, 6), (31, 23)], [], (53, 133)),
    # exactly 3/4
    'E3': (1116, 4, [(64, 30), (27, 6), (5, 3), (7, 6), (11, 3), (19, 12), (29, 4), (31, 12), (47, 2)], [], (62, 114)),
    # proportion 2820/3735 = 188/249 > 3/4
    'E4': (705, 4, [(4, 0), (9, 0), (7, 2), (13, 9), (17, 3), (31, 25), (47, 34)], [], (15, 43)),
}


def kf_exact(H):
    """Kf(H) = N tr(M) - 1^T M 1, M = inverse of the grounded Laplacian (exact); H connected."""
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
    Mone = M*fmpq_mat(n, 1, [1]*n)
    tr = fmpq(0)
    s = fmpq(0)
    for i in range(n):
        tr += M[i, i]
        s += Mone[i, 0]
    return (n+1)*tr-s


def build(a, b, gadgets, fillers, tree):
    G = hub_graph([b]*a, 1)
    w = ('w', 0)
    info = [('block', ('u', 1, 0)), ('hub', ('h', 1)), ('root', w)]
    cnt = 0
    for (c, y) in fillers:
        for i in range(y):
            vs = add_clique(G, w, c, ('F', cnt))
            if i == 0:
                info.append(('filler K_%d' % c, vs[0]))
            cnt += 1
    for gi, (c, j) in enumerate(gadgets):
        vs = add_gadget(G, w, c, j, ('Z', gi))
        info.append(('gadget K_%d j=%d y' % (c, j), vs[0]))
        if c >= 3:
            info.append(('gadget K_%d j=%d clique' % (c, j), vs[1]))
        if j:
            info.append(('gadget K_%d j=%d leaf' % (c, j), (('Z', gi), 'leaf', 0)))
    new = add_tree(G, w, tree[0], tree[1])
    for z in new:
        if G.degree(z) == 1:
            info.append(('tree leaf', z))
            break
    return G, info


def test_classes(C, name, spec, only=None):
    a, b = spec[0], spec[1]
    G, info = build(*spec)
    n = G.number_of_nodes()
    assert nx.is_connected(G)
    t0 = time.time()
    K0 = kf_exact(G)
    print('%s: n=%d edges=%d, %d block vertices, proportion %d/%d, Kf(G) = %s  (%.1f s)' % (
        name, n, G.number_of_edges(), a*b, a*b, n, K0, time.time()-t0), flush=True)
    for lab, v in info:
        if only is not None and lab not in only:
            continue
        Gv = G.copy()
        Gv.remove_node(v)
        if not nx.is_connected(Gv):
            C.check(lab != 'block', '%s %s: cut vertex, not good' % (name, lab))
            continue
        t1 = time.time()
        good = kf_exact(Gv) == K0
        C.check(good == (lab == 'block'), '%s %s: %s (%.1f s)' % (name, lab, 'good' if good else 'not good', time.time()-t1))


if __name__ == '__main__':
    C = Checker()
    mode = sys.argv[1] if len(sys.argv) > 1 else 'small'
    if mode == 'small':
        for spec in SMALL:
            a, b = spec[0], spec[1]
            G, info = build(*spec)
            K0 = kf_exact(G)
            got = []
            for v in list(G.nodes()):
                Gv = G.copy()
                Gv.remove_node(v)
                if nx.is_connected(Gv) and kf_exact(Gv) == K0:
                    got.append(v)
            want = {('u', i, j) for i in range(a) for j in range(b)}
            C.check(set(got) == want, 'a=%d b=%d n=%d: K(G) = the %d block vertices (computed for all n vertices)' % (
                a, b, G.number_of_nodes(), a*b))
    elif mode == 'family':
        c = construct(1, 100, 2, 2)
        spec = (c['a'], c['b'], [(2*q, j) for q, j in zip(c['Q'], c['j'])], [(2*c['B'], c['B']*c['k'])], (c['l'], c['t']))
        test_classes(C, 'family beta=1/100 b=2 s=2', spec)
    elif mode in EXAMPLES:
        test_classes(C, mode, EXAMPLES[mode])
    elif mode == 'E3quick':
        test_classes(C, 'E3', EXAMPLES['E3'], only=('block', 'hub', 'root'))
    C.finish()
