#!/usr/bin/env python3
"""check_family.py -- exact certificate for the helical family H_L, L = 20..LMAX (default 130, n = 60..280).

For every L the script checks with integer arithmetic (PASS/FAIL, exit 1 on failure):
  C1  H_L is a simple plane triangulation (planar, 3V-6 edges, closed surface), degrees 4 and 6 only, six vertices of
      degree 4, pairwise non-adjacent;
  C2  N_b^T N_b = N_w^T N_w = A(H) + D(H)/2 (exact integer matrices);
  C3  N_c N_c^T = 3I + A(T_c), c = b, w;
  C4  T_b and T_w are fullerene duals (planar, simple, 3m-6 edges, twelve vertices of degree 5, the others of degree 6,
      minimum degree 5, 3-connected), n = 2L + 20 = 2(m-2);
  C5  charpoly A(T_b) = charpoly A(T_w) (exact, python-flint);
  C6  the numbers of hexagons with exactly j pentagon neighbours differ; for L >= 25 they are the constants
      (black) 11, 6, 5, 1 and (white) 8, 9, 4, 1 for j = 1, 2, 3, 4;
  C7  nauty canonical forms of T_b, T_w differ.
Theorem B of PROGRESS.md needs C1-C7 only at L = 20..25 (the rest follows from the arguments given there); the larger
range is a cross-check.

Usage: python check_family.py [LMAX]
Environment: python 3.12, python-flint 0.9, networkx, nauty labelg.  Runtime about 1 min for LMAX = 130.
"""
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family import *

FAILED = []


def report(name, ok, detail=''):
    print('%s  %s%s' % ('PASS' if ok else 'FAIL', name, ('  ' + detail) if detail else ''), flush=True)
    if not ok:
        FAILED.append(name)


def fm(M):
    return flint.fmpz_mat(len(M), len(M[0]), [x for r in M for x in r])


def main():
    lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 130
    ok = {k: True for k in ('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7')}
    hist = {}
    for L in range(20, lmax + 1):
        faces, col = colour_faces(L)
        V = L + 14
        h = adjacency_H(L)
        G = nx.Graph()
        G.add_nodes_from(range(V))
        G.add_edges_from((v, w) for v in range(V) for w in h[v] if v < w)
        planar, _ = nx.check_planarity(G)
        deg = collections.Counter(len(a) for a in h)
        four = [v for v in range(V) if len(h[v]) == 4]
        ok['C1'] &= planar and G.number_of_edges() == 3 * V - 6 and set(deg) == {4, 6} and deg[4] == 6 and \
            all(w not in h[v] for v in four for w in four)
        # every vertex link is one cycle, Euler characteristic 2
        link = collections.defaultdict(list)
        for f in faces:
            for v in f:
                link[v].append(tuple(x for x in f if x != v))
        for v, es in link.items():
            LG = nx.Graph(es)
            ok['C1'] &= all(d == 2 for _, d in LG.degree()) and nx.is_connected(LG)
        ok['C1'] &= V - G.number_of_edges() + len(faces) == 2
        Nb, rb = incidence(faces, col, 0, V)
        Nw, rw = incidence(faces, col, 1, V)
        Nb, Nw = fm(Nb), fm(Nw)
        lhs = Nb.transpose() * Nb
        Aplus = flint.fmpz_mat(V, V, [(1 if j in h[i] else 0) + (len(h[i]) // 2 if i == j else 0)
                                       for i in range(V) for j in range(V)])
        ok['C2'] &= lhs == Nw.transpose() * Nw == Aplus and lhs == Aplus
        tb, _ = colour_graph(faces, col, 0)
        tw, _ = colour_graph(faces, col, 1)
        m = len(tb)
        for t, N in ((tb, Nb), (tw, Nw)):
            Mt = flint.fmpz_mat(m, m, [(1 if j in t[i] else 0) + (3 if i == j else 0) for i in range(m) for j in range(m)])
            ok['C3'] &= (N * N.transpose()) == Mt
            GT = nx.Graph()
            GT.add_nodes_from(range(m))
            GT.add_edges_from((v, w) for v in range(m) for w in t[v] if v < w)
            pl, _ = nx.check_planarity(GT)
            d = collections.Counter(len(a) for a in t)
            ok['C4'] &= pl and GT.number_of_edges() == 3 * m - 6 and set(d) == {5, 6} and d[5] == 12 and \
                nx.node_connectivity(GT) >= 3 and 2 * (m - 2) == 2 * L + 20
        ok['C5'] &= charpoly(tb) == charpoly(tw)
        hb, hw = hex_histogram(tb), hex_histogram(tw)
        hist[L] = (hb, hw)
        ok['C6'] &= hb != hw
        if L >= 25:
            ok['C6'] &= hb[1:5] == (11, 6, 5, 1) and hw[1:5] == (8, 9, 4, 1)
        can = labelg_canon([tb, tw])
        ok['C7'] &= can[0] != can[1]
        if L % 20 == 0 or L < 26:
            print('L=%d n=%d hist black %s white %s' % (L, 2 * L + 20, hb, hw), flush=True)
    names = {'C1': 'H_L is a triangulation with six non-adjacent 4-vertices, the others of degree 6',
             'C2': 'N_b^T N_b = N_w^T N_w = A(H) + D(H)/2', 'C3': 'N_c N_c^T = 3I + A(T_c)',
             'C4': 'T_b, T_w are fullerene duals with n = 2L + 20', 'C5': 'equal integer characteristic polynomials',
             'C6': 'different hexagon histograms; constants for L >= 25', 'C7': 'nauty: T_b and T_w are not isomorphic'}
    for k in sorted(ok):
        report('%s %s (L = 20..%d)' % (k, names[k], lmax), ok[k])
    if FAILED:
        sys.exit(1)


if __name__ == '__main__':
    main()
