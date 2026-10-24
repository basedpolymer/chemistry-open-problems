#!/usr/bin/env python3
"""describe_pairs.py -- human-readable certificates for the headline pairs, computed from graph6 strings only.

For each pair: Atlas-style canonical face spirals, graph6 of the duals T_n AND of the fullerene graphs C_n
themselves (cubic graphs, obtained as duals of the planar triangulations), the integer characteristic
polynomials of A(T_n), A(T_n^5), A(T_n^6), A(T_n)+D, of the fullerene graph A(C_n), automorphism group orders,
and a simple non-isomorphism certificate (an isomorphism invariant that differs).
Writes certificates/headline_pairs.json and prints a summary.

Usage: python3 -I describe_pairs.py
Environment: python3.13, python-flint 0.9, networkx.  Runtime < 1 min.
"""
import json
import os
import sys

import flint
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from planarcode import to_graph6  # noqa: E402
from spiral import canonical_spiral  # noqa: E402
from verify_pairs import adj_from_graph6_embedded  # noqa: E402

CERT = os.path.join(ROOT, 'certificates')


def faces(rot):
    pos = [{w: i for i, w in enumerate(rot[v])} for v in range(len(rot))]
    seen, out = set(), []
    for v in range(len(rot)):
        for w in rot[v]:
            if (v, w) in seen:
                continue
            f = []
            a, b = v, w
            while (a, b) not in seen:
                seen.add((a, b))
                f.append(a)
                c = rot[b][(pos[b][a] + 1) % len(rot[b])]
                a, b = b, c
            out.append(f)
    return out


def fullerene_graph(rot):
    """the cubic plane graph dual to the triangulation: one vertex per triangle, adjacent if sharing an edge."""
    F = faces(rot)
    edge2f = {}
    for i, f in enumerate(F):
        for j in range(len(f)):
            e = frozenset((f[j], f[(j + 1) % len(f)]))
            edge2f.setdefault(e, []).append(i)
    adj = [[] for _ in F]
    for e, fs in edge2f.items():
        assert len(fs) == 2
        a, b = fs
        adj[a].append(b)
        adj[b].append(a)
    assert all(len(a) == 3 for a in adj)
    return adj


def cp_mat(adj, verts=None, diag=False):
    vs = list(range(len(adj))) if verts is None else verts
    loc = {v: i for i, v in enumerate(vs)}
    k = len(vs)
    M = [[0] * k for _ in range(k)]
    for v in vs:
        for w in adj[v]:
            if w in loc:
                M[loc[v]][loc[w]] = 1
        if diag:
            M[loc[v]][loc[v]] = len(adj[v])
    return str(flint.fmpz_mat(M).charpoly()) if k else '1'


def nxg(adj):
    G = nx.Graph()
    G.add_nodes_from(range(len(adj)))
    G.add_edges_from((v, w) for v in range(len(adj)) for w in adj[v] if v < w)
    return G


def describe(g6):
    rot = adj_from_graph6_embedded(g6)
    m = len(rot)
    deg = [len(a) for a in rot]
    V5 = [v for v in range(m) if deg[v] == 5]
    V6 = [v for v in range(m) if deg[v] == 6]
    C = fullerene_graph(rot)
    naut = sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(nxg(rot), nxg(rot)).isomorphisms_iter())
    return {'n': 2 * (m - 2), 'canonical_spiral': list(canonical_spiral(rot)), 'graph6_dual': g6,
            'graph6_fullerene': to_graph6(C),
            'charpoly_A_dual': cp_mat(rot), 'charpoly_A_T5': cp_mat(rot, V5), 'charpoly_A_T6': cp_mat(rot, V6),
            'charpoly_A_plus_D_dual': cp_mat(rot, diag=True), 'charpoly_A_fullerene': cp_mat(C),
            'automorphism_group_order': naut}


def main():
    pairsT = json.load(open(os.path.join(CERT, 'cospectral_pairs_T.json')))
    pairsQ = json.load(open(os.path.join(CERT, 'cospectral_pairs_signless_laplacian.json')))
    sel = [('smallest A(T)-cospectral pair (refutes arXiv Conjecture 1(a)); Atlas C42:6 and C42:22',
            [c for c in pairsT if c['n'] == 42][0]),
           ('pair of Bille et al. (C44:37, C44:38), cap-rotation type', [c for c in pairsT if c['n'] == 44][0]),
           ('smallest signless-Laplacian (A+D) cospectral pair: refutes Conjecture 3 of arXiv:2410.19322 for alpha = beta',
            [c for c in pairsQ if c['n'] == 44][0])]
    out = []
    for title, c in sel:
        d = [describe(g) for g in c['graph6_dual']]
        inv = [k for k in ('charpoly_A_T5', 'charpoly_A_T6', 'charpoly_A_plus_D_dual', 'charpoly_A_fullerene',
                           'canonical_spiral') if d[0][k] != d[1][k]]
        out.append({'title': title, 'atlas_numbers': c.get('atlas_numbers'), 'members': d,
                    'invariants_that_differ (non-isomorphism certificates)': inv,
                    'A_dual_charpolys_equal': d[0]['charpoly_A_dual'] == d[1]['charpoly_A_dual'],
                    'A_plus_D_charpolys_equal': d[0]['charpoly_A_plus_D_dual'] == d[1]['charpoly_A_plus_D_dual']})
        print('==', title, c.get('atlas_numbers'))
        for x in d:
            print('   spiral', x['canonical_spiral'], '|Aut|', x['automorphism_group_order'])
        print('   A(T) equal:', out[-1]['A_dual_charpolys_equal'], ' A+D equal:', out[-1]['A_plus_D_charpolys_equal'],
              ' differing invariants:', inv)
    with open(os.path.join(CERT, 'headline_pairs.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
