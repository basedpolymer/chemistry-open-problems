#!/usr/bin/env python3
"""find_caps.py -- caps of the thinnest "6-row" nanotube, harvested from small fullerene duals.

A ROW of a fullerene dual T is an induced 6-cycle c_0..c_5 of degree-6 vertices (hexagons) such that T - C has
exactly two components and every c_j has exactly two neighbours in each component (the local structure of a
row of the triangulated cylinder in which row r vertex j is adjacent to row r+1 vertices j, j+1).  Each
component is a CAP.  For every cap Z found (n <= NMAX) we record the patch H = T[Z u C] (graph6, nauty canonical
form), |Z|, |Aut(H)| (networkx VF2), and the connected components of the pentagon graph T5(Z).

Usage: python3 find_caps.py NMAX [out.json]
Environment: python3.13, networkx, nauty labelg, buckygen.  Runtime ~1 min for NMAX = 70.
"""
import json
import os
import subprocess
import sys

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from planarcode import read_planar_code, to_graph6  # noqa: E402
from scan_driver import BUCKYGEN, LABELG  # noqa: E402


def rows_of(adj):
    m = len(adj)
    deg = [len(a) for a in adj]
    six = [v for v in range(m) if deg[v] == 6]
    S6 = set(six)
    found = set()
    for s in six:
        # paths s = p0 - p1 - ... - p5 - s through degree-6 vertices, p0 minimal
        stack = [[s]]
        while stack:
            path = stack.pop()
            if len(path) == 6:
                if s in adj[path[-1]]:
                    cyc = path
                    key = frozenset(cyc)
                    if key in found:
                        continue
                    # induced?
                    ok = all(sum(1 for w in adj[v] if w in key) == 2 for v in cyc)
                    if ok:
                        found.add(key)
                        yield cyc
                continue
            for w in adj[path[-1]]:
                if w in S6 and w not in path and w > s:
                    stack.append(path + [w])


def caps_of(adj):
    m = len(adj)
    G = nx.Graph()
    G.add_nodes_from(range(m))
    G.add_edges_from((v, w) for v in range(m) for w in adj[v] if v < w)
    for cyc in rows_of(adj):
        C = set(cyc)
        H = G.subgraph([v for v in range(m) if v not in C])
        comps = [set(c) for c in nx.connected_components(H)]
        if len(comps) != 2:
            continue
        if not all(all(sum(1 for w in adj[v] if w in comp) == 2 for v in cyc) for comp in comps):
            continue
        for comp in comps:
            yield cyc, comp


def main():
    nmax = int(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else None
    caps = {}
    for n in [24] + list(range(26, nmax + 1, 2)):
        data = subprocess.run([BUCKYGEN, '-q', '%dd' % n], capture_output=True, check=True).stdout
        for adj in read_planar_code(data):
            for cyc, Z in caps_of(adj):
                verts = sorted(Z) + list(cyc)
                loc = {v: i for i, v in enumerate(verts)}
                h = [[loc[w] for w in adj[v] if w in loc] for v in verts]
                g6 = to_graph6(h)
                can = subprocess.run([LABELG, '-q'], input=(g6 + '\n').encode(), capture_output=True,
                                     check=True).stdout.decode().split()[0]
                if can in caps:
                    continue
                Hg = nx.Graph()
                Hg.add_nodes_from(range(len(h)))
                Hg.add_edges_from((i, j) for i in range(len(h)) for j in h[i] if i < j)
                naut = sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(Hg, Hg).isomorphisms_iter())
                P = [v for v in Z if len(adj[v]) == 5]
                T5 = nx.Graph()
                T5.add_nodes_from(P)
                T5.add_edges_from((v, w) for v in P for w in adj[v] if w in P and v < w)
                comps = sorted((len(c), T5.subgraph(c).number_of_edges()) for c in nx.connected_components(T5))
                caps[can] = {'size': len(Z), 'n_found': n, 'aut': naut, 'pentagons': len(P),
                             'T5_components_(vertices,edges)': comps, 'patch_graph6': g6,
                             'patch_order': 'cap vertices first, then the row c_0..c_5 in cyclic order'}
        print('n=%d: %d distinct caps so far' % (n, len(caps)), flush=True)
    res = sorted(caps.values(), key=lambda c: (c['size'], c['aut']))
    for c in res:
        print(c['size'], c['size'] % 6, 'aut', c['aut'], 'T5', c['T5_components_(vertices,edges)'], 'first n', c['n_found'])
    if out:
        with open(out, 'w') as f:
            json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
