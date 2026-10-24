#!/usr/bin/env python3
"""search_rotation_caps.py -- caps of zigzag (L,0) tubes whose rotation by s steps preserves the spectrum of ANY
fullerene dual containing them (strong cap condition), harvested from all fullerene duals with n <= NMAX.

For a cap Y attached to a zigzag row C = c_0..c_{L-1} (c_j ~ y_j, y_{j+1}, y_j = common Y-neighbour of c_{j-1},
c_j), let B_s be the Y-C incidence after rotating Y by s steps and W_s(x) = B_s^T adj(xI - A_Y) B_s.
If W_s = W_0 then, for every graph R glued to C on the other side, the two fullerene duals Y u R and
rot_s(Y) u R have the same characteristic polynomial (Schur complement).  We report every such (cap, s) for
which the rotation is NOT induced by a symmetry of the cap (i.e. the rotated attachment differs).

Usage: python3 search_rotation_caps.py NMAX [LMIN LMAX] [out.json]
"""
import json
import os
import subprocess
import sys

import flint
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from planarcode import read_planar_code, to_graph6  # noqa: E402
from scan_driver import BUCKYGEN, LABELG  # noqa: E402


def canon1(adj):
    r = subprocess.run([LABELG, '-q'], input=(to_graph6(adj) + '\n').encode(), capture_output=True, check=True)
    return r.stdout.decode().split()[0]


def zigzag_rows(adj, Lmin, Lmax):
    m = len(adj)
    six = set(v for v in range(m) if len(adj[v]) == 6)
    seen = set()
    for s in sorted(six):
        stack = [[s]]
        while stack:
            path = stack.pop()
            last = path[-1]
            if len(path) >= Lmin and s in adj[last]:
                key = frozenset(path)
                if key not in seen:
                    seen.add(key)
                    yield list(path)
            if len(path) < Lmax:
                for w in adj[last]:
                    if w in six and w > s and w not in path and not any(u in adj[w] for u in path[1:-1]):
                        stack.append(path + [w])


def cap_sides(adj, cyc):
    C = set(cyc)
    m = len(adj)
    G = nx.Graph()
    G.add_nodes_from(v for v in range(m) if v not in C)
    G.add_edges_from((v, w) for v in range(m) if v not in C for w in adj[v] if w not in C)
    comps = [set(c) for c in nx.connected_components(G)]
    if len(comps) != 2:
        return []
    out = []
    for Y in comps:
        if not all(sum(1 for w in adj[v] if w in Y) == 2 for v in cyc):
            return []
        L = len(cyc)
        y = [None] * L
        ok = True
        for j in range(L):
            com = [w for w in adj[cyc[j]] if w in Y and w in adj[cyc[(j + 1) % L]]]
            if len(com) != 1:
                ok = False
                break
            y[(j + 1) % L] = com[0]
        if not ok:
            return []
        if not all(set(w for w in adj[cyc[j]] if w in Y) == {y[j], y[(j + 1) % L]} for j in range(L)):
            return []
        out.append((Y, y))
    return out


def patch(adj, cyc, Y, y, s):
    """the cap Y glued, rotated by s, to an abstract row 0..L-1 (row vertices get labels after Y)."""
    Yl = sorted(Y)
    loc = {v: i for i, v in enumerate(Yl)}
    k, L = len(Yl), len(cyc)
    h = [[loc[w] for w in adj[v] if w in loc] for v in Yl] + [[] for _ in range(L)]
    for j in range(L):
        h[k + j] += [k + (j + 1) % L, k + (j - 1) % L]
        for w in (y[(j - s) % L], y[(j - s + 1) % L]):
            h[k + j].append(loc[w])
            h[loc[w]].append(k + j)
    return h, Yl, loc


def W_equal(adj, Y, y, L, s):
    Yl = sorted(Y)
    loc = {v: i for i, v in enumerate(Yl)}
    k = len(Yl)

    def W_at(x0, shift):
        M = flint.fmpq_mat(k, k, [(-1 if Yl[j] in adj[Yl[i]] else 0) + (x0 if i == j else 0)
                                   for i in range(k) for j in range(k)])
        d = M.det()
        B = flint.fmpq_mat(k, L, [0] * (k * L))
        for j in range(L):
            for w in (y[(j - shift) % L], y[(j - shift + 1) % L]):
                B[loc[w], j] = 1
        return (B.transpose() * M.inv() * B) * d
    return all(W_at(1000 + 7 * t, s) == W_at(1000 + 7 * t, 0) for t in range(k + 3))


def labelled_equal(h0, h1):
    """are the two attachments identical up to an automorphism of the cap fixing the row pointwise?  We test the
    weaker, sufficient-for-'different' criterion: is there an isomorphism h0 -> h1 that is the identity on the row?
    (row vertices are the last L labels).  Implemented with VF2 + node attributes."""
    G0, G1 = nx.Graph(), nx.Graph()
    n = len(h0)
    for G, h in ((G0, h0), (G1, h1)):
        G.add_nodes_from(range(n))
        G.add_edges_from((i, j) for i in range(n) for j in h[i] if i < j)
    return G0, G1


def main():
    nmax = int(sys.argv[1])
    Lmin = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    Lmax = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    out = sys.argv[4] if len(sys.argv) > 4 else None
    seen = {}
    hits = []
    for n in [24] + list(range(26, nmax + 1, 2)):
        data = subprocess.run([BUCKYGEN, '-q', '%dd' % n], capture_output=True, check=True).stdout
        for adj in read_planar_code(data):
            for cyc in zigzag_rows(adj, Lmin, Lmax):
                for Y, y in cap_sides(adj, cyc):
                    L = len(cyc)
                    h0, Yl, loc = patch(adj, cyc, Y, y, 0)
                    key = canon1(h0)
                    if key in seen:
                        continue
                    seen[key] = (n, len(Y), L)
                    k = len(Yl)
                    for s in range(1, L):
                        h1, _, _ = patch(adj, cyc, Y, y, s)
                        # is the rotated attachment the same labelled patch (row fixed pointwise)?
                        G0, G1 = labelled_equal(h0, h1)
                        nm = lambda a, b: a['c'] == b['c']
                        for G in (G0, G1):
                            for v in G.nodes:
                                G.nodes[v]['c'] = v if v >= k else -1
                        same = nx.is_isomorphic(G0, G1, node_match=nm)
                        if same:
                            continue
                        if W_equal(adj, Y, y, L, s):
                            P = sum(1 for v in Y if len(adj[v]) == 5)
                            item = {'n_found': n, 'L': L, 'cap_size': len(Y), 'shift': s, 'pentagons': P,
                                    'patch_graph6_shift0': to_graph6(h0), 'patch_graph6_shift': to_graph6(h1)}
                            hits.append(item)
                            print('HIT', item, flush=True)
        print('n=%d caps=%d hits=%d' % (n, len(seen), len(hits)), flush=True)
    if out:
        with open(out, 'w') as f:
            json.dump({'caps_by_L': {str(L): sum(1 for v in seen.values() if v[2] == L) for L in range(Lmin, Lmax + 1)},
                       'hits': hits}, f, indent=1)


if __name__ == '__main__':
    main()
