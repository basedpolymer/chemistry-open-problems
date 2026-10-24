#!/usr/bin/env python3
"""classify_pairs_general.py -- like classify_pairs.py, but for every separating induced cycle C (length <= LMAX,
vertices of any degree) and every shift s that preserves the cyclic sequence (deg(c_j), #Y-neighbours of c_j):
re-glue the side Y rotated by s along C and compare with the partner (nauty).  Reports the mechanism with the
smallest rotated side and whether the strong cap condition W_s = W_0 holds.

Usage: python3 classify_pairs_general.py NMIN NMAX [LMAX] [out.json]
"""
import json
import os
import sys

import flint
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from planarcode import from_graph6  # noqa: E402
from classify_pairs import canon, sides  # noqa: E402
from verify_pairs import adj_from_graph6_embedded  # noqa: E402

CERT = os.path.join(ROOT, 'certificates')


def induced_cycles(adj, Lmax):
    m = len(adj)
    seen = set()
    for s in range(m):
        stack = [[s]]
        while stack:
            path = stack.pop()
            last = path[-1]
            if len(path) >= 4 and s in adj[last]:
                key = frozenset(path)
                if key not in seen:
                    seen.add(key)
                    yield list(path)
            if len(path) < Lmax:
                for w in adj[last]:
                    if w > s and w not in path and not any(u in adj[w] for u in path[1:-1]):
                        if len(path) >= 2 and len(path) < 3 and False:
                            continue
                        stack.append(path + [w])


def blocks(emb, cyc, Y):
    """for each c_j, its Y-neighbours in rotation order, from the one shared with c_{j-1} to the one shared
    with c_{j+1}; None if the structure is not a clean boundary."""
    L = len(cyc)
    out = []
    for j in range(L):
        v = cyc[j]
        rot = emb[v]
        k = len(rot)
        isY = [w in Y for w in rot]
        if not any(isY) or all(isY):
            return None
        # contiguous run of Y
        start = [i for i in range(k) if isY[i] and not isY[(i - 1) % k]]
        if len(start) != 1:
            return None
        i0 = start[0]
        run = []
        i = i0
        while isY[i % k]:
            run.append(rot[i % k])
            i += 1
        prev, nxt = cyc[(j - 1) % L], cyc[(j + 1) % L]
        if run[0] in emb[prev] and run[-1] in emb[nxt]:
            pass
        elif run[-1] in emb[prev] and run[0] in emb[nxt]:
            run = run[::-1]
        else:
            return None
        out.append(run)
    for j in range(L):
        if out[j][-1] != out[(j + 1) % L][0]:
            return None
    return out


def rotate_general(adj, cyc, bl, s):
    L = len(cyc)
    new = [set(a) for a in adj]
    for j in range(L):
        for w in bl[j]:
            new[cyc[j]].discard(w)
            new[w].discard(cyc[j])
    for j in range(L):
        for w in bl[(j - s) % L]:
            new[cyc[j]].add(w)
            new[w].add(cyc[j])
    return [sorted(a) for a in new]


def strong_general(adj, cyc, Y, bl, s):
    L = len(cyc)
    Yl = sorted(Y)
    loc = {v: i for i, v in enumerate(Yl)}
    k = len(Yl)

    def W_at(x0, shift):
        M = flint.fmpq_mat(k, k, [(-1 if Yl[j] in adj[Yl[i]] else 0) + (x0 if i == j else 0)
                                   for i in range(k) for j in range(k)])
        d = M.det()
        Minv = M.inv()
        B = flint.fmpq_mat(k, L, [0] * (k * L))
        for j in range(L):
            for w in bl[(j - shift) % L]:
                B[loc[w], j] = 1
        return (B.transpose() * Minv * B) * d
    pts = [1000 + 7 * t for t in range(k + 3)]
    return all(W_at(x0, s) == W_at(x0, 0) for x0 in pts)


def main():
    nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    Lmax = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    out = sys.argv[4] if len(sys.argv) > 4 else None
    pairs = json.load(open(os.path.join(CERT, 'cospectral_pairs_T.json')))
    res = []
    for c in pairs:
        if not (nmin <= c['n'] <= nmax) or len(c['graph6_dual']) != 2:
            continue
        g1 = from_graph6(c['graph6_dual'][0])
        emb = adj_from_graph6_embedded(c['graph6_dual'][0])
        target = canon([from_graph6(c['graph6_dual'][1])])[0]
        deg = [len(a) for a in g1]
        hits = []
        for cyc in induced_cycles(g1, Lmax):
            sd = sides(g1, cyc) if False else None
            C = set(cyc)
            G = nx.Graph()
            G.add_nodes_from(v for v in range(len(g1)) if v not in C)
            G.add_edges_from((v, w) for v in range(len(g1)) if v not in C for w in g1[v] if w not in C)
            comps = [set(x) for x in nx.connected_components(G)]
            if len(comps) != 2:
                continue
            L = len(cyc)
            for Y in comps:
                bl = blocks(emb, cyc, Y)
                if bl is None:
                    continue
                pat = [(deg[cyc[j]], len(bl[j])) for j in range(L)]
                for s in range(1, L):
                    if any(pat[j] != pat[(j - s) % L] for j in range(L)):
                        continue
                    g1r = rotate_general(g1, cyc, bl, s)
                    if canon([g1r])[0] == target:
                        hits.append({'L': L, 'pattern_Y_nbrs': [p[1] for p in pat], 'cycle_degrees': [p[0] for p in pat],
                                     'side_size': len(Y), 'other_side_size': len(g1) - len(Y) - L, 'shift': s,
                                     'strong': strong_general(g1, cyc, Y, bl, s)})
        best = None
        if hits:
            strong = [h for h in hits if h['strong']]
            best = min(strong or hits, key=lambda h: (h['side_size'], h['L']))
        res.append({'n': c['n'], 'atlas': c.get('atlas_numbers'), 'graph6_dual': c['graph6_dual'],
                    'mechanism': best, 'n_hits': len(hits), 'n_strong_hits': sum(1 for h in hits if h['strong'])})
        print(c['n'], c.get('atlas_numbers'), best, len(hits), flush=True)
    if out:
        with open(out, 'w') as f:
            json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
