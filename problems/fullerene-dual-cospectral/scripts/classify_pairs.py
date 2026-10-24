#!/usr/bin/env python3
"""classify_pairs.py -- which A(T)-cospectral pairs are "cap rotations" along a zigzag row?

A zigzag row of length L of a fullerene dual T is an induced L-cycle c_0..c_{L-1} of degree-6 vertices that
separates T into two components Y, Z with every c_j having exactly two neighbours in each component; then
consecutive c_j, c_{j+1} have exactly one common neighbour y_{j+1} in Y and c_j ~ y_j, y_{j+1}.  ROTATING Y by s
steps means replacing the edges c_j y_j, c_j y_{j+1} by c_j y_{j-s}, c_j y_{j-s+1}; the result is again a
fullerene dual.  For each certified pair (G1, G2) we search a row of G1, a side Y and a shift s such that the
rotated graph is isomorphic to G2 (nauty canonical forms), and we test the STRONG CAP CONDITION for Y:
equal det(xI - A_Y) for both gluings is automatic (same Y); the boundary matrix
W_s(x) = B_s^T adj(xI - A_Y) B_s (L x L, polynomial entries) must satisfy W_s = W_0, where B_s is the Y-C
adjacency after rotation by s.  If it holds, rotating Y by s preserves the characteristic polynomial for ANY
graph on the other side of the row (Schur complement lemma, see result.md).

Usage: python3 classify_pairs.py NMAX [out.json]
Environment: python3.13, python-flint 0.9, networkx, nauty labelg.
"""
import json
import os
import subprocess
import sys

import flint
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from planarcode import from_graph6, to_graph6  # noqa: E402
from scan_driver import LABELG  # noqa: E402

CERT = os.path.join(ROOT, 'certificates')


def canon(adjs):
    g6 = [to_graph6(a) for a in adjs]
    r = subprocess.run([LABELG, '-q'], input=('\n'.join(g6) + '\n').encode(), capture_output=True, check=True)
    return r.stdout.decode().split()


def zigzag_rows(adj, Lmax=12):
    m = len(adj)
    deg = [len(a) for a in adj]
    six = set(v for v in range(m) if deg[v] == 6)
    seen = set()
    for s in sorted(six):
        stack = [[s]]
        while stack:
            path = stack.pop()
            last = path[-1]
            if len(path) >= 3 and s in adj[last]:
                key = frozenset(path)
                if key not in seen:
                    seen.add(key)
                    C = key
                    if all(sum(1 for w in adj[v] if w in C) == 2 for v in path):
                        yield list(path)
            if len(path) < Lmax:
                for w in adj[last]:
                    if w in six and w > s and w not in path:
                        # keep the path induced: w may touch only `last` (and s, which closes the cycle)
                        if any(u in adj[w] for u in path[1:-1]):
                            continue
                        stack.append(path + [w])


def sides(adj, cyc):
    C = set(cyc)
    m = len(adj)
    G = nx.Graph()
    G.add_nodes_from(v for v in range(m) if v not in C)
    G.add_edges_from((v, w) for v in range(m) if v not in C for w in adj[v] if w not in C)
    comps = [set(c) for c in nx.connected_components(G)]
    if len(comps) != 2:
        return None
    for comp in comps:
        if not all(sum(1 for w in adj[v] if w in comp) == 2 for v in cyc):
            return None
    return comps


def rotate(adj, cyc, Y, s):
    L = len(cyc)
    y = [None] * L
    for j in range(L):
        a, b = cyc[j], cyc[(j + 1) % L]
        com = [w for w in adj[a] if w in Y and w in adj[b]]
        if len(com) != 1:
            return None, None
        y[(j + 1) % L] = com[0]
    for j in range(L):
        if set(w for w in adj[cyc[j]] if w in Y) != {y[j], y[(j + 1) % L]}:
            return None, None
    new = [set(a) for a in adj]
    for j in range(L):
        for w in (y[j], y[(j + 1) % L]):
            new[cyc[j]].discard(w)
            new[w].discard(cyc[j])
    for j in range(L):
        for w in (y[(j - s) % L], y[(j - s + 1) % L]):
            new[cyc[j]].add(w)
            new[w].add(cyc[j])
    return [sorted(a) for a in new], y


def strong_condition(adj, cyc, Y, y, s):
    """W_s == W_0 with W = B^T adj(xI - A_Y) B over Z[x] (exact, python-flint)."""
    L = len(cyc)
    Yl = sorted(Y)
    loc = {v: i for i, v in enumerate(Yl)}
    k = len(Yl)
    # adj(xI - A_Y) = det * inverse; compute via fmpz_poly matrix by evaluating at enough integer points?
    # simpler and exact: work over Q(x) with fmpq_mat at many rational points and compare -- here we use
    # polynomial identity testing at k+L+2 distinct integers x0, which is EXACT for polynomials of degree < k+L+2.
    def W_at(x0, shift):
        M = flint.fmpq_mat(k, k, [(-1 if Yl[j] in adj[Yl[i]] else 0) + (x0 if i == j else 0)
                                   for i in range(k) for j in range(k)])
        d = M.det()
        Minv = M.inv()
        B = flint.fmpq_mat(k, L, [0] * (k * L))
        for j in range(L):
            for w in (y[(j - shift) % L], y[(j - shift + 1) % L]):
                B[loc[w], j] = 1
        return (B.transpose() * Minv * B) * d
    pts = [1000 + 7 * t for t in range(k + 3)]   # adj entries have degree <= k-1; x0 large avoids eigenvalues
    return all(W_at(x0, s) == W_at(x0, 0) for x0 in pts)


def boundary_y(adj, cyc, Y):
    L = len(cyc)
    y = [None] * L
    for j in range(L):
        a, b = cyc[j], cyc[(j + 1) % L]
        com = [w for w in adj[a] if w in Y and w in adj[b]]
        if len(com) != 1:
            return None
        y[(j + 1) % L] = com[0]
    for j in range(L):
        if set(w for w in adj[cyc[j]] if w in Y) != {y[j], y[(j + 1) % L]}:
            return None
    return y


def attach_map(kind, L, par):
    """y-index pair attached to c_j: rotation by s: (j-s, j-s+1); reflection with constant c: (c-j, c-j+1)."""
    if kind == 'rot':
        return lambda j: ((j - par) % L, (j - par + 1) % L)
    return lambda j: ((par - j) % L, (par - j + 1) % L)


def reglue(adj, cyc, y, f):
    L = len(cyc)
    new = [set(a) for a in adj]
    for j in range(L):
        for w in (y[j], y[(j + 1) % L]):
            new[cyc[j]].discard(w)
            new[w].discard(cyc[j])
    for j in range(L):
        for t in f(j):
            new[cyc[j]].add(y[t])
            new[y[t]].add(cyc[j])
    return [sorted(a) for a in new]


def strong_map(adj, cyc, Y, y, f):
    L = len(cyc)
    Yl = sorted(Y)
    loc = {v: i for i, v in enumerate(Yl)}
    k = len(Yl)
    f0 = attach_map('rot', L, 0)

    def W_at(x0, g):
        M = flint.fmpq_mat(k, k, [(-1 if Yl[j] in adj[Yl[i]] else 0) + (x0 if i == j else 0)
                                   for i in range(k) for j in range(k)])
        d = M.det()
        B = flint.fmpq_mat(k, L, [0] * (k * L))
        for j in range(L):
            for t in g(j):
                B[loc[y[t]], j] = 1
        return (B.transpose() * M.inv() * B) * d
    return all(W_at(1000 + 7 * t, f) == W_at(1000 + 7 * t, f0) for t in range(k + 3))


def main():
    nmax = int(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else None
    pairs = json.load(open(os.path.join(CERT, 'cospectral_pairs_T.json')))
    res = []
    for c in pairs:
        if c['n'] > nmax or len(c['graph6_dual']) != 2:
            continue
        g1 = from_graph6(c['graph6_dual'][0])
        g2 = from_graph6(c['graph6_dual'][1])
        target = canon([g2])[0]
        hits = []
        for cyc in zigzag_rows(g1):
            sd = sides(g1, cyc)
            if not sd:
                continue
            for Y in sd:
                y = boundary_y(g1, cyc, Y)
                if y is None:
                    continue
                L = len(cyc)
                for kind, rng in (('rot', range(1, L)), ('refl', range(L))):
                    for par in rng:
                        f = attach_map(kind, L, par)
                        g1r = reglue(g1, cyc, y, f)
                        if canon([g1r])[0] == target:
                            hits.append({'L': L, 'side_size': len(Y), 'other_side_size': len(g1) - len(Y) - L,
                                         'kind': kind, 'param': par, 'strong': strong_map(g1, cyc, Y, y, f)})
        best = None
        if hits:
            strong = [h for h in hits if h['strong']]
            best = min(strong or hits, key=lambda h: h['side_size'])
        res.append({'n': c['n'], 'atlas': c.get('atlas_numbers'), 'rotation_mechanism': best,
                    'n_hits': len(hits)})
        print(c['n'], c.get('atlas_numbers'), best, flush=True)
    if out:
        with open(out, 'w') as f:
            json.dump(res, f, indent=1)


if __name__ == '__main__':
    main()
