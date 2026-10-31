"""fdc.py -- small library for the all-n work on Conjecture 1 of Bille-Buchstaber-Spodarev.
Plane triangulations are given by rotation systems: rot[v] = neighbours of v in cyclic order.
All exact arithmetic is python-flint.
"""
import os
import subprocess
import sys

import flint
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from planarcode import from_graph6, to_graph6, read_planar_code  # noqa: E402,F401


def rotation_system(adj):
    """rotation system (from the unique planar embedding of a 3-connected planar graph)."""
    G = nx.Graph()
    G.add_nodes_from(range(len(adj)))
    G.add_edges_from((v, w) for v in range(len(adj)) for w in adj[v] if v < w)
    ok, emb = nx.check_planarity(G)
    assert ok
    return [list(emb.neighbors_cw_order(v)) for v in range(len(adj))]


def adjset(rot):
    return [set(r) for r in rot]


def charpoly(adj, verts=None):
    """integer characteristic polynomial of the subgraph induced on verts (fmpz_poly)."""
    if verts is None:
        verts = list(range(len(adj)))
    verts = list(verts)
    idx = {v: i for i, v in enumerate(verts)}
    k = len(verts)
    M = [0] * (k * k)
    for v in verts:
        for w in adj[v]:
            if w in idx:
                M[idx[v] * k + idx[w]] = 1
    return flint.fmpz_mat(k, k, M).charpoly() if k else flint.fmpz_poly([1])


def W_matrix_at(adj, I, Bdry, x0):
    """value at x0 of p_I(x) and W(x) = B^T adj(xI - A_I) B, B = incidence I x Bdry."""
    I = list(I)
    loc = {v: i for i, v in enumerate(I)}
    k = len(I)
    M = [0] * (k * k)
    for v in I:
        i = loc[v]
        M[i * k + i] = x0
        for w in adj[v]:
            if w in loc:
                M[i * k + loc[w]] = -1
    M = flint.fmpq_mat(k, k, M)
    d = M.det()
    B = [0] * (k * len(Bdry))
    for j, b in enumerate(Bdry):
        for w in adj[b]:
            if w in loc:
                B[loc[w] * len(Bdry) + j] = 1
    B = flint.fmpq_mat(k, len(Bdry), B)
    return d, (B.transpose() * M.inv() * B) * d


def W_identity(adj1, I1, Bd1, adj2, I2, Bd2):
    """exact test: same p_I and same W(x) on the matched boundaries Bd1[j] <-> Bd2[j].
    Entries of W have degree <= |I| - 1, p_I degree |I|, evaluated at |I| + 2 points beyond the spectrum."""
    if len(I1) != len(I2):
        return False
    k = len(I1)
    for t in range(k + 2):
        x0 = 1000 + 7 * t
        d1, W1 = W_matrix_at(adj1, I1, Bd1, x0)
        d2, W2 = W_matrix_at(adj2, I2, Bd2, x0)
        if d1 != d2 or W1 != W2:
            return False
    return True


def walk(rot1, rot2, u, v, a, b, sgn):
    """grow the isomorphism u->a, v->b (orientation sign sgn: +1 same, -1 reversed).  Returns the dict of the
    mapped vertices and the set of vertices whose whole star is matched."""
    deg1 = [len(r) for r in rot1]
    deg2 = [len(r) for r in rot2]
    if deg1[u] != deg2[a] or deg1[v] != deg2[b]:
        return {}, set()
    mp = {u: a}
    inv = {a: u}
    done = set()
    bad = set()
    stack = [(u, v, a, b)]
    while stack:
        u, v, a, b = stack.pop()
        if u in done or u in bad:
            continue
        k = deg1[u]
        if k != deg2[a]:
            bad.add(u)
            continue
        i = rot1[u].index(v)
        j = rot2[a].index(b)
        pairs = [(rot1[u][(i + t) % k], rot2[a][(j + sgn * t) % k]) for t in range(k)]
        ok = True
        for x, y in pairs:
            if (x in mp and mp[x] != y) or (y in inv and inv[y] != x):
                ok = False
                break
        if not ok:
            bad.add(u)
            continue
        for x, y in pairs:
            mp[x] = y
            inv[y] = x
        done.add(u)
        for x, y in pairs:
            if x not in done and x not in bad and deg1[x] == deg2[y]:
                stack.append((x, u, y, a))
    return mp, done


def best_common_patch(rot1, rot2):
    """largest fully matched region between two triangulations, anchored at a pentagon.
    Returns (mp, core), core = vertices of rot1 with fully matched star."""
    pent1 = [v for v in range(len(rot1)) if len(rot1[v]) == 5]
    pent2 = [v for v in range(len(rot2)) if len(rot2[v]) == 5]
    best = ({}, set())
    for p in pent1:
        q0 = rot1[p][0]
        for r in pent2:
            for b in rot2[r]:
                for sgn in (1, -1):
                    mp, core = walk(rot1, rot2, p, q0, r, b, sgn)
                    if len(core) > len(best[1]):
                        best = (mp, core)
    return best


def labelg_canon(adjs):
    """nauty canonical graph6 strings."""
    inp = ''.join(to_graph6(a) + '\n' for a in adjs).encode()
    r = subprocess.run(['labelg', '-q'], input=inp, capture_output=True, check=True)
    return r.stdout.decode().split()
