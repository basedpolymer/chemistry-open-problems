"""family.py -- the helical family H_L of even triangulations (tube (7,1)) and its two colour twins.

H_L has the vertices 0,1,...,L (the helix) and thirteen further vertices
    p0 p1 p2 p3 p4 q1 q2 (first cap), r0 r1 r2 r3 s1 s2 (second cap).
Faces (e = L):
    tube:   B_t = {t, t+1, t+8} and W_t = {t, t+7, t+8}, 0 <= t <= L-8;
    cap 1:  the 20 triangles of CAP1 below, a disk bounded by the 8-cycle 0,1,...,7;
    cap 2:  the 18 triangles of CAP2 below, a disk bounded by the 8-cycle e-7,...,e.
Black faces are those of the proper 2-colouring in which every B_t is black.
"""
import collections
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdc import *  # noqa: F401,F403

CAP1 = [('p0', 0, 1), ('p0', 1, 'p4'), ('p0', 'p4', 'p1'), ('p0', 'p1', 6), ('p0', 6, 7), ('p0', 7, 0),
        ('p1', 'p4', 'q2'), ('p1', 'q2', 'p2'), ('p1', 'p2', 5), ('p1', 5, 6),
        ('p2', 'q2', 'p3'), ('p2', 'p3', 'q1'), ('p2', 'q1', 4), ('p2', 4, 5),
        ('p3', 3, 'q1'), ('p3', 'q2', 'p4'), ('p3', 'p4', 2), ('p3', 2, 3),
        ('p4', 1, 2), ('q1', 3, 4)]
# in CAP2 an integer -k stands for the helix vertex e - k, 0 stands for e
CAP2 = [('r0', -7, -6), ('r0', -6, 'r1'), ('r0', 'r1', 'r3'), ('r0', 'r3', -1), ('r0', -1, 0), ('r0', 0, -7),
        ('r1', 'r3', 'r2'), ('r1', 'r2', 's1'), ('r1', 's1', -5), ('r1', -5, -6),
        ('s1', -5, -4), ('s1', -4, 'r2'),
        ('r2', -4, -3), ('r2', -3, 's2'), ('r2', 's2', 'r3'),
        ('s2', -3, -2), ('s2', -2, 'r3'),
        ('r3', -2, -1)]
NAMES = ['p0', 'p1', 'p2', 'p3', 'p4', 'q1', 'q2', 'r0', 'r1', 'r2', 'r3', 's1', 's2']


def faces_H(L):
    """list of (triple of vertex labels, colour hint) with labels integers: helix t = t, named vertices = L+1+i"""
    idx = {nm: L + 1 + i for i, nm in enumerate(NAMES)}
    F = []
    for t in range(L - 8 + 1):
        F.append(((t, t + 1, t + 8), 'B'))
        F.append(((t, t + 7, t + 8), 'W'))
    for f in CAP1:
        F.append((tuple(idx.get(x, x) for x in f), None))
    for f in CAP2:
        F.append((tuple(idx[x] if isinstance(x, str) else L + x for x in f), None))
    return F


def colour_faces(L):
    """proper 2-colouring of the faces of H_L with all B_t black: returns (faces, colours), colour 0 = black."""
    F = faces_H(L)
    faces = [tuple(sorted(f)) for f, _ in F]
    assert len(set(faces)) == len(faces), 'repeated face'
    ef = collections.defaultdict(list)
    for i, f in enumerate(faces):
        for e in itertools.combinations(f, 2):
            ef[e].append(i)
    assert all(len(v) == 2 for v in ef.values()), 'not a closed surface'
    col = [None] * len(faces)
    i0 = [i for i, (f, c) in enumerate(F) if c == 'B'][0]
    col[i0] = 0
    stack = [i0]
    while stack:
        i = stack.pop()
        for e in itertools.combinations(faces[i], 2):
            j = [x for x in ef[e] if x != i][0]
            if col[j] is None:
                col[j] = 1 - col[i]
                stack.append(j)
            elif col[j] == col[i]:
                raise ValueError('faces of H_L are not 2-colourable')
    assert all(c is not None for c in col)
    for (f, c), cc in zip(F, col):
        if c == 'B':
            assert cc == 0
        if c == 'W':
            assert cc == 1
    return faces, col


def adjacency_H(L):
    faces, _ = colour_faces(L)
    V = L + 1 + len(NAMES)
    adj = [set() for _ in range(V)]
    for f in faces:
        for a, b in itertools.combinations(f, 2):
            adj[a].add(b)
            adj[b].add(a)
    return [sorted(a) for a in adj]


def incidence(faces, col, c, V):
    """N_c: rows = faces of colour c, columns = vertices of H (python lists)."""
    rows = [f for f, k in zip(faces, col) if k == c]
    return [[1 if v in f else 0 for v in range(V)] for f in rows], rows


def colour_graph(faces, col, c):
    """the graph T_c: adjacency lists on the faces of colour c."""
    rows = [f for f, k in zip(faces, col) if k == c]
    at = collections.defaultdict(list)
    for i, f in enumerate(rows):
        for v in f:
            at[v].append(i)
    adj = [set() for _ in rows]
    for v, lst in at.items():
        for a, b in itertools.combinations(lst, 2):
            adj[a].add(b)
            adj[b].add(a)
    return [sorted(a) for a in adj], rows


def hex_histogram(t):
    """number of hexagons (degree 6 vertices of T) having exactly j pentagon neighbours"""
    pent = {v for v in range(len(t)) if len(t[v]) == 5}
    h = collections.Counter(sum(1 for w in t[v] if w in pent) for v in range(len(t)) if v not in pent)
    return tuple(h[j] for j in range(7))
