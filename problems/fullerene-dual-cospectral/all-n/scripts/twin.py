"""twin.py -- even triangulations H (degrees 4 and 6) and their two colour-class fullerene duals T_b, T_w.

H: plane triangulation, all degrees even, exactly six vertices of degree 4, the others of degree 6.
Faces of H are properly 2-coloured.  T_c (c = black, white) has the faces of colour c as vertices, two of them
adjacent iff they share a vertex of H.
"""
import collections
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdc import *  # noqa: F401,F403


def faces_of(rot):
    """faces of a triangulation as sorted vertex triples."""
    F = set()
    for v, r in enumerate(rot):
        k = len(r)
        for i in range(k):
            F.add(tuple(sorted((v, r[i], r[(i + 1) % k]))))
    return sorted(F)


def two_colour_faces(F):
    """proper 2-colouring of the faces of an even triangulation (faces adjacent across an edge differ)."""
    ef = collections.defaultdict(list)
    for i, f in enumerate(F):
        for e in itertools.combinations(f, 2):
            ef[e].append(i)
    col = [None] * len(F)
    col[0] = 0
    stack = [0]
    while stack:
        i = stack.pop()
        for e in itertools.combinations(F[i], 2):
            j = [x for x in ef[e] if x != i][0]
            if col[j] is None:
                col[j] = 1 - col[i]
                stack.append(j)
            elif col[j] == col[i]:
                raise ValueError('not 2-colourable')
    return col


def colour_class_graph(F, col, c):
    """adjacency lists (sets) of T_c; vertices numbered in the order of the faces of colour c."""
    idx = [i for i in range(len(F)) if col[i] == c]
    pos = {i: k for k, i in enumerate(idx)}
    at = collections.defaultdict(list)
    for i in idx:
        for v in F[i]:
            at[v].append(pos[i])
    adj = [set() for _ in idx]
    for v, lst in at.items():
        for a, b in itertools.combinations(lst, 2):
            adj[a].add(b)
            adj[b].add(a)
    return [sorted(a) for a in adj]


def twins(rot):
    """the two colour-class graphs of the even triangulation with rotation system rot."""
    F = faces_of(rot)
    col = two_colour_faces(F)
    return colour_class_graph(F, col, 0), colour_class_graph(F, col, 1)


def perfect_matchings(pent, adj):
    pent = sorted(pent)
    if not pent:
        yield []
        return
    u = pent[0]
    for v in adj[u]:
        if v in pent and v != u:
            rest = [w for w in pent if w not in (u, v)]
            for M in perfect_matchings(rest, adj):
                yield [(u, v)] + M


def lines_of(rot, M):
    """for a perfect matching M of the pentagons by edges of T: the partition of E(T) into triangular faces and the
    edges of M, or None if there is none."""
    F = faces_of(rot)
    Mset = {frozenset(e) for e in M}
    ef = collections.defaultdict(list)
    for i, f in enumerate(F):
        for e in itertools.combinations(f, 2):
            ef[frozenset(e)].append(i)
    col = [None] * len(F)
    col[0] = 0
    stack = [0]
    while stack:
        i = stack.pop()
        for e in itertools.combinations(F[i], 2):
            fe = frozenset(e)
            j = [x for x in ef[fe] if x != i][0]
            want = col[i] if fe in Mset else 1 - col[i]
            if col[j] is None:
                col[j] = want
                stack.append(j)
            elif col[j] != want:
                return None
    qc = {col[ef[e][0]] for e in Mset}
    if len(qc) != 1:
        return None
    down = qc.pop()
    return [list(f) for f, c in zip(F, col) if c != down] + [list(e) for e in M]


def clique_partitions(rot):
    m = len(rot)
    adj = [set(r) for r in rot]
    pent = [v for v in range(m) if len(rot[v]) == 5]
    res = []
    for M in perfect_matchings(pent, adj):
        L = lines_of(rot, M)
        if L is not None:
            res.append((M, L))
    return res


def line_graph(lines):
    """the even triangulation H as a graph: lines of the partition, adjacent iff they meet."""
    L = len(lines)
    adj = [set() for _ in range(L)]
    at = collections.defaultdict(list)
    for i, l in enumerate(lines):
        for v in l:
            at[v].append(i)
    for v, lst in at.items():
        for a, b in itertools.combinations(lst, 2):
            adj[a].add(b)
            adj[b].add(a)
    return [sorted(a) for a in adj]
