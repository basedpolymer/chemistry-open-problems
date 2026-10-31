"""Exact linear algebra and network helpers for the verifier's reimplementation.

Written from the mathematical specification, without reading the solver's scripts.
Python 3.12, standard library only. Deterministic.
"""
from fractions import Fraction as Fr
from itertools import product
from math import gcd
from functools import reduce


def rref(rows, ncols):
    """Reduced row echelon form over Q. Returns (nonzero rows as Fractions, pivot columns)."""
    M = [[Fr(x) for x in r] for r in rows]
    piv, r = [], 0
    for c in range(ncols):
        s = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if s is None:
            continue
        M[r], M[s] = M[s], M[r]
        h = M[r][c]
        M[r] = [x / h for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    return M[:r], piv


def primitive(v):
    """Scale a rational vector to a primitive integer one, first nonzero entry positive."""
    den = reduce(lambda a, b: a * b // gcd(a, b), [Fr(x).denominator for x in v], 1)
    w = [int(Fr(x) * den) for x in v]
    g = reduce(gcd, [abs(x) for x in w], 0)
    if g:
        w = [x // g for x in w]
    f = next((x for x in w if x != 0), 1)
    return [(-x if f < 0 else x) for x in w]


def nullspace(rows, ncols):
    """Integer basis of {v : rows * v = 0}."""
    R, piv = rref(rows, ncols)
    free = [c for c in range(ncols) if c not in piv]
    out = []
    for f in free:
        v = [Fr(0)] * ncols
        v[f] = Fr(1)
        for i, p in enumerate(piv):
            v[p] = -R[i][f]
        out.append(primitive(v))
    return out


def det(M):
    """Exact determinant by fraction-free (Bareiss) elimination."""
    n = len(M)
    if n == 0:
        return 1
    A = [list(map(int, r)) for r in M]
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            s = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if s is None:
                return 0
            A[k], A[s] = A[s], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def complexes(n):
    """All complexes of molecularity at most two on n species."""
    return [y for y in product(range(3), repeat=n) if sum(y) <= 2]


def diff(u, v):
    return tuple(a - b for a, b in zip(u, v))


def named(names, text):
    """Parse a complex like 'A+B', '2E', '0' into an exponent tuple."""
    y = [0] * len(names)
    if text == '0':
        return tuple(y)
    for tok in text.split('+'):
        k, s = (int(tok[0]), tok[1:]) if tok[0].isdigit() else (1, tok)
        y[names.index(s)] += k
    return tuple(y)


def add_n4_types():
    """Register the 26 four-species types found by enum_sources.py as T0, ..., T25."""
    import json
    import os
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sources-verifier.json')
    data = json.load(open(path))
    for i, rec in enumerate(data['n4']):
        PAPER_TYPES[f'T{i}'] = ('ABCD', rec['sources'], ''.join('ABCD'[j] for j in rec['acr']))


# Source table of the manuscript, transcribed by hand, used only for comparison.
PAPER_TYPES = {
    'S1': ('ABCDE', ['0', 'E', 'C+D', 'A+B'], 'E'),
    'S2': ('ABCDE', ['0', '2E', 'C+D', 'A+B'], 'E'),
    'S3': ('ABCDE', ['E', '2E', 'C+D', 'A+B'], 'E'),
    'S4': ('ABCDE', ['E', 'D', 'C+E', 'A+B'], 'C'),
    'S5': ('ABCDE', ['E', 'D+E', 'C+E', 'A+B'], 'CD'),
    'S6': ('ABCDE', ['E', 'D+E', 'C+D', 'A+B'], 'D'),
    'S7': ('ABCDE', ['E', 'D+E', '2C', 'A+B'], 'D'),
    'S8': ('ABCDE', ['E', 'D+E', 'B+C', 'A+C'], 'D'),
    'S9': ('ABCDEF', ['F', 'E+F', 'C+D', 'A+B'], 'E'),
}
