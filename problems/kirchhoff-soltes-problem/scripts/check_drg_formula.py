#!/usr/bin/env python3
"""check_drg_formula.py -- validates drg_formula.drg_kf_pair against direct exact computation (flint) on
explicit distance-regular graphs: cycles C_4..C_9, K_{4,4}, K_6 (a_1 > 0), Petersen, Q_3, Q_4, Heawood,
Clebsch (16,5,0,2), odd graph O_4 (35 vertices), Hamming H(2,3), H(3,3) (a_1 > 0, local graph d K_{q-1}),
Paley(13) (local graph spectrum computed directly).  PASS/FAIL; exit 1 on failure.  Runtime ~1 min.
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr  # noqa: E402
from kflib import kf_inv, delete  # noqa: E402
from drg_formula import drg_kf_pair  # noqa: E402
import numpy as np  # noqa: E402


def tofr(q):
    return Fr(int(q.p), int(q.q))


def from_edges(n, E):
    adj = {i: set() for i in range(n)}
    for u, v in E:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def cycle(n):
    return from_edges(n, [(i, (i + 1) % n) for i in range(n)])


def petersen():
    return from_edges(10, [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)]
                      + [(5 + i, 5 + (i + 2) % 5) for i in range(5)])


def cube(d):
    return from_edges(2 ** d, [(u, u ^ (1 << i)) for u in range(2 ** d) for i in range(d) if u < u ^ (1 << i)])


def heawood():
    lines = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2)]
    return from_edges(14, [(p, 7 + j) for j, L in enumerate(lines) for p in L])


def clebsch():
    vs = list(range(16))
    E = [(u, w) for u in vs for w in vs if u < w and bin(u ^ w).count('1') in (1, 4)]
    return from_edges(16, E)   # folded 5-cube = Clebsch graph


def odd(kk):
    subs = list(itertools.combinations(range(2 * kk - 1), kk - 1))
    idx = {s: i for i, s in enumerate(subs)}
    E = [(idx[a], idx[b]) for a in subs for b in subs if idx[a] < idx[b] and not set(a) & set(b)]
    return from_edges(len(subs), E)


def hamming(d, q):
    vs = list(itertools.product(range(q), repeat=d))
    idx = {s: i for i, s in enumerate(vs)}
    E = [(idx[a], idx[b]) for a in vs for b in vs if idx[a] < idx[b] and sum(x != y for x, y in zip(a, b)) == 1]
    return from_edges(len(vs), E)


def paley(q):
    sq = {(x * x) % q for x in range(1, q)}
    return from_edges(q, [(a, b) for a in range(q) for b in range(a + 1, q) if (b - a) % q in sq])


def complete(n):
    return from_edges(n, [(i, j) for i in range(n) for j in range(i + 1, n)])


def kmm(m):
    return from_edges(2 * m, [(i, m + j) for i in range(m) for j in range(m)])


def local_spectrum(adj, v):
    N = sorted(adj[v])
    A = np.array([[1.0 if b in adj[a] else 0.0 for b in N] for a in N])
    ev = sorted(np.linalg.eigvalsh(A))
    return ev


cases = []
for nn in range(4, 10):
    D = nn // 2
    b = [2] + [1] * (D - 1)
    c = [1] * (D - 1) + ([2] if nn % 2 == 0 else [1])
    cases.append((f'C{nn}', cycle(nn), b, c, None))
cases += [
    ('K44', kmm(4), [4, 3], [1, 4], None),
    ('K6', complete(6), [5], [1], [-1] * 4),
    ('Petersen', petersen(), [3, 2], [1, 1], None),
    ('Q3', cube(3), [3, 2, 1], [1, 2, 3], None),
    ('Q4', cube(4), [4, 3, 2, 1], [1, 2, 3, 4], None),
    ('Heawood', heawood(), [3, 2, 2], [1, 1, 3], None),
    ('Clebsch', clebsch(), [5, 4], [1, 2], None),
    ('O4', odd(4), [4, 3, 3], [1, 1, 2], None),
    ('H(2,3)', hamming(2, 3), [4, 2], [1, 2], [1] + [-1] * 2),       # local 2K_2: eigs 1,1,-1,-1; minus a_1=1
    ('H(3,3)', hamming(3, 3), [6, 4, 2], [1, 2, 3], [1, 1] + [-1] * 3),
    ('Paley13', paley(13), [6, 3], [1, 3], 'direct'),
]
ok = True
for name, adj, b, c, theta in cases:
    if theta == 'direct':
        ev = local_spectrum(adj, 0)
        a1 = b[0] - b[1] - 1
        # remove one copy of a1, round to exact algebraic? Paley(13) local graph is C6 -> exact cos values
        # use exact values for C6: 2cos(2 pi j/6): 2,1,1,-1,-1,-2  (a1 = 3? no: Paley(13) has lambda=2)
        theta = None
    if name == 'Paley13':
        # Paley(13): srg(13,6,2,3); local graph of 0 is the 2-regular graph on the squares {1,3,4,9,10,12}
        N = sorted(adj[0])
        A = [[1 if y in adj[x] else 0 for y in N] for x in N]
        # 2-regular on 6 vertices: check it is C6 or 2C3 and use exact spectrum
        ev = sorted(np.linalg.eigvalsh(np.array(A, float)))
        if np.allclose(ev, [-2, -1, -1, 1, 1, 2]):
            theta = [Fr(-2), Fr(-1), Fr(-1), Fr(1), Fr(1)]     # C6, remove eigenvalue 2 = a_1
        elif np.allclose(ev, [-1, -1, -1, -1, 2, 2]):
            theta = [Fr(-1)] * 4 + [Fr(2)]                    # 2C3, remove one eigenvalue 2
        else:
            print('unexpected local graph', ev)
    n, k, kfG, kfH = drg_kf_pair(b, c, theta)
    exG = tofr(kf_inv(adj))
    exH = tofr(kf_inv(delete(adj, 0)))
    good = (n == len(adj) and kfG == exG and kfH == exH)
    ok &= good
    print(('PASS ' if good else 'FAIL ') + f'{name}: n={n} Kf(G)={kfG} Kf(G-v)={kfH} direct=({exG}, {exH})')
print('PASS DRG closed form' if ok else 'FAIL DRG closed form')
sys.exit(0 if ok else 1)
