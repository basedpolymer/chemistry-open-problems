#!/usr/bin/env python3
"""drg_families.py -- exact Delta = Kf(G) - Kf(G - v) for distance-regular families via drg_formula.py
(intersection array + local-graph spectrum).  For a DRG the formula holds at EVERY vertex (no
vertex-transitivity needed): Delta > 0 means no vertex is good.  Prints one line per family and
PASS if Delta != 0 for every member except C5 (Delta = 0); the sign is reported
(Delta > 0: no vertex good because Kf drops; Delta < 0: Kf rises).  Exit 1 otherwise.
Families: cycles C_n (4..200), hypercubes Q_d (2..120), odd graphs O_{k+1} (k = 1..80), incidence graphs
of generalized triangles/quadrangles/hexagons of order (q,q) (q = 2..150, any integer q: the formula
does not need existence), Hamming H(d,q) (2..30 x 2..30), Johnson J(m,k) (m <= 40), the 13 cubic DRGs,
Moore graphs {k,k-1;1,k} (k = 2..300), triangle-free SRG arrays {k,k-1;1,mu} (k <= 120, all mu | k(k-1)
with mu <= k).  Runtime: a few minutes.
"""
import os
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drg_formula import drg_kf_pair  # noqa: E402

ok = True
log = []


def run(name, b, c, theta=None, expect_zero=False):
    global ok
    n, k, kfG, kfH = drg_kf_pair(b, c, theta)
    d = kfG - kfH
    good = (d == 0) if expect_zero else (d != 0)
    if not good:
        ok = False
        print('FAIL', name, b, c, d)
    return n, d


def family(label, items):
    mn = None
    cnt = pos = neg = zero = 0
    for name, b, c, th, z in items:
        n, d = run(name, b, c, th, z)
        cnt += 1
        pos += d > 0
        neg += d < 0
        zero += d == 0
        if not z:
            if mn is None or abs(d) < abs(mn[0]):
                mn = (d, name, n)
    line = (f'{label}: {cnt} graphs, Delta>0: {pos}, Delta<0: {neg}, Delta=0: {zero}'
            + (f'; smallest |Delta| (Delta != 0) = {float(mn[0]):.6g} at {mn[1]} (n={mn[2]})' if mn else ''))
    print(line)
    log.append(line)


def cyc(nn):
    D = nn // 2
    return [2] + [1] * (D - 1), [1] * (D - 1) + ([2] if nn % 2 == 0 else [1])


# cycles: C5 is the only zero
items = []
for nn in range(4, 201):
    b, c = cyc(nn)
    items.append((f'C{nn}', b, c, None, nn == 5))
family('cycles C_4..C_200 (C_5: Delta = 0)', [x for x in items])

family('hypercubes Q_2..Q_120', [(f'Q{d}', [d - i for i in range(d)], [i + 1 for i in range(d)], None, False)
                                 for d in range(2, 121)])


def odd_array(kk):   # O_{kk+1}: degree kk+1, diameter kk
    b, c = [], []
    for i in range(kk):
        j = i // 2
        b.append(kk + 1 - j if i % 2 == 0 else kk - j)
    for i in range(1, kk + 1):
        c.append((i + 1) // 2)
    return b, c


family('odd graphs O_3..O_81', [(f'O{kk+1}', *odd_array(kk), None, False) for kk in range(2, 81)])


def gen_polygon(q, m):   # incidence graph of a generalized m-gon of order (q,q): diameter m, bipartite
    return [q + 1] + [q] * (m - 1), [1] * (m - 1) + [q + 1]


for m, lab in ((3, 'generalized triangles (projective planes)'), (4, 'generalized quadrangles'), (6, 'generalized hexagons')):
    family(f'incidence graphs of {lab} of order (q,q), q=2..150', [(f'GP{m}({q})', *gen_polygon(q, m), None, False)
                                                                   for q in range(2, 151)])

items = []
for d in range(2, 31):
    for q in range(2, 31):
        b = [(d - i) * (q - 1) for i in range(d)]
        c = [i + 1 for i in range(d)]
        theta = None if q == 2 else [Fr(q - 2)] * (d - 1) + [Fr(-1)] * (d * (q - 2))
        items.append((f'H({d},{q})', b, c, theta, False))
family('Hamming graphs H(d,q), 2<=d,q<=30', items)

items = []
for m in range(4, 41):
    for kk in range(2, m // 2 + 1):
        D = kk
        b = [(kk - i) * (m - kk - i) for i in range(D)]
        c = [(i + 1) ** 2 for i in range(D)]
        # local graph K_kk x K_{m-kk} (rook graph): eigenvalues m-2 (=a_1), kk-2 (m-kk-1), m-kk-2 (kk-1), -2 ((kk-1)(m-kk-1))
        theta = [Fr(kk - 2)] * (m - kk - 1) + [Fr(m - kk - 2)] * (kk - 1) + [Fr(-2)] * ((kk - 1) * (m - kk - 1))
        items.append((f'J({m},{kk})', b, c, theta, False))
family('Johnson graphs J(m,k), m<=40', items)

cubic = [
    ('K4', [3], [1], [Fr(-1)] * 2), ('K33', [3, 2], [1, 3], None), ('Petersen', [3, 2], [1, 1], None),
    ('Q3', [3, 2, 1], [1, 2, 3], None), ('Heawood', [3, 2, 2], [1, 1, 3], None), ('Pappus', [3, 2, 2, 1], [1, 1, 2, 3], None),
    ('Coxeter', [3, 2, 2, 1], [1, 1, 1, 2], None), ('Tutte8cage', [3, 2, 2, 2], [1, 1, 1, 3], None),
    ('Dodecahedron', [3, 2, 1, 1, 1], [1, 1, 1, 2, 3], None), ('Desargues', [3, 2, 2, 1, 1], [1, 1, 2, 2, 3], None),
    ('Tutte12cage', [3, 2, 2, 2, 2, 2], [1, 1, 1, 1, 1, 3], None), ('BiggsSmith', [3, 2, 2, 2, 1, 1, 1], [1, 1, 1, 1, 1, 1, 3], None),
    ('Foster', [3, 2, 2, 2, 2, 1, 1, 1], [1, 1, 1, 1, 2, 2, 2, 3], None)]
family('the 13 cubic distance-regular graphs', [(nm, b, c, th, False) for nm, b, c, th in cubic])

family('Moore graphs {k,k-1;1,1}, k=3..300 (formula; exist only for k=3,7,(57))',
       [(f'Moore{kk}', [kk, kk - 1], [1, 1], None, False) for kk in range(3, 301)])
items = []
for kk in range(2, 121):
    for mu in range(1, kk + 1):
        if (kk * (kk - 1)) % mu == 0:
            items.append((f'tfSRG({kk},{mu})', [kk, kk - 1], [1, mu], None, (kk, mu) == (2, 1)))
family('triangle-free SRG arrays {k,k-1;1,mu}, k<=120 (k=2,mu=1 is C5: Delta=0)', items)

print('PASS all DRG families: Delta != 0 except C5' if ok else 'FAIL')
sys.exit(0 if ok else 1)
