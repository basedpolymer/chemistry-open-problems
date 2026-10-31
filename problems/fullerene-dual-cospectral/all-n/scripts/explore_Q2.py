#!/usr/bin/env python3
"""exploration: structure of the orthogonal intertwiner Q (Q A2 = A1 Q) relative to the best partial isomorphism"""
import json, os, sys, collections, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdc import *
CERT = os.path.join(ROOT, 'certificates')
pairs = json.load(open(os.path.join(CERT, 'cospectral_pairs_T.json')))
idx = int(sys.argv[1])
c = pairs[idx]
a1 = from_graph6(c['graph6_dual'][0]); a2 = from_graph6(c['graph6_dual'][1])
r1 = rotation_system(a1); r2 = rotation_system(a2)
mp, core = best_common_patch(r1, r2)
m = len(a1)


def A_of(adj):
    return flint.fmpz_mat(m, m, [1 if w in adj[v] else 0 for v in range(m) for w in range(m)])


def walkm(A):
    v = flint.fmpz_mat(m, 1, [1] * m)
    cols = []
    for k in range(m):
        cols.append([int(v[i, 0]) for i in range(m)])
        v = A * v
    return flint.fmpq_mat(flint.fmpz_mat(m, m, [cols[j][i] for i in range(m) for j in range(m)]))


A1 = A_of(a1); A2 = A_of(a2)
W1 = walkm(A1); W2 = walkm(A2)
Q = W1 * W2.inv()
print('check Q A2 = A1 Q', Q * flint.fmpq_mat(A2) == flint.fmpq_mat(A1) * Q)
print('check Q Q^T = I', (Q * Q.transpose()) == flint.fmpq_mat(m, m, [1 if i == j else 0 for i in range(m) for j in range(m)]))
# permutation aligning core
perm2 = {}
for v in core:
    perm2[v] = mp[v]
rest1 = [v for v in range(m) if v not in core]
rest2 = [v for v in range(m) if v not in {mp[u] for u in core}]
for a, b in zip(rest1, rest2):
    perm2[a] = b
P = flint.fmpq_mat(m, m, [1 if perm2[j] == i else 0 for i in range(m) for j in range(m)])  # P maps G1 labels to G2
Qp = Q * P      # now Q' relates A1 to P^T A2 P with labels of G1
Iden = flint.fmpq_mat(m, m, [1 if i == j else 0 for i in range(m) for j in range(m)])
print('rank(Q P - I)', (Qp - Iden).rank(), 'm', m, '|core|', len(core))
print('diag entries on core:', [str(Qp[v, v]) for v in sorted(core)][:6])
den = 1
for i in range(m):
    for j in range(m):
        den = den * int(Qp[i, j].q) // math.gcd(den, int(Qp[i, j].q))
print('level', den, flint.fmpz(den).factor())
print('Q symmetric after P?', Qp == Qp.transpose(), ' Qp^2=I?', Qp * Qp == Iden)
zeros = sum(1 for i in range(m) for j in range(m) if Q[i, j] == 0)
print('zeros in Q', zeros, 'of', m * m)
sc = flint.fmpq(den)
rows = [[int(Qp[i, j] * sc) for j in range(m)] for i in range(m)]
order = sorted(core) + rest1
print('order', order, 'degrees', [len(a1[v]) for v in order])
vals = collections.Counter(abs(rows[i][j]) for i in range(m) for j in range(m))
print('distinct abs values', len(vals), vals.most_common(10))
# sum of squares structure: image of Qp - I
E = Qp - Iden
rr, nullity = E.rref()[0], None
print('rank', E.rank())
import numpy as np
Qf = np.array([[float(Qp[i, j].p) / float(Qp[i, j].q) for j in order_] for i in order_ for order_ in [order]]) if False else None
Qf = np.array([[float(int(Qp[i, j].p)) / float(int(Qp[i, j].q)) for j in order] for i in order])
np.set_printoptions(linewidth=250, precision=2, suppress=True)
print(Qf)
w = np.linalg.eigvals(Qf)
print('eig angles', np.sort(np.round(np.angle(w), 3)))
# sparse basis of the fixed space U = ker(Qp - I)
Em = flint.fmpq_mat(E)
den2 = 1
for i in range(m):
    for j in range(m):
        den2 = den2 * int(E[i, j].q) // math.gcd(den2, int(E[i, j].q))
Ez = flint.fmpz_mat(m, m, [int(E[i, j] * flint.fmpq(den2)) for i in range(m) for j in range(m)])
N, nullity = Ez.nullspace()
print('nullity', nullity)
basis = [[int(N[i, j]) for i in range(m)] for j in range(nullity)]
Bm = flint.fmpz_mat(nullity, m, [x for r in basis for x in r])
# clear content and saturate: use LLL on row lattice
L = Bm.lll()
for i in range(nullity):
    r = [int(L[i, j]) for j in range(m)]
    from functools import reduce
    g = reduce(math.gcd, r)
    r = [x // g for x in r]
    print({('%d(%d)' % (j, len(a1[j]))): r[j] for j in range(m) if r[j]})
print('--- U perp')
# U perp = nullspace of Bm (row lattice)
Up, nul2 = Bm.nullspace()
print('nullity', nul2)
b2 = [[int(Up[i, j]) for i in range(m)] for j in range(nul2)]
L2 = flint.fmpz_mat(nul2, m, [x for r in b2 for x in r]).lll()
for i in range(nul2):
    r = [int(L2[i, j]) for j in range(m)]
    g = reduce(math.gcd, r)
    r = [x // g for x in r]
    print({('%d(%d)' % (j, len(a1[j]))): r[j] for j in range(m) if r[j]})
print('edges among pentagons', [(u, v) for u in range(m) for v in a1[u] if u < v and len(a1[u]) == 5 and len(a1[v]) == 5])
print('adj', {v: a1[v] for v in range(m)})
print('core', sorted(core), 'non core', rest1, 'mp', mp)
