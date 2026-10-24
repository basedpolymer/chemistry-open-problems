#!/usr/bin/env python3
"""check_thread_lemma.py -- exact check of the thread formula (Lemma 2 of result.md).

G = R + thread x - v_1 - ... - v_{l-1} - y  (R connected, x != y, internal vertices of degree 2).
Removing v_i (a = i-1, b = l-1-i) leaves R with pendant paths of a vertices at x and b at y and
  Kf(G - v_i) = Kf(R) + n_R [a(a+1)/2 + b(b+1)/2] + a Rt_R(x) + b Rt_R(y)
               + (a^3-a)/6 + (b^3-b)/6 + ab (a+b+2)/2 + ab r_R(x,y),
a quadratic in a (a+b fixed) with leading coefficient n_R - 1 - r_R(x,y) > 0 unless R is an x-y path.
The script verifies the formula exactly (python-flint) on deterministic random instances, checks the
sign of the leading coefficient, and checks Corollary 3 (two-internal-vertex thread: the two internal
vertices give equal values iff Rt_R(x) = Rt_R(y)).  PASS/FAIL; exit 1 on failure.  Runtime < 1 min.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import fmpq  # noqa: E402
from kflib import kf_inv, delete, is_connected, fmpq_mat, laplacian_rows  # noqa: E402


def resistances(adj):
    """exact resistance matrix via grounded inverse"""
    vs = sorted(adj)
    n = len(vs)
    L = laplacian_rows(adj)
    A = fmpq_mat([row[:-1] for row in L[:-1]]).inv()   # grounded at vs[-1]
    M = [[A[i, j] if i < n - 1 and j < n - 1 else fmpq(0) for j in range(n)] for i in range(n)]
    return {(vs[i], vs[j]): M[i][i] + M[j][j] - 2 * M[i][j] for i in range(n) for j in range(n)}


def rand_connected(n, p, rng):
    while True:
        adj = {i: set() for i in range(n)}
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < p:
                    adj[i].add(j)
                    adj[j].add(i)
        if is_connected(adj):
            return adj


rng = random.Random(7)
ok = True
ntests = 0
for trial in range(60):
    nR = rng.randint(3, 9)
    R = rand_connected(nR, rng.choice([0.3, 0.5, 0.8]), rng)
    x, y = rng.sample(range(nR), 2)
    ell = rng.randint(2, 6)          # number of edges of the thread; ell-1 internal vertices
    G = {u: set(R[u]) for u in R}
    path = [x] + [nR + k for k in range(ell - 1)] + [y]
    for k in range(ell - 1):
        G[nR + k] = set()
    for a_, b_ in zip(path, path[1:]):
        G[a_].add(b_)
        G[b_].add(a_)
    r = resistances(R)
    RtX = sum(r[(x, u)] for u in R)
    RtY = sum(r[(y, u)] for u in R)
    kR = kf_inv(R)
    rxy = r[(x, y)]
    vals = []
    for i in range(1, ell):
        a, b = i - 1, ell - 1 - i
        H = delete(G, path[i])
        exact = kf_inv(H)
        formula = (kR + nR * fmpq(a * (a + 1) + b * (b + 1), 2) + a * RtX + b * RtY
                   + fmpq(a ** 3 - a + b ** 3 - b, 6) + fmpq(a * b * (a + b + 2), 2) + a * b * rxy)
        ntests += 1
        if exact != formula:
            ok = False
            print('FAIL formula', trial, i, exact, formula)
        vals.append(exact)
    # leading coefficient: second difference of vals (in a) equals 2*(nR - 1 - rxy)
    c2 = nR - 1 - rxy
    for k in range(len(vals) - 2):
        if vals[k + 2] - 2 * vals[k + 1] + vals[k] != 2 * c2:
            ok = False
            print('FAIL second difference', trial)
    isP = (sum(len(R[u]) for u in R) // 2 == nR - 1) and all(len(R[u]) <= 2 for u in R) and len(R[x]) == 1 and len(R[y]) == 1
    if (c2 > 0) != (not isP):
        ok = False
        print('FAIL sign of leading coefficient', trial)
    if ell == 3:
        if (vals[0] == vals[1]) != (RtX == RtY):
            ok = False
            print('FAIL two-vertex thread criterion', trial)
print(f'{ntests} vertex deletions checked exactly')
print('PASS thread formula, leading coefficient n_R-1-r_R(x,y) > 0, two-vertex criterion' if ok else 'FAIL')
sys.exit(0 if ok else 1)
