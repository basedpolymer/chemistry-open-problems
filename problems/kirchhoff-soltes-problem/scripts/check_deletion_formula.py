#!/usr/bin/env python3
"""check_deletion_formula.py -- exact check of the local vertex-deletion formula (deletion_formula.py)
against direct exact computation of Kf(G - v), on deterministic random connected graphs (n = 4..16).
Exact arithmetic: python-flint fmpq_mat.  PASS/FAIL, exit 1 on failure.  Runtime < 1 min.
"""
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import fmpq, fmpq_mat  # noqa: E402
from kflib import kf_inv, delete, is_connected, laplacian_rows  # noqa: E402
from deletion_formula import deletion_quantities  # noqa: E402


def pinv_exact(adj):
    n = len(adj)
    L = fmpq_mat(laplacian_rows(adj))
    J = fmpq_mat([[fmpq(1, n)] * n for _ in range(n)])
    G = (L + J).inv() - J       # L^+ = (L + J/n)^{-1} - J/n
    return G


def tr_inv_mul(Z, W):
    Zm = fmpq_mat(Z)
    Wm = fmpq_mat(W)
    P = Zm.inv() * Wm
    return sum((P[i, i] for i in range(P.nrows())), fmpq(0))


rng = random.Random(11)
ok = True
cnt = 0
for trial in range(40):
    n = rng.randint(4, 16)
    while True:
        adj = {i: set() for i in range(n)}
        p = rng.choice([0.25, 0.4, 0.7])
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < p:
                    adj[i].add(j)
                    adj[j].add(i)
        if is_connected(adj):
            break
    G = pinv_exact(adj)
    G2 = G * G
    trG = sum((G[i, i] for i in range(n)), fmpq(0))
    for v in range(n):
        H = delete(adj, v)
        if not is_connected(H):
            continue
        N = sorted(adj[v])
        kfG, kfH = deletion_quantities(n, v, N, lambda x, y: G[x, y], lambda x, y: G2[x, y], trG,
                                       fmpq(1), fmpq(0), tr_inv_mul)
        cnt += 1
        if kfG != kf_inv(adj) or kfH != kf_inv(H):
            ok = False
            print('FAIL', trial, v)
print(f'{cnt} vertex deletions checked exactly')
print('PASS local deletion formula' if ok else 'FAIL local deletion formula')
sys.exit(0 if ok else 1)
