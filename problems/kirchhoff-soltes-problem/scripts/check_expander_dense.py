#!/usr/bin/env python3
"""check_expander_dense.py -- checks for Lemma 4 and Theorems 5, 6 of result.md (expander and dense bounds).

Exact part (python-flint rationals), on deterministic random 2-connected graphs n = 4..11:
  (E1) sum_v n sigma_v = 2 Kf(G), where n sigma_v = sum_{u<w} E_v(u,w) = n tr(Gamma^2 L*_v);
  (E2) for every v and every pair u,w != v:  r_{G-v}(u,w) - r_G(u,w) >= E_v(u,w)   (lower half of Lemma 4);
  (E3) eps_v = sum_w E_v(v,w) >= (n-1)/d_v.
Floating-point sanity part (numpy, tolerance 1e-9; the inequalities involve irrational eigenvalues):
  (F1) r_{G-v}(u,w) - r_G(u,w) <= (1 + 1/a(G-v)) E_v(u,w)        (upper half of Lemma 4);
  (F2) sum_v [Kf(G)-Kf(G-v)] >= (1+1/a_*) sum_v eps_v - 2Kf/a_*    (Theorem 5, sharp form);
  (F3) a(G) > 1  =>  Kf(G)-Kf(G-v) >= 1 - sum_i (n-mu_i)/(mu_i(mu_i-1)) for every v   (Theorem 6).
PASS/FAIL, exit 1 on failure.  Runtime ~1-2 min.
"""
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import fmpq, fmpq_mat  # noqa: E402
from kflib import kf_inv, delete, is_connected, laplacian_rows  # noqa: E402

rng = random.Random(5)
ok = True
tol = 1e-9


def pinv_exact(adj):
    n = len(adj)
    L = fmpq_mat(laplacian_rows(adj))
    J = fmpq_mat([[fmpq(1, n)] * n for _ in range(n)])
    return (L + J).inv() - J


def res_matrix(G, n):
    return [[G[i, i] + G[j, j] - 2 * G[i, j] for j in range(n)] for i in range(n)]


ntest = 0
while ntest < 25:
    n = rng.randint(4, 11)
    p = rng.choice([0.3, 0.5, 0.7, 0.9])
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p:
                adj[i].add(j)
                adj[j].add(i)
    if not is_connected(adj) or any(not is_connected(delete(adj, v)) for v in range(n)):
        continue
    ntest += 1
    Gm = pinv_exact(adj)
    R = res_matrix(Gm, n)
    kG = kf_inv(adj)
    Lf = np.array(laplacian_rows(adj), float)
    evG = np.linalg.eigvalsh(Lf)
    aG = evG[1]
    tot_sigma = fmpq(0)
    sumdelta = 0.0
    sumeps = 0.0
    astar = 1e18
    for v in range(n):
        H = delete(adj, v)
        idx = [u for u in range(n) if u != v]
        GH = pinv_exact({i: {idx.index(w) for w in H[u]} for i, u in enumerate(idx)})
        RH = res_matrix(GH, n - 1)
        LH = np.array(laplacian_rows({i: {idx.index(w) for w in H[u]} for i, u in enumerate(idx)}), float)
        aH = np.linalg.eigvalsh(LH)[1]
        astar = min(astar, aH)

        def E(u, w):   # energy in the star of v for unit current u -> w in G
            phi = [Gm[x, u] - Gm[x, w] for x in range(n)]
            return sum(((phi[v] - phi[x]) ** 2 for x in adj[v]), fmpq(0))
        sig = fmpq(0)
        for u in range(n):
            for w in range(u + 1, n):
                e = E(u, w)
                sig += e
                if u != v and w != v:
                    inc = RH[idx.index(u)][idx.index(w)] - R[u][w]
                    if inc < e:
                        ok = False
                        print('FAIL (E2)', n, v, u, w)
                    if float(inc) > (1 + 1 / aH) * float(e) + tol:
                        ok = False
                        print('FAIL (F1)', n, v, u, w)
        tot_sigma += sig
        eps = sum((E(v, w) for w in range(n) if w != v), fmpq(0))
        if eps < fmpq(n - 1, len(adj[v])):
            ok = False
            print('FAIL (E3)')
        sumeps += float(eps)
        d = kG - kf_inv(H)
        sumdelta += float(d)
        if aG > 1 + 1e-12:
            bound = 1 - sum((n - m) / (m * (m - 1)) for m in evG[1:])
            if float(d) < bound - tol:
                ok = False
                print('FAIL (F3)', n, v)
    if tot_sigma != 2 * kG:
        ok = False
        print('FAIL (E1)', n)
    if sumdelta < (1 + 1 / astar) * sumeps - 2 * float(kG) / astar - tol:
        ok = False
        print('FAIL (F2)', n)
print(f'{ntest} random 2-connected graphs')
print('PASS exact identities (E1)-(E3) and floating-point checks (F1)-(F3)' if ok else 'FAIL')
sys.exit(0 if ok else 1)
