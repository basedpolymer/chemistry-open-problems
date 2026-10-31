#!/usr/bin/env python3
"""Exact maximum of W over trees on n = 2d+1 vertices with diameter at most d, for 3 <= d <= 400.

A tree of diameter <= d has a centre: a vertex with all depths <= d/2 (d even) or an edge with both
sides of height <= (d-1)/2 (d odd).  With W(T) = sum_e a_e (n - a_e), a knapsack over the children of
the root gives the exact optimum.  Integer arithmetic only.  Prints PASS if the maximum is <= W(C_{2d+1})
for every d, and FAIL otherwise, with the expected equalities at d = 3, 4.
Environment: Python 3 with numpy and numba. Runtime: under one minute. Deterministic.
Log: certificates/tree-dp.txt.
"""
import sys
import numba as nb
import numpy as np


@nb.njit
def solve(d):
    n = 2 * d + 1
    R = d // 2
    M = np.full(n + 2, -1, dtype=np.int64)
    M[1] = 0
    for _ in range(1, R + 1):
        val = np.full(n + 2, -1, dtype=np.int64)
        for a in range(1, n + 1):
            if M[a] >= 0:
                val[a] = a * (n - a) + M[a]
        P = np.full(n + 2, -1, dtype=np.int64)
        P[0] = 0
        for s in range(1, n + 1):
            best = -1
            for a in range(1, s + 1):
                if val[a] >= 0 and P[s - a] >= 0:
                    c = val[a] + P[s - a]
                    if c > best:
                        best = c
            P[s] = best
        Mn = np.full(n + 2, -1, dtype=np.int64)
        for m in range(1, n + 1):
            Mn[m] = P[m - 1]
        M = Mn
    if d % 2 == 0:
        return M[n]
    best = -1
    for a in range(1, n):
        if M[a] >= 0 and M[n - a] >= 0:
            c = a * (n - a) + M[a] + M[n - a]
            if c > best:
                best = c
    return best


ok = True
print("d  max W(tree)  W(C_{2d+1})  gap  ratio")
for d in list(range(3, 41)) + [50, 60, 80, 100, 150, 200, 300, 400]:
    v = int(solve(d))
    Wc = (2 * d + 1) * d * (d + 1) // 2
    print(d, v, Wc, Wc - v, round(v / Wc, 5))
    ok &= (v == Wc) if d in (3, 4) else (v < Wc)
print("PASS" if ok else "FAIL", "tree maximum below W(C_{2d+1}) for 5<=d, equal for d=3,4")
sys.exit(0 if ok else 1)
