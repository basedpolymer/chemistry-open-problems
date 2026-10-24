#!/usr/bin/env python3
"""exact_small_cases.py -- exact Kirchhoff indices (two independent methods: grounded inverse and
characteristic polynomial, python-flint) for: C5 (all deletions), the seven theta graphs theta(a,b,c) with
path lengths in {1,2,3} (at most one of length 1) -- by Corollary 3 these are the only bicyclic candidates --
and the paper's K_{4,3} (smaller block good).  Writes certificates/small_cases.txt; PASS/FAIL; exit 1 on failure.
Runtime: seconds.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kflib import kf_inv, kf_charpoly, delete, is_connected, g6_encode  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(here, '..', 'certificates', 'small_cases.txt'), 'w')
ok = True


def report(name, adj, expect_ks):
    global ok
    k1, k2 = kf_inv(adj), kf_charpoly(adj)
    if k1 != k2:
        ok = False
    lines = [f'{name}  graph6={g6_encode(adj)}  n={len(adj)}  Kf(G)={k1}']
    allgood = True
    for v in sorted(adj):
        H = delete(adj, v)
        if not is_connected(H):
            lines.append(f'   v={v}: G-v disconnected')
            allgood = False
            continue
        h1, h2 = kf_inv(H), kf_charpoly(H)
        if h1 != h2:
            ok = False
        lines.append(f'   v={v} (deg {len(adj[v])}): Kf(G-v)={h1}  {"GOOD" if h1 == k1 else "Kf(G)-Kf(G-v)=" + str(k1 - h1)}')
        allgood &= (h1 == k1)
    res = (allgood == expect_ks)
    ok &= res
    lines.append(('PASS ' if res else 'FAIL ') + f'{name}: Kirchhoff Soltes = {allgood} (expected {expect_ks})')
    out.write('\n'.join(lines) + '\n')
    print(lines[-1])


def theta(a, b, c):
    adj = {0: set(), 1: set()}
    nxt = 2
    for L in (a, b, c):
        prev = 0
        for _ in range(L - 1):
            adj[nxt] = set()
            adj[prev].add(nxt)
            adj[nxt].add(prev)
            prev = nxt
            nxt += 1
        adj[prev].add(1)
        adj[1].add(prev)
    return adj


report('C5', {i: {(i - 1) % 5, (i + 1) % 5} for i in range(5)}, True)
for t in [(1, 2, 2), (1, 2, 3), (1, 3, 3), (2, 2, 2), (2, 2, 3), (2, 3, 3), (3, 3, 3)]:
    report(f'theta{t}', theta(*t), False)
k43 = {i: (set(range(4, 7)) if i < 4 else set(range(4))) for i in range(7)}
report('K_{4,3}', k43, False)
out.close()
print('PASS exact small cases (written to certificates/small_cases.txt)' if ok else 'FAIL exact small cases')
sys.exit(0 if ok else 1)
