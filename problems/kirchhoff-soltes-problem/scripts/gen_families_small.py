#!/usr/bin/env python3
"""gen_families_small.py -- graph6 generator for the exact (kscheck) part of the family scans:
  gp     : generalized Petersen graphs GP(n,k), 3 <= n <= 32, 1 <= k < n/2       (2n <= 64 vertices)
  tori4  : Cay(Z^2/Lambda, {+-e1,+-e2}) for all sublattices Lambda of index 5..64 (simple 4-regular)
  tri6   : Cay(Z^2/Lambda, {+-e1,+-e2,+-(e1+e2)}) for all Lambda of index 7..64  (simple 6-regular)
Usage: gen_families_small.py FAMILY | kscheck      (kscheck tests every vertex)
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kflib import g6_encode  # noqa: E402

fam = sys.argv[1]
cnt = 0
if fam == 'gp':
    for n in range(3, 33):
        for k in range(1, (n - 1) // 2 + 1):
            adj = {i: set() for i in range(2 * n)}
            for i in range(n):
                for (u, w) in [(i, (i + 1) % n), (i, n + i), (n + i, n + (i + k) % n)]:
                    adj[u].add(w)
                    adj[w].add(u)
            print(g6_encode(adj))
            cnt += 1
else:
    gens = [(1, 0), (-1, 0), (0, 1), (0, -1)] + ([(1, 1), (-1, -1)] if fam == 'tri6' else [])
    for n in range(5, 65):
        for a in range(1, n + 1):
            if n % a:
                continue
            d = n // a
            for b in range(d):
                def red(x, y):
                    p = x // a
                    return (x - p * a, (y - p * b) % d)
                elems = sorted({red(x, y) for x in range(a) for y in range(d)})
                assert len(elems) == n
                idx = {e: i for i, e in enumerate(elems)}
                adj = {i: set() for i in range(n)}
                simple = True
                for e in elems:
                    imgs = [red(e[0] + s[0], e[1] + s[1]) for s in gens]
                    if len(set(imgs)) < len(gens) or e in imgs:
                        simple = False
                        break
                    for t in imgs:
                        adj[idx[e]].add(idx[t])
                if simple:
                    print(g6_encode(adj))
                    cnt += 1
print(f'{fam}: {cnt} graphs generated', file=sys.stderr)
