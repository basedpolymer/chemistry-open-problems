#!/usr/bin/env python3
"""modular_check.py -- exact (modular) test Kf(G) =? Kf(G - v) for one graph of any size, using
python-flint nmod_mat (same criterion as kscheck.c: k_p(H) = N tr(A^{-1}) - 1^T A^{-1} 1 mod p, A the
grounded Laplacian; k_p(G) != k_p(G-v) proves Kf(G) != Kf(G-v)).  Used for near-misses of the floating
point scans with more than 64 vertices.
Usage: modular_check.py FILE.g6 [v ...]   (default v = 0).  Prints REJECTED (rigorous) or PASSES-MOD-P.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import nmod_mat  # noqa: E402
from kflib import g6_decode, delete, is_connected, laplacian_rows  # noqa: E402

PRIMES = [1073741789, 1073741783, 1073741741]


def kf_mod(adj, p):
    n = len(adj)
    L = laplacian_rows(adj)
    A = nmod_mat([row[:-1] for row in L[:-1]], p)
    if int(A.det()) == 0:
        return None
    Ai = A.inv()
    N = n - 1
    tr = sum(int(Ai[i, i]) for i in range(N)) % p
    tot = sum(int(Ai[i, j]) for i in range(N) for j in range(N)) % p
    return (n * tr - tot) % p


def test(adj, v):
    H = delete(adj, v)
    if not is_connected(H):
        return 'REJECTED (cut vertex)'
    for p in PRIMES:
        a, b = kf_mod(adj, p), kf_mod(H, p)
        if a is None or b is None:
            continue
        if a != b:
            return f'REJECTED (Kf(G) != Kf(G-v) mod {p})'
    return 'PASSES-MOD-P (needs exact rational check)'


if __name__ == '__main__':
    for line in open(sys.argv[1]):
        if not line.strip():
            continue
        adj = g6_decode(line)
        vs = [int(x) for x in sys.argv[2:]] or [0]
        for v in vs:
            print(line.split()[0][:30], 'n =', len(adj), 'v =', v, test(adj, v))
