#!/usr/bin/env python3
"""test_kscheck_modular.py -- validates the C filter kscheck against exact rational arithmetic.

For a deterministic sample of random connected graphs (n = 3..40) and a list of named graphs,
compares the residues printed by `kscheck -k` (Kf(G) mod p1 and Kf(G-v) mod p1, p1 = 1073741789)
with exact rational values (python-flint) reduced mod p1, and checks that the two exact methods
(grounded inverse, characteristic polynomial) agree.  Prints PASS/FAIL, exits non-zero on failure.
Runtime: about 1-2 minutes.  Environment: python 3.13, python-flint 0.9, gcc.
Usage: python3 test_kscheck_modular.py [path/to/kscheck]
"""
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kflib import g6_encode, kf_inv, kf_charpoly, delete, is_connected  # noqa: E402

P1 = 1073741789
here = os.path.dirname(os.path.abspath(__file__))
exe = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, 'kscheck')


def modp(q):
    num, den = int(q.p), int(q.q)
    if den % P1 == 0:
        return None
    return num * pow(den, -1, P1) % P1


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


rng = random.Random(20261024)
graphs = []
for n in range(3, 41):
    for _ in range(4 if n <= 20 else 2):
        graphs.append(rand_connected(n, rng.choice([0.15, 0.3, 0.5, 0.8]) if n > 8 else 0.6, rng))
# named graphs: cycles, complete graphs, Petersen
for n in (3, 5, 8, 11):
    graphs.append({i: {(i - 1) % n, (i + 1) % n} for i in range(n)})
for n in (4, 7):
    graphs.append({i: set(range(n)) - {i} for i in range(n)})
pet = {i: set() for i in range(10)}
for i in range(5):
    for a, b in ((i, (i + 1) % 5), (i, i + 5), (i + 5, (i + 2) % 5 + 5)):
        pet[a].add(b)
        pet[b].add(a)
graphs.append(pet)

inp = ''.join(g6_encode(g) + '\n' for g in graphs)
out = subprocess.run([exe, '-k'], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
out = [l for l in out if l.strip()]
ok = True
nchk = 0
for g, line in zip(graphs, out):
    vals = [int(x) for x in line.split()]
    n = vals[0]
    assert n == len(g)
    kG = kf_inv(g)
    if len(g) <= 22:
        if kf_charpoly(g) != kG:
            print('FAIL: kf_inv and kf_charpoly disagree', g6_encode(g))
            ok = False
    if vals[1] != -1 and vals[1] != modp(kG):
        print('FAIL: Kf(G) mod p mismatch', g6_encode(g))
        ok = False
    nchk += 1
    for v in range(n):
        H = delete(g, v)
        if not is_connected(H):
            if vals[2 + v] != -1:
                print('FAIL: cut vertex not flagged', g6_encode(g), v)
                ok = False
            continue
        kH = kf_inv(H)
        if vals[2 + v] != -1 and vals[2 + v] != modp(kH):
            print('FAIL: Kf(G-v) mod p mismatch', g6_encode(g), v)
            ok = False
        nchk += 1
print(f'{len(graphs)} graphs, {nchk} Kirchhoff indices compared mod p1')
# known values
c5 = {i: {(i - 1) % 5, (i + 1) % 5} for i in range(5)}
checks = [
    ('Kf(C5) = 10', kf_inv(c5) == 10),
    ('Kf(P4) = 10', kf_inv(delete(c5, 0)) == 10),
    ('Kf(Petersen) = 135/... cross-check charpoly', kf_inv(pet) == kf_charpoly(pet)),
    ('Kf(K7) = 6', kf_inv({i: set(range(7)) - {i} for i in range(7)}) == 6),
]
for name, c in checks:
    print(('PASS ' if c else 'FAIL ') + name)
    ok &= c
print('PASS kscheck residues agree with exact arithmetic' if ok else 'FAIL kscheck residues')
sys.exit(0 if ok else 1)
