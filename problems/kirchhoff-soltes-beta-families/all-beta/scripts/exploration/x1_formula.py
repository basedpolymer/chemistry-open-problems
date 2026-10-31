#!/usr/bin/env python3
"""Exploration: closed form of Kf for K_N + hubs on a partition; compare with brute force (exact)."""
import sys, random
sys.path.insert(0, '../../../scripts')
sys.path.insert(0, '/mnt/c/Users/demor/chemistry-open-problems/problems/kirchhoff-soltes-beta-families/scripts')
from fractions import Fraction as Fr
import networkx as nx
import kfexact as kx
from flint import fmpq

def build(blocks, hubless=0):
    G = nx.Graph()
    base = []
    for i, b in enumerate(blocks):
        for j in range(b):
            base.append(('u', i, j))
    for j in range(hubless):
        base.append(('w', j))
    for p in range(len(base)):
        for q in range(p+1, len(base)):
            G.add_edge(base[p], base[q])
    for i, b in enumerate(blocks):
        for j in range(b):
            G.add_edge(('h', i), ('u', i, j))
    return G

def Kf_formula(N, a, ap, Y):
    # Y = sum over hubs of 1/b_i, a hubs, ap = blocks + hubless
    return Y*Fr((N+1)*(N+a-1), N) + (N+a-2) + Fr((1+ap)*(1-a), N+1) + Fr(a*(ap-2), N)

if __name__ != "__main__": raise_skip = True
random.seed(1)
bad = 0
for trial in range(40):
    r = random.randint(1, 5)
    blocks = [random.randint(1, 4) for _ in range(r)]
    hl = random.randint(0, 3)
    N = sum(blocks)+hl
    if N < 2: continue
    G = build(blocks, hl)
    a = len(blocks); ap = a+hl
    Y = sum(Fr(1, b) for b in blocks)
    f = Kf_formula(N, a, ap, Y)
    k = kx.kf(G)
    k = Fr(int(k.p), int(k.q))
    if f != k:
        bad += 1
        print('MISMATCH', blocks, hl, f, k)
print('mismatches', bad)
