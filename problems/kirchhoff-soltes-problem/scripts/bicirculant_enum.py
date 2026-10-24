#!/usr/bin/env python3
"""bicirculant_enum.py -- graph6 stream of all bicirculants on 2m vertices (graphs with a semiregular
automorphism of order m with two orbits): vertices u_x, w_x (x in Z_m), edges u_x u_{x+s} (s in S0),
w_x w_{x+s} (s in S1), u_x w_{x+t} (t in T), S0, S1 symmetric subsets of Z_m \\ {0}, T a non-empty subset of
Z_m.  Translating the w-orbit makes 0 in T without loss of generality (T empty gives a disconnected
graph).  No isomorph rejection (duplicates are harmless).  Every bicirculant on 2m vertices is
isomorphic to one of these.  Usage: bicirculant_enum.py m  |  kscheck_fast   (counts on stderr)
"""
import itertools
import sys

m = int(sys.argv[1])
n = 2 * m
half = list(range(1, m // 2 + 1))


def symsets():
    for r in range(len(half) + 1):
        for comb in itertools.combinations(half, r):
            s = set()
            for a in comb:
                s.add(a % m)
                s.add((-a) % m)
            yield sorted(s)


S = list(symsets())
nbits = n * (n - 1) // 2
pos = {}
for j in range(1, n):
    for i in range(j):
        pos[(i, j)] = j * (j - 1) // 2 + i


def bitindex(a, b):
    if a > b:
        a, b = b, a
    return pos[(a, b)]


# precompute edge bit lists for each component
bits0 = {tuple(s): [bitindex(x, (x + d) % m) for x in range(m) for d in s if x < (x + d) % m] for s in S}
bits1 = {tuple(s): [bitindex(m + x, m + (x + d) % m) for x in range(m) for d in s if x < (x + d) % m] for s in S}
bitsT = {t: [bitindex(x, m + (x + t) % m) for x in range(m)] for t in range(m)}
nchar = (nbits + 5) // 6
head = chr(n + 63)
out = sys.stdout
cnt = 0
for s0 in S:
    b0 = bits0[tuple(s0)]
    for s1 in S:
        b01 = b0 + bits1[tuple(s1)]
        for r in range(m):
            for rest in itertools.combinations(range(1, m), r):
                T = (0,) + rest
                arr = bytearray(nchar)
                for b in b01:
                    arr[b // 6] |= 1 << (5 - b % 6)
                for t in T:
                    for b in bitsT[t]:
                        arr[b // 6] |= 1 << (5 - b % 6)
                out.write(head + ''.join(chr(c + 63) for c in arr) + '\n')
                cnt += 1
print(f'm={m} n={n} bicirculants generated {cnt}', file=sys.stderr)
