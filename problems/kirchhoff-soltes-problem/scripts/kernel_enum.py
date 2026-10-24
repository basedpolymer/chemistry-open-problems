#!/usr/bin/env python3
"""kernel_enum.py -- generate every 2-connected simple graph with cyclomatic number c that is not a
cycle and in which every thread (maximal path of degree-2 vertices) has at most 2 internal vertices.

By Corollary 3 of result.md every Kirchhoff Soltes graph other than a cycle is of this kind.
Such a graph G is a subdivision of its kernel K (suppress degree-2 vertices): K is a loopless
multigraph without cut vertex, minimum degree >= 3, with m_K - n_K = c - 1, hence n_K <= 2(c-1),
m_K <= 3(c-1); each kernel edge is subdivided t in {0,1,2} times, at most one edge of each
parallel class with t = 0 (simplicity).  Kernels come from nauty: geng -C (underlying simple graph,
2-connected; for n_K = 2 the single edge K2) and multig -e<m_K> (multiplicities); kernels with a
vertex of degree < 3 are discarded.  Output: graph6 on stdout (feed to kscheck), counts on stderr.
Usage: kernel_enum.py c [NAUTY_DIR]
Each isomorphism class of G may appear several times (no isomorph rejection); this is harmless.
"""
import itertools
import subprocess
import sys
from collections import Counter

c = int(sys.argv[1])
NAUTY = sys.argv[2] if len(sys.argv) > 2 else 'certificates/sources/nauty/nauty2_8_9'


def kernels(c):
    out = []
    for nk in range(2, 2 * (c - 1) + 1):
        mk = nk + c - 1
        if 2 * mk < 3 * nk:
            continue
        if nk == 2:
            simple = ['A_']     # K2 in graph6
        else:
            if mk < nk:
                continue
            emax = min(mk, nk * (nk - 1) // 2)
            res = subprocess.run([f'{NAUTY}/geng', '-Cq', str(nk), f'{nk}:{emax}'],
                                 capture_output=True, text=True, check=True)
            simple = res.stdout.split()
        if not simple:
            continue
        res = subprocess.run([f'{NAUTY}/multig', '-q', '-T', f'-e{mk}'], input='\n'.join(simple) + '\n',
                             capture_output=True, text=True, check=True)
        for line in res.stdout.strip().split('\n'):
            if not line.strip():
                continue
            a = list(map(int, line.split()))
            nv, ne = a[0], a[1]
            edges = [(a[2 + 3 * k], a[3 + 3 * k], a[4 + 3 * k]) for k in range(ne)]
            deg = Counter()
            for u, v, mlt in edges:
                deg[u] += mlt
                deg[v] += mlt
            if nv != nk or min(deg[u] for u in range(nv)) < 3:
                continue
            assert sum(m for _, _, m in edges) == mk
            out.append((nk, edges))
    return out


def g6(n, edges):
    bits = [0] * (n * (n - 1) // 2)
    for u, v in edges:
        if u > v:
            u, v = v, u
        bits[v * (v - 1) // 2 + u] = 1
    while len(bits) % 6:
        bits.append(0)
    s = [chr(n + 63)]
    for k in range(0, len(bits), 6):
        x = 0
        for t in range(6):
            x = (x << 1) | bits[k + t]
        s.append(chr(x + 63))
    return ''.join(s)


ks = kernels(c)
print(f'c={c} kernels {len(ks)}', file=sys.stderr)
count = Counter()
w = sys.stdout.write
for nk, edges in ks:
    # per parallel class: multisets of subdivision numbers with at most one zero
    choices = []
    for (u, v, mlt) in edges:
        opts = [ms for ms in itertools.combinations_with_replacement((0, 1, 2), mlt) if ms.count(0) <= 1]
        choices.append(opts)
    for combo in itertools.product(*choices):
        n = nk
        E = []
        for (u, v, mlt), ms in zip(edges, combo):
            for t in ms:
                if t == 0:
                    E.append((u, v))
                else:
                    prev = u
                    for _ in range(t):
                        E.append((prev, n))
                        prev = n
                        n += 1
                    E.append((prev, v))
        count[n] += 1
        w(g6(n, E) + '\n')
for n in sorted(count):
    print(f'c={c} n={n} generated {count[n]}', file=sys.stderr)
print(f'c={c} total generated {sum(count.values())}', file=sys.stderr)
