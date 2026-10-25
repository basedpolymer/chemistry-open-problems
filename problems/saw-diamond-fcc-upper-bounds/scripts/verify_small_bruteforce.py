#!/usr/bin/env python3
"""Brute-force checker for small LATCERT1 certificates (diamond and fcc lattices), pure Python, exact integers.

Environment: Python 3.10+, standard library only.  Runtime: about 1 min for 20 000 states (diamond K = 16) and
about 2 min for 3 000 states (fcc K = 7); use it for certificates up to ~10^5 states.

It re-implements the checks of src/latsaw_check.cpp in a different way:
  * the lattice is given by its vertex test and its neighbour function only;
  * the automorphisms gamma(x) = e + g x with gamma(o) = e are found by trying all 48 signed permutation matrices g
    and keeping those for which gamma maps the neighbours of every vertex of a box onto the neighbours of its image
    (a local test of the automorphism property; the global statement is proved in the paper);
  * containment gamma(F(child)) <= F(s) u {o} is tested with Python sets;
  * the ratio max_s sum_e V(child)/V(s) is computed with fractions.Fraction and compared with P/Q.
Prints PASS/FAIL per claim; exit code 0 iff everything passes.

Usage: python3 verify_small_bruteforce.py file.lcert
"""
import itertools
import struct
import sys
from fractions import Fraction

FCC = [(-1, -1, 0), (-1, 1, 0), (1, -1, 0), (1, 1, 0), (-1, 0, -1), (-1, 0, 1), (1, 0, -1), (1, 0, 1),
       (0, -1, -1), (0, -1, 1), (0, 1, -1), (0, 1, 1)]
DIA = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]


def main():
    path = sys.argv[1]
    data = open(path, 'rb').read()
    ok_all = True

    def report(ok, msg):
        nonlocal ok_all
        print(('PASS ' if ok else 'FAIL ') + msg)
        ok_all = ok_all and ok

    if data[:8] != b'LATCERT1':
        report(False, 'magic')
        return 1
    lat, K = struct.unpack_from('<II', data, 8)
    N, P, Q = struct.unpack_from('<QQQ', data, 16)
    if lat == 1:
        steps, bits, name = FCC, 4, 'fcc'
    else:
        steps, bits, name = DIA, 2, 'diamond'
    nd = len(steps)
    print(f'{path}: lattice {name}, K = {K}, N = {N}, claimed mu <= {P}/{Q}')
    report(len(data) == 40 + 16 * N + 4 * nd * N, 'file size')
    codes = struct.unpack_from(f'<{N}Q', data, 40)
    V = struct.unpack_from(f'<{N}Q', data, 40 + 8 * N)
    child = struct.unpack_from(f'<{N * nd}I', data, 40 + 16 * N)

    def is_vertex(x):
        if lat == 1:
            return sum(x) % 2 == 0
        if all(c % 2 == 0 for c in x):
            return sum(x) % 4 == 0
        if all(c % 2 == 1 for c in x):
            return sum(x) % 4 == 3
        return False

    def neighbours(x):
        if lat == 1:
            return {tuple(a + b for a, b in zip(x, u)) for u in steps}
        sgn = 1 if x[0] % 2 == 0 else -1
        return {tuple(a + sgn * b for a, b in zip(x, u)) for u in steps}

    origin = (0, 0, 0)
    report(neighbours(origin) == set(steps) and all(is_vertex(u) for u in steps), 'neighbours of the origin')

    # automorphisms gamma(x) = e + g x with gamma(o) = e, for every neighbour e of o
    mats = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            mats.append((perm, signs))

    def gmap(g, e, x):
        perm, signs = g
        y = [0, 0, 0]
        for c in range(3):
            y[perm[c]] = signs[c] * x[c]
        return (e[0] + y[0], e[1] + y[1], e[2] + y[2])

    box = [x for x in itertools.product(range(-4, 5), repeat=3) if is_vertex(x)]
    gam = {}
    for e in steps:
        gam[e] = []
        for g in mats:
            good = gmap(g, e, origin) == e
            for x in box:
                if not good:
                    break
                y = gmap(g, e, x)
                if not is_vertex(y) or {gmap(g, e, z) for z in neighbours(x)} != neighbours(y):
                    good = False
            if good:
                gam[e].append(g)
    expected = 48 if lat == 1 else 24
    report(all(len(gam[e]) == expected for e in steps), f'{expected} automorphisms fixing the frame for every neighbour e')

    def positions(c):
        assert c >> 63 == 1
        l = (c >> 57) & 63
        assert l * bits <= 57 and ((c & ((1 << 57) - 1)) >> (l * bits)) == 0
        p = [origin]
        for i in range(l):
            t = (c >> (bits * i)) & ((1 << bits) - 1)
            assert t < nd
            u = steps[t]
            s = 1 if lat == 1 else (1 if i % 2 == 1 else -1)
            # p_{i+1} = p_i - s_i, s_i = displacement from v_{i+1} to v_i
            p.append(tuple(a - s * b for a, b in zip(p[-1], u)))
        return p

    report(all(codes[i] < codes[i + 1] for i in range(N - 1)), 'codes strictly increasing')
    report(all(1 <= v <= 1 << 60 for v in V), 'weights in [1, 2^60]')
    report((1 << 63) in codes, 'empty walk listed')
    pos = [positions(c) for c in codes]
    report(all(len(set(p)) == len(p) and all(is_vertex(x) for x in p) for p in pos),
           'every listed walk is a self-avoiding walk of the lattice')
    fails_child = fails_cont = fails_ineq = 0
    best = Fraction(0)
    for s in range(N):
        Fs = set(pos[s][1:])
        allowed_set = Fs | {origin}
        total = 0
        for j, e in enumerate(steps):
            c = child[s * nd + j]
            if e in Fs:
                continue
            if c >= N:
                fails_child += 1
                continue
            Fc = pos[c][1:]
            if not any(all(gmap(g, e, x) in allowed_set for x in Fc) for g in gam[e]):
                fails_cont += 1
                continue
            total += V[c]
        if Q * total > P * V[s]:
            fails_ineq += 1
        best = max(best, Fraction(total, V[s]))
    report(fails_child == 0, f'every allowed step has a listed child (failures {fails_child})')
    report(fails_cont == 0, f'containment gamma(F(child)) <= F(s) u {{o}} (failures {fails_cont})')
    report(fails_ineq == 0, f'Q*sum <= P*V(s) for all states (failures {fails_ineq})')
    report(best <= Fraction(P, Q), f'exact max ratio {best.numerator}/{best.denominator} = {float(best):.12f} <= P/Q')
    print(('PASS' if ok_all else 'FAIL') + f' => mu({name}) <= {P}/{Q}' if ok_all else 'FAIL')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
