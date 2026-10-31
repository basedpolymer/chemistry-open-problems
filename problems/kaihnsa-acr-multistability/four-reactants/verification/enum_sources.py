#!/usr/bin/env python3
"""Verifier's source enumeration (runtime about 10 s, Python 3.12, standard library).

For 3 <= n <= 8 enumerate all sets of four distinct at-most-bimolecular complexes that involve all n
species, keep those with rank A = 3 and a coordinate i with e_i in row(A) (computed from an exact
nullspace, not from minors), and count orbits under S_n by union-find on the retained family.
Prints PASS/FAIL. Writes verification/sources-verifier.json.
"""
import json
import os
import sys
from itertools import combinations
from common import rref, nullspace, complexes, diff, named, PAPER_TYPES

HERE = os.path.dirname(os.path.abspath(__file__))


def analyse(S, n):
    S = sorted(S)
    A = [diff(S[i], S[3]) for i in range(3)]
    if len(rref(A, n)[1]) != 3:
        return None
    V = nullspace(A, n)
    acr = [i for i in range(n) if all(v[i] == 0 for v in V)]
    return (V, acr) if acr else None


def orbits(family, n):
    parent = {s: s for s in family}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for s in family:
        for i in range(n - 1):
            p = list(range(n))
            p[i], p[i + 1] = p[i + 1], p[i]
            t = frozenset(tuple(y[p[j]] for j in range(n)) for y in s)
            parent[find(s)] = find(t)
    groups = {}
    for s in family:
        groups.setdefault(find(s), []).append(s)
    return list(groups.values())


def fmt(y):
    names = 'ABCDEFGH'
    t = []
    for i, c in enumerate(y):
        t += [names[i]] * c
    return '+'.join(t) if t else '0'


def main():
    ok = True

    def report(name, cond, extra=''):
        nonlocal ok
        print(('PASS' if cond else 'FAIL'), name, extra)
        ok = ok and bool(cond)

    result = {}
    all_orbits = {}
    for n in range(3, 9):
        C = complexes(n)
        fam = {}
        for S in combinations(C, 4):
            if any(all(y[i] == 0 for y in S) for i in range(n)):
                continue
            a = analyse(S, n)
            if a:
                fam[frozenset(S)] = a
        orbs = orbits(list(fam), n) if fam else []
        all_orbits[n] = (orbs, fam)
        result[n] = {'sets': len(fam), 'orbits': len(orbs)}
        print(f'n={n}: {len(fam)} labelled source sets, {len(orbs)} orbits under S_n')

    report('no type with 7 species', result[7]['orbits'] == 0)
    report('no type with 8 species', result[8]['orbits'] == 0)
    report('one type with 6 species', result[6]['orbits'] == 1)
    report('eight types with 5 species', result[5]['orbits'] == 8)
    report('26 types with 4 species', result[4]['orbits'] == 26, f"found {result[4]['orbits']}")

    # n = 4: distinct kernel entries, including zero, and a listing
    listing = []
    for orb in all_orbits[4][0]:
        s = min(orb, key=lambda z: sorted(z))
        V, acr = all_orbits[4][1][s]
        listing.append({'sources': [fmt(y) for y in sorted(s)], 'v': V[0], 'acr': acr,
                        'distinct': len(set(V[0]))})
    # checked on every labelled set, not only on orbit representatives
    worst = max(len(set(V[0])) for V, _ in all_orbits[4][1].values())
    report('n=4: at most 3 distinct kernel entries, incl. 0',
           worst <= 3 and all(0 in x['v'] for x in listing), f'max {worst}')
    report('n=4: kernel is one-dimensional', all(len(V) == 1 for V, _ in all_orbits[4][1].values()))
    report('n=3: toric set is a singleton', all(len(V) == 0 for V, _ in all_orbits[3][1].values())
           if all_orbits[3][1] else True, f"{result[3]['sets']} sets with an ACR coordinate")

    # the paper's S-table against the orbits
    orbs_named = {}
    for name, (names, src, acr) in PAPER_TYPES.items():
        n = len(names)
        S = frozenset(named(names, t) for t in src)
        hit = S in all_orbits[n][1]
        got = ''.join(names[i] for i in all_orbits[n][1][S][1]) if hit else None
        report(f'{name} is a retained type with ACR {acr}', hit and got == acr, f'got {got}')
        for i, orb in enumerate(all_orbits[n][0]):
            if S in orb:
                orbs_named[(n, i)] = name
    report('S-table hits every orbit exactly once (n=5,6)', len(orbs_named) == 9)

    with open(os.path.join(HERE, 'sources-verifier.json'), 'w') as f:
        json.dump({'counts': result, 'n4': listing}, f, indent=1)
    print('RESULT', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
