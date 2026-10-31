#!/usr/bin/env python3
"""Span-by-span comparison of the verifier's classification with the solver's certificate.

Reads spans-*.json (verifier) and ../certificates/subspaces.json (solver) and compares, for every span
that reaches the injectivity stage, the exact RREF of S^perp and the label. Runtime about 20 s.
"""
import json
import os
import sys
from fractions import Fraction as Fr

from common import rref

HERE = os.path.dirname(os.path.abspath(__file__))


def key(W, n):
    R = rref([[Fr(x) for x in r] for r in W], n)[0]
    return tuple(tuple(r) for r in R)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'certificates', 'subspaces.json')
    solver = json.load(open(path))
    ok = True
    for name in ['S1', 'S2', 'S3', 'S7', 'S9']:
        mine = json.load(open(os.path.join(HERE, f'spans-{name}.json')))
        n = mine['n']
        a = {key(W, n): lab for W, lab in mine['labels']}
        b = {key(r['W'], n): r['classification'] for r in solver[name]['records']}
        same = a == b
        counts_same = solver[name]['spans'] == mine['counts']['spans'] and \
            solver[name]['counts'].get('missing-source', 0) == mine['counts']['missing'] and \
            solver[name]['counts'].get('no-positive-flux', 0) == mine['counts']['noflux']
        print('PASS' if same and counts_same else 'FAIL', name,
              f'{len(a)} spans labelled identically; totals of the two early filters equal')
        ok = ok and same and counts_same
    print('RESULT', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
