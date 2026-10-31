#!/usr/bin/env python3
"""Why does the kinetic-rank exclusion hold? Exact structural audit (runtime about 20 s).

The spans excluded only by the kinetic-rank test matter: on some of them the toric system has four positive
solutions in one class (see stress_count.py), so the theorem depends on this exclusion. The test is the vanishing
of the cubic det[L1 L2 L3] on ker G, an algebraic statement. This script reproves it without the cubic.

Let V_a be the span of the allowed reaction vectors of source a that lie in S. Each column N_a lies in V_a. A
nondegenerate equilibrium needs rank N = 3, so any three columns of N are independent. By Rado's theorem an
independent choice from V_a, V_b, V_c exists only if dim(sum of V_j, j in J) >= |J| for every subset J of the
triple. So a failure of this condition for one triple excludes the span, whatever the rate constants.

All ranks are computed in exact rational arithmetic. PASS if every kinetic-degenerate span fails the condition.
"""
import json
import os
import sys
from collections import Counter
from itertools import combinations

import numpy as np

from common import rref, named, diff, complexes, PAPER_TYPES

HERE = os.path.dirname(os.path.abspath(__file__))


def vectors_by_source(name, W):
    names, src, _ = PAPER_TYPES[name]
    n = len(names)
    sources = [named(names, t) for t in src]
    Wm = np.array(W, dtype=np.int64).reshape(-1, n)
    out = [[] for _ in sources]
    for a, y in enumerate(sources):
        for p in complexes(n):
            v = diff(p, y)
            if any(v) and not np.any(Wm @ np.array(v)):
                out[a].append(v)
    return out, n


def rado_failure(vs, n):
    for tr in combinations(range(4), 3):
        for size in (1, 2, 3):
            for J in combinations(tr, size):
                rows = [v for a in J for v in vs[a]]
                if len(rref(rows, n)[1] if rows else []) < size:
                    return tr, J
    return None


def main():
    ok = True
    for name in ['S1', 'S2', 'S3', 'S7', 'S9']:
        data = json.load(open(os.path.join(HERE, f'spans-{name}.json')))
        tally = Counter()
        for W, lab in data['labels']:
            vs, n = vectors_by_source(name, W)
            if lab == 'kinetic-degenerate':
                f = rado_failure(vs, n)
                tally['excluded by Rado' if f else 'NOT explained'] += 1
                if f:
                    tally[f'failing subset of size {len(f[1])}'] += 1
            elif lab == 'normal-form':
                # the exclusion must not fire on the spans the theorem keeps, or it would prove too much
                tally['surviving span also fails Rado' if rado_failure(vs, n) else 'surviving span satisfies Rado'] += 1
        print(name, dict(sorted(tally.items())), flush=True)
        ok = ok and tally['NOT explained'] == 0
    print('RESULT', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
