#!/usr/bin/env python3
"""Why does the mutant "dual cone test disabled" survive? (runtime about 1 min)

Applies the injectivity, class-degeneracy and kinetic-rank tests to the spans that the positive-flux test
rejects, and reports how they are labelled. If every rejected span lands in one of the three exclusions or
in one of the 19 normal forms, the positive-flux exclusion is redundant for the proof.
"""
import json
import os
import sys
from collections import Counter

import numpy as np

from common import rref, nullspace, named, diff, PAPER_TYPES, complexes
from enum_spans import all_spans, perp_basis, dual_cone_nonzero, cauchy_binet, cubic_vanishes, automorphisms

HERE = os.path.dirname(os.path.abspath(__file__))


def run(name):
    names, src_txt, acr_txt = PAPER_TYPES[name]
    n = len(names)
    sources = [named(names, t) for t in src_txt]
    i0 = names.index(acr_txt[0])
    A = [diff(sources[i], sources[3]) for i in range(3)]
    V = nullspace(A, n)
    reactions = [(a, diff(p, y)) for a, y in enumerate(sources) for p in complexes(n) if p != y]
    vecs = sorted({v for _, v in reactions})
    spans = all_spans(vecs, n)
    allv = np.array([v for _, v in reactions], dtype=np.int64)
    alla = np.array([a for a, _ in reactions])
    auto = automorphisms(sources, n)
    known = set()
    fname = os.path.join(HERE, f'spans-{name}.json')
    data = json.load(open(fname))
    known = {tuple(tuple(r) for r in x['Wprime']) for x in data['forms']}
    tally = Counter()
    extra_forms = Counter()
    for key, piv in spans.items():
        W = perp_basis([list(r) for r in key], piv, n)
        mask = np.all(np.array(W, dtype=np.int64) @ allv.T == 0, axis=0)
        src_here = alla[mask]
        if len(set(src_here.tolist())) < 4:
            continue
        g = allv[mask][:, piv]
        if not dual_cone_nonzero(g):
            continue
        coef = cauchy_binet(W, V, n)
        nz = [c for c in coef if c]
        if not nz:
            tally['class-degenerate'] += 1
        elif all(c > 0 for c in nz) or all(c < 0 for c in nz):
            tally['injective'] += 1
        elif cubic_vanishes(g.tolist(), src_here.tolist(), nullspace(g.T.tolist(), len(g))):
            tally['kinetic-degenerate'] += 1
        else:
            best = None
            for p in auto:
                Wp = [[row[p[i]] for i in range(n)] for row in W]
                R = rref([r[:i0] + r[i0 + 1:] for r in Wp], n - 1)[0]
                k = tuple(tuple(str(x) for x in r) for r in R)
                if best is None or k < best:
                    best = k
            if best in known:
                tally['in a known normal form'] += 1
            else:
                tally['NEW FORM'] += 1
                extra_forms[best] += 1
    return tally, extra_forms


def main():
    ok = True
    for name in ['S1', 'S2', 'S3', 'S7', 'S9']:
        tally, extra = run(name)
        print(name, 'spans rejected by the flux test, relabelled:', dict(tally), flush=True)
        for k, v in extra.items():
            print('   new form', v, k)
        ok = ok and tally.get('NEW FORM', 0) == 0
    print('RESULT', 'PASS: the positive-flux exclusion is redundant for the proof' if ok else
          'INFO: the positive-flux exclusion is needed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
