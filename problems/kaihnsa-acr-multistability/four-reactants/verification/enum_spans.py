#!/usr/bin/env python3
"""Verifier's stoichiometric classification (runtime about 5 min on one core per type group).

Written from the mathematical specification of the manuscript, before reading the solver's scripts.

For each source type:
  * all triples of distinct allowed reaction vectors p - y_a; a rank-3 span is identified by its exact
    RREF over Q (the solver uses primitive Pluecker coordinates);
  * W = integer basis of the orthogonal complement; all allowed reactions (source, vector) with W v = 0;
  * "missing a source": a source with no allowed reaction in the span;
  * positive maximal flux: exact dual-cone test by cross products, compared on every span with the
    HiGHS LP "maximise t, Gx = 0, t <= x <= 1" (floating point, used only as a second opinion);
  * D_W coefficients det W_I det V_I (Cauchy-Binet): injective if nonzero and of one sign,
    class-degenerate if all zero;
  * kinetic rank: integer basis K of ker G, columns L_a = sum over reactions of source a; the cubic
    det[L1 L2 L3] on ker G is tested by exact symmetrisation of its coefficient tensor;
  * the rest: orbit of the row-reduced W without the ACR column under the full group of species
    permutations preserving the source set, computed by brute force over S_n.

Usage: enum_spans.py [S1 S2 ...]. Writes verification/spans-<type>.json. Prints PASS/FAIL per type.
"""
import json
import os
import sys
from fractions import Fraction as Fr
from itertools import combinations, permutations

import numpy as np
from scipy.optimize import linprog

from common import rref, nullspace, det, complexes, diff, named, PAPER_TYPES, add_n4_types

HERE = os.path.dirname(os.path.abspath(__file__))

# Manuscript's table, only for comparison at the end:
# spans, missing, no flux, injective, class-degenerate, kinetically degenerate, remaining
PAPER_COUNTS = {
    'S1': (4333, 2965, 1126, 169, 34, 39, 0),
    'S2': (4287, 3161, 782, 301, 6, 37, 0),
    'S3': (4333, 2965, 1126, 129, 6, 43, 64),
    'S7': (4397, 3050, 1042, 231, 32, 20, 22),
    'S9': (23035, 17734, 4442, 569, 144, 92, 54),
}
PAPER_FORMS = {'S3': 11, 'S7': 5, 'S9': 3}


def all_spans(vecs, n):
    """Distinct rank-3 spans of triples of vectors, keyed by exact RREF."""
    spans = {}
    # incremental elimination to save Fractions: reduce pair, then add a third vector
    for i, j in combinations(range(len(vecs)), 2):
        R2, p2 = rref([vecs[i], vecs[j]], n)
        if len(p2) < 2:
            continue
        for k in range(j + 1, len(vecs)):
            v = [Fr(x) for x in vecs[k]]
            for r, p in zip(R2, p2):
                if v[p] != 0:
                    f = v[p]
                    v = [a - f * b for a, b in zip(v, r)]
            if not any(v):
                continue
            R, piv = rref(R2 + [v], n)
            key = tuple(tuple(r) for r in R)
            if key not in spans:
                spans[key] = piv
    return spans


def perp_basis(R, piv, n):
    return nullspace([list(r) for r in R], n)


def dual_cone_nonzero(g):
    """True iff some nonzero y has g_i . y >= 0 for all i (columns g span R^3). Exact on int64."""
    m = len(g)
    I, J = np.triu_indices(m, 1)
    c = np.cross(g[I], g[J])
    c = c[np.any(c != 0, axis=1)]
    ip = c @ g.T
    return bool(np.any(np.all(ip >= 0, axis=1)) or np.any(np.all(ip <= 0, axis=1)))


def lp_positive(g):
    m = len(g)
    res = linprog(c=[0.0] * m + [-1.0], A_eq=np.hstack([g.T.astype(float), np.zeros((3, 1))]),
                  b_eq=np.zeros(3),
                  A_ub=np.hstack([-np.eye(m), np.ones((m, 1))]), b_ub=np.zeros(m),
                  bounds=[(0, 1)] * m + [(None, 1)], method='highs')
    return res.status == 0 and -res.fun > 1e-7


def cauchy_binet(W, V, n):
    k = n - 3
    out = []
    for I in combinations(range(n), k):
        a = det([[row[i] for i in I] for row in W])
        b = det([[v[i] for v in V] for i in I])  # V has basis vectors as list entries
        out.append(a * b)
    return out


def cubic_vanishes(g, src, K):
    """det[L1 L2 L3] on ker G, identically zero? g: m x 3 ints, src: source index per reaction,
    K: list of d integer kernel vectors of length m."""
    d = len(K)
    if d == 0:
        return True
    K = np.array(K, dtype=object).reshape(d, len(g))
    L = []
    for a in range(3):
        rows = [r for r in range(len(g)) if src[r] == a]
        M = np.zeros((3, d), dtype=object)
        for r in rows:
            M += np.outer(np.array(g[r], dtype=object), K[:, r])
        L.append(M)
    A, B, C = L
    # exterior coordinates: cross product of columns of A and B, dotted with columns of C
    T = {}
    for i in range(d):
        for j in range(d):
            u = A[:, i]
            w = B[:, j]
            x = (u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0])
            for k in range(d):
                c = x[0] * C[0, k] + x[1] * C[1, k] + x[2] * C[2, k]
                if c:
                    key = tuple(sorted((i, j, k)))
                    T[key] = T.get(key, 0) + c
    return all(v == 0 for v in T.values())


def automorphisms(sources, n):
    S = set(sources)
    out = []
    for p in permutations(range(n)):
        if {tuple(y[p[i]] for i in range(n)) for y in sources} == S:
            out.append(p)
    return out


def classify(name):
    names, src_txt, acr_txt = PAPER_TYPES[name]
    n = len(names)
    sources = [named(names, t) for t in src_txt]
    acr = [names.index(c) for c in acr_txt]
    A = [diff(sources[i], sources[3]) for i in range(3)]
    V = nullspace(A, n)
    assert len(rref(A, n)[1]) == 3 and all(all(v[i] == 0 for v in V) for i in acr)

    reactions = []
    for a, y in enumerate(sources):
        for p in complexes(n):
            if p != y:
                reactions.append((a, diff(p, y)))
    vecs = sorted({v for _, v in reactions})
    spans = all_spans(vecs, n)
    counts = dict(spans=len(spans), missing=0, noflux=0, injective=0, classdeg=0, kindeg=0, remaining=0)
    remaining, labels, lp_disagree = [], [], 0
    auto = automorphisms(sources, n)

    allv = np.array([v for _, v in reactions], dtype=np.int64)
    alla = np.array([a for a, _ in reactions])

    for key, piv in spans.items():
        R = [list(r) for r in key]
        W = perp_basis(R, piv, n)
        Wm = np.array(W, dtype=np.int64)
        mask = np.all(Wm @ allv.T == 0, axis=0)
        src_here = alla[mask]
        if len(set(src_here.tolist())) < 4:
            counts['missing'] += 1
            continue
        g = allv[mask][:, piv]             # S is isomorphic to R^3 by the pivot coordinates
        assert np.linalg.matrix_rank(g) == 3
        cone_fails = dual_cone_nonzero(g)
        if cone_fails == lp_positive(g):   # exact: positive iff not cone_fails
            lp_disagree += 1
        if cone_fails:
            counts['noflux'] += 1
            continue
        coef = cauchy_binet(W, V, n)
        nz = [c for c in coef if c]
        if not nz:
            counts['classdeg'] += 1
            labels.append((W, 'class-degenerate'))
            continue
        if all(c > 0 for c in nz) or all(c < 0 for c in nz):
            counts['injective'] += 1
            labels.append((W, 'injective'))
            continue
        K = nullspace(g.T.tolist(), len(g))
        if cubic_vanishes(g.tolist(), src_here.tolist(), K):
            counts['kindeg'] += 1
            labels.append((W, 'kinetic-degenerate'))
            continue
        counts['remaining'] += 1
        labels.append((W, 'normal-form'))
        remaining.append((W, coef))

    # Orbits of remaining spans. W without the (first) ACR column, row reduced, modulo automorphisms.
    i0 = acr[0]
    forms = {}
    records = []
    for W, coef in remaining:
        best = None
        for p in auto:
            # permuting species by p sends the column of p[i] to position i
            Wp = [[row[p[i]] for i in range(n)] for row in W]
            # the ACR coordinate is unique for the types that reach this stage
            Rp = rref([r[:i0] + r[i0 + 1:] for r in Wp], n - 1)[0]
            k = tuple(tuple(str(x) for x in r) for r in Rp)
            if best is None or k < best:
                best = k
        forms.setdefault(best, 0)
        forms[best] += 1
        records.append({'W': W, 'coef': coef, 'form': [list(r) for r in best]})

    out = {'type': name, 'n': n, 'names': names, 'sources': src_txt, 'counts': counts,
           'lp_disagreements': lp_disagree, 'automorphisms': len(auto),
           'forms': [{'Wprime': [list(r) for r in k], 'spans': v} for k, v in sorted(forms.items())],
           'remaining': records, 'labels': labels}
    with open(os.path.join(HERE, f'spans-{name}.json'), 'w') as f:
        json.dump(out, f)
    return out


def main():
    add_n4_types()
    if sys.argv[1:] == ['n4']:
        todo = [k for k in PAPER_TYPES if k.startswith('T')]
    else:
        todo = sys.argv[1:] or [k for k in PAPER_TYPES if not k.startswith('T')]
    ok = True
    for name in todo:
        r = classify(name)
        c = r['counts']
        tup = (c['spans'], c['missing'], c['noflux'], c['injective'], c['classdeg'], c['kindeg'],
               c['remaining'])
        print(name, 'counts', tup, 'forms', len(r['forms']), 'LP disagreements', r['lp_disagreements'],
              'automorphisms', r['automorphisms'], flush=True)
        good = sum(tup[1:]) == tup[0] and r['lp_disagreements'] == 0
        if name in PAPER_COUNTS:
            good = good and tup == PAPER_COUNTS[name]
            print(('PASS' if tup == PAPER_COUNTS[name] else 'FAIL'), name, 'matches manuscript table', flush=True)
        if name in PAPER_FORMS:
            print(('PASS' if len(r['forms']) == PAPER_FORMS[name] else 'FAIL'), name, 'number of normal forms',
                  flush=True)
            good = good and len(r['forms']) == PAPER_FORMS[name]
        print(('PASS' if r['lp_disagreements'] == 0 else 'FAIL'), name, 'exact cone test agrees with LP', flush=True)
        ok = ok and good
    print('RESULT', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
