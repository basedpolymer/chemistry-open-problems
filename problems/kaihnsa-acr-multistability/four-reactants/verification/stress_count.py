#!/usr/bin/env python3
"""Solves the full equilibrium system of random classes without any scalar reduction (runtime 3-6 min, 3 cores).

For each span labelled normal-form in spans-*.json (types S3, S7, S9, also S4, S6, S8 which the manuscript treats
by the group argument), draw a random positive rational point x0, take the toric constants and the class of x0,
and solve   W'x = W'x0, binomials   with msolve (certified real root isolation in exact arithmetic). Count the
positive solutions and evaluate D = det(W' diag(x) V') at them.

PASS requires at most two positive solutions, and opposite signs of D when there are two.
Control: the same harness on mixed-sign spans labelled kinetic-degenerate, where the theorem says nothing.
Finding three or more solutions there shows that the harness can see a violation.

This is a test, not a proof. Sampling is deterministic (seeded).
"""
import json
import multiprocessing as mp
import os
import random
import re
import subprocess
import sys
import tempfile
from collections import Counter
from fractions import Fraction as Fr
from math import gcd
from functools import reduce

import numpy as np

from common import rref, nullspace, named, diff, PAPER_TYPES

HERE = os.path.dirname(os.path.abspath(__file__))
NUM = r'(-?\d+(?: / \d+(?:\^\d+)?)?)'
PAIR = re.compile(r'\[' + NUM + r', ' + NUM + r'\]')


def num(s):
    if '/' in s:
        a, b = s.split(' / ')
        base, _, e = b.partition('^')
        return Fr(int(a), int(base)**int(e or 1))
    return Fr(int(s))


def setup(name, W):
    names, src, acr = PAPER_TYPES[name]
    n = len(names)
    idx = [i for i in range(n) if names[i] not in acr]
    Wr = rref([[Fr(r[i]) for i in idx] for r in W], len(idx))[0]
    Y = [named(names, s) for s in src]
    A = [diff(Y[i], Y[3]) for i in range(3)]
    V = nullspace(A, n)
    Vr = [[v[i] for i in idx] for v in V]
    b = nullspace(Vr, len(idx))
    return Wr, Vr, b


def lcm_den(xs):
    return reduce(lambda a, b: a * b // gcd(a, b), [Fr(x).denominator for x in xs], 1)


def poly_text(Wr, b, x0, xclass=None):
    """Binomial constants from x0; the class from xclass (default x0)."""
    m = len(x0)
    var = [f'x{i}' for i in range(m)]
    polys = []
    for row in Wr:
        T = sum(r * x for r, x in zip(row, xclass if xclass is not None else x0))
        L = lcm_den(list(row) + [T])
        terms = [f'{int(r * L)}*{v}' for r, v in zip(row, var) if r != 0]
        polys.append('+'.join(terms + [str(int(-T * L))]).replace('+-', '-'))
    for row in b:
        c = Fr(1)
        for e, x in zip(row, x0):
            c *= x**e
        pos = '*'.join(f'{v}^{e}' for v, e in zip(var, row) if e > 0) or '1'
        neg = '*'.join(f'{v}^{-e}' for v, e in zip(var, row) if e < 0) or '1'
        polys.append(f'{c.denominator}*{pos}-{c.numerator}*{neg}')
    return ','.join(var) + '\n0\n' + ',\n'.join(polys) + '\n'


def msolve(text, m):
    with tempfile.TemporaryDirectory() as d:
        fi, fo = os.path.join(d, 'i'), os.path.join(d, 'o')
        open(fi, 'w').write(text)
        subprocess.run(['msolve', '-t', '1', '-f', fi, '-o', fo], capture_output=True, timeout=60)
        out = open(fo).read()
    head = out.strip()
    if head.startswith('[1,') or head.startswith('[1 ,'):
        return None                     # positive-dimensional
    if head.startswith('[-1'):
        return []                       # no complex solution
    p = [(num(a), num(b)) for a, b in PAIR.findall(out)]
    assert len(p) % m == 0, out[:300]
    return [p[i:i + m] for i in range(0, len(p), m)]


def job(arg):
    name, W, seed, nsamples, x0pool = arg
    rng = random.Random(seed)
    Wr, Vr, b = setup(name, W)
    m = len(Vr[0])
    hist = Counter()
    bad = []
    Wn = np.array([[float(x) for x in r] for r in Wr])
    Vn = np.array(Vr, dtype=float).T
    for _ in range(nsamples):
        x0 = [Fr(rng.randint(1, x0pool), rng.randint(1, x0pool)) for _ in range(m)]
        sols = msolve(poly_text(Wr, b, x0), m)
        if sols is None:
            hist['positive-dimensional'] += 1
            continue
        pos = [s for s in sols if all(lo > 0 for lo, hi in s)]
        hist[len(pos)] += 1
        if len(pos) == 0:
            bad.append(('x0 lost', x0))
        if len(pos) >= 2:
            D = [np.linalg.det(Wn @ np.diag([float((lo + hi) / 2) for lo, hi in s]) @ Vn) for s in pos]
            if len(pos) == 2 and not (D[0] * D[1] < 0):
                bad.append(('same sign', x0, D))
            if len(pos) > 2:
                bad.append(('more than two', x0, len(pos)))
    return name, hist, bad


def main():
    nsamples = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    tasks, control = [], []
    seed = 20261031
    for name in ['S3', 'S7', 'S9', 'S4', 'S6', 'S8']:
        data = json.load(open(os.path.join(HERE, f'spans-{name}.json')))
        for W, lab in data['labels']:
            seed += 1
            if lab == 'normal-form':
                tasks.append((name, W, seed, nsamples, 9))
            elif lab == 'kinetic-degenerate' and name in ('S3', 'S7', 'S9'):
                control.append((name, W, seed, nsamples // 2, 9))
    print(f'{len(tasks)} spans, {nsamples} random classes each; control {len(control)} spans', flush=True)
    ok = True
    with mp.Pool(3) as pool:
        tot = {}
        badall = []
        for name, hist, bad in pool.imap_unordered(job, tasks, chunksize=4):
            tot.setdefault(name, Counter()).update(hist)
            badall += [(name, b) for b in bad]
        for name in sorted(tot):
            print(name, 'positive solutions per class:', dict(sorted(tot[name].items(), key=str)), flush=True)
            if any(isinstance(k, int) and k > 2 for k in tot[name]):
                ok = False
        print('PASS' if ok else 'FAIL', 'no class with three or more positive solutions in surviving spans')
        same = [b for b in badall if b[1][0] == 'same sign']
        print('PASS' if not same else 'FAIL', 'two positive solutions always have opposite D signs', len(same))
        lost = [b for b in badall if b[1][0] == 'x0 lost']
        print('PASS' if not lost else 'FAIL', 'the seed point is always recovered', len(lost))
        ok = ok and not same and not lost
        ctot = Counter()
        for name, hist, bad in pool.imap_unordered(job, control, chunksize=4):
            ctot.update(hist)
        print('control, kinetic-degenerate mixed-sign spans:', dict(sorted(ctot.items(), key=str)))
        seen3 = any(isinstance(k, int) and k >= 3 for k in ctot)
        print('PASS' if seen3 else 'FAIL', 'harness detects classes with three or more solutions in the control')
        ok = ok and seen3
    print('RESULT', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
