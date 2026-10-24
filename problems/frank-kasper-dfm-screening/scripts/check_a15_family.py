#!/usr/bin/env python3
"""Certify a lower bound for Reddy's A15 formula: SI Eq. S16 / arXiv Eq. A16.

Python 3.12+, python-flint; runtime < 1 second. No optimization is used.
The proof covers the polynomial expression for every c >= 0, hence every
admissible member of this one-parameter flat-face family. It does NOT bound
arbitrary A15 curved partitions or arbitrary weighted Voronoi diagrams.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json

from flint import arb, ctx

ctx.prec = 192
ROOT = Path(__file__).resolve().parents[1]


def ball(q):
    q = Q(q)
    return arb(q.numerator)/q.denominator


def multiply(p, q):
    r = [arb(0) for _ in range(len(p)+len(q)-1)]
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i+j] += a*b
    return r


def bernstein(p, lo, hi):
    """Power coefficients -> Bernstein coefficients on the exact interval."""
    n = len(p)-1
    x, h = ball(lo), ball(hi-lo)
    powers = [h**i*sum((p[j]*comb(j, i)*x**(j-i)
                       for j in range(i, n+1)), arb(0)) for i in range(n+1)]
    return [sum((powers[i]*comb(k, i)/comb(n, i)
                 for i in range(k+1)), arb(0)) for k in range(n+1)]


def main():
    a = 12+12*arb(6).sqrt()
    b = 6*arb(5).sqrt()-4*arb(6).sqrt()-3
    if not (a > 0 and b > 0):
        raise ValueError('coefficient signs not certified')
    # F(c)^3 = 5/55296 (a+b c^2)^2 (3c^4-5c^3+10).
    lower = Q(10699, 10000)
    p = [5*x for x in multiply(multiply([a, arb(0), b], [a, arb(0), b]),
                              [arb(10), arb(0), arb(0), arb(-5), arb(3)])]
    p[0] -= 55296*ball(lower)**3
    pending, leaves = [(Q(0), Q(2), 0)], []
    while pending:
        lo, hi, depth = pending.pop()
        cs = bernstein(p, lo, hi)
        if all(x > 0 for x in cs):
            leaves.append([str(lo), str(hi)])
        else:
            if depth >= 24:
                raise ValueError('could not certify polynomial positivity')
            mid = (lo+hi)/2
            pending.extend(((lo, mid, depth+1), (mid, hi, depth+1)))
    # For c>=2, derivative sign is that of
    # q(c)=b*c^3*(24*c-35)+3*a*c*(4*c-5)+40*b > 0.
    # This proves the expression increasing there. Positivity on [0,2]
    # therefore completes the proof for the whole half-line.
    limit_cube = 50*a*a/55296
    limit = limit_cube.root(3)
    if not (limit > ball(Q(1157, 1000)) and limit < ball(Q(1158, 1000))):
        raise ValueError('unexpected c=0 limit')
    at_sample = lambda c: (5*(a+b*c*c)**2*(3*c**4-5*c**3+10)/55296).root(3)
    sample = at_sample(ball(Q(122, 100)))
    if not sample < ball(Q(10701, 10000)):
        raise ValueError('sample upper bound failed')
    report = {'status': 'PASS', 'source': 'PNAS DOI 10.1073/pnas.1809655115, SI Eq. (S16); arXiv:1806.01726v1, Eq. (A16)',
              'scope': 'one-parameter flat A15 formula, c >= 0; not general curved A15',
              'strict_lower_bound': str(lower), 'upper_at_c_1_22': str(sample),
              'bernstein_cover_of_0_to_2': sorted(leaves, key=lambda row: Q(row[0])),
              'c_ge_2': 'derivative has strictly positive grouped terms',
              'c_zero_limit': str(limit), 'precision_bits': ctx.prec}
    (ROOT/'certificates'/'a15_family_report.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS flat A15 family: F(c) > 1.0699 for all c >= 0')
    print('PASS F(1.22) < 1.0701; this is only a feasible sample, not a claimed minimizer')
    print('PASS finite degeneration limit in (1.157, 1.158):', limit)
    print('PASS Bernstein cover:', len(leaves), 'intervals; 192-bit outward enclosures')


if __name__ == '__main__':
    main()
