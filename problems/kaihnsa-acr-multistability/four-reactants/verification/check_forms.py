#!/usr/bin/env python3
"""Verifier's check of the nineteen normal forms (runtime about 1 min, sympy 1.14).

Part 1. The conservation tables of the manuscript, typed in by hand from paper.tex, are compared with the
orbits computed by enum_spans.py (files spans-S3/S7/S9.json).

Part 2. For each form the scalar reduction is re-derived with the quantity that the determinant really
equals up to a fixed factor, namely D = det(W' diag(x') V') (see the orientation formula: the augmented
determinant is c * D_W(x) / prod(x_i) with a fixed c != 0, proved in verification.md). The solver works
with ordinary binomial rows; here V' is the kernel basis from the source exponents.

Along a curve x(t) satisfying the class equations and all but one log-binomial equation L = 0,
D(t) = mu(t) * L'(t) holds identically. The script computes mu symbolically, prints it, and checks
  (a) the curve identities,
  (b) that L' * mu - D vanishes identically,
  (c) on 2000 exact rational sample points of the positive domain, that mu has one sign.
The factorisations of mu printed here were read by hand: each is a constant times positive quantities.

Part 3. Root-count lemmas for the scalar equations, exact.
"""
import json
import os
import random
import sys
from fractions import Fraction as Fr

import sympy as sp

from common import rref, nullspace, named, diff, PAPER_TYPES

HERE = os.path.dirname(os.path.abspath(__file__))
ok = True


def report(name, cond, extra=''):
    global ok
    print(('PASS' if cond else 'FAIL'), name, extra, flush=True)
    ok = ok and bool(cond)


H = Fr(1, 2)
# Hand transcription of the three tables, rows in the coordinates of the table header.
TABLE = {
    'S3': [
        [[0, 1, 0, 0], [0, 0, 1, 1]],
        [[1, 0, -2, 2], [0, 1, 1, 0]],
        [[1, 0, -1, 1], [0, 1, H, H]],
        [[1, 0, -1, 1], [0, 1, 1, 0]],
        [[1, 0, -1, 1], [0, 1, 1, 1]],
        [[1, 0, -1, 1], [0, 1, 2, 0]],
        [[1, 0, -H, 1], [0, 1, H, 0]],
        [[1, 0, 0, 0], [0, 1, 1, 1]],
        [[1, 0, 0, 1], [0, 1, 1, 0]],
        [[1, 0, H, H], [0, 1, H, H]],
        [[1, 0, 1, 1], [0, 1, 1, 1]],
    ],
    'S7': [  # columns a,b,c,e
        [[0, 1, 0, 1], [0, 0, 1, -1]],
        [[1, 0, -1, 2], [0, 1, 1, 0]],
        [[1, 0, -H, H], [0, 1, H, H]],
        [[1, 0, 0, -1], [0, 1, 0, 0]],
        [[1, 1, 0, 1], [0, 0, 1, 0]],
    ],
    'S9': [  # columns a,b,c,d,f
        [[0, 1, 0, 0, 1], [0, 0, 1, 0, -1], [0, 0, 0, 1, 0]],
        [[1, 0, -H, 0, H], [0, 1, H, 0, H], [0, 0, 0, 1, 0]],
        [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 1, 1]],
    ],
}
ACR_INDEX = {'S3': 4, 'S7': 3, 'S9': 4}


def canon(name, Rows):
    names, src, _ = PAPER_TYPES[name]
    n = len(names)
    sources = [named(names, t) for t in src]
    from itertools import permutations
    auto = [p for p in permutations(range(n))
            if {tuple(y[p[i]] for i in range(n)) for y in sources} == set(sources)]
    i0 = ACR_INDEX[name]
    best = None
    for p in auto:
        assert p[i0] == i0
        # Rows has the ACR column deleted, restore a zero column to reuse the permutation
        full = [r[:i0] + [Fr(0)] + r[i0:] for r in [[Fr(x) for x in r] for r in Rows]]
        Wp = [[r[p[i]] for i in range(n)] for r in full]
        R = rref([r[:i0] + r[i0 + 1:] for r in Wp], n - 1)[0]
        k = tuple(tuple(str(x) for x in r) for r in R)
        if best is None or k < best:
            best = k
    return best


def part1():
    for name, forms in TABLE.items():
        with open(os.path.join(HERE, f'spans-{name}.json')) as f:
            data = json.load(f)
        mine = {tuple(tuple(r) for r in x['Wprime']) for x in data['forms']}
        theirs = {canon(name, F) for F in forms}
        report(f'{name}: computed orbits equal the table in the manuscript',
               mine == theirs and len(theirs) == len(forms), f'{len(mine)} forms')
        print('   spans per form', sorted(x['spans'] for x in data['forms']))


t = sp.Symbol('t', positive=True)
p, q, e0 = sp.symbols('p q e0', positive=True)
T, U = sp.symbols('T U', real=True)
Vs = sp.Symbol('V', positive=True)


def S3_curves():
    """(a,b,c,d) with ab=p, cd=q. The third entry names the residual exponent row."""
    yield 0, (p / T, T, t, U - t), (0, 0, 1, 1), sp.log(q)
    for i, (al, be, ga, de) in zip(range(1, 7),
                                   [(2, 2, 1, 0), (1, 1, sp.Rational(1, 2), sp.Rational(1, 2)), (1, 1, 1, 0),
                                    (1, 1, 1, 1), (1, 1, 2, 0), (sp.Rational(1, 2), 1, sp.Rational(1, 2), 0)]):
        yield i, (T + al * t - be * q / t, U - ga * t - de * q / t, t, q / t), (1, 1, 0, 0), sp.log(p)
    yield 7, (T, p / T, t, U - p / T - t), (0, 0, 1, 1), sp.log(q)
    yield 8, (t, p / t, U - p / t, T - t), (0, 0, 1, 1), sp.log(q)
    a0 = sp.Symbol('a0', positive=True)
    # forms 9, 10: a is the unique positive solution of a - p/a = T - U, s = c + d is then fixed
    yield 9, (a0, p / a0, t, Vs - t), (0, 0, 1, 1), sp.log(q)
    yield 10, (a0, p / a0, t, Vs - t), (0, 0, 1, 1), sp.log(q)


def S7_curves():
    """(a,b,c,e) with ab=pe, c^2=qe."""
    yield 0, (p * (t - U) / (T - t + U), T - t + U, t, t - U), (0, 0, 2, -1), sp.log(q)
    # forms 1, 2: drop ab=pe; parametrise by c=t, e=t^2/q
    yield 1, (T + t - 2 * t**2 / q, U - t, t, t**2 / q), (1, 1, 0, -1), sp.log(p)
    yield 2, (T + t / 2 - t**2 / (2 * q), U - t / 2 - t**2 / (2 * q), t, t**2 / q), (1, 1, 0, -1), sp.log(p)
    yield 3, (T + t**2 / q, U, t, t**2 / q), (1, 1, 0, -1), sp.log(p)
    yield 4, (t, T - U**2 / q - t, U, U**2 / q), (1, 1, 0, -1), sp.log(p)


def S9_curves():
    """(a,b,c,d,f) with ab=pf, cd=qf."""
    yield 0, (p * t / (T - t), T - t, t + U, Vs, t), (0, 0, 1, 1, -1), sp.log(q)
    yield 1, (T + (q / Vs - 1) * t / 2, U - (q / Vs + 1) * t / 2, q * t / Vs, Vs, t), (1, 1, 0, 0, -1), sp.log(p)
    yield 2, (T, U, t, Vs - T * U / p - t, T * U / p), (0, 0, 1, 1, -1), sp.log(q)


OTHER = {
    'S3': {(0, 0, 1, 1): (1, 1, 0, 0), (1, 1, 0, 0): (0, 0, 1, 1)},
    'S7': {(0, 0, 2, -1): (1, 1, 0, -1), (1, 1, 0, -1): (0, 0, 2, -1)},
    'S9': {(0, 0, 1, 1, -1): (1, 1, 0, 0, -1), (1, 1, 0, 0, -1): (0, 0, 1, 1, -1)},
}


def kernel_basis(name):
    names, src, _ = PAPER_TYPES[name]
    n = len(names)
    Y = [named(names, s) for s in src]
    A = [diff(Y[i], Y[3]) for i in range(3)]
    V = nullspace(A, n)
    i0 = ACR_INDEX[name]
    assert all(v[i0] == 0 for v in V)
    return [v[:i0] + v[i0 + 1:] for v in V]


def rational_point(sym_list, rng):
    return {s: sp.Rational(rng.randint(1, 60), rng.randint(1, 12)) * (rng.choice([1, -1]) if not s.is_positive else 1)
            for s in sym_list}


def part2():
    rng = random.Random(20261031)
    curves = {'S3': S3_curves, 'S7': S7_curves, 'S9': S9_curves}
    for name, gen in curves.items():
        V = kernel_basis(name)
        for idx, x, resid, logc in gen():
            R = sp.Matrix([[sp.Rational(str(v)) for v in row] for row in TABLE[name][idx]])
            X = sp.Matrix(x)
            # (a) class equations along the curve
            cons = [sp.simplify(sp.diff(e, t)) for e in R * X]
            good_a = all(c == 0 for c in cons)
            # the other binomial holds identically: its monomial ratio is constant along the curve
            r0 = OTHER[name][tuple(resid)]
            assert all(sum(a * b for a, b in zip(r0, v)) == 0 for v in V)
            assert all(sum(a * b for a, b in zip(resid, v)) == 0 for v in V)
            mono = sp.prod([xi**e for xi, e in zip(x, r0)])
            good_c = sp.simplify(sp.diff(sp.simplify(mono), t)) == 0
            # L(t) = sum resid_i log x_i - log const; L'(t) = sum resid_i x_i'/x_i
            Lp = sp.together(sum(e * sp.diff(xi, t) / xi for e, xi in zip(resid, x)))
            D = (R * sp.diag(*x) * sp.Matrix(V).T).det()
            mu = sp.factor(sp.cancel(sp.together(D / Lp)))
            good_d = sp.simplify(mu * Lp - D) == 0
            # (c) sign of mu on exact points of the positive domain
            free = sorted((mu.free_symbols | X.free_symbols) - {t}, key=str)
            signs = set()
            tries = 0
            while tries < 2000 and len(signs) < 3:
                tries += 1
                sub = rational_point(free, rng)
                sub[t] = sp.Rational(rng.randint(1, 80), rng.randint(1, 12))
                vals = [xi.subs(sub) if hasattr(xi, 'subs') else xi for xi in x]
                if any(not v.is_real for v in vals) or not all(v > 0 for v in vals):
                    continue
                m = mu.subs(sub)
                signs.add(sp.sign(m))
            report(f'{name}/{idx}: curve, D = mu L\', mu has one sign on sampled domain',
                   good_a and good_c and good_d and len(signs) == 1 and 0 not in signs,
                   f'mu = {mu}, signs {signs}')


def part3():
    # exact root-count facts used in the scalar reductions, on a sample of exact parameters
    u = sp.Symbol('u', positive=True)
    # strictly log-concave product for forms 1-6 of S3 and forms 1,2 of S7 (symbolic second derivative)
    for (al, be, ga, de) in [(2, 2, 1, 0), (1, 1, sp.Rational(1, 2), sp.Rational(1, 2)), (1, 1, 1, 0), (1, 1, 1, 1),
                             (1, 1, 2, 0), (sp.Rational(1, 2), 1, sp.Rational(1, 2), 0)]:
        a = T + al * u - be * q / u
        b = U - ga * u - de * q / u
        report('S3: a strictly concave, b concave', sp.simplify(sp.diff(a, u, 2) + 2 * be * q / u**3) == 0
               and sp.simplify(sp.diff(b, u, 2) + 2 * de * q / u**3) == 0 and be > 0 and de >= 0)
    for (al, be, ga, de) in [(1, 2, 1, 0), (sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2))]:
        a = T * u + al - be / (q * u)
        b = U * u - ga - de / (q * u)
        report('S7: A strictly concave, B concave', sp.simplify(sp.diff(a, u, 2) + 2 * be / (q * u**3)) == 0
               and sp.simplify(sp.diff(b, u, 2) + 2 * de / (q * u**3)) == 0 and be > 0 and de >= 0)


if __name__ == '__main__':
    part1()
    part2()
    part3()
    print('RESULT', 'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)
