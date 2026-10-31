#!/usr/bin/env python3
"""Verifier's re-derivation of the five-source network (runtime about 10 s).

Everything is rebuilt from the reaction list with this file's own mass-action code. No coordinates of the
manuscript are used for stability: the Jacobian is restricted to S in a basis of S computed from the
reaction vectors, the Routh-Hurwitz conditions are evaluated exactly on the 3x3 matrix, and the real roots
of the characteristic polynomial are isolated exactly. The list of positive equilibria in the class a+b=9
is obtained with a lex Groebner basis (SymPy) and a second time with msolve.
Also checks: r = 1 at every positive equilibrium of every class, existence of a positive equilibrium in every
class, molecularity, the variant coupled to C, and the behaviour for nearby rate vectors.
"""
import random
import sys
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

from stress_count import msolve

ok = True


def report(name, cond, extra=''):
    global ok
    print(('PASS' if cond else 'FAIL'), name, extra, flush=True)
    ok = ok and bool(cond)


NAMES = 'ABCR'
BASE = [('0', 'C', 64), ('A+B', '2B', 3), ('A+C', 'B+C', 5), ('B+C', 'A', 64), ('0', 'R', 1), ('R', '0', 1)]
COUPLED = [('0', 'C', 62), ('A+B', '2B', 3), ('A+C', 'B+C', 5), ('B+C', 'A', 64), ('0', 'C+R', 1), ('R', 'C', 1)]
x = sp.symbols('a b c r', positive=True)


def cplx(text):
    y = [0, 0, 0, 0]
    if text == '0':
        return tuple(y)
    for tok in text.split('+'):
        k, s = (int(tok[0]), tok[1:]) if tok[0].isdigit() else (1, tok)
        y[NAMES.index(s)] += k
    return tuple(y)


def field(net):
    f = [sp.Integer(0)] * 4
    for src, prod, k in net:
        y, z = cplx(src), cplx(prod)
        rate = k * sp.prod([xi**e for xi, e in zip(x, y)])
        f = [fi + (zj - yj) * rate for fi, yj, zj in zip(f, y, z)]
    return sp.Matrix([sp.expand(e) for e in f])


def stoich(net):
    return sp.Matrix([[b - a for a, b in zip(cplx(s), cplx(p))] for s, p, _ in net]).T


def routh_hurwitz3(M):
    lam = sp.Symbol('lam')
    cp = sp.Poly(M.charpoly(lam).as_expr(), lam)
    c = cp.all_coeffs()                         # 1, a2, a1, a0
    a2, a1, a0 = c[1], c[2], c[3]
    return a2 > 0 and a0 > 0 and a2 * a1 - a0 > 0, (a2, a1, a0), cp


def restrict(J, Sb):
    """Matrix of J|_S in the basis Sb (columns), exact."""
    C = (Sb.T * Sb).inv() * Sb.T * J * Sb
    assert sp.simplify(J * Sb - Sb * C) == sp.zeros(*(J * Sb).shape), 'S is not J-invariant'
    return C


def analyse(net, label):
    srcs = {s for s, _, _ in net}
    report(f'{label}: exactly five distinct reactant complexes', len(srcs) == 5, sorted(srcs))
    report(f'{label}: all complexes at most bimolecular',
           all(sum(cplx(s)) <= 2 and sum(cplx(p)) <= 2 for s, p, _ in net))
    N = stoich(net)
    report(f'{label}: dim S = 3', N.rank() == 3)
    Sb = sp.Matrix.hstack(*N.columnspace())
    report(f'{label}: S-perp spanned by (1,1,0,0)', [list(v) for v in N.T.nullspace()] == [[1, 1, 0, 0]])
    f = field(net)
    a, b, c, r = x
    # ACR for all compatibility classes: r is forced by a polynomial identity; existence for each T
    report(f'{label}: r-dot = 1 - r', sp.expand(f[3] - (1 - r)) == 0)
    # every class a+b=T has a positive equilibrium? (existence in at least one class is what ACR needs)
    report(f'{label}: a-dot + b-dot = 0', sp.expand(f[0] + f[1]) == 0)
    # positive equilibria in the class: Groebner basis, then exact real root isolation
    A, B, C, R = sp.symbols('A B C R')
    eqs = [e.subs(dict(zip(x, (A, B, C, R)))) for e in f] + [A + B - 9]
    G = sp.groebner(eqs, A, B, C, R, order='lex', domain='QQ')
    uni = [g for g in G.exprs if g.free_symbols == {R}]
    last = [g for g in G.exprs if g.free_symbols == {C}]
    # Eliminate directly: list of solutions
    sols = sp.solve_poly_system(list(G.exprs), A, B, C, R)
    sols = [s for s in sols if all(v.is_real and v > 0 for v in s)]
    sols = sorted(sols)
    print('   positive equilibria', sols)
    report(f'{label}: exact equilibria in a+b=9 by Groebner basis',
           sols == [(4, 5, sp.Rational(1, 5), 1), (6, 3, sp.Rational(1, 3), 1), (8, 1, 1, 1)])
    J = f.jacobian(x)
    verdicts = []
    for s in sols:
        Js = J.subs(dict(zip(x, s)))
        M = restrict(Js, Sb)
        stable, coeffs, cp = routh_hurwitz3(M)
        verdicts.append(stable)
        roots = sp.real_roots(cp)
        print('   ', s, 'char poly coefficients', coeffs, 'real roots', [sp.N(z, 8) for z in roots])
        if not stable:
            report(f'{label}: {s} has a positive real eigenvalue', any(z > 0 for z in roots))
            report(f'{label}: {s} nondegenerate', M.det() != 0)
    report(f'{label}: stable, unstable, stable', verdicts == [True, False, True])
    return sols


def msolve_check():
    """Certified real root isolation of the whole system; a second route to the same three positive points."""
    text = 'a,b,c,r\n0\n-3*a*b-5*a*c+64*b*c,\n64-64*b*c,\n1-r,\na+b-9\n'
    sols = msolve(text, 4)
    pos = [s for s in sols if all(lo > 0 for lo, hi in s)]
    want = [[Fr(4), Fr(5), Fr(1, 5), Fr(1)], [Fr(6), Fr(3), Fr(1, 3), Fr(1)], [Fr(8), Fr(1), Fr(1), Fr(1)]]
    hit = sum(any(all(lo <= w <= hi for (lo, hi), w in zip(s, p)) for s in pos) for p in want)
    return len(sols), len(pos), hit


def persistence():
    """Robustness of the behaviour: 200 random rate vectors within 2 percent. Floating point."""
    rng = random.Random(7)
    good = 0
    for _ in range(200):
        k = [64, 3, 5, 64, 1, 1]
        k = [Fr(v) * (1 + Fr(rng.randint(-200, 200), 10000)) for v in k]
        # in the class a+b=9: bc = k0/k3, so c = (k0/k3)/b; g(a) = k0 - k1 a b - k2 a c
        kk = [float(v) for v in k]
        c0 = kk[0] / kk[3]
        a = np.polynomial.Polynomial([0, 1])
        b = 9 - a
        g = kk[0] * b - kk[1] * a * b * b - kk[2] * c0 * a
        rts = [z.real for z in g.roots() if abs(z.imag) < 1e-9 and 0 < z.real < 9]
        stable = 0
        for ai in rts:
            bi = 9 - ai
            ci = c0 / bi
            # (a,c) block
            J = np.array([[-kk[1] * bi + kk[1] * ai - kk[2] * ci - kk[3] * ci, -kk[2] * ai + kk[3] * bi],
                          [kk[3] * ci, -kk[3] * bi]])
            stable += np.linalg.eigvals(J).real.max() < 0
        good += stable == 2 and len(rts) == 3
    return good


def acr_all_classes():
    """Existence of a positive equilibrium in every class a+b=T, T>0, for every T (not needed for ACR)."""
    a, T = sp.symbols('a T', positive=True)
    b = T - a
    g = 64 * b - 3 * a * b**2 - 5 * a
    # g(0+) = 64 T > 0 and g(T-) = -5 T < 0: a root in (0,T) by the intermediate value theorem
    return sp.simplify(g.subs(a, 0) - 64 * T) == 0 and sp.simplify(g.subs(a, T) + 5 * T) == 0


if __name__ == '__main__':
    analyse(BASE, 'base')
    analyse(COUPLED, 'coupled')
    report('positive equilibrium in every class a+b=T', acr_all_classes())
    total, pos, hit = msolve_check()
    report('msolve: exactly three positive solutions, boxes contain P1, P2, P3', pos == 3 and hit == 3,
           f'{total} real solutions in all')
    good = persistence()
    report('200 of 200 nearby rate vectors keep 3 equilibria, 2 stable (floating point)', good == 200, good)
    print('RESULT', 'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)
