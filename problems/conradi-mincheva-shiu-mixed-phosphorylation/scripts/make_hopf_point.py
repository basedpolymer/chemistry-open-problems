#!/usr/bin/env python3
"""Generate certificates/hopf_point.json: rational isolating interval for the Hopf parameter mu0 and the
normalisation indices used by certify_subcritical_hopf.py. (Not itself a certificate: the certificate
script re-checks everything from the stored data.)

Environment: Python 3.12, sympy 1.14, python-flint 0.9, numpy. Runtime: a few seconds.
Usage: python make_hopf_point.py "k1,...,k10" x1 x6 mu_guess out.json
"""
import json
import sys

import numpy as np
import sympy as sp
from flint import fmpq, fmpq_poly

kstr = sys.argv[1].split(','); x1 = sp.Rational(sys.argv[2]); x6 = sp.Rational(sys.argv[3])
mu_guess = sp.Rational(sys.argv[4]); out = sys.argv[5]
k = [sp.Rational(s) for s in kstr]
# model definitions (same conventions as certify_subcritical_hopf.py)
REACTIONS = [((0, 1), (2,)), ((2,), (0, 1)), ((2,), (3,)), ((3,), (4, 1)), ((4, 5), (6,)),
             ((6,), (4, 5)), ((6,), (7, 5)), ((7, 5), (8,)), ((8,), (7, 5)), ((8,), (0, 5))]
Nm = sp.zeros(9, 10)
for j, (re, pr) in enumerate(REACTIONS):
    for i in re:
        Nm[i, j] -= 1
    for i in pr:
        Nm[i, j] += 1
ZIDX = [2, 3, 4, 6, 7, 8]
E = sp.zeros(9, 6)
for j, i in enumerate(ZIDX):
    E[i, j] = 1
E[0, :] = sp.Matrix([[-1] * 6]); E[1, :] = sp.Matrix([[-1, -1, 0, 0, 0, 0]]); E[5, :] = sp.Matrix([[0, 0, 0, -1, 0, -1]])
Psel = sp.zeros(6, 9)
for j, i in enumerate(ZIDX):
    Psel[j, i] = 1
Wm = sp.Matrix([[0, 1, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 1], [1, 0, 1, 1, 1, 0, 1, 1, 1]])
mu, lam = sp.symbols('mu lam')
k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = k
xs = sp.Matrix([x1, mu, k1 / (k2 + k3) * x1 * mu, k1 * k3 / ((k2 + k3) * k4) * x1 * mu,
                k1 * k3 * (k6 + k7) / ((k2 + k3) * k5 * k7) * x1 * mu / x6, x6,
                k1 * k3 / ((k2 + k3) * k7) * x1 * mu, k1 * k3 * (k9 + k10) / ((k2 + k3) * k8 * k10) * x1 * mu / x6,
                k1 * k3 / ((k2 + k3) * k10) * x1 * mu])
D = sp.zeros(10, 9)
for j, (re, pr) in enumerate(REACTIONS):
    for a_ in re:
        v = k[j]
        for b_ in re:
            if b_ != a_:
                v = v * xs[b_]
        D[j, a_] += v
A = (Psel * Nm * D * E).applyfunc(sp.expand)
b = [sp.expand(c) for c in sp.Poly((lam * sp.eye(6) - A).det(method='berkowitz'), lam).all_coeffs()]
H5 = sp.Matrix(5, 5, lambda r, s: b[2 * (r + 1) - (s + 1)] if 0 <= 2 * (r + 1) - (s + 1) <= 6 else 0).det(method='berkowitz')
cs = sp.Poly(sp.expand(H5), mu).all_coeffs()[::-1]
h = fmpq_poly([fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in cs])
a = fmpq(int((mu_guess * sp.Rational(999, 1000)).p), int((mu_guess * sp.Rational(999, 1000)).q))
bb = fmpq(int((mu_guess * sp.Rational(1001, 1000)).p), int((mu_guess * sp.Rational(1001, 1000)).q))
assert h(a) * h(bb) < 0
for _ in range(200):
    m = (a + bb) / 2
    if h(a) * h(m) <= 0:
        bb = m
    else:
        a = m
mu0 = float(a)
An = np.array([[float(A[i, j].subs(mu, sp.Rational(a.p, a.q))) for j in range(6)] for i in range(6)])
ev, R = np.linalg.eig(An); i = np.argmax(ev.imag); om = ev[i].imag
Mq = An - 1j * om * np.eye(6)


def choose(Mat):
    U, S, Vh = np.linalg.svd(Mat)
    v = Vh[-1].conj(); j = int(np.argmax(np.abs(v)))
    best = None
    for r in range(6):
        sub = np.delete(np.delete(Mat, r, 0), j, 1)
        c = np.linalg.cond(sub)
        if best is None or c < best[0]:
            best = (c, r)
    return j, best[1]

jq, rq = choose(Mq); jw, rw = choose(Mq.T)
T = (Wm * xs).applyfunc(sp.expand)
json.dump({'k': kstr, 'x1': str(x1), 'x6': str(x6),
           'mu_interval': [str(sp.Rational(int(a.p), int(a.q))), str(sp.Rational(int(bb.p), int(bb.q)))],
           'mu0_approx': mu0, 'omega0_approx': om,
           'normalization': [jq, rq, jw, rw],
           'totals_of_mu': [str(t) for t in T],
           'note': 'mu = x2 (free kinase at the steady state); steady state phi(x1, mu, x6); totals (Ktot, Ptot, Stot)'},
          open(out, 'w'), indent=1)
print('mu0 ~ %.15g, omega0 ~ %.10g, interval width %.3e, normalization %s' % (mu0, om, float(bb - a), (jq, rq, jw, rw)))
print('totals(mu):', [str(t) for t in T])
