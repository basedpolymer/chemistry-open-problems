#!/usr/bin/env python3
"""
check_gale_formula.py -- sanity check (exact rational arithmetic) of Proposition 2.1(2) of the brief-E write-up:
for a 3x5 matrix N of rank 3 with kernel basis [r1, r2] and a reactant matrix Y with rank [Y; 1^T] = 4,
    det(N diag(lambda) Y^T) = gamma*eta * alpha^3 * prod(mu) * sum_i c_i v_i / mu_i
for lambda = alpha * mu(s), mu(s) = (1-s) r1 + s r2, v = r2 - r1, with gamma*eta independent of s and alpha.
Checks the identity on random integer instances at several rational s and alpha (the constant gamma*eta is
fitted on the first point and must reproduce all others).  Runtime: about 2 seconds.
"""
import random, sys
import sympy as sp

random.seed(7)
COMP = [(a, b, c) for a in range(3) for b in range(3) for c in range(3) if a + b + c <= 2]
ok_all = True
trials = 0
while trials < 40:
    Ys = random.sample(COMP, 5)
    Y = sp.Matrix([[y[k] for y in Ys] for k in range(3)])
    if sp.Matrix.vstack(Y, sp.ones(1, 5)).rank() < 4:
        continue
    N = sp.Matrix(3, 5, lambda i, j: random.randint(-3, 3))
    if N.rank() < 3:
        continue
    K = N.nullspace()
    r1, r2 = K[0], K[1]
    c = sp.Matrix.vstack(Y, sp.ones(1, 5)).nullspace()[0]
    v = r2 - r1
    vals = []
    for s in [sp.Rational(1, 3), sp.Rational(2, 7), sp.Rational(5, 11), sp.Rational(-3, 4)]:
        mu = (1 - s) * r1 + s * r2
        if any(m == 0 for m in mu):
            continue
        for alpha in [1, sp.Rational(3, 2)]:
            lam = alpha * mu
            lhs = (N * sp.diag(*lam) * Y.T).det()
            gp = sum(c[i] * v[i] / mu[i] for i in range(5))
            core = alpha ** 3 * sp.prod(list(mu)) * gp
            vals.append((lhs, core))
    if len(vals) < 3 or vals[0][1] == 0:
        continue
    const = vals[0][0] / vals[0][1]
    ok = all(sp.simplify(l - const * cst) == 0 for l, cst in vals)
    ok_all = ok_all and ok
    trials += 1
print(("PASS" if ok_all else "FAIL"), "det(N diag(lambda) Y^T) = const * alpha^3 * prod(mu) * g'(s) on", trials,
      "random instances (exact)")
sys.exit(0 if ok_all else 1)
