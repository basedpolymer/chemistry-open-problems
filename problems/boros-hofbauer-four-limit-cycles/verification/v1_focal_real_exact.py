"""
v1_focal_real_exact.py -- VERIFIER's independent exact computation of the focal values of the
Boros-Hofbauer system on the trace-zero hypersurface, by a DIFFERENT algorithm from the solver's:
real-coordinate Lyapunov-function method (homogeneous linear algebra degree by degree, real monomial basis
xi^i eta^(n-i), explicit circle-average normalisation), with the vector field derived independently below.

Written by the verifier (Claude Opus 5.5), without reading the solver's notes.
Environment: Python 3.12, python-flint 0.9.0 (fmpq_mpoly, exact), sympy 1.14 (exact matrix inverses and a
symbolic identity). Runtime: a few seconds.

Derivation (verifier's own).
  BH: x' = 1/y - x^a y^b, y' = K(-(d-1) - 1/y + d x^a y^b), K = a/(1+bd) (tr J = 0).
  Put c = 1+b, s = (d-1)/(1+bd) = w^2 (w > 0), u = x-1, v = y-1, I = 1/(1+v), E = (1+u)^a (1+v)^b.
  One checks symbolically (sympy, below):  a u' + c v' = a s [ (c-1) I + E - c ].
  Linear change (u, v) = D * (c xi, a (w eta - xi)) with D = a c w, time tau = omega t, omega = a w:
      xi'  = u' / D^2,        eta' = w [ (c-1) I + E - c ] / D^2.
  In the unscaled variables (xi0, eta0) = D (xi, eta) the degree-n homogeneous part of any series picks up
  D^n, so the degree-n coefficient of the new field is D^(n-2) * (degree-n part computed with
  u = c xi, v = a (w eta - xi)) -- a polynomial in (a, c, w) for n >= 2, and the linear part is (-eta, xi)
  (checked below).
  Lyapunov function F = xi^2 + eta^2 + sum_{n>=3} F_n, F_n homogeneous; degree by degree solve
      Lambda F_n + T_n = V * (xi^2+eta^2)^(n/2)   (V = 0 if n odd),  Lambda = -eta d/dxi + xi d/deta,
      T_n = sum_{m=2}^{n-1} [ dF_m/dxi * P_{n-m+1} + dF_m/deta * Q_{n-m+1} ],
  with the normalisation  mean_{theta} F_n(cos theta, sin theta) = 0  for even n.
  Then dF/dtau = sum_k Vhat_k (xi^2+eta^2)^(k+1) + O(|.|^(2j+3)).
  This normalisation is the same as "coefficient of (z zbar)^(n/2) equal to 0" in complex notation, so the
  solver's V~_k should coincide EXACTLY with Vhat_k (and V~_k = w G_k(a, c, w^2)).

Checks (PASS/FAIL):
  E1 the symbolic identity a u' + c v' = a s [(c-1) I + E - c] (sympy);
  E2 linear part of the transformed field is (-eta, xi);
  E3 Vhat_k == w * G_k(a, c, w^2) for k = 1..4, G_k read from ../certificates/G_k.json (exact identity);
  E4 G_1 equals the closed form of L_1 printed by Boros-Hofbauer (sympy, exact).
"""
import os
import sys
import json
from fractions import Fraction
import sympy as sp
from flint import fmpq_mpoly_ctx, fmpq, Ordering

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "..", "certificates")
CTX = fmpq_mpoly_ctx.get(('a', 'c', 'w'), Ordering.lex)
A, C, W = CTX.gens()
ONE = CTX.from_dict({(0, 0, 0): 1})
ZERO = CTX.from_dict({})
results = []


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


# ---------------------------------------------------------------- E1: symbolic identity
def check_identity():
    a, b, d, I, E = sp.symbols('a b d I E')
    K = a / (1 + b * d)
    s = (d - 1) / (1 + b * d)
    c = 1 + b
    udot = I - E
    vdot = K * (-(d - 1) - I + d * E)
    expr = sp.simplify(a * udot + c * vdot - a * s * ((c - 1) * I + E - c))
    report("E1 a u' + c v' == a s [(c-1) I + E - c] on tr J = 0 (sympy)", expr == 0)
    # also: tr J = 0 for this K, det J = a^2 s (= omega^2)
    x, y = sp.symbols('x y', positive=True)
    f = 1 / y - x**a * y**b
    g = K * (-(d - 1) - 1 / y + d * x**a * y**b)
    J = sp.Matrix([[sp.diff(f, x), sp.diff(f, y)], [sp.diff(g, x), sp.diff(g, y)]]).subs({x: 1, y: 1})
    report("E1' tr J(1,1) = 0 and det J(1,1) = a^2 s", sp.simplify(J.trace()) == 0
           and sp.simplify(J.det() - a**2 * s) == 0)


# ---------------------------------------------------------------- bivariate series in (xi, eta)
# series: dict {(i, j): coeff} meaning coeff * xi^i * eta^j, truncated at total degree NT
NT = 9


def smul(P, Q):
    out = {}
    for (i1, j1), c1 in P.items():
        for (i2, j2), c2 in Q.items():
            if i1 + j1 + i2 + j2 <= NT:
                k = (i1 + i2, j1 + j2)
                out[k] = out[k] + c1 * c2 if k in out else c1 * c2
    return out


def sadd(P, Q, cq=1):
    out = dict(P)
    for k, v in Q.items():
        out[k] = out[k] + v * cq if k in out else v * cq
    return out


def binom_poly(alpha, i):
    r = ONE
    for k in range(i):
        r = r * (alpha - k) * fmpq(1, k + 1)
    return r


def build_field():
    u = {(1, 0): C}
    v = {(1, 0): -A, (0, 1): A * W}
    upow = [{(0, 0): ONE}]
    vpow = [{(0, 0): ONE}]
    for _ in range(NT):
        upow.append(smul(upow[-1], u))
        vpow.append(smul(vpow[-1], v))
    I = {}
    for j in range(NT + 1):
        I = sadd(I, vpow[j], (-1) ** j)
    Ua = {}
    for i in range(NT + 1):
        Ua = sadd(Ua, upow[i], binom_poly(A, i))
    Vb = {}
    for j in range(NT + 1):
        Vb = sadd(Vb, vpow[j], binom_poly(C - 1, j))
    E = smul(Ua, Vb)
    udot = sadd(I, E, -1)
    Wf = {k: val * (C - 1) for k, val in I.items()}
    Wf = sadd(Wf, E)
    Wf = sadd(Wf, {(0, 0): C}, -1)
    D = A * C * W
    # degree-n homogeneous parts, rescaled by D^(n-2)
    P = {}
    Q = {}
    lin = {}
    for (i, j), val in udot.items():
        n = i + j
        if n == 0:
            if val != 0:
                raise RuntimeError("udot has a constant term")
            continue
        if n == 1:
            lin[('P', i, j)] = val
            continue
        P[(i, j)] = val * D ** (n - 2)
    for (i, j), val in Wf.items():
        n = i + j
        if n == 0:
            if val != 0:
                raise RuntimeError("W has a constant term")
            continue
        if n == 1:
            lin[('Q', i, j)] = val * W
            continue
        Q[(i, j)] = val * W * D ** (n - 2)
    # linear part: xi' = lin_P / D, eta' = lin_Q / D ; must be (-eta, xi)
    okP = (lin.get(('P', 1, 0), ZERO) == 0) and (lin.get(('P', 0, 1), ZERO) == -D)
    okQ = (lin.get(('Q', 1, 0), ZERO) == D) and (lin.get(('Q', 0, 1), ZERO) == 0)
    report("E2 linear part of the transformed field is (-eta, xi)", okP and okQ)
    return P, Q


# ---------------------------------------------------------------- homogeneous linear algebra
def hom_parts(S):
    """series dict -> {n: list of coefficients h[i] of xi^i eta^(n-i)}"""
    out = {}
    for (i, j), val in S.items():
        n = i + j
        out.setdefault(n, [ZERO] * (n + 1))
        out[n][i] = out[n][i] + val
    return out


def hmul(p, q):
    n1, n2 = len(p) - 1, len(q) - 1
    out = [ZERO] * (n1 + n2 + 1)
    for i1, c1 in enumerate(p):
        if c1 == 0:
            continue
        for i2, c2 in enumerate(q):
            if c2 == 0:
                continue
            out[i1 + i2] = out[i1 + i2] + c1 * c2
    return out


def dxi(p):
    n = len(p) - 1
    return [p[i] * i for i in range(1, n + 1)] if n >= 1 else [ZERO]


def deta(p):
    n = len(p) - 1
    return [p[i] * (n - i) for i in range(0, n)] if n >= 1 else [ZERO]


def circle_mean(i, j):
    """(1/2pi) int_0^{2pi} cos^i sin^j: exact rational."""
    if i % 2 or j % 2:
        return Fraction(0)
    p, q = i // 2, j // 2
    from math import factorial
    return Fraction(factorial(2 * p) * factorial(2 * q), 4 ** (p + q) * factorial(p) * factorial(q) * factorial(p + q))


def solver_matrix(n):
    """Inverse of the matrix of the linear system for degree n (exact, sympy)."""
    if n % 2 == 1:
        M = sp.zeros(n + 1, n + 1)
        for i in range(n + 1):
            if i >= 1:
                M[i - 1, i] += -i
            if i + 1 <= n:
                M[i + 1, i] += (n - i)
        return M.inv()
    M = sp.zeros(n + 2, n + 2)
    for i in range(n + 1):
        if i >= 1:
            M[i - 1, i] += -i
        if i + 1 <= n:
            M[i + 1, i] += (n - i)
    k = n // 2
    for i in range(0, n + 1, 2):        # (xi^2 + eta^2)^k = sum binom(k, i/2) xi^i eta^(n-i)
        M[i, n + 1] = -sp.binomial(k, i // 2)
    for i in range(n + 1):
        cm = circle_mean(i, n - i)
        M[n + 1, i] = sp.Rational(cm.numerator, cm.denominator)
    return M.inv()


def lyapunov(P, Q, nmax=10):
    Ph, Qh = hom_parts(P), hom_parts(Q)
    F = {2: [ONE, ZERO, ONE]}             # xi^2 + eta^2 : coefficients of eta^2, xi eta, xi^2
    Vs = []
    for n in range(3, nmax + 1):
        T = [ZERO] * (n + 1)
        for m in range(2, n):
            k = n - m + 1
            if k not in Ph:
                continue
            t1 = hmul(dxi(F[m]), Ph[k])
            t2 = hmul(deta(F[m]), Qh[k])
            for i in range(n + 1):
                T[i] = T[i] + t1[i] + t2[i]
        Minv = solver_matrix(n)
        rhs = [-t for t in T] + ([ZERO] if n % 2 == 0 else [])
        sol = []
        for r in range(Minv.rows):
            acc = ZERO
            for cidx in range(Minv.cols):
                e = Minv[r, cidx]
                if e != 0:
                    acc = acc + rhs[cidx] * fmpq(int(e.p), int(e.q))
            sol.append(acc)
        F[n] = sol[:n + 1]
        if n % 2 == 0:
            Vs.append(sol[n + 1])
            # sanity: residual check  Lambda F_n + T_n - V e_n == 0
            kk = n // 2
            LF = [ZERO] * (n + 1)
            for i in range(n + 1):
                if i >= 1:
                    LF[i - 1] = LF[i - 1] - F[n][i] * i
                if i + 1 <= n:
                    LF[i + 1] = LF[i + 1] + F[n][i] * (n - i)
            for i in range(n + 1):
                e = int(sp.binomial(kk, i // 2)) if i % 2 == 0 else 0
                if LF[i] + T[i] - sol[n + 1] * e != 0:
                    raise RuntimeError("residual check failed at degree %d" % n)
    return Vs


def load_G(k):
    with open(os.path.join(CERT, "G_%d.json" % k)) as fh:
        data = json.load(fh)
    d = {}
    for ea, ec, es, cf in data["terms"]:
        num, _, den = cf.partition('/')
        d[(ea, ec, 2 * es + 1)] = fmpq(int(num), int(den) if den else 1)   # w * G(a, c, w^2)
    return CTX.from_dict(d)


def main():
    check_identity()
    P, Q = build_field()
    Vs = lyapunov(P, Q, 10)
    for k, V in enumerate(Vs, 1):
        G = load_G(k)
        report("E3 Vhat_%d (verifier, real method) == w*G_%d(a,c,w^2) (solver's certificate)" % (k, k), V == G,
               "%d terms" % len(V))
    # E4: closed form of the paper for L_1
    a, b, d = sp.symbols('a b d', positive=True)
    c = 1 + b
    s = (d - 1) / (1 + b * d)
    with open(os.path.join(CERT, "G_1.json")) as fh:
        g1 = json.load(fh)
    G1 = sum(sp.Rational(t[3]) * a**t[0] * c**t[1] * s**t[2] for t in g1["terms"])
    # L_1 = pi G_1/(a^2 c^4 s^(1/2)); compare squares to avoid sqrt branch issues (both have the same sign
    # factor check separately)
    L1_ours_over_pi = G1 / (a**2 * c**4 * sp.sqrt(s))
    L1_paper_over_pi = a * (2 * a**2 * d + a * (1 + b * d) - (1 + b * d)**2) / (8 * a * sp.sqrt((d - 1) * (1 + b * d)**3))
    ratio = sp.simplify(sp.powsimp(sp.radsimp(L1_ours_over_pi / L1_paper_over_pi), force=True))
    # numeric confirmation at random points as well
    import random
    random.seed(1)
    worst = 0
    for _ in range(20):
        av = sp.Rational(random.randint(1, 300), 100)
        bv = sp.Rational(random.randint(-90, 300), 100)
        dv = sp.Rational(random.randint(101, 600), 100)
        if 1 + bv * dv <= 0:
            continue
        x1 = L1_ours_over_pi.subs({a: av, b: bv, d: dv})
        x2 = L1_paper_over_pi.subs({a: av, b: bv, d: dv})
        worst = max(worst, abs(sp.N(x1 - x2, 50)))
    report("E4 G_1 gives exactly the paper's closed form of L_1", ratio == 1 and worst < 1e-45,
           "simplified ratio = %s, max |diff| at 20 random points = %s" % (ratio, sp.N(worst, 3)))
    print("ALL PASS" if all(results) else "SOME CHECKS FAILED")
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
