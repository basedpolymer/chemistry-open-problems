#!/usr/bin/env python3
"""Verifier's own exact checks (sympy) of the analytic steps of result.md / paper.tex, Sections 2-3.

Written independently of scripts/certify_lienard.py: the families are parametrised by the raw reaction-vector
entries (c1, c3, d1, d3 and d2 left symbolic where possible), the Liénard identity is checked on the actual
rational functions A, B, C, D of each family, and the key inequality of Lemma 3.5 is proved through a different
factorisation of P'(a).

V1  normalisation (F3): x = xbar X, y = ybar Y, t = (xbar/mu) T gives X' = sum c_i u_i m_i, Y' = rho sum d_i u_i m_i.
V2  Dulac 1/(xy): divergence formulas for Lemma 2.1(iv), (v).
V3  for c2 = -1 and d2 symbolic: alpha2 = alpha1 + alpha3, beta2 = beta1 + beta3, BC - AD = u1u2u3 x^2 (1-x) (case 9)
    and u1u2u3 x (1 - x^2) (case 10, d2 = -1), D(x0) = -u1 u3 (case 9), N(1) = 1, N'(1).
V4  Liénard identity on the families' A, B, C, D: with W = y - N(x) and f = N' - rho D/B,
    (W + F)' = -B g where g = -rho (BC - AD)/B^2, and xdot/B = W; f(1) = -tr J(1,1)/B(1).
V5  eta substitution: N = (eta+k)^2/((1+k)^2 eta), (g/rho) dx/deta = K (eta+k)^2 (eta-1)/eta^2, K > 0.
V6  H_k primitive; P(a) = H_k(k^2/a) - H_k(a) has P'(a) = -(a+k)^2 ((a-1)^2 + k^2 - 1)/a^3 (< 0 on (0,1) for k >= 1);
    P(1) = H_k(k^2); N(eta_a) - N(eta_b) factorisation.
V7  case 10: g/rho = (u1 u3/u2) h'(x), h = x^2/2 - log x; d/dt [h(1/t) - h(t)] = -(1-t^2)^2/t^3;
    N(a) - N(b) factorisation; Remark 2: F(a) - F(b) = (b-a)[rho - (alpha1 - alpha3/(ab))/alpha2] and
    tr J(1,1) = alpha2 (N'(1) - rho).
V8  exact spot checks of the sign claims on random rational data (illustration only).
Runtime: under a minute. Environment: WSL Fedora, Python 3, sympy.
"""
import random
import sys

import sympy as sp

fails = 0


def report(ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + msg, flush=True)
    if not ok:
        fails += 1


def zero(e):
    e = sp.cancel(sp.together(e))
    if e == 0:
        return True
    return sp.simplify(e) == 0


x, y, X, Y = sp.symbols('x y X Y', positive=True)
rho, mu, xb, yb = sp.symbols('rho mu xbar ybar', positive=True)
c1, c2, c3, d1, d2, d3 = sp.symbols('c1 c2 c3 d1 d2 d3')
SRC = {9: [(2, 0), (1, 1), (0, 1)], 10: [(2, 0), (1, 1), (0, 0)]}


def cross(c, d):
    return (c[1] * d[2] - c[2] * d[1], c[2] * d[0] - c[0] * d[2], c[0] * d[1] - c[1] * d[0])


def m(case, i, X_, Y_):
    a, b = SRC[case][i]
    return X_ ** a * Y_ ** b


# V1 -----------------------------------------------------------------------------------------------
for case in (9, 10):
    c, d = (c1, c2, c3), (d1, d2, d3)
    u = cross(c, d)
    kap = [mu * u[i] / m(case, i, xb, yb) for i in range(3)]
    xdot = sum(c[i] * kap[i] * m(case, i, x, y) for i in range(3))
    ydot = sum(d[i] * kap[i] * m(case, i, x, y) for i in range(3))
    Xp = (xdot / xb * xb / mu).subs({x: xb * X, y: yb * Y})
    Yp = (ydot / yb * xb / mu).subs({x: xb * X, y: yb * Y})
    ok = zero(Xp - sum(c[i] * u[i] * m(case, i, X, Y) for i in range(3))) and \
        zero(Yp - (xb / yb) * sum(d[i] * u[i] * m(case, i, X, Y) for i in range(3)))
    report(ok, f"V1 case {case}: normalised system X' = sum c_i u_i m_i, Y' = rho sum d_i u_i m_i, rho = xbar/ybar")

# V2 -----------------------------------------------------------------------------------------------
k1, k2, k3 = sp.symbols('kappa1 kappa2 kappa3', positive=True)
for case in (9, 10):
    f = sum([c1, c2, c3][i] * [k1, k2, k3][i] * m(case, i, x, y) for i in range(3))
    g = sum([d1, d2, d3][i] * [k1, k2, k3][i] * m(case, i, x, y) for i in range(3))
    div = sp.diff(f / (x * y), x) + sp.diff(g / (x * y), y)
    num = sp.expand(sp.simplify(div * x ** 2 * y ** 2))
    tgt = -d1 * k1 * x ** 3 + c1 * k1 * x ** 2 * y - c3 * k3 * y ** 2 if case == 9 else \
        -d1 * k1 * x ** 3 + c1 * k1 * x ** 2 * y - c3 * k3 * y - d3 * k3 * x
    report(zero(num - tgt), f"V2 case {case}: x^2 y^2 div(F/(xy)) = {tgt}")
# generalised LVA type
f = c1 * k1 * x ** 2 - k2 * x * y
g = d2 * k2 * x * y - k3 * y
report(zero(sp.diff(f / (x * y), x) + sp.diff(g / (x * y), y) - c1 * k1 / y), "V2 LVA type: div(F/(xy)) = c1 kappa1 / y")

# V3, V4 -------------------------------------------------------------------------------------------
fams = {
    '9a': (9, (c1, -1, c3), (d1, -1, d3)),
    '9b': (9, (c1, -1, c3), (d1, 0, -1)),
    '9(c2=-1, d2 free)': (9, (c1, -1, c3), (d1, d2, d3)),
    '10a': (10, (c1, -1, c3), (d1, -1, d3)),
}
data = {}
for name, (case, c, d) in fams.items():
    u1, u2, u3 = [sp.expand(t) for t in cross(c, d)]
    a1, a2, a3 = c[0] * u1, -c[1] * u2, c[2] * u3
    b1, b2, b3 = d[0] * u1, -d[1] * u2, d[2] * u3
    report(zero(a2 - a1 - a3) and zero(b2 - b1 - b3), f"V3 {name}: alpha2 = alpha1 + alpha3, beta2 = beta1 + beta3")
    if case == 9:
        A_, B_, C_, D_ = a1 * x ** 2, a3 - a2 * x, b1 * x ** 2, b3 - b2 * x
        tgt = u1 * u2 * u3 * x ** 2 * (1 - x)
    else:
        A_, B_, C_, D_ = a1 * x ** 2 + a3, -a2 * x, b1 * x ** 2 + b3, -b2 * x
        tgt = u1 * u2 * u3 * x * (1 - x ** 2)
    Xn = sum(c[i] * [u1, u2, u3][i] * m(case, i, x, y) for i in range(3))
    Yn = sum(d[i] * [u1, u2, u3][i] * m(case, i, x, y) for i in range(3))
    report(zero(Xn - A_ - B_ * y) and zero(Yn - C_ - D_ * y), f"V3 {name}: xdot = A + B y, ydot/rho = C + D y")
    report(zero(B_ * C_ - A_ * D_ - tgt), f"V3 {name}: BC - AD = u1u2u3 * " + ("x^2(1-x)" if case == 9 else "x(1-x^2)"))
    N = -A_ / B_
    report(zero(N.subs(x, 1) - 1), f"V3 {name}: N(1) = 1")
    Np1 = sp.diff(N, x).subs(x, 1)
    report(zero(Np1 - ((a1 - a3) / a1 if case == 9 else (a1 - a3) / a2)),
           f"V3 {name}: N'(1) = " + ("(alpha1-alpha3)/alpha1" if case == 9 else "(alpha1-alpha3)/alpha2"))
    if case == 9:
        x0 = a3 / a2
        report(zero(D_.subs(x, x0) + u1 * u3), f"V3 {name}: D(x0) = -u1 u3 with x0 = alpha3/alpha2")
    else:
        report(zero(D_ / B_ - b2 / a2), f"V3 {name}: D/B = beta2/alpha2 (= 1 when d2 = -1)")
    # V4 Liénard identity on the actual functions
    fL = sp.diff(N, x) - rho * D_ / B_
    gL = -rho * (B_ * C_ - A_ * D_) / B_ ** 2
    xdot = A_ + B_ * y
    ydot = rho * (C_ + D_ * y)
    W = y - N
    Wdot = sp.diff(W, x) * xdot + sp.diff(W, y) * ydot
    Ydot = Wdot + fL * xdot           # d/dt (W + F(x)) with F' = f
    report(zero(xdot / B_ - W) and zero(Ydot / B_ + gL),
           f"V4 {name}: xdot/B = Y - F, Ydot/B = -g  (Y = y - N + F, f = N' - rho D/B, g = -rho(BC-AD)/B^2)")
    J = sp.Matrix([xdot, ydot]).jacobian([x, y]).subs({x: 1, y: 1})
    report(zero(fL.subs(x, 1) + J.trace() / B_.subs(x, 1)), f"V4 {name}: f(1) = -tr J(1,1)/B(1)")
    data[name] = dict(case=case, u=(u1, u2, u3), al=(a1, a2, a3), be=(b1, b2, b3), A=A_, B=B_, C=C_, D=D_, N=N)

# V5 ----------------------------------------------------------------------------------------------
eta, k, K = sp.symbols('eta k K', positive=True)
for name in ('9a', '9b'):
    D_ = data[name]
    u1, u2, u3 = D_['u']
    a1, a2, a3 = D_['al']
    kk = a3 / a1
    xe = (eta + kk) / (1 + kk)
    report(zero(xe.subs(eta, 0) - a3 / a2) and zero(xe.subs(eta, 1) - 1),
           f"V5 {name}: x = (eta+k)/(1+k), k = alpha3/alpha1, maps eta = 0 to x0 and eta = 1 to x = 1")
    report(zero(D_['N'].subs(x, xe) - (eta + kk) ** 2 / ((1 + kk) ** 2 * eta)), f"V5 {name}: N = (eta+k)^2/((1+k)^2 eta)")
    gr = (u1 * u2 * u3 * x ** 2 * (x - 1) / D_['B'] ** 2).subs(x, xe) * sp.diff(xe, eta)
    Kval = sp.cancel(sp.together(gr / ((eta + kk) ** 2 * (eta - 1) / eta ** 2)))
    report(zero(sp.diff(Kval, eta)) and zero(Kval - u1 * u2 * u3 / (a1 ** 2 * (1 + kk) ** 4)),
           f"V5 {name}: (g/rho) dx = K (eta+k)^2 (eta-1)/eta^2 deta, K = u1u2u3/(alpha1^2 (1+k)^4) > 0")

# V6 ----------------------------------------------------------------------------------------------
a, b, t = sp.symbols('a b t', positive=True)
Hk = eta ** 2 / 2 + (2 * k - 1) * eta + (k ** 2 - 2 * k) * sp.log(eta) + k ** 2 / eta - (k ** 2 + 2 * k - sp.Rational(1, 2))
report(zero(sp.diff(Hk, eta) - (eta + k) ** 2 * (eta - 1) / eta ** 2) and zero(Hk.subs(eta, 1)),
       "V6 H_k' = (eta+k)^2(eta-1)/eta^2, H_k(1) = 0 (primitive as printed)")
P = Hk.subs(eta, k ** 2 / a) - Hk.subs(eta, a)
report(zero(sp.diff(P, a) + (a + k) ** 2 * ((a - 1) ** 2 + k ** 2 - 1) / a ** 3),
       "V6 P'(a) = -(a+k)^2((a-1)^2 + k^2 - 1)/a^3 < 0 for 0 < a < 1, k >= 1")
report(zero(P.subs(a, 1) - Hk.subs(eta, k ** 2)), "V6 P(1) = H_k(k^2) (>= 0 as k^2 >= 1 and H_k increasing on [1, oo))")
Ne = lambda e: (e + k) ** 2 / ((1 + k) ** 2 * e)
report(zero(Ne(a) - Ne(b) - (a - b) * (1 - k ** 2 / (a * b)) / (1 + k) ** 2), "V6 N(a) - N(b) = (a-b)(1-k^2/(ab))/(1+k)^2")

# V7 ----------------------------------------------------------------------------------------------
D_ = data['10a']
u1, u2, u3 = D_['u']
a1, a2, a3 = D_['al']
h = x ** 2 / 2 - sp.log(x)
gr = u1 * u2 * u3 * x * (x ** 2 - 1) / D_['B'] ** 2
report(zero(gr - (u1 * u3 / u2) * sp.diff(h, x)), "V7 10a: g/rho = (u1u3/u2) h'(x), h = x^2/2 - log x")
ht = h.subs(x, 1 / t) - h.subs(x, t)
report(zero(sp.diff(ht, t) + (1 - t ** 2) ** 2 / t ** 3) and zero(ht.subs(t, 1)),
       "V7 h(1/t) - h(t): zero at t = 1, derivative -(1-t^2)^2/t^3")
Nx = D_['N']
report(zero(Nx.subs(x, a) - Nx.subs(x, b) - (a - b) * (a1 - a3 / (a * b)) / a2), "V7 N(a) - N(b) = (a-b)(alpha1 - alpha3/(ab))/alpha2")
# Remark 2 (any case-10 network with c2 = d2 = -1): psi = x - 1, F = N - 1 - rho (x - 1)
Fx = Nx - 1 - rho * (x - 1)
report(zero(Fx.subs(x, a) - Fx.subs(x, b) - (b - a) * (rho - (a1 - a3 / (a * b)) / a2)),
       "V7 Remark 2: F(a) - F(b) = (b-a)[rho - (alpha1 - alpha3/(ab))/alpha2]")
xdot = D_['A'] + D_['B'] * y
ydot = rho * (D_['C'] + D_['D'] * y)
J = sp.Matrix([xdot, ydot]).jacobian([x, y]).subs({x: 1, y: 1})
report(zero(J.trace() - a2 * (sp.diff(Nx, x).subs(x, 1) - rho)), "V7 Remark 2: tr J(1,1) = alpha2 (N'(1) - rho)")

# V8 ----------------------------------------------------------------------------------------------
random.seed(20261025)
ok = True
for _ in range(400):
    kv = 1 + sp.Rational(random.randint(0, 3000), random.randint(1, 300))
    av = sp.Rational(random.randint(1, 9999), 10000)
    ok &= bool(sp.N(P.subs({k: kv, a: av}), 50) > 0)
report(ok, "V8 P(a) > 0 at 400 random rational (k >= 1, 0 < a < 1), 50 digits (illustration)")
# k = 1 limit: P(a) -> 0 only at a = 1
report(zero(P.subs({k: 1, a: 1})), "V8 k = 1: P(1) = 0 (the inequality is tight only at a = 1)")

print("ALL PASS" if fails == 0 else f"{fails} FAIL")
sys.exit(1 if fails else 0)
