#!/usr/bin/env python3
"""Exact symbolic certificate (sympy) for every identity and sign claim used in result.md.

Environment: Python 3.12/3.13 + sympy 1.14. Single thread. Runtime: < 1 minute.
Prints PASS/FAIL for each claim and exits with status 1 on any failure.

Notation (result.md): reaction i has source (a_i, b_i) and reaction vector (c_i, d_i); u = c x d;
case 9 monomials x^2, xy, y; case 10 monomials x^2, xy, 1.
Sections:
  S1  normalisation: eq. at (1,1), the free parameter is rho = xbar/ybar; positive equilibrium is unique
  S2  the general Liénard identity for xdot = A + B y, ydot = rho (C + D y)
  S3  the three families 9a, 9b, 10a: A, B, C, D, BC - AD = u1u2u3 x^2(1-x) (case 9) / x(1-x^2) (case 10),
      sign of D/B, f(1) < 0, equivalence k >= 1 <=> BBH's inequality fails (no Andronov-Hopf)
  S4  case 9 reduction to the eta variable: N and Gtilde in terms of k only
  S5  the key inequality (Lemma 3.3): H_k(k^2/a) > H_k(a) for 0 < a < 1, k >= 1
  S6  case 10: Gtilde ~ x^2/2 - log x, and the inequality h(1/a) > h(a)
  S7  the remaining classes: Dulac 1/(xy) (c1 <= 0, LVA type), boundary behaviour
"""
import sys
import random
import sympy as sp

fails = 0


def check(cond, msg):
    global fails
    cond = bool(cond)
    print(("PASS " if cond else "FAIL ") + msg, flush=True)
    fails += (not cond)


def is_zero(e):
    return sp.simplify(sp.expand(e)) == 0


def all_pos_coeffs(expr, gens):
    """expr is a polynomial in gens with only positive coefficients (and is not 0)."""
    P = sp.Poly(sp.expand(expr), *gens)
    return P.is_zero is False and all(c > 0 for c in P.coeffs())


x, y, X, Y, t = sp.symbols("x y X Y t", positive=True)
rho, mu, xb, yb = sp.symbols("rho mu xbar ybar", positive=True)
k1, k2, k3 = sp.symbols("kappa1 kappa2 kappa3", positive=True)
C1, C2, C3, D1, D2, D3 = sp.symbols("c1 c2 c3 d1 d2 d3")
SRC = {9: [(2, 0), (1, 1), (0, 1)], 10: [(2, 0), (1, 1), (0, 0)]}


def mon(src, X, Y):
    return X ** src[0] * Y ** src[1]


def mass_action(case, c, d, kap, X, Y):
    m = [mon(s, X, Y) for s in SRC[case]]
    return (sum(c[i] * kap[i] * m[i] for i in range(3)), sum(d[i] * kap[i] * m[i] for i in range(3)))


def uvec(c, d):
    return (c[1] * d[2] - c[2] * d[1], c[2] * d[0] - c[0] * d[2], c[0] * d[1] - c[1] * d[0])


# ---------------------------------------------------------------- S1 normalisation
print("== S1 normalisation")
c = (C1, C2, C3); d = (D1, D2, D3)
u = uvec(c, d)
for case in (9, 10):
    # rate constants in terms of the equilibrium (xb, yb) and mu: kappa_i xb^a_i yb^b_i = mu u_i
    kap = [mu * u[i] / mon(SRC[case][i], xb, yb) for i in range(3)]
    fx, fy = mass_action(case, c, d, kap, x, y)
    # x = xb X, y = yb Y, t = (xb/mu) T:  dX/dT = (xb/mu) (1/xb) xdot,  dY/dT = (xb/mu)(1/yb) ydot
    FX = sp.simplify((xb / mu) / xb * fx.subs({x: xb * X, y: yb * Y}))
    FY = sp.simplify((xb / mu) / yb * fy.subs({x: xb * X, y: yb * Y}))
    m = [mon(s, X, Y) for s in SRC[case]]
    NX = sum(c[i] * u[i] * m[i] for i in range(3))
    NY = (xb / yb) * sum(d[i] * u[i] * m[i] for i in range(3))
    check(is_zero(FX - NX) and is_zero(FY - NY),
          f"S1 case {case}: normalised system is X' = sum c_i u_i M_i, Y' = (xbar/ybar) sum d_i u_i M_i")
    check(is_zero(NX.subs({X: 1, Y: 1}) - sum(c[i] * u[i] for i in range(3))) and
          is_zero(sum(c[i] * u[i] for i in range(3))) and is_zero(sum(d[i] * u[i] for i in range(3))),
          f"S1 case {case}: (1,1) is an equilibrium (sum c_i u_i = sum d_i u_i = 0 identically)")
# Jacobian of the normalised system at (1,1): det and trace
for case in (9, 10):
    m = [mon(s_, X, Y) for s_ in SRC[case]]
    NX = sum(c[i] * u[i] * m[i] for i in range(3))
    NY = rho * sum(d[i] * u[i] * m[i] for i in range(3))
    J = sp.Matrix([NX, NY]).jacobian([X, Y]).subs({X: 1, Y: 1})
    fac = 1 if case == 9 else 2
    Abbh = C1 * u[0] - C3 * u[2]
    trt = Abbh - rho * D1 * u[0] if case == 9 else Abbh + rho * D2 * u[1]
    check(is_zero(J.det() - fac * rho * u[0] * u[1] * u[2]) and is_zero(J.trace() - trt),
          f"S1 case {case}: at (1,1) det J = {fac} rho u1u2u3 and tr J = (c1u1 - c3u3) "
          + ("- rho d1 u1" if case == 9 else "+ rho d2 u2"))
# uniqueness of the positive equilibrium (u > 0): explicit solution of kappa_i M_i(x,y) = mu u_i
U1, U2, U3 = sp.symbols("u1 u2 u3", positive=True)
xs9 = k3 * U2 / (k2 * U3); ys9 = k1 * k3 * U2 ** 2 / (k2 ** 2 * U1 * U3); mus9 = k3 * ys9 / U3
check(is_zero(k1 * xs9 ** 2 - mus9 * U1) and is_zero(k2 * xs9 * ys9 - mus9 * U2) and is_zero(k3 * ys9 - mus9 * U3),
      "S1 case 9: kappa1 x^2 = mu u1, kappa2 x y = mu u2, kappa3 y = mu u3 solved by x = k3u2/(k2u3), "
      "y = k1k3u2^2/(k2^2u1u3) (unique: x is forced by eqs 2/3, then y by eqs 1/2)")
mus10 = k3 / U3; xs10 = sp.sqrt(mus10 * U1 / k1); ys10 = mus10 * U2 / (k2 * xs10)
check(is_zero(k1 * xs10 ** 2 - mus10 * U1) and is_zero(k2 * xs10 * ys10 - mus10 * U2),
      "S1 case 10: mu = k3/u3, x = sqrt(mu u1/k1), y = mu u2/(k2 x) (unique positive solution)")

# ---------------------------------------------------------------- S2 general Liénard identity
print("== S2 Liénard identity")
Af, Bf, Cf, Df, Ff = [sp.Function(n)(x) for n in ("A", "B", "C", "D", "F")]
xdot = Af + Bf * y
ydot = rho * (Cf + Df * y)
N = -Af / Bf
f = sp.diff(N, x) - rho * Df / Bf
g = -rho * (Bf * Cf - Af * Df) / Bf ** 2
Yl = y - N + Ff                       # Liénard variable Y = w + F(x), w = y - N(x)
Ydot = sp.diff(Yl, x) * xdot + sp.diff(Yl, y) * ydot
Ydot = Ydot.subs(sp.Derivative(Ff, x), f)
check(is_zero(xdot / Bf - (Yl - Ff)), "S2 (1/B) xdot = Y - F(x)  with Y = y - N(x) + F(x), N = -A/B")
check(is_zero(sp.together(Ydot / Bf + g)), "S2 (1/B) dY/dt = -g(x)  with F' = f = N' - rho D/B, g = -rho (BC-AD)/B^2")
# trace relation f(1) = -tr J(1,1) / B(1) when (1,1) is an equilibrium (A(1) = -B(1))
a0, a1, b0, b1, cc0, dd0 = sp.symbols("A1 Ap1 B1 Bp1 C1v D1v")
# tr J at (1,1): d/dx (A + B y) + d/dy rho(C + D y) = A'(1) + B'(1) + rho D(1); N'(1) = -(A'B - AB')/B^2 at 1
Np1 = -(a1 * b0 - a0 * b1) / b0 ** 2
check(is_zero((Np1 - rho * dd0 / b0).subs(a0, -b0) + (a1 + b1 + rho * dd0) / b0),
      "S2 f(1) = -trJ(1,1)/B(1) (uses A(1) = -B(1), i.e. (1,1) on the x-nullcline)")

# ---------------------------------------------------------------- S3 the three families
print("== S3 families 9a, 9b, 10a")
c1, c3, p, q = sp.symbols("c1 c3 p q", positive=True)   # integers >= 1 in the networks
dd1 = sp.Symbol("d1", positive=True)
fam = {
    # name: (case, c, d, substitution description)
    "9a": (9, (c1, -1, c3), (c1 + p, -1, c3 - q)),    # d1 = c1 + p (p>=1), d3 = c3 - q (q>=1)
    "9b": (9, (c1, -1, c3), (dd1, 0, -1)),
    "10a": (10, (c1, -1, c3), (c1 + p, -1, c3 - q)),  # additionally d3 = c3 - q >= 0
}
store = {}
for name, (case, cv, dv) in fam.items():
    u1, u2, u3 = [sp.expand(e) for e in uvec(cv, dv)]
    gens = (c1, c3, p, q) if name != "9b" else (c1, c3, dd1)
    check(all_pos_coeffs(u1, gens) and all_pos_coeffs(u2, gens) and all_pos_coeffs(u3, gens),
          f"S3 {name}: u = ({u1}, {u2}, {u3}) > 0 (polynomials with positive coefficients)")
    m = [mon(s, x, y) for s in SRC[case]]
    fx = sp.expand(sum(cv[i] * [u1, u2, u3][i] * m[i] for i in range(3)))
    fy = sp.expand(sum(dv[i] * [u1, u2, u3][i] * m[i] for i in range(3)))
    A = fx.subs(y, 0); B = sp.diff(fx, y); C = fy.subs(y, 0); D = sp.diff(fy, y)
    check(is_zero(A + B * y - fx) and is_zero(C + D * y - fy), f"S3 {name}: system is linear in y: A={A}, B={B}")
    al1, al2, al3 = cv[0] * u1, u2, cv[2] * u3
    check(is_zero(al2 - al1 - al3), f"S3 {name}: alpha2 = alpha1 + alpha3 (alpha1 = c1u1, alpha2 = u2, alpha3 = c3u3)")
    kk = sp.simplify(al3 / al1)
    Abbh = sp.expand(cv[0] * u1 - cv[2] * u3)
    # BBH inequality (1/2)(d3/c3 + d1/c1) < d2/c2  <=>  c1 d3 + c3 d1 < -2 c1 c3 d2  (c2 = -1)
    bbh_gap = sp.expand(-2 * cv[0] * cv[2] * dv[1] - cv[0] * dv[2] - cv[2] * dv[0])  # > 0 iff left BBH inequality holds
    check(is_zero(Abbh - bbh_gap) and is_zero(Abbh - (al1 - al3)),
          f"S3 {name}: c1u1 - c3u3 = alpha1 - alpha3 = 2c1c3(-d2) - c1d3 - c3d1 (so BBH-Hopf fails <=> k = alpha3/alpha1 >= 1)")
    BCAD = sp.factor(B * C - A * D)
    u123 = u1 * u2 * u3
    target = u123 * x ** 2 * (1 - x) if case == 9 else u123 * x * (1 - x ** 2)
    check(is_zero(BCAD - target), f"S3 {name}: BC - AD = u1u2u3 * {'x^2(1-x)' if case == 9 else 'x(1-x^2)'}")
    # sign of D/B on the Liénard domain I
    if case == 9:
        x0 = al3 / al2
        check(is_zero(B - (al3 - al2 * x)), f"S3 {name}: B = alpha3 - alpha2 x, negative exactly for x > x0 = alpha3/alpha2")
        check(is_zero(A.subs(x, x0) - al1 * x0 ** 2), f"S3 {name}: on x = x0, xdot = A(x0) = alpha1 x0^2 > 0")
        Dx0 = sp.factor(D.subs(x, x0))
        slope = sp.diff(D, x)
        check(sp.simplify(-Dx0).is_positive and (is_zero(slope) or sp.simplify(-slope).is_positive),
              f"S3 {name}: D(x0) = {Dx0} < 0 and D' = {slope} <= 0, hence D < 0 < -B on (x0, oo): D/B > 0")
    else:
        check(is_zero(B + u2 * x) and is_zero(D / B - 1), f"S3 {name}: B = -u2 x < 0 on (0, oo) and D/B = 1")
        check(all_pos_coeffs(u3, gens) and all_pos_coeffs(sp.expand(cv[2] * u3), gens),
              f"S3 {name}: on x = 0, xdot = alpha3 = c3 u3 > 0")
    Nf = sp.simplify(-A / B)
    Np = sp.simplify(sp.diff(Nf, x).subs(x, 1))
    target_Np = (al1 - al3) / al1 if case == 9 else (al1 - al3) / al2
    check(is_zero(Np - target_Np), f"S3 {name}: N'(1) = {'(alpha1-alpha3)/alpha1 = 1-k' if case==9 else '(alpha1-alpha3)/alpha2'} <= 0 iff k >= 1")
    DB1 = sp.simplify((D / B).subs(x, 1))
    check(sp.simplify(DB1).is_positive, f"S3 {name}: D(1)/B(1) = {sp.factor(DB1)} > 0, so f(1) = N'(1) - rho D(1)/B(1) < 0 when k >= 1")
    store[name] = dict(A=A, B=B, C=C, D=D, N=Nf, al=(al1, al2, al3), u=(u1, u2, u3), k=kk, case=case)

# ---------------------------------------------------------------- S4 case 9: eta variable
print("== S4 case 9 in the variable eta")
eta, k, s = sp.symbols("eta k s", positive=True)
for name in ("9a", "9b"):
    S = store[name]
    al1, al2, al3 = S["al"]; u1, u2, u3 = S["u"]
    kk = S["k"]
    xeta = (eta + kk) / (1 + kk)
    check(is_zero(xeta.subs(eta, 0) - al3 / al2) and is_zero(xeta.subs(eta, 1) - 1) and is_zero(sp.diff(xeta, eta) - 1 / (1 + kk)),
          f"S4 {name}: x = (eta+k)/(1+k) maps eta in (0,oo) increasingly onto (x0,oo), eta = 1 <-> x = 1")
    Neta = sp.simplify(S["N"].subs(x, xeta))
    check(is_zero(Neta - (eta + kk) ** 2 / ((1 + kk) ** 2 * eta)), f"S4 {name}: N(x(eta)) = (eta+k)^2 / ((1+k)^2 eta)")
    gt = u1 * u2 * u3 * x ** 2 * (x - 1) / S["B"] ** 2      # g / rho
    lhs = sp.simplify(gt.subs(x, xeta) * sp.diff(xeta, eta))
    const = sp.simplify(lhs / ((eta + kk) ** 2 * (eta - 1) / eta ** 2))
    check(sp.simplify(sp.diff(const, eta)) == 0 and sp.simplify(const).is_positive,
          f"S4 {name}: (g/rho) dx = K (eta+k)^2 (eta-1)/eta^2 d eta with K = {sp.factor(const)} > 0")
    gens4 = (c1, c3, p, q) if name == "9a" else (c1, c3, dd1)
    check(is_zero(gt * S["B"] ** 2 - u1 * u2 * u3 * x ** 2 * (x - 1)) and all_pos_coeffs(sp.expand(u1 * u2 * u3), gens4),
          f"S4 {name}: g/rho = u1u2u3 x^2 (x-1)/B^2 with u1u2u3 > 0, so (x-1) g(x) > 0 on I minus {{1}}")

# ---------------------------------------------------------------- S5 the key inequality
print("== S5 Lemma: H_k(k^2/a) > H_k(a) for 0<a<1, k>=1")
a = sp.Symbol("a", positive=True)
Hk = eta ** 2 / 2 + (2 * k - 1) * eta + (k ** 2 - 2 * k) * sp.log(eta) + k ** 2 / eta - (k ** 2 + 2 * k - sp.Rational(1, 2))
check(is_zero(sp.diff(Hk, eta) - (eta + k) ** 2 * (eta - 1) / eta ** 2) and is_zero(Hk.subs(eta, 1)),
      "S5 H_k(eta) := eta^2/2 + (2k-1)eta + (k^2-2k)log eta + k^2/eta - (k^2+2k-1/2) satisfies H_k(1)=0, H_k' = (eta+k)^2(eta-1)/eta^2")
P = Hk.subs(eta, k ** 2 / a) - Hk.subs(eta, a)
Pp = sp.simplify(sp.diff(P, a))
Q = a ** 4 + (2 * k - 2) * a ** 3 + 2 * k * (k - 2) * a ** 2 + (2 * k - 2) * k ** 2 * a + k ** 4
check(is_zero(Pp + Q / a ** 3), "S5 P(a) := H_k(k^2/a) - H_k(a) has P'(a) = -Q(a)/a^3, Q = a^4+(2k-2)a^3+2k(k-2)a^2+(2k-2)k^2 a+k^4")
sv = (a - k) ** 2 / a
check(is_zero(Q / a ** 2 - (8 * k * (k - 1) + (6 * k - 2) * sv + sv ** 2)),
      "S5 Q(a)/a^2 = 8k(k-1) + (6k-2)s + s^2 with s = (a-k)^2/a >= 0; hence Q > 0 for k >= 1 and a != k")
check(is_zero(P.subs(a, 1) - Hk.subs(eta, k ** 2)), "S5 P(1) = H_k(k^2) >= 0 (H_k >= 0 on [1,oo) since H_k' >= 0 there)")
# exact sanity: random rational k >= 1, a in (0,1): P(a) > 0 evaluated with 60-digit arithmetic
random.seed(7)
okP = True
for _ in range(300):
    kv = 1 + sp.Rational(random.randint(0, 4000), random.randint(1, 400))
    av = sp.Rational(random.randint(1, 9999), 10000)
    val = sp.N(P.subs({k: kv, a: av}), 60)
    okP &= bool(val > 0)
check(okP, "S5 spot check: P(a) > 0 at 300 random rational (k >= 1, 0 < a < 1) (60-digit evaluation; illustration only)")
# N(a) - N(b) in eta: (a-b)(1 - k^2/(ab))/(1+k)^2
bb = sp.Symbol("b", positive=True)
Ne = lambda e: (e + k) ** 2 / ((1 + k) ** 2 * e)
check(is_zero(Ne(a) - Ne(bb) - (a - bb) * (1 - k ** 2 / (a * bb)) / (1 + k) ** 2),
      "S5 N(a) - N(b) = (a-b)(1 - k^2/(ab))/(1+k)^2 (eta variable); > 0 when a < b and ab < k^2")

# ---------------------------------------------------------------- S6 case 10
print("== S6 case 10")
S = store["10a"]
al1, al2, al3 = S["al"]; u1, u2, u3 = S["u"]
gt = u1 * u2 * u3 * x * (x ** 2 - 1) / S["B"] ** 2
hx = x ** 2 / 2 - sp.log(x)
check(is_zero(gt - (u1 * u3 / u2) * sp.diff(hx, x)), "S6 g/rho = (u1u3/u2) h'(x) with h(x) = x^2/2 - log x; (x-1)g > 0 on (0,oo) minus {1}")
check(is_zero(S["N"] - (al1 * x ** 2 + al3) / (al2 * x)), "S6 N(x) = (alpha1 x^2 + alpha3)/(alpha2 x)")
Nc = lambda e: (al1 * e ** 2 + al3) / (al2 * e)
check(is_zero(Nc(a) - Nc(bb) - (a - bb) * (al1 - al3 / (a * bb)) / al2),
      "S6 N(a) - N(b) = (a-b)(alpha1 - alpha3/(ab))/alpha2; > 0 when a < b and ab < k = alpha3/alpha1")
tt = sp.Symbol("t", positive=True)
kt = hx.subs(x, 1 / tt) - hx.subs(x, tt)
check(is_zero(sp.diff(kt, tt) + (1 - tt ** 2) ** 2 / tt ** 3) and is_zero(kt.subs(tt, 1)),
      "S6 h(1/t) - h(t) vanishes at t = 1 and has derivative -(1-t^2)^2/t^3 < 0: so h(1/a) > h(a) for 0 < a < 1")

# ---------------------------------------------------------------- S7 other classes, boundary
print("== S7 Dulac 1/(xy) and boundary")
for case in (9, 10):
    fx, fy = mass_action(case, (C1, C2, C3), (D1, D2, D3), (k1, k2, k3), x, y)
    div = sp.diff(fx / (x * y), x) + sp.diff(fy / (x * y), y)
    if case == 9:
        tgt = (-D1 * k1 * x ** 3 + C1 * k1 * x ** 2 * y - C3 * k3 * y ** 2) / (x ** 2 * y ** 2)
        txt = "(-d1 k1 x^3 + c1 k1 x^2 y - c3 k3 y^2)/(x^2 y^2)"
    else:
        tgt = (-D1 * k1 * x ** 3 + C1 * k1 * x ** 2 * y - C3 * k3 * y - D3 * k3 * x) / (x ** 2 * y ** 2)
        txt = "(-d1 k1 x^3 + c1 k1 x^2 y - c3 k3 y - d3 k3 x)/(x^2 y^2)"
    check(is_zero(div - tgt), f"S7 case {case}: div(F/(xy)) = {txt}")
    print("     -> with d1 >= 0, c3 >= 0 (and d3 >= 0 in case 10): < 0 on the open quadrant if c1 < 0, or c1 = 0 (then d1 > 0)")
# LVA type: 2X -> (2+c1)X, X+Y -> (1+d2)Y, Y -> 0
fx, fy = mass_action(9, (C1, -1, 0), (0, D2, -1), (k1, k2, k3), x, y)
div = sp.simplify(sp.diff(fx / (x * y), x) + sp.diff(fy / (x * y), y))
check(is_zero(div - C1 * k1 / y), f"S7 LVA type: div(F/(xy)) = c1 k1 / y > 0 for c1 >= 1")
# boundary: for the three families the axes are crossed inwards
r = sp.Symbol("r", nonnegative=True)   # family 10a: d3 = r >= 0, c3 = r + q (q >= 1)
for name in ("9a", "9b", "10a"):
    case, cv, dv = fam[name]
    fx, fy = mass_action(case, cv, dv, (k1, k2, k3), x, y)
    ex = sp.factor(fx.subs(x, 0)); ey = sp.factor(fy.subs(y, 0))
    if name == "10a":
        ex = sp.expand(ex.subs(c3, r + q)); ey = sp.expand(ey.subs(c3, r + q))
    print(f"     {name}: xdot|x=0 = {ex},  ydot|y=0 = {ey}")
    okb = sp.simplify(ex).is_positive and sp.simplify(ey).is_positive
    check(okb, f"S7 {name}: xdot > 0 on {{x=0, y>0}} and ydot > 0 on {{y=0, x>0}}")

print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
