#!/usr/bin/env python3
"""
lt_uniform.py -- the uniform-in-n part of the proof that the long-time skewness coefficient is positive
for every regular n-gon with n >= N1 (default N1 = 30), using RADIAL approximations around the disk.

Setting (result.md, Section 5): P = regular n-gon with apothem 1, t = tan(pi/n), R = sqrt(1+t^2).
Radial approximations (exact, independent of n except through t):
    ut_N = (1 - r^2)/2 + t^2/6,         u_N = ut_N - <ut_N> = 1/4 + t^2/12 - r^2/2,
    g_N  = r^4/32 - (1/16 + t^2/48) r^2 + b_0(t)   (so -Lap g_N = u_N, <g_N> = 0).
For every radial polynomial phi(r^2) the average over P is a polynomial in t only:
    <r^(2k)>_P = (1/(k+1)) sum_j C(k,j) t^(2j)/(2j+1)      (Lemma 6),
and all a-posteriori quantities of lt_certify.py, divided by |P| = n t, are explicit functions of t:
    ||ut_N||^2_{L2(dP)} = 2 n t^5/45, ||e||^2_{L2(dP)} = n t^5/360, osc_dP(ut_N) = t^2/2,
    sup_P |u_N| = 1/4 + 5 t^2/12.
The claim  I_n/(n t) - Err/(n t) =: F(t) > 0 for all 0 <= t <= tan(pi/N1)  is verified with sympy
(exact rational polynomial identities) and arb interval arithmetic on a subdivision of [0, tan(pi/N1)].

Usage: python3 lt_uniform.py [N1] [pieces]     (defaults 30, 400)
Runtime: ~10 s (python3.13, sympy 1.14, python-flint 0.9.0). Prints PASS/FAIL, exit 1 on failure.
"""
import sys
import sympy as sp
from flint import arb, fmpq, ctx

ctx.prec = 200
N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 30
PIECES = int(sys.argv[2]) if len(sys.argv) > 2 else 400

t, s, tau = sp.symbols('t s tau', nonnegative=True)   # s = r^2
ok_all = True


def check(name, cond):
    global ok_all
    print(("PASS" if cond else "FAIL") + ": " + name)
    ok_all &= bool(cond)


def avg(phi):
    """<phi(r^2)>_P for a polynomial phi in s = r^2 (Lemma 6)."""
    P = sp.Poly(sp.expand(phi), s)
    out = 0
    for (k,), c in P.terms():
        out += c * sp.Rational(1, k + 1) * sum(sp.binomial(k, j) * t ** (2 * j) / (2 * j + 1) for j in range(k + 1))
    return sp.expand(out)


def lap(phi):
    """Laplacian of a radial function phi(s), s = r^2:  4 s phi'' + 4 phi'."""
    return sp.expand(4 * s * sp.diff(phi, s, 2) + 4 * sp.diff(phi, s))


# --- Lemma 6 sanity check by direct integration over the half-sector T ---------------------------
x, y = sp.symbols('x y', real=True)
for k in range(0, 4):
    direct = sp.integrate(sp.integrate((x ** 2 + y ** 2) ** k, (y, 0, t * x)), (x, 0, 1)) / (t / 2)
    check(f"Lemma 6 formula for <r^{2*k}>_P", sp.simplify(direct - avg(s ** k)) == 0)

utN = (1 - s) / 2 + t ** 2 / 6
uN = sp.expand(utN - avg(utN))
check("u_N = 1/4 + t^2/12 - r^2/2", sp.expand(uN - (sp.Rational(1, 4) + t ** 2 / 12 - s / 2)) == 0)
gN = s ** 2 / 32 - (sp.Rational(1, 16) + t ** 2 / 48) * s
gN = sp.expand(gN - avg(gN))
check("-Lap g_N = u_N", sp.expand(-lap(gN) - uN) == 0)
check("<g_N> = 0", sp.expand(avg(gN)) == 0)
# boundary data on side x = 1: r^2 = 1 + tau^2, 0 <= tau <= t
ut_b = sp.expand(utN.subs(s, 1 + tau ** 2))
check("ut_N = t^2/6 - tau^2/2 on the side", sp.expand(ut_b - (t ** 2 / 6 - tau ** 2 / 2)) == 0)
e_b = sp.expand((2 * sp.diff(gN, s)).subs(s, 1 + tau ** 2))   # d/dx g(r^2) at x=1 = 2 g'(s)
check("e = d_n g_N = (3 tau^2 - t^2)/24 on the side", sp.expand(e_b - (3 * tau ** 2 - t ** 2) / 24) == 0)
ut_b2 = sp.expand(2 * sp.integrate(ut_b ** 2, (tau, 0, t)) / t)    # ||ut_N||^2_dP / (n t)
check("||ut_N||^2_{L2(dP)} = 2 n t^5/45", sp.expand(ut_b2 - 2 * t ** 4 / 45) == 0)
e_b2 = sp.expand(2 * sp.integrate(e_b ** 2, (tau, 0, t)) / t)
check("||e||^2_{L2(dP)} = n t^5/360", sp.expand(e_b2 - t ** 4 / 360) == 0)

I_avg = avg(uN * gN ** 2)
G2 = avg(gN ** 2)
var_G2 = sp.expand(avg(gN ** 4) - G2 ** 2)
UG = avg(uN * gN)
var_UG = sp.expand(avg(uN ** 2 * gN ** 2) - UG ** 2)
check("disk limit <u g^2> = 1/184320", I_avg.subs(t, 0) == sp.Rational(1, 184320))
print("   <u_N g_N^2> =", sp.series(I_avg, t, 0, 6))


def to_arb_poly_eval(expr):
    P = sp.Poly(sp.expand(expr), t)
    coeffs = [(m[0], c) for m, c in P.terms()]

    def f(tb):
        v = arb(0)
        for k, c in coeffs:
            c = sp.Rational(c)
            v += arb(fmpq(int(c.p), int(c.q))) * tb ** k
        return v
    return f


fI, fG2, fvG2, fvUG = (to_arb_poly_eval(e) for e in (I_avg, G2, var_G2, var_UG))


def F_lower(tb):
    """Rigorous lower bound of F(t) for t in the ball tb (Lemma 5 with the radial data)."""
    R = (1 + tb * tb).sqrt()
    p = 2 * R / arb.pi()
    tau_tr = (2 * p * p + 2 * R * p).sqrt()
    epsE = tb * tb * ((1 + R) / 45).sqrt()
    eN = tb * tb / arb(360).sqrt()
    eta = p * epsE + tau_tr * eN
    om = p * eta
    osc = tb * tb / 2
    T1 = epsE * fvG2(tb).sqrt()
    T2 = osc * (2 * fG2(tb).sqrt() * om + om * om)
    T3 = 2 * fvUG(tb).sqrt() * om
    T4 = (arb(1) / 4 + 5 * tb * tb / 12) * om * om
    return fI(tb) - (T1 + T2 + T3 + T4)


t0 = (arb.pi() / N1).tan()
t0u = arb(t0.upper())
worst = None
allpos = True
for i in range(PIECES):
    a = t0u * i / PIECES
    b = t0u * (i + 1) / PIECES
    tb = arb((a + b) / 2, ((b - a) / 2).upper())
    tb = tb.union(a).union(b)
    v = F_lower(tb)
    lo = v.lower()
    if worst is None or lo < worst:
        worst = lo
    if not (v > 0):
        allpos = False
        print("   FAIL on piece", i, tb, v)
print(f"   min over pieces of the lower bound of F(t) on [0, tan(pi/{N1})]: {arb(worst).str(5)}")
check(f"F(t) > 0 for all 0 <= t <= tan(pi/{N1}), i.e. long-time coefficient > 0 for all n >= {N1}", allpos)
print("PASS (all claims)" if ok_all else "FAIL (some claim)")
sys.exit(0 if ok_all else 1)
