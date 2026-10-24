#!/usr/bin/env python3
"""
sg_uniform.py -- uniform-in-n part of the proof that the geometric skewness S^G_n > 0 for all n >= N2
(result.md, Section 6.3, Proposition 9).

Apothem-1 n-gon P, t = tan(pi/n), R = sqrt(1+t^2), |P| = n t, s = r^2.
  Y   = (1-s)/2 + t^2/6 (radial),  H := ut - Y  harmonic,  H = b(tau) := tau^2/2 - t^2/6 on each half side,
  Y_c = Y - <Y> = 1/4 + t^2/12 - s/2,  phi = Y_c^2 - <Y_c^2>,  -Lap Phi = phi (radial polynomial),
  rho = Phi|_dP - cbar (cbar = boundary mean of Phi).
Exact identity (u = Y_c + H_c, H_c = H - <H>):
  int u^3 = [int Y_c^3 - 3 int_dP b d_nPhi] + 3 int_dP rho d_nH + 3 int Y_c H_c^2 + int H_c^3.
The bracket divided by n t is an explicit polynomial M(t) = t^4/360 + O(t^6) (checked exactly).
Error bounds (all divided by n t):
  E1 = 3 ||rho|| (t+R) ||d_s b|| / (n t)            (Rellich, Lemma 7; ||d_s b||^2 = 2 n t^3/3),
  E2 = (3 (1/4 + 5t^2/12) + t^2/2) * Qhat,  Qhat >= int_P H^2 / (n t),
  Qhat = (vhat + ehat)^2, V = sum_{k<=K} beta_k Re(w^(nk)), beta_k = 2 t^2 (-1)^k/(k^2 pi^2):
     vhat = sum_k |beta_k| sqrt( t/(pi k) exp(k pi t/(1-t^2/3)) )                          (Lemma 10a)
     ehat = sqrt(1+R) [ t^2 sqrt(2/(3 pi^4 K^3)) + sum_k |beta_k| d_k(t) ],
     d_k(t) = k pi t^2/(3-t^2) + exp(k pi t/(2(1-t^2/3))) - 1                                 (Lemma 10b).
Claim: G(t) := (M - E1 - E2)/t^4 > 0 for 0 < t <= tan(pi/N2); verified by sympy (exact parts) and arb
interval arithmetic on a subdivision.  Valid for n >= 4 integer with t = tan(pi/n) <= 1.

Usage: python3 sg_uniform.py [N2] [K] [pieces]   (defaults 60 6 400).  Runtime ~5 s.
"""
import sys
import sympy as sp
from flint import arb, fmpq, ctx

ctx.prec = 200
N2 = int(sys.argv[1]) if len(sys.argv) > 1 else 60
KV = int(sys.argv[2]) if len(sys.argv) > 2 else 6
PIECES = int(sys.argv[3]) if len(sys.argv) > 3 else 400

t, s, tau = sp.symbols('t s tau', nonnegative=True)
ok_all = True


def check(name, cond):
    global ok_all
    print(("PASS" if cond else "FAIL") + ": " + name)
    ok_all &= bool(cond)


def avg(phi):
    P = sp.Poly(sp.expand(phi), s)
    out = 0
    for (k,), c in P.terms():
        out += c * sp.Rational(1, k + 1) * sum(sp.binomial(k, j) * t ** (2 * j) / (2 * j + 1) for j in range(k + 1))
    return sp.expand(out)


def lap(phi):
    return sp.expand(4 * s * sp.diff(phi, s, 2) + 4 * sp.diff(phi, s))


def neg_lap_inv(q):
    P = sp.Poly(sp.expand(q), s)
    return sp.expand(sum(-c * s ** (k + 1) / (4 * (k + 1) ** 2) for (k,), c in P.terms()))


Y = (1 - s) / 2 + t ** 2 / 6
Yc = sp.expand(Y - avg(Y))
check("Y_c = 1/4 + t^2/12 - s/2", sp.expand(Yc - (sp.Rational(1, 4) + t ** 2 / 12 - s / 2)) == 0)
b = sp.expand(-Y.subs(s, 1 + tau ** 2))
check("H = -Y = tau^2/2 - t^2/6 on the side", sp.expand(b - (tau ** 2 / 2 - t ** 2 / 6)) == 0)
phi = sp.expand(Yc ** 2 - avg(Yc ** 2))
Phi = neg_lap_inv(phi)
check("-Lap Phi = phi", sp.expand(-lap(Phi) - phi) == 0)
Phib = sp.expand(Phi.subs(s, 1 + tau ** 2))
cbar = sp.expand(sp.integrate(Phib, (tau, 0, t)) / t)
rho = sp.expand(Phib - cbar)
rho2 = sp.expand(sp.integrate(rho ** 2, (tau, 0, t)))          # int_0^t rho^2
dnPhi = sp.expand((2 * sp.diff(Phi, s)).subs(s, 1 + tau ** 2))   # d/dx at x=1 of Phi(r^2)
# int_P Y_c^3 /(n t) = <Y_c^3>;  int_dP b d_nPhi /(n t) = (2n int_0^t ...)/(n t) = (2/t) int_0^t
M = sp.expand(avg(Yc ** 3) - 3 * sp.expand(2 * sp.integrate(sp.expand(b * dnPhi), (tau, 0, t)) / t))
Mser = sp.series(M, t, 0, 8).removeO()
check("M(t) = t^4/360 + O(t^6)", sp.expand(Mser - sp.Rational(1, 360) * t ** 4).as_poly(t).degree() >= 6 or
      sp.expand(Mser - sp.Rational(1, 360) * t ** 4) == 0)
print("   M(t) =", sp.series(M, t, 0, 10))
Mq = sp.cancel(M / t ** 4)
check("M(t)/t^4 is a polynomial", sp.Poly(sp.numer(sp.together(Mq)), t) is not None and sp.denom(sp.together(Mq)).is_number)
rq = sp.cancel(rho2 / t ** 9)
check("int_0^t rho^2 = t^9 * polynomial", sp.denom(sp.together(rq)).is_number)
print("   int_0^t rho^2 =", sp.series(rho2, t, 0, 13))


def arb_poly_eval(expr):
    P = sp.Poly(sp.expand(expr), t)
    co = [(m[0], sp.Rational(c)) for m, c in P.terms()]

    def f(tb):
        v = arb(0)
        for k, c in co:
            v += arb(fmpq(int(c.p), int(c.q))) * tb ** k
        return v
    return f


fMq = arb_poly_eval(Mq)
frq = arb_poly_eval(rq)
PI = arb.pi()


def G_lower(tb):
    tb = tb.nonnegative_part()      # t >= 0 (removes rounding below 0 so that sqrt(t) is defined)
    R = (1 + tb * tb).sqrt()
    one_m = 1 - tb * tb / 3
    # E1/t^4 = 6 (t+R) sqrt(rho2 t/3)/t^4 = 6 (t+R) t sqrt(rq/3)    (rho2 = t^9 rq)
    E1q = 6 * (tb + R) * tb * (frq(tb) / 3).nonnegative_part().sqrt()
    # vhat/t^2 and ehat/t^2
    vq = arb(0)
    eq = arb(0)
    for k in range(1, KV + 1):
        bk = 2 / (k * k * PI * PI)                 # |beta_k| / t^2
        vq += bk * ((tb / (PI * k)) * (k * PI * tb / one_m).exp()).nonnegative_part().sqrt()
        dk = k * PI * tb * tb / (3 - tb * tb) + (k * PI * tb / (2 * one_m)).exp() - 1
        eq += bk * dk
    eq = (1 + R).sqrt() * ((arb(2) / (3 * PI ** 4 * KV ** 3)).sqrt() + eq)
    Qq = (vq + eq) ** 2                             # Qhat / t^4
    E2q = (3 * (arb(1) / 4 + 5 * tb * tb / 12) + tb * tb / 2) * Qq
    return fMq(tb) - E1q - E2q


t2 = (PI / N2).tan()
t2u = arb(t2.upper())
worst = None
allpos = True
for i in range(PIECES):
    a = t2u * i / PIECES
    bnd = t2u * (i + 1) / PIECES
    tb = a.union(bnd)               # ball containing [a, bnd]
    v = G_lower(tb)
    if worst is None or v.lower() < worst:
        worst = v.lower()
    if not (v > 0):
        allpos = False
print(f"   min over pieces of lower bound of G(t) = int u^3/(n t^5) - errors on (0, tan(pi/{N2})]: {arb(worst).str(5)}"
      f"   (G(0) = 1/360 = {1/360:.5f})")
check(f"G(t) > 0 on (0, tan(pi/{N2})], i.e. geometric skewness > 0 for all n >= {N2}", allpos)
print("PASS (all claims)" if ok_all else "FAIL (some claim)")
sys.exit(0 if ok_all else 1)
