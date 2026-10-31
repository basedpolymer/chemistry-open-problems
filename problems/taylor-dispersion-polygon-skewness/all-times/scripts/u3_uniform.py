#!/usr/bin/env python3
"""
u3_uniform.py -- short proof that <u^3> > 0 for every n >= N1 (default N1 = 40), independent of the older script
sg_uniform.py: the L2 norm of the harmonic part H is bounded by the elementary Fourier argument of the paper
(Lemma "H small") instead of an explicit harmonic approximant.

Apothem-1 n-gon, t = tan(pi/n), R = sqrt(1+t^2), Y = (1-r^2)/2 + t^2/6, H = ut - Y harmonic, H = tau^2/2 - t^2/6 on the side.
With Y_c = Y - <Y>, phi = Y_c^2 - <Y_c^2>, -Lap Phi = phi (radial polynomial), rho = Phi|_dP - boundary mean:
   <u^3> = M(t) + 3 (1/|P|) int_dP rho d_nH + 3 <Y_c H_c^2> + <H_c^3>,  M(t) = t^4/360 + 19 t^6/7560,
   |(3/|P|) int_dP rho d_nH| <= 6 (t+R) t^5 sqrt(rq(t)/3),  int_0^t rho^2 = t^9 rq(t),
   <H_c^2> <= (t^4/16) (1/(n+1) + 4 t^2/3)  and  1/(n+1) <= t/pi,
   |3 <Y_c H_c^2> + <H_c^3>| <= (3 (1/4 + 5t^2/12) + t^2/2) <H_c^2>.
Claim: F(t) = M/t^4 - 6(t+R) t sqrt(rq/3) - (3(1/4+5t^2/12)+t^2/2)(t/pi + 4t^2/3)/16 > 0 on [0, tan(pi/N1)].
Runtime: seconds.  PASS/FAIL, exit 1 on failure.
"""
import sys
import sympy as sp
from flint import arb, fmpq, ctx

ctx.prec = 200
N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 40
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


def neg_lap_inv(q):
    P = sp.Poly(sp.expand(q), s)
    return sp.expand(sum(-c * s ** (k + 1) / (4 * (k + 1) ** 2) for (k,), c in P.terms()))


Y = (1 - s) / 2 + t ** 2 / 6
Yc = sp.expand(Y - avg(Y))
b = sp.expand(-Y.subs(s, 1 + tau ** 2))
phi = sp.expand(Yc ** 2 - avg(Yc ** 2))
Phi = neg_lap_inv(phi)
Phib = sp.expand(Phi.subs(s, 1 + tau ** 2))
cbar = sp.expand(sp.integrate(Phib, (tau, 0, t)) / t)
rho = sp.expand(Phib - cbar)
rho2 = sp.expand(sp.integrate(rho ** 2, (tau, 0, t)))
dnPhi = sp.expand((2 * sp.diff(Phi, s)).subs(s, 1 + tau ** 2))
M = sp.expand(avg(Yc ** 3) - 3 * sp.expand(2 * sp.integrate(sp.expand(b * dnPhi), (tau, 0, t)) / t))
check("M(t) = t^4/360 + 19 t^6/7560 exactly", sp.expand(M - t ** 4 / 360 - 19 * t ** 6 / 7560) == 0)
rq = sp.cancel(rho2 / t ** 9)
Prq = sp.Poly(rq, t)
co = [(m[0], sp.Rational(c)) for m, c in Prq.terms()]


def frq(tb):
    v = arb(0)
    for k, c in co:
        v += arb(fmpq(int(c.p), int(c.q))) * tb ** k
    return v


PI = arb.pi()


def F(tb):
    tb = tb.nonnegative_part()
    R = (1 + tb * tb).sqrt()
    E1 = 6 * (tb + R) * tb * (frq(tb) / 3).nonnegative_part().sqrt()
    E2 = (3 * (arb(1) / 4 + 5 * tb * tb / 12) + tb * tb / 2) * (tb / PI + 4 * tb * tb / 3) / 16
    return arb(1) / 360 + 19 * tb * tb / 7560 - E1 - E2


tN = arb(((PI / N1).tan()).upper())
pieces = 400
allpos = True
worst = None
for i in range(pieces):
    tb = (tN * i / pieces).union(tN * (i + 1) / pieces)
    v = F(tb)
    worst = v.lower() if worst is None or v.lower() < worst else worst
    allpos &= bool(v > 0)
print(f"   min F on (0, tan(pi/{N1})] >= {float(worst):.4e}  (F(0) = {1/360:.5f})")
check(f"<u^3> > 0 for all n >= {N1}", allpos)
print("PASS (all claims)" if ok_all else "FAIL (some claim)")
sys.exit(0 if ok_all else 1)
