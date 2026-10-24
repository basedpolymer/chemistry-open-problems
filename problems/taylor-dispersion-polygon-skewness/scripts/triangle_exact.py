#!/usr/bin/env python3
"""
triangle_exact.py -- exact (symbolic, sympy) treatment of the equilateral triangle, n = 3.

Triangle (inradius 1, centroid 0, as in Aminian-Camassa-McLaughlin eq. (4) of the published version):
    sides  y = -1,  y = 2 + sqrt(3) z,  y = 2 - sqrt(3) z.
Claims checked (exact arithmetic):
  1. ut = (1+y)((2-y)^2 - 3 z^2)/6 solves -Lap ut = 2, ut = 0 on the boundary.
  2. g1 below (a polynomial of degree 6) solves -Lap g1 = u := ut - <ut>, d_n g1 = 0 on all three sides,
     <g1> = 0.
  3. <u g1> = 9/1540 and <u g1^2> = -7101/476476000 < 0, hence <u g2> = 2<u g1^2> < 0:
     the long-time skewness is NEGATIVE (for every Peclet number), and the coefficient
     3<u g2>/(2<u g1>)^(3/2) = -789 sqrt(770)/309400 = -0.0707623... (paper: -0.07076).
  4. (<u g2> = 2 <u g1^2> by Lemma 1 of result.md; g2 itself is not a polynomial.)
  5. Geometric (short-time) skewness <u^3>/<u^2>^(3/2) = 4 sqrt(21)/99 = sqrt(112/3267) = 0.18515 > 0
     (paper eq. (12): 0.185), with <u^2> = 27/700, <u^3> = 27/19250 (paper eq. (14)).
Runtime: ~1-2 min (sympy 1.14). Prints PASS/FAIL per claim; exit 1 on failure.
"""
import sys
import sympy as sp

y, z = sp.symbols('y z', real=True)
s3 = sp.sqrt(3)
ok_all = True


def check(name, cond):
    global ok_all
    print(("PASS" if cond else "FAIL") + ": " + name)
    ok_all &= bool(cond)


def lap(f):
    return sp.expand(sp.diff(f, y, 2) + sp.diff(f, z, 2))


def integ(f):
    return sp.integrate(sp.integrate(sp.expand(f), (z, -(2 - y) / s3, (2 - y) / s3)), (y, -1, 2))


AREA = sp.simplify(integ(1))


def avg(f):
    return sp.nsimplify(sp.simplify(integ(f) / AREA))


ut = (1 + y) * ((2 - y) ** 2 - 3 * z ** 2) / 6
check("-Lap ut = 2", sp.simplify(lap(ut) + 2) == 0)
check("ut = 0 on y=-1", sp.expand(ut.subs(y, -1)) == 0)
check("ut = 0 on y=2+sqrt3 z", sp.expand(ut.subs(y, 2 + s3 * z)) == 0)
check("ut = 0 on y=2-sqrt3 z", sp.expand(ut.subs(y, 2 - s3 * z)) == 0)
mu = avg(ut)
check("<ut> = 3/10", mu == sp.Rational(3, 10))
u = sp.expand(ut - mu)

R = sp.Rational
g1 = (y ** 6 / 960 - y ** 5 / 96 - y ** 4 * z ** 2 / 64 + y ** 4 / 32 + y ** 3 * z ** 2 / 48 + y ** 2 * z ** 4 / 64
      + y ** 2 * z ** 2 / 16 - R(11, 120) * y ** 2 + y * z ** 4 / 32 - z ** 6 / 960 + z ** 4 / 32 - R(11, 120) * z ** 2)
g1 = sp.expand(g1 - avg(g1))
check("-Lap g1 = u", sp.expand(-lap(g1) - u) == 0)


def neumann_ok(g):
    c1 = sp.expand(sp.diff(g, y).subs(y, -1)) == 0                    # outward normal (-1,0): -g_y
    # side y = 2 - sqrt3 z (z > 0 part): outward normal proportional to (1, sqrt3)
    c2 = sp.expand((sp.diff(g, y) + s3 * sp.diff(g, z)).subs(y, 2 - s3 * z)) == 0
    c3 = sp.expand((sp.diff(g, y) - s3 * sp.diff(g, z)).subs(y, 2 + s3 * z)) == 0
    return c1 and c2 and c3


check("d_n g1 = 0 on all three sides", neumann_ok(g1))
check("<g1> = 0", avg(g1) == 0)
ug = avg(u * g1)
ug2 = avg(u * g1 ** 2)
check("<u g1> = 9/1540", ug == R(9, 1540))
check("<u g1^2> = -7101/476476000 (negative)", ug2 == R(-7101, 476476000) and ug2 < 0)
kappa = sp.nsimplify(sp.simplify(6 * ug2 / (2 * ug) ** R(3, 2)))
check("LT coefficient 3<u g2>/(2<u g1>)^(3/2) = -789 sqrt(770)/309400",
      sp.simplify(kappa + 789 * sp.sqrt(770) / 309400) == 0)
print("   kappa_3 (inradius 1) =", sp.N(kappa, 15), "; circumradius 1:", sp.N(kappa / 2, 15))

# Note: <u g2> = 2<u g1^2> is Lemma 1 of result.md (Green's identities, valid on any Lipschitz domain).
# (A polynomial ansatz for g2 of degree <= 12 has no solution, so g2 is not a polynomial; it is not needed.)

u2 = avg(u ** 2)
u3 = avg(u ** 3)
SG = sp.nsimplify(sp.simplify(u3 / u2 ** R(3, 2)))
check("<u^2> = 27/700, <u^3> = 27/19250", u2 == R(27, 700) and u3 == R(27, 19250))
check("geometric skewness = 4 sqrt(21)/99 > 0", sp.simplify(SG - 4 * sp.sqrt(21) / 99) == 0)
print("   S^G =", sp.N(SG, 15))
print("PASS (all claims)" if ok_all else "FAIL (some claim)")
sys.exit(0 if ok_all else 1)
