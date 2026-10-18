#!/usr/bin/env python3
"""Symbolic certificate: no periodic orbit for the Banaji-Boros-Hofbauer example, and for a subfamily.

Example (BBH 2024, Section 6):  2X -> 3X + 2Y (k1),  X + Y -> 0 (k2),  Y -> X + Y (k3)
    xdot = k1 x^2 - k2 x y + k3 y,   ydot = 2 k1 x^2 - k2 x y.
Subfamily (case 9 with c2 = d2 = -1, d3 = 0):
    2X -> (2+c1) X + d1 Y,  X + Y -> 0,  Y -> (1+c3) X + Y,   c1, c3 >= 1, d1 >= 1 integers,
    xdot = c1 k1 x^2 - k2 x y + c3 k3 y,   ydot = d1 k1 x^2 - k2 x y.
Claim: with the Dulac function B = y^(-b), b = 1 + 2 k1 (d1 - c1)/k2, the divergence of B F times y^(b+1) is the
quadratic form Q(x,y) = A x^2 + Bc x y + C y^2 with C = -k2 < 0 and discriminant 4 d1 k1 (k1 (2 c1 - d1) - k2),
negative for all k > 0 whenever d1 >= 2 c1; hence Q < 0 on the open quadrant. The example is c1 = c3 = 1, d1 = 2.
All checks are exact symbolic identities (sympy). Runtime < 10 s. Prints PASS/FAIL, exits 1 on failure.
"""
import sys
import sympy as sp

x, y = sp.symbols("x y", positive=True)
k1, k2, k3 = sp.symbols("k1 k2 k3", positive=True)
c1, c3, d1 = sp.symbols("c1 c3 d1", positive=True)
fails = 0


def check(cond, msg):
    global fails
    print(("PASS " if cond else "FAIL ") + msg)
    fails += (not cond)


def analyse(f, g, label, d1v, c1v):
    b = 1 + 2 * k1 * (d1v - c1v) / k2
    B = y ** (-b)
    div = sp.diff(B * f, x) + sp.diff(B * g, y)
    Q = sp.expand(sp.simplify(div * y ** (b + 1)))
    P = sp.Poly(Q, x, y)
    check(all(sum(m) == 2 for m in P.monoms()), f"{label}: y^(b+1) div(B F) is a quadratic form in (x, y)")
    A = sp.simplify(P.coeff_monomial(x ** 2))
    Bc = sp.simplify(P.coeff_monomial(x * y))
    C = sp.simplify(P.coeff_monomial(y ** 2))
    disc = sp.factor(sp.simplify(Bc ** 2 - 4 * A * C))
    print(f"   coefficients: x^2: {A},  xy: {Bc},  y^2: {C}")
    print(f"   discriminant: {disc}")
    check(sp.simplify(C + k2) == 0, f"{label}: coefficient of y^2 is -k2 < 0")
    target = 4 * d1v * k1 * (k1 * (2 * c1v - d1v) - k2)
    check(sp.simplify(disc - target) == 0, f"{label}: discriminant = 4 d1 k1 (k1 (2 c1 - d1) - k2)")
    return disc


# --- the example
f = k1 * x ** 2 - k2 * x * y + k3 * y
g = 2 * k1 * x ** 2 - k2 * x * y
disc = analyse(f, g, "example", 2, 1)
check(sp.simplify(disc + 8 * k1 * k2) == 0, "example: discriminant = -8 k1 k2 < 0 for all k > 0")
check(sp.simplify(g.subs(y, 0) - 2 * k1 * x ** 2) == 0, "example: on y = 0, ydot = 2 k1 x^2 > 0 (axis not invariant)")
check(sp.simplify(f.subs(x, 0) - k3 * y) == 0, "example: on x = 0, xdot = k3 y > 0 (axis not invariant)")
eq = [s for s in sp.solve([f, g], [x, y], dict=True) if s[x] != 0]
check(len(eq) == 1, f"example: unique positive equilibrium {eq}")
J = sp.Matrix([f, g]).jacobian([x, y]).subs(eq[0])
check(sp.simplify(J.det() - 4 * k1 * k3 ** 2 / k2) == 0 and sp.simplify(J.trace() + 2 * k3) == 0,
      f"example: det J = {sp.simplify(J.det())} > 0, tr J = {sp.simplify(J.trace())} < 0")

# --- the subfamily
f = c1 * k1 * x ** 2 - k2 * x * y + c3 * k3 * y
g = d1 * k1 * x ** 2 - k2 * x * y
disc = analyse(f, g, "subfamily", d1, c1)
# Theorem 9 (case 9, c3 > 0): Hopf iff (1/2)(d3/c3 + d1/c1) < d2/c2 < d1/c1, here d2/c2 = 1, d3 = 0: c1 < d1 < 2 c1
# positive kernel u = c x d: u1 = c2 d3 - c3 d2 = c3, u2 = c3 d1 - c1 d3 = c3 d1, u3 = c1 d2 - c2 d1 = d1 - c1
cv = sp.Matrix([c1, -1, c3]); dv = sp.Matrix([d1, -1, 0])
u = cv.cross(dv)
check(list(u) == [c3, c3 * d1, d1 - c1], f"subfamily: u = c x d = {list(u)}")
print("   => periodic orbits require u > 0, i.e. d1 > c1; Hopf region is c1 < d1 < 2 c1;")
print("      for d1 >= 2 c1 the discriminant is < 0 for every k > 0 (2 c1 - d1 <= 0), so Q is negative definite.")
check(sp.simplify(g.subs(y, 0) - d1 * k1 * x ** 2) == 0 and sp.simplify(f.subs(x, 0) - c3 * k3 * y) == 0,
      "subfamily: both axes are crossed inwards (ydot > 0 on y = 0, xdot > 0 on x = 0)")
print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
