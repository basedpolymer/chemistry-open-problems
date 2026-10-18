#!/usr/bin/env python3
"""Independent symbolic verification (verifier's own code, written without reading the solver's reasoning).

Environment: WSL Fedora, Python 3 + sympy. Runtime: a few seconds. Single thread.
Prints PASS/FAIL for each claim and exits 1 on any failure.

What is checked:
 V1  The ODE of the BBH example is the mass-action ODE of the network 2X -> 3X+2Y, X+Y -> 0, Y -> X+Y
     (built from source/target complexes by a generic mass-action constructor, not typed in by hand).
 V2  Case-9 source configuration (2X, X+Y, Y) is positively oriented (BBH (15)); the subfamily is
     c = (c1, -1, c3), d = (d1, -1, 0); its *targets* are (2+c1)X + d1Y, 0, c3X + Y (NOT (1+c3)X + Y).
 V3  For a general Dulac function B = x^alpha y^beta the polynomial x^(1-alpha) y^(1-beta) div(BF) is computed;
     with alpha = 0, beta = -b it reduces to the quadratic form Q_b = -b d1 k1 x^2 + (2 c1 k1 + (b-1) k2) x y - k2 y^2.
 V4  Delta(b) = disc(Q_b) is a convex quadratic in b, minimised at b* = 1 + 2 k1 (d1-c1)/k2,
     Delta(b*) = 4 d1 k1 (k1 (2c1-d1) - k2).
 V5  For d1 = 2 c1 + e (e >= 0) the quantity -Delta(b*)/4 = d1 k1 (k2 + e k1) is a polynomial with nonnegative
     coefficients and a positive term, hence > 0; and the y^2 coefficient is -k2 < 0 => Q negative definite.
 V6  Example (c1=c3=1, d1=2): Delta = -8 k1 k2.
 V7  Kernel of the stoichiometric matrix: u = (c3, c3 d1, d1 - c1); positive equilibrium iff d1 > c1;
     equilibrium (x*, y*) solves k1 x^2 = mu c3, k2 x y = mu c3 d1, k3 y = mu (d1-c1).
 V8  Sign of tr J at the positive equilibrium equals sign of k1 (2c1 - d1) - k2, i.e. the Hopf curve tr J = 0
     is exactly the curve Delta(b*) = 0 (so the Dulac exclusion is sharp: it covers exactly tr J <= 0);
     det J > 0 at the positive equilibrium.
 V9  At Delta(b*) = 0 the form is -k2 (y - r x)^2 with r > 0: semidefinite, zero only on a ray (measure zero).
 V10 Example: equilibrium (2k3/k2, 4k1k3/k2^2), tr J = -2k3, det J = 4 k1 k3^2/k2.
 V11 Axes: on y = 0, ydot = d1 k1 x^2; on x = 0, xdot = c3 k3 y.
 V12 BBH Theorem 9 (arXiv v1 numbering; Theorem 3 in the journal version), case 9 with c3 > 0, specialised
     to c2 = d2 = -1, d3 = 0, reads c1 < d1 < 2 c1.
 V13 The Dulac case also covers the admissible c1 in {-2,-1,0} (with d1 >= 1): Delta(b*) < 0.
"""
import sys
import random
import sympy as sp

x, y = sp.symbols("x y", positive=True)
k1, k2, k3 = sp.symbols("k1 k2 k3", positive=True)
c1, c3, d1, e = sp.symbols("c1 c3 d1 e", positive=True)
b, al, be = sp.symbols("b alpha beta", real=True)
fails = 0


def check(cond, msg):
    global fails
    cond = bool(cond)
    print(("PASS " if cond else "FAIL ") + msg)
    fails += (not cond)


def mass_action(reactions, rates):
    """reactions: list of (source, target) with complexes as (nX, nY)."""
    fx, fy = 0, 0
    for (src, tgt), k in zip(reactions, rates):
        rate = k * x ** src[0] * y ** src[1]
        fx += (tgt[0] - src[0]) * rate
        fy += (tgt[1] - src[1]) * rate
    return sp.expand(fx), sp.expand(fy)


# ---------- V1: the example ODE from the network
ex_net = [((2, 0), (3, 2)), ((1, 1), (0, 0)), ((0, 1), (1, 1))]
fx, fy = mass_action(ex_net, [k1, k2, k3])
check(sp.expand(fx - (k1 * x**2 - k2 * x * y + k3 * y)) == 0 and sp.expand(fy - (2 * k1 * x**2 - k2 * x * y)) == 0,
      f"V1 example ODE from network: xdot = {fx}, ydot = {fy}")

# ---------- V2: case 9 orientation and subfamily targets
s1, s2, s3 = sp.Matrix([2, 0]), sp.Matrix([1, 1]), sp.Matrix([0, 1])
orient = (s2 - s1)[0] * (s3 - s1)[1] - (s2 - s1)[1] * (s3 - s1)[0]
check(orient > 0, f"V2 sources (2X, X+Y, Y) positively oriented (det = {orient})")
# subfamily with targets source + reaction vector
cvec, dvec = [c1, -1, c3], [d1, -1, 0]
targets = [(s[0] + cvec[i], s[1] + dvec[i]) for i, s in enumerate([(2, 0), (1, 1), (0, 1)])]
print("   subfamily targets:", targets)
check(targets[2] == (c3, 1), "V2 third target is c3 X + Y (so the example c3 = 1 gives Y -> X + Y)")
check(targets[0] == (2 + c1, d1) and targets[1] == (0, 0), "V2 first two targets (2+c1)X + d1Y and 0")
# ODE of the subfamily from the generic constructor, for concrete integers, against the claimed ODE
for (C1, C3, D1) in [(1, 1, 2), (2, 3, 5), (3, 1, 1), (1, 4, 7)]:
    net = [((2, 0), (2 + C1, D1)), ((1, 1), (0, 0)), ((0, 1), (C3, 1))]
    gx, gy = mass_action(net, [k1, k2, k3])
    ok = (sp.expand(gx - (C1 * k1 * x**2 - k2 * x * y + C3 * k3 * y)) == 0
          and sp.expand(gy - (D1 * k1 * x**2 - k2 * x * y)) == 0)
    check(ok, f"V2 subfamily ODE from network (c1,c3,d1)=({C1},{C3},{D1}): xdot={gx}, ydot={gy}")
# with the misprinted target (1+c3)X + Y the ODE would have (1+c3) k3 y instead
net_bad = [((2, 0), (3, 2)), ((1, 1), (0, 0)), ((0, 1), (2, 1))]
gx, _ = mass_action(net_bad, [k1, k2, k3])
check(sp.expand(gx - fx) != 0, "V2 target (1+c3)X+Y with c3=1 (i.e. Y -> 2X+Y) does NOT give the BBH example")

# ---------- V3: general Dulac computation
f = c1 * k1 * x**2 - k2 * x * y + c3 * k3 * y
g = d1 * k1 * x**2 - k2 * x * y
Bgen = x**al * y**be
divgen = sp.diff(Bgen * f, x) + sp.diff(Bgen * g, y)
Pgen = sp.expand(sp.simplify(divgen / (x**(al - 1) * y**(be - 1))))
print("   x^(1-a) y^(1-b) div(x^a y^b F) =", sp.collect(Pgen, [x, y]))
Q = sp.expand(Pgen.subs({al: 0, be: -b}) / x)  # alpha = 0: overall factor x^(-1) y^(-b-1) * x
Qclaim = -b * d1 * k1 * x**2 + (2 * c1 * k1 + (b - 1) * k2) * x * y - k2 * y**2
check(sp.expand(Q - Qclaim) == 0, f"V3 y^(b+1) div(y^-b F) = {sp.collect(Q, [x, y])}")
# direct recomputation, no general formula
Bd = y**(-b)
Qd = sp.expand(sp.simplify((sp.diff(Bd * f, x) + sp.diff(Bd * g, y)) * y**(b + 1)))
check(sp.expand(Qd - Qclaim) == 0, "V3 direct computation agrees")

# ---------- V4: optimisation in b
A = -b * d1 * k1
Bc = 2 * c1 * k1 + (b - 1) * k2
C = -k2
Delta = sp.expand(Bc**2 - 4 * A * C)
check(sp.Poly(Delta, b).degree() == 2 and sp.Poly(Delta, b).LC() == k2**2, "V4 Delta(b) is a convex quadratic in b")
bstar = sp.solve(sp.diff(Delta, b), b)
check(len(bstar) == 1 and sp.simplify(bstar[0] - (1 + 2 * k1 * (d1 - c1) / k2)) == 0, f"V4 argmin b* = {bstar}")
Dmin = sp.factor(Delta.subs(b, bstar[0]))
check(sp.simplify(Dmin - 4 * d1 * k1 * (k1 * (2 * c1 - d1) - k2)) == 0, f"V4 Delta(b*) = {Dmin}")

# ---------- V5: sign for d1 >= 2 c1
negq = sp.expand(-Dmin.subs(d1, 2 * c1 + e) / 4)
P = sp.Poly(negq, c1, e, k1, k2)
check(all(co > 0 for co in P.coeffs()), f"V5 -Delta/4 with d1 = 2c1+e: {negq} has positive coefficients")
negq0 = sp.expand(-Dmin.subs(d1, 2 * c1) / 4)
check(sp.expand(negq0 - 2 * c1 * k1 * k2) == 0, "V5 at e = 0: -Delta/4 = 2 c1 k1 k2 > 0")
check(sp.simplify(bstar[0].subs(d1, 2 * c1 + e)) .is_positive, "V5 b* > 0 for d1 >= 2c1 (x^2 coefficient -b d1 k1 < 0)")

# randomized exact check of negativity of Q at b* on the open quadrant (rationals)
random.seed(1)
bad = 0
for _ in range(400):
    C1v = random.randint(1, 6); D1v = 2 * C1v + random.randint(0, 6)
    K1 = sp.Rational(random.randint(1, 10**4), random.randint(1, 10**4))
    K2 = sp.Rational(random.randint(1, 10**4), random.randint(1, 10**4))
    X = sp.Rational(random.randint(1, 10**4), random.randint(1, 10**4))
    Y = sp.Rational(random.randint(1, 10**4), random.randint(1, 10**4))
    val = Qclaim.subs({b: 1 + 2 * K1 * (D1v - C1v) / K2, d1: D1v, c1: C1v, k1: K1, k2: K2, x: X, y: Y})
    bad += int(bool(val >= 0))
check(bad == 0, "V5 400 random exact rational samples with d1 >= 2c1: Q(x,y) < 0")

# ---------- V6: example
check(sp.simplify(Dmin.subs({c1: 1, d1: 2}) + 8 * k1 * k2) == 0, "V6 example: Delta = -8 k1 k2")

# ---------- V7: kernel and equilibrium
Gam = sp.Matrix([[c1, -1, c3], [d1, -1, 0]])
u = sp.Matrix(cvec).cross(sp.Matrix(dvec))
check(sp.simplify(Gam * u) == sp.zeros(2, 1) and list(u) == [c3, c3 * d1, d1 - c1], f"V7 u = c x d = {list(u)} in ker Gamma")
check(Gam.subs({c3: 1}).rank() == 2, "V7 rank Gamma = 2 (generic c1,d1 with c3 = 1)")
mu = sp.symbols("mu", positive=True)
# equilibrium: (k1 x^2, k2 x y, k3 y) = mu u, needs d1 > c1
sol = sp.solve([k1 * x**2 - mu * c3, k2 * x * y - mu * c3 * d1, k3 * y - mu * (d1 - c1)], [x, y, mu], dict=True)
check(len(sol) == 1, f"V7 unique positive equilibrium (for d1 > c1): {sol}")
xs, ys = sol[0][x], sol[0][y]
check(sp.simplify(f.subs({x: xs, y: ys})) == 0 and sp.simplify(g.subs({x: xs, y: ys})) == 0, "V7 it is an equilibrium")

# ---------- V8: trace / Hopf curve
J = sp.Matrix([f, g]).jacobian([x, y])
trJ = sp.simplify(J.trace().subs({x: xs, y: ys}))
detJ = sp.simplify(J.det().subs({x: xs, y: ys}))
# claim: tr J = (positive factor) * (k1 (2c1 - d1) - k2)
ratio = sp.simplify(trJ / (k1 * (2 * c1 - d1) - k2))
print("   tr J =", sp.factor(trJ), "; tr J / (k1(2c1-d1)-k2) =", ratio)
check(ratio.is_positive or sp.simplify(ratio.subs(d1, c1 + e)).is_positive,
      "V8 sign tr J = sign(k1 (2c1 - d1) - k2) at the positive equilibrium (d1 = c1 + e > c1)")
detf = sp.factor(detJ.subs(d1, c1 + e))
print("   det J (d1 = c1+e) =", detf)
check(sp.simplify(detJ.subs(d1, c1 + e)).is_positive, "V8 det J > 0 at the positive equilibrium")

# ---------- V9: Delta = 0 case
k2h = k1 * (2 * c1 - d1)
Qh = sp.expand(Qclaim.subs(b, bstar[0]).subs(k2, k2h))
r = sp.symbols("r", positive=True)
rr = sp.simplify((2 * c1 * k1 + (bstar[0] - 1) * k2) / (2 * k2)).subs(k2, k2h)
check(sp.simplify(Qh + k2h * (y - rr * x)**2) == 0, f"V9 at k2 = k1(2c1-d1): Q = -k2 (y - r x)^2 with r = {sp.simplify(rr)}")

# ---------- V10: example equilibrium and Jacobian
fe, ge = fx, fy
Je = sp.Matrix([fe, ge]).jacobian([x, y]).subs({x: 2 * k3 / k2, y: 4 * k1 * k3 / k2**2})
check(sp.simplify(fe.subs({x: 2 * k3 / k2, y: 4 * k1 * k3 / k2**2})) == 0 and sp.simplify(ge.subs({x: 2 * k3 / k2, y: 4 * k1 * k3 / k2**2})) == 0,
      "V10 example equilibrium (2k3/k2, 4k1k3/k2^2)")
check(sp.simplify(Je.trace() + 2 * k3) == 0 and sp.simplify(Je.det() - 4 * k1 * k3**2 / k2) == 0, "V10 tr J = -2k3, det J = 4k1k3^2/k2")
exsol = sp.solve([fe, ge], [x, y], dict=True)
check(len(exsol) == 1, f"V10 only positive equilibrium of the example: {exsol}")

# ---------- V11: axes
check(sp.expand(g.subs(y, 0) - d1 * k1 * x**2) == 0 and sp.expand(f.subs(x, 0) - c3 * k3 * y) == 0,
      "V11 ydot|_{y=0} = d1 k1 x^2 > 0 (x>0), xdot|_{x=0} = c3 k3 y > 0 (y>0)")

# ---------- V12: BBH Theorem 9, case 9, c3 > 0 branch specialised
C2, D2, D3 = -1, -1, 0
lhs = sp.Rational(1, 2) * (D3 / c3 + d1 / c1)
mid = sp.Rational(D2, C2)
cond1 = sp.simplify((mid - lhs) * 2 * c1)      # >0  <=>  d1 < 2c1
cond2 = sp.simplify((d1 / c1 - mid) * c1)      # >0  <=>  d1 > c1
check(sp.expand(cond1 - (2 * c1 - d1)) == 0 and sp.expand(cond2 - (d1 - c1)) == 0,
      "V12 Hopf condition 1/2(d3/c3 + d1/c1) < d2/c2 < d1/c1 becomes c1 < d1 < 2c1")

# ---------- V13: the hypothesis c1 >= 1 is not needed for the Dulac case (admissible c1 in {-2,-1,0}, d1 >= 1)
for C1v in (-2, -1, 0):
    for D1v in range(1, 6):
        Dv = sp.factor(Dmin.subs({c1: C1v, d1: D1v}))
        negpart = sp.expand(-Dv)
        Pn = sp.Poly(negpart, k1, k2)
        ok = all(co > 0 for co in Pn.coeffs())
        if not ok:
            check(False, f"V13 c1={C1v}, d1={D1v}: Delta = {Dv}")
check(True, "V13 for c1 in {-2,-1,0}, d1 in 1..5: -Delta(b*) has positive coefficients in (k1,k2) (Delta < 0)")
# general statement: c1 = -m (m >= 0), d1 >= 1: -Delta/4 = d1 k1 (k2 + k1 (d1 + 2m)) > 0
m = sp.symbols("m", nonnegative=True)
check(sp.expand(-Dmin.subs(c1, -m) / 4 - d1 * k1 * (k2 + k1 * (d1 + 2 * m))) == 0,
      "V13 c1 = -m <= 0: -Delta(b*)/4 = d1 k1 (k2 + k1 (d1 + 2m)) > 0")
# targets stay admissible: 2X -> (2+c1)X + d1Y needs 2 + c1 >= 0
check(all(2 + C1v >= 0 for C1v in (-2, -1, 0)), "V13 targets (2+c1)X + d1Y are admissible complexes for c1 >= -2")

print("ALL PASS" if not fails else f"{fails} FAILURES")
sys.exit(1 if fails else 0)
