"""Why the global ACR step of the dimension-two proof cannot be dropped: an exact example.

Runtime: under 2 s (WSL, Python 3.12, sympy 1.14).   Usage: python global_step_example.py

The bimolecular network on X1, X2, X3 (dimension 2, S = {x1 + x2 + x3 = 0})
    X2+X3 -> 2X1 (3/2),  X1+X2 -> 2X2 (2),  X1+X3 -> 2X3 (2),  X1 -> X2 (1/2),  X1 -> X3 (1/2),
    X2 -> X3 (1),  X3 -> X2 (1)
has the normal form of the manuscript with x1 as "ACR candidate" and c = 1:
    f = psi * alpha + (x1 - 1) * P * ell,   alpha = (1, -1/2, -1/2),  ell = (0, 1, -1),
    psi = 3 x2 x3 - 2 x1 x2 - 2 x1 x3 - x1,   P = x2 - x3,
and it satisfies both coefficient-cone conditions (it is a genuine bimolecular network).  In the class
x1 + x2 + x3 = 5 the points p = (1, 1, 3) and q = (1, 3, 1) lie on H = {x1 = 1}, on a common line p + R ell, and are
BOTH exponentially stable, with P(p) P(q) < 0: every local hypothesis of the lemmas on the ACR hyperplane holds.
The network does not have ACR: on the plane P = 0 the equation psi = 0 has positive solutions with x1 != 1 (a third
steady state of the same class is exhibited).  So the global ACR argument (lemmas "global ACR" and "no transversal
sign change") is what excludes this configuration; a purely local argument cannot.
"""
import sys
import sympy as sp

x1, x2, x3, lam = sp.symbols('x1 x2 x3 lam')
X = (x1, x2, x3)
FAIL = []


def check(c, m):
    print(("PASS " if c else "FAIL ") + m)
    if not c:
        FAIL.append(m)


def cx(s):
    v = [0, 0, 0]
    for part in s.split('+'):
        c = 1
        if part[0].isdigit():
            c, part = int(part[0]), part[1:]
        v[int(part[1]) - 1] += c
    return v


R = [("X2+X3", "2X1", sp.Rational(3, 2)), ("X1+X2", "2X2", 2), ("X1+X3", "2X3", 2), ("X1", "X2", sp.Rational(1, 2)),
     ("X1", "X3", sp.Rational(1, 2)), ("X2", "X3", 1), ("X3", "X2", 1)]
f = sp.zeros(3, 1)
for a, b, k in R:
    ya, yb = cx(a), cx(b)
    mono = x1**ya[0] * x2**ya[1] * x3**ya[2]
    for i in range(3):
        f[i] += k * mono * (yb[i] - ya[i])
f = f.applyfunc(sp.expand)
check(all(sum(cx(b)) <= 2 and sum(cx(a)) <= 2 for a, b, k in R), "at most bimolecular")
V = sp.Matrix([[cx(b)[i] - cx(a)[i] for i in range(3)] for a, b, k in R])
check(V.rank() == 2 and all(sum(V.row(r)) == 0 for r in range(V.rows)), "dimension 2, S = {x1 + x2 + x3 = 0}")
alpha = sp.Matrix([1, sp.Rational(-1, 2), sp.Rational(-1, 2)])
ell = sp.Matrix([0, 1, -1])
psi = 3 * x2 * x3 - 2 * x1 * x2 - 2 * x1 * x3 - x1
P = x2 - x3
check((f - psi * alpha - (x1 - 1) * P * ell).applyfunc(sp.expand) == sp.zeros(3, 1),
      "normal form f = psi*alpha + (x1 - 1) P ell with c = 1")
B = sp.Matrix([[1, 0], [-1, 1], [0, -1]])          # basis of S
Bp = (B.T * B).inv() * B.T
J = f.jacobian(X)
p = {x1: 1, x2: 1, x3: 3}
q = {x1: 1, x2: 3, x3: 1}
for name, pt in (("p", p), ("q", q)):
    check(f.subs(pt) == sp.zeros(3, 1), "%s = %s is a steady state" % (name, tuple(pt.values())))
    JS = Bp * J.subs(pt) * B
    tr, det = JS.trace(), JS.det()
    print("   %s: trace %s, det %s, P = %s" % (name, tr, det, P.subs(pt)))
    check(tr < 0 and det > 0, "%s exponentially stable relative to S" % name)
check(P.subs(p) * P.subs(q) < 0, "P(p) P(q) < 0")
check(sp.simplify((sp.Matrix([q[v] - p[v] for v in X])).cross(ell)) == sp.zeros(3, 1), "q - p is parallel to ell")
# third steady state of the class x1+x2+x3 = 5 on the plane P = 0 (x2 = x3)
sol = sp.solve([psi.subs(x3, x2), x1 + 2 * x2 - 5], [x1, x2], dict=True)
pos = [s for s in sol if s[x1].is_positive and s[x2].is_positive]
for s in pos:
    pt = {x1: s[x1], x2: s[x2], x3: s[x2]}
    print("   third steady state:", tuple(sp.nsimplify(v) for v in pt.values()), "~",
          tuple(round(float(v), 6) for v in pt.values()))
    check(f.subs(pt).applyfunc(sp.simplify) == sp.zeros(3, 1) and sp.simplify(pt[x1] - 1) != 0,
          "positive steady state with x1 != 1 in the same class: no ACR in X1")
check(len(pos) >= 1, "global ACR fails")
print("ALL PASS" if not FAIL else "FAILURES")
sys.exit(1 if FAIL else 0)
