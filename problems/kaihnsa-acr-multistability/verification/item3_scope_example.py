"""The original wording of item 3 of Section 6 was too broad: an exact counterexample of dimension 2.

Runtime: under 1 s (Python 3.12, sympy 1.14).   Usage: python item3_scope_example.py

The original text said: "No bimolecular network with 3 species and at most 5 reactions has three positive steady
states in a compatibility class".  The enumeration behind it only covers full-dimensional networks.  The
two-dimensional network with 4 reactions
    0 -> Z (k1 = 10),  Y+Z -> X (k2 = 10),  X+Z -> Y+Z (k3 = 1),  X+Y -> 2Y (k4 = 1)
has, in the class x + y = 6, the three positive steady states (5, 1, 1), (4, 2, 1/2), (3, 3, 1/3), all nondegenerate,
and the outer two are exponentially stable.  (It has no ACR, so nothing in the paper's main results is affected;
only the scope of item 3 had to be stated correctly.)
"""
import sys
import sympy as sp

x, y, z, lam = sp.symbols('x y z lam')
k1, k2, k3, k4 = 10, 10, 1, 1
f = sp.Matrix([k2 * y * z - k3 * x * z - k4 * x * y,
               -k2 * y * z + k3 * x * z + k4 * x * y,
               k1 - k2 * y * z])
# reaction vectors: (0,0,1), (1,-1,-1), (-1,1,0), (-1,1,0); S = {x + y = 0 component}: span{(0,0,1), (1,-1,0)}
V = sp.Matrix([[0, 0, 1], [1, -1, -1], [-1, 1, 0], [-1, 1, 0]])
FAIL = []


def check(c, m):
    print(("PASS " if c else "FAIL ") + m)
    if not c:
        FAIL.append(m)


check(V.rank() == 2, "dimension 2 (conservation law x + y)")
B = sp.Matrix([[0, 1], [0, -1], [1, 0]])
Bp = (B.T * B).inv() * B.T
J = f.jacobian((x, y, z))
pts = [(5, 1, 1), (4, 2, sp.Rational(1, 2)), (3, 3, sp.Rational(1, 3))]
nstable = 0
for p in pts:
    sub = dict(zip((x, y, z), p))
    check(f.subs(sub) == sp.zeros(3, 1) and p[0] + p[1] == 6, "%s is a positive steady state in x + y = 6" % (p,))
    JS = Bp * J.subs(sub) * B
    tr, det = JS.trace(), JS.det()
    st = bool(tr < 0) and bool(det > 0)
    nstable += int(st)
    print("   trace %s, det %s -> %s" % (tr, det, "exponentially stable" if st else "unstable"))
    check(det != 0, "nondegenerate at %s" % (p,))
# all positive steady states of the class: z = 1/y, cubic in y
cub = sp.factor(sp.numer(sp.together(f[0].subs({z: sp.Rational(k1, k2) / y, x: 6 - y}))))
print("   class equation:", cub)
check(nstable == 2, "two exponentially stable positive steady states in one class (no ACR)")
print("ALL PASS" if not FAIL else "FAILURES")
sys.exit(1 if FAIL else 0)
