#!/usr/bin/env python3
"""Small explicit instance (n = 3) of the family N_n, checked with sympy in exact arithmetic.

Network N_3:  X1 -> 2X1 + X3,  2X1 + X2 -> X1 + 2X2,  2X2 + X3 -> X2
Rates kappa = (1/25, 1/22500, 1/1296); class x1 + x2 - x3 = 60.
Claim: exactly three positive equilibria on this class, (25,36,1), (50,18,8), (75,12,27), all nondegenerate.
Runtime < 5 s. Prints PASS/FAIL, exits 1 on failure.
"""
import sys
import sympy as sp
x1, x2, x3 = sp.symbols("x1 x2 x3", positive=True)
k = [sp.Rational(1, 25), sp.Rational(1, 22500), sp.Rational(1, 1296)]
v = [k[0]*x1, k[1]*x1**2*x2, k[2]*x2**2*x3]
f = sp.Matrix([v[0] - v[1], v[1] - v[2], v[0] - v[2]])          # Gamma v, Gamma columns (1,0,1), (-1,1,0), (0,-1,-1)
T = 60
fails = 0
def check(c, m):
    global fails
    print(("PASS " if c else "FAIL ") + m); fails += (not c)
# all positive solutions of f = 0 on the class: solve exactly
sols = sp.solve([v[0] - v[1], v[1] - v[2], x1 + x2 - x3 - T], [x1, x2, x3], dict=True)
pos = sorted([tuple(s[x] for x in (x1, x2, x3)) for s in sols
              if all(s[x].is_real and s[x] > 0 for x in (x1, x2, x3))])
check(pos == [(25, 36, 1), (50, 18, 8), (75, 12, 27)], f"positive equilibria on x1+x2-x3=60: {pos}")
J = f.jacobian([x1, x2, x3])
B = sp.Matrix([[1, -1], [0, 1], [1, 0]])            # basis of im Gamma = {x1 + x2 - x3 = 0}^perp ... columns (1,0,1),(-1,1,0)
for p in pos:
    Jp = J.subs({x1: p[0], x2: p[1], x3: p[2]})
    M = (B.T*B).inv()*B.T*Jp*B
    check(sp.simplify(B*M - Jp*B) == sp.zeros(3, 2), f"J maps im Gamma into itself at {p}")
    ev = [sp.N(e, 6) for e in M.eigenvals()]
    check(M.det() != 0, f"restricted det at {p} = {M.det()} (eigenvalues {ev})")
print("ALL PASS" if not fails else f"{fails} FAILURES"); sys.exit(1 if fails else 0)
