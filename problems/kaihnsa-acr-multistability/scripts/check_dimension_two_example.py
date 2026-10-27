"""Illustration of the dimension-two theorem (Theorem 15 of paper.pdf) on the KNS network of Example 4.10.

Network {X1+X2 -> 2X3 (k1), X3 -> X1 (k2), 2X3 -> 2X2 (k3)}: bimolecular, dimension 2, ACR in X3 (x3 = k2/(2 k3)).
For random rational rate constants and totals T, the class x1+x2+x3 = T contains two positive steady states; the
script checks with exact arithmetic (sympy) that the determinant of the Jacobian restricted to S has opposite signs
at the two, i.e. one of them is a saddle, as that theorem predicts. It also checks the determinant formula of
Lemma 10 (determinant on the ACR hyperplane), det = -P(x) * (grad psi . ell), in the normal form
f = psi*alpha + (x3 - c) P ell.
This is a sanity check of the statements, not part of their proof.
Runtime: a few seconds (Python 3.12, sympy 1.14).
"""
import random, sys
import sympy as sp

x1, x2, x3 = sp.symbols('x1 x2 x3', positive=True)
ok = True
random.seed(1)
for trial in range(20):
    k1, k2, k3 = [sp.Rational(random.randint(1, 40), random.randint(1, 10)) for _ in range(3)]
    f = sp.Matrix([-k1*x1*x2 + k2*x3, -k1*x1*x2 + 2*k3*x3**2, 2*k1*x1*x2 - k2*x3 - 2*k3*x3**2])
    c = k2/(2*k3)
    K = k2*c/k1                      # x1*x2 at steady state
    T = c + 2*sp.sqrt(K) + sp.Rational(random.randint(1, 20), 3)   # large enough for two positive roots
    # steady states: x3 = c, x1 + x2 = T - c, x1 x2 = K
    s = T - c
    roots = [(s + sp.sqrt(s**2 - 4*K))/2, (s - sp.sqrt(s**2 - 4*K))/2]
    # basis of S (x1+x2+x3 conserved): sigma1 = (1,0,-1), sigma2 = (0,1,-1)
    B = sp.Matrix([[1, 0], [0, 1], [-1, -1]])
    J = f.jacobian([x1, x2, x3])
    dets = []
    for r1 in roots:
        pt = {x1: r1, x2: s - r1, x3: c}
        assert all(sp.simplify(e.subs(pt)) == 0 for e in f)
        M = (B.T*B).inv()*B.T*J.subs(pt)*B          # matrix of Df|_S in the basis of S
        dets.append(sp.nsimplify(sp.simplify(M.det())))
    opposite = sp.N(dets[0]*dets[1]) < 0
    # normal form: ACR species X3 plays the role of X1; alpha = f/h on H, ell spans S ∩ e3^perp = (1,-1,0)
    fH = f.subs(x3, c)
    h = sp.simplify(fH[2])                          # psi = f_3 (alpha_3 = 1 after normalization)
    alpha = sp.Matrix([sp.simplify(fH[i]/h) for i in range(3)])
    ell = sp.Matrix([1, -1, 0])
    psi = f[2]/alpha[2]
    phi_vec = sp.simplify(f - psi*alpha)            # = (x3 - c) * P * ell
    P = sp.simplify(phi_vec[0]/((x3 - c)*ell[0]))
    formula_ok = True
    for r1, d in zip(roots, dets):
        pt = {x1: r1, x2: s - r1, x3: c}
        grad = sp.Matrix([sp.diff(psi, v) for v in (x1, x2, x3)]).subs(pt)
        pred = -P.subs(pt)*(grad.dot(ell))
        # det in basis (alpha, ell) vs basis (sigma1, sigma2): same value (basis independent)
        formula_ok &= sp.simplify(pred - d) == 0
    ok &= opposite and formula_ok
    print("k=(%s,%s,%s) T=%s  det signs %s  formula %s" % (k1, k2, k3, T, [sp.sign(sp.N(d)) for d in dets],
          "PASS" if formula_ok else "FAIL"))
print("opposite determinant signs and determinant formula (Lemma 10) on 20 random instances:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
