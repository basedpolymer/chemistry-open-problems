"""Verifier V2: Theorem T*' (round2/ec/notes.md §4.1): for WF-EC N and every kappa > 0, N* has a positive steady
state.  Exact for m = 2 (the root pi_1 of H_1(pi_1, 1 - pi_1) is an algebraic number; Gamma^T omega(k) is reduced
modulo its minimal polynomial and must vanish); numerical for m = 3, 4 (the cooperative flow pi' = H(pi) on the
simplex, then a check of Gamma^T omega).  Own WF-EC generator (impure/unbalanced allowed).

Environment: WSL Fedora 44, Python 3.12, sympy 1.14, numpy, scipy.  Usage: python tstar_v2.py SEED COUNT
"""
import sys, os, time, random
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa


T0 = time.time()
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
COUNT = int(sys.argv[2]) if len(sys.argv) > 2 else 20
rng = random.Random(1000 + seed)
fails = 0
exact_ok = exact_tot = 0
num_ok = num_tot = 0
unb = 0
while exact_tot + num_tot < COUNT:
    n = rng.randint(3, 5)
    m = rng.randint(2, min(4, n))
    N = gen_wfec(rng, n, m)
    if N is None:
        continue
    srcs, ct = ec_structure(N)
    m = len(srcs)
    r = len(N)
    unb += not is_balanced(N)
    rays = flux_rays(N)
    v = [sum(ray[t] for ray in rays) for t in range(r)]
    A = A_matrix(N, v)
    cp = [next(c for c in range(m) if srcs[c][p]) for p in range(n)]
    ap = [srcs[cp[p]][p] for p in range(n)]
    kap = [sp.Rational(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(n)]
    G = Gamma(N)
    if m == 2:
        exact_tot += 1
        z = sp.symbols("z")
        pi = [z, 1 - z]
        om = [kap[p] * pi[cp[p]] ** ap[p] for p in range(n)]
        H1 = sp.expand(sum(om[p] * sp.Rational(A[p][0].numerator, A[p][0].denominator) for p in range(n)))
        roots = [rt for rt in sp.Poly(H1, z).real_roots() if 0 < rt < 1]
        ok = False
        if roots:
            rt = roots[0]
            mp = sp.minimal_polynomial(rt, z)
            # omega as polynomials in z; Gamma^T omega reduced mod minimal polynomial
            res = [sp.rem(sp.expand(sum(G[p][t] * om[p] for p in range(n))), mp, z) for t in range(r)]
            ok = all(sp.simplify(e) == 0 for e in res)
        exact_ok += ok
        if not ok:
            fails += 1
            print("FAIL exact", N, kap)
    else:
        num_tot += 1
        Af = np.array([[float(x) for x in row] for row in A])
        kf = np.array([float(x) for x in kap])
        apf = np.array(ap, dtype=float)

        def H(t, pi):
            pi = np.maximum(pi, 0)
            om = kf * pi[cp] ** apf
            return om @ Af
        sol = solve_ivp(H, (0, 2000), np.full(m, 1.0 / m), method="LSODA", rtol=1e-12, atol=1e-15)
        pis = sol.y[:, -1]
        om = kf * pis[cp] ** apf
        res = np.abs(np.array(G, dtype=float).T @ om).max() / om.max()
        ok = (pis > 1e-9).all() and res < 1e-8
        num_ok += ok
        if not ok:
            fails += 1
            print("FAIL numeric", N, [str(x) for x in kap], pis, res)
print(f"T*': exact (m = 2) {exact_ok}/{exact_tot} PASS; numerical (m = 3, 4) {num_ok}/{num_tot} PASS; "
      f"unbalanced among all: {unb}")
print(("PASS" if fails == 0 else "FAIL") + f" T*' on {exact_tot + num_tot} WF-EC networks; runtime {time.time() - T0:.1f} s")
sys.exit(1 if fails else 0)
