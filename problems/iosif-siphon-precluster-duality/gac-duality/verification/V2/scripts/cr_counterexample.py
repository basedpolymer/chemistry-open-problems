"""Verifier V2: counterexample to the uniqueness/stability part of Conjecture CR (round2/ec/notes.md §7.3).

NC = {C -> B, 2A -> D, 2B+D -> 2A+2B, 2B+2D -> 2A+2C+D}  (species A, B, C, D; obtained from the network found by
random_wfec_v2.py, seed 3, by deleting the pure catalyst X1 and the reaction 2C -> 2B).

Claims (exact unless stated):
 1. NC is connected, (H1), (H2), (H3); s = 2; Petri clusters {C}|{R1}, {A}|{R2}, {B,D}|{R3,R4}: q = 3, so the
    cluster rank equation s = q - 1 holds; NC is TEFC (2B+D and 2B+2D have the same support) but not PC
    (they are not proportional), and not EC.
 2. NC and NC* have no critical siphon (Theorem C^gen applies: persistence, no boundary steady state).
 3. There are rational x > 0 and k > 0 with f_k(x) = 0 and det(J(x)|_S) < 0 (= (-1)^(s+1) sign).
 4. Hence (degree argument: the class is a compact polytope with no boundary steady state, total index (-1)^s = +1,
    and x is a nondegenerate zero of index -1) the class of x contains at least two more positive steady states:
    NC is multistationary.  Conjecture CR (uniqueness and exponential stability in every class) is false.
 5. Illustration: all positive steady states in that class, by an exact resultant (real roots isolated exactly,
    then evaluated to 30 digits), and the eigenvalues of J|_S at each.

Environment: WSL Fedora 44, Python 3.12, sympy 1.14, mpmath.  PASS/FAIL, exit 1 on FAIL.
"""
import sys, os, time, random
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

T0 = time.time()
FAILS = []


def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + name + ((" -- " + info) if info else ""))
    if not cond:
        FAILS.append(name)


NC = parse("C->B, 2A->D, 2B+D->2A+2B, 2B+2D->2A+2C+D", ["A", "B", "C", "D"])
n, r = 4, 4
s = rank(Gamma(NC))
cl = clusters(NC)
check("1. connected, (H1), (H2), (H3)", connected(NC) and conservative(NC) and consistent(NC) and rows_distinct(NC))
check("1. s = 2, q = 3 clusters {C},{A},{B,D}: cluster rank equation s = q - 1", s == 2 and len(cl) == 3,
      str([(sorted(S), sorted(T)) for S, T in cl]))
check("1. not PC, not EC", not is_PC(NC) and not is_EC(NC))
check("2. no critical siphon in NC nor NC*", not critical_siphons(NC) and not critical_siphons(dual(NC)))
print("   conservation laws:", [[str(t) for t in z] for z in cons_rays(NC)], " flux rays:",
      [[str(t) for t in v] for v in flux_rays(NC)])
print("   NC* =", dual(NC))

# 3. search a rational steady state with det(J|_S) < 0
rng = random.Random(5)
rays = flux_rays(NC)
best = None
for trial in range(20000):
    x = [F(rng.randint(1, 12), rng.randint(1, 12)) for _ in range(n)]
    wts = [F(rng.randint(1, 12), rng.randint(1, 12)) for _ in rays]
    v = [sum(w * ray[t] for w, ray in zip(wts, rays)) for t in range(r)]
    k = [v[t] / monomial(x, NC[t][0]) for t in range(r)]
    cp = charpoly(jacobian(NC, k, x))
    q = strip_zero_roots(cp, n - s)
    detJS = q[-1]  # (-1)^2 q(0)
    if detJS < 0:
        size = sum(t.numerator + t.denominator for t in x + k)
        if best is None or size < best[0]:
            best = (size, x, k, detJS)
check("3. a rational steady state with det(J|_S) < 0 exists", best is not None)
_, x, k, detJS = best
check("3. f_k(x) = 0 exactly", all(e == 0 for e in field(NC, k, x)))
print("   x =", [str(t) for t in x], " k =", [str(t) for t in k], " det(J|_S) =", str(detJS))
T1 = x[1] + x[2]
T2 = x[0] + 2 * x[3]
print(f"   class: B + C = {T1}, A + 2D = {T2}")
check("4. degree argument: >= 2 further positive steady states in the class (index -1 at x, total +1)", detJS < 0)

# 5. all steady states in the class: P1 = k1 C - 2 k4 B^2 D^2 (B-balance), P2 = k2 A^2 - k3 B^2 D - k4 B^2 D^2
B, D = sp.symbols("B D")
kk = [sp.Rational(t.numerator, t.denominator) for t in k]
T1s, T2s = sp.Rational(T1.numerator, T1.denominator), sp.Rational(T2.numerator, T2.denominator)
C_ = T1s - B
A_ = T2s - 2 * D
P1 = sp.expand(kk[0] * C_ - 2 * kk[3] * B ** 2 * D ** 2)
P2 = sp.expand(kk[1] * A_ ** 2 - kk[2] * B ** 2 * D - kk[3] * B ** 2 * D ** 2)
Rz = sp.Poly(sp.resultant(P1, P2, B), D)
Droots = [rt for rt in Rz.real_roots() if 0 < rt < T2s / 2]
sols = []
for rt in Droots:
    Dv = sp.N(rt, 40)
    # P1 is quadratic in B with exactly one positive root, which lies in (0, T1)
    a2 = 2 * kk[3] * Dv ** 2
    Bv = (-kk[0] + sp.sqrt(kk[0] ** 2 + 4 * a2 * kk[0] * T1s)) / (2 * a2)
    res = abs(sp.N(P2.subs({B: Bv, D: Dv}), 40))
    if res < sp.Float("1e-25"):
        sols.append((Bv, Dv))
check("5. exactly three positive steady states in the class (resultant roots, verified to 25 digits)",
      len(sols) == 3, f"{len(Droots)} resultant roots in (0, T2/2)")
for Bv, Dv in sols:
    xv = [T2s - 2 * Dv, Bv, T1s - Bv, Dv]
    J = sp.Matrix(n, n, lambda i, j: 0)
    for t, (a, b) in enumerate(NC):
        w = kk[t] * sp.prod([xv[i] ** a[i] for i in range(n)])
        for p in range(n):
            if a[p]:
                for i in range(n):
                    J[i, p] += (b[i] - a[i]) * w * a[p] / xv[p]
    import numpy as np
    ev = list(np.linalg.eigvals(np.array(J.evalf(30).tolist(), dtype=float)))
    print("   steady state (A,B,C,D) =", [f"{float(t):.10g}" for t in xv], " eigenvalues:",
          sorted([f"{e.real:.4g}{e.imag:+.4g}j" for e in ev]))
print(f"runtime {time.time() - T0:.1f} s; {len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
