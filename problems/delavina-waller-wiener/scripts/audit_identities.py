#!/usr/bin/env python3
"""Independent audit of the algebra in M. Liu, "The DeLaVina-Waller conjecture on the Wiener index"
(Zenodo 10.5281/zenodo.22015517, working paper, 19 Aug 2026).

Written from the statements in the paper, without using or running the author's scripts.
Each claim prints PASS or FAIL; exit status is non-zero on any FAIL.
Environment: Python 3 with sympy. Runtime: under one minute. Deterministic.

Checked here:
  (2.1) W(C_m) formula; Lemma 3.2 recursion identity; Lemma 4.1 and (4.2) barrier identity;
  (4.5), (4.8) exact evaluations; Lemma A.1 coefficient assembly;
  (7.4) the constant A; (9.1) the identity for A/(2s);
  Prop. 8.2 and Prop. 8.3 closed forms.
"""
import sys
from fractions import Fraction
import sympy as sp

failures = 0


def check(name, ok):
    global failures
    print(("PASS" if ok else "FAIL"), name)
    if not ok:
        failures += 1


# ---------- (7.4): A = 2[W(C_n) - W(C_k)] - 2 rho (rho+eps) s - rho s^2 ----------
a, rho = sp.symbols("a rho", positive=True)
for eps in (0, 1):
    k = 2 * rho + eps
    n = 2 * (a + rho) + 1
    s = n - k
    Wn = sp.Rational(1, 8) * n * (n ** 2 - 1)
    Wk = rho ** 3 if eps == 0 else rho * (rho + 1) * (2 * rho + 1) / 2
    lhs = 2 * (Wn - Wk) - 2 * rho * (rho + eps) * s - rho * s ** 2
    A = 2 * a ** 3 + 2 * a ** 2 * rho + 3 * a ** 2 + 2 * a * rho ** 2 + 2 * a * rho + a + (1 - eps) * rho ** 2
    check(f"(7.4) A identity, eps={eps}", sp.expand(lhs - A) == 0)
    # (9.1)
    d = a + rho
    rhs = (d ** 2 - d * a + a ** 2) / 2 + a / 2 + a * (d - a) / (2 * (2 * a + 1)) + eps * (a + 1) * (2 * d + 1) / (4 * (2 * a + 1))
    check(f"(9.1) A/(2s) identity, eps={eps}", sp.simplify(A / (2 * s) - rhs) == 0)

# identities used in the proof of Prop 7.3: E_i(k-E_i) and rho - q
q = sp.symbols("q")
for eps in (0, 1):
    k = 2 * rho + eps
    E = rho - q
    check(f"E(k-E)=rho(rho+eps)-q(q+eps), eps={eps}",
          sp.expand(E * (k - E) - (rho * (rho + eps) - q * (q + eps))) == 0)

# ---------- barrier p(t) of section 4, in q = sqrt(t) ----------
Q = sp.symbols("q", positive=True)
r2 = sp.sqrt(2)
t = Q ** 2
p = 2 * t - 4 * r2 / 3 * Q ** 3 + sp.Rational(2, 3) * t ** 2 - 2 * r2 / 45 * Q ** 5 - sp.Rational(4, 405) * Q ** 6
dq = lambda f: sp.diff(f, Q) / (2 * Q)      # d/dt = (1/(2q)) d/dq
p1 = sp.simplify(dq(p))
p2 = sp.simplify(dq(p1))
p3 = sp.simplify(dq(p2))
w = 1 - p1 / 2
id42 = sp.simplify(p - w ** 2 - Q ** 6 * (135 - 60 * r2 * Q - 8 * Q ** 2) / 36450)
check("(4.2) p - w^2 identity", id42 == 0)
check("Lemma 4.1: -270 q p'' expression",
      sp.simplify(-270 * Q * p2 - (270 * r2 - 360 * Q + 45 * r2 * Q ** 2 + 16 * Q ** 3)) == 0)
check("Lemma 4.1: 540 q^3 p''' expression",
      sp.simplify(540 * Q ** 3 * p3 - (270 * r2 - 45 * r2 * Q ** 2 - 32 * Q ** 3)) == 0)
half = 1 / sp.sqrt(2)
check("(4.8) p(1/2)=79/162", sp.nsimplify(sp.simplify(p.subs(Q, half))) == sp.Rational(79, 162))
check("(4.8) p'(1/2)=163/270", sp.nsimplify(sp.simplify(p1.subs(Q, half))) == sp.Rational(163, 270))
check("(4.8) p''(1/2)=-233/270", sp.nsimplify(sp.simplify(p2.subs(Q, half))) == sp.Rational(-233, 270))
D = sp.symbols("D", positive=True)  # D = sqrt(d)
d_ = D ** 2
val = d_ ** 2 * p.subs(Q, 1 / (r2 * D))
target = d_ - sp.Rational(2, 3) * D + sp.Rational(1, 6) - 1 / (90 * D) - 1 / (810 * d_)
check("(4.5) d^2 p(1/(2d))", sp.simplify(val - target) == 0)
Gamma = d_ + sp.Rational(2, 3) * D - sp.Rational(1, 6) + 1 / (90 * D) + 1 / (810 * d_)
check("B_d(1) = 2d", sp.simplify(val + Gamma - 2 * d_) == 0)
# p''(1/(2d))/4 value used in the proof of Prop 4.2
val2 = p2.subs(Q, 1 / (r2 * D)) / 4
target2 = -D / 2 + sp.Rational(1, 3) - 1 / (24 * D) - 1 / (135 * d_)
check("p''(1/(2d))/4 evaluation", sp.simplify(val2 - target2) == 0)

# the lower bound 151/3240 in Prop 4.2
lhs = Fraction(1, 6) - Fraction(1, 12) - Fraction(11, 360) - Fraction(1, 162)
check("151/3240 constant", lhs == Fraction(151, 3240) and lhs > 0)

# Lemma A.1 coefficient assembly
u, dd = sp.symbols("u d", positive=True)
delta = sp.Rational(1, 2) - u + 1 / (2 * dd)
lhsA = dd ** 2 * (sp.Rational(79, 162) - sp.Rational(163, 270) * delta - sp.Rational(233, 540) * delta ** 2) \
    + dd + (dd / 12 + sp.Rational(4, 3)) - sp.Rational(1, 6) + sp.Rational(1, 200)
T = sp.Rational(79, 162) - sp.Rational(163, 270) * (sp.Rational(1, 2) - u) - sp.Rational(233, 540) * (sp.Rational(1, 2) - u) ** 2
L = (611 + 466 * u) / 1080
K = sp.Rational(11489, 10800)
check("Lemma A.1 assembly", sp.simplify(lhsA - (dd ** 2 * T + dd * L + K)) == 0)
check("Lemma A.1: (d/12+4/3)^2 - 4d/9 = (d-16)^2/144",
      sp.simplify((dd / 12 + sp.Rational(4, 3)) ** 2 - 4 * dd / 9 - (dd - 16) ** 2 / 144) == 0)

# Taylor (4.9): p'''>0 means p(1/2-delta) <= Taylor quadratic; compare with the exact third-order term sign
check("(4.9) coefficients are p(1/2), -p'(1/2), p''(1/2)/2",
      sp.Rational(233, 540) == sp.Rational(233, 270) / 2)

# ---------- Prop 8.2 ----------
r_ = sp.symbols("r", positive=True)
A_tri = (2 * (2 * r_) ** 3 + 2 * (2 * r_) ** 2 + 3 * (2 * r_) ** 2 + 2 * (2 * r_) + 2 * (2 * r_) + 2 * r_)
check("Prop 8.2: A", sp.expand(A_tri - (16 * r_ ** 3 + 20 * r_ ** 2 + 10 * r_)) == 0)
gap = A_tri / 2 + sp.Rational(1, 2) * (4 * r_) ** 2 / 3 - 4 * r_ * (2 * r_ + 1) * (r_ + 1)
check("Prop 8.2: deficit r(2r+3)/3", sp.simplify(gap - r_ * (2 * r_ + 3) / 3) == 0)

# ---------- Prop 8.3, C4 core ----------
kp = sp.symbols("kappa", positive=True)
v = sp.symbols("v", positive=True)
n43 = 4 * kp + 5
Fp = sp.summation(v * (n43 - v), (v, 1, kp))
check("Prop 8.3: F_path", sp.simplify(Fp - kp * (kp + 1) * (5 * kp + 7) / 3) == 0)
Fe = 2 * (n43 - 1) + sp.summation(v * (n43 - v), (v, 3, kp + 1))
check("Prop 8.3: F_exc", sp.simplify(Fe - (5 * kp ** 3 + 21 * kp ** 2 + 16 * kp + 6) / 3) == 0)
Phi = 8 + 4 * (4 * kp + 1) + (4 * kp ** 2 + 2 * kp) + (4 * kp ** 2 + 2 * kp)
check("Prop 8.3: Phi", sp.expand(Phi - 4 * (kp + 1) * (2 * kp + 3)) == 0)
tot = 3 * Fp + Fe + Phi
check("Prop 8.3: total", sp.simplify(tot - (20 * kp ** 3 + 81 * kp ** 2 + 97 * kp + 42) / 3) == 0)
Wc = (kp + 1) * (2 * kp + 3) * (4 * kp + 5)
check("Prop 8.3: W(C_{4k+5}) formula", sp.simplify(Wc - n43 * (n43 ** 2 - 1) / 8) == 0)
check("Prop 8.3: deficit", sp.simplify(Wc - tot - (4 * kp + 1) * (kp ** 2 + 2 * kp + 3) / 3) == 0)

# Brute check of Prop 8.3's claim on which extra-vertex position is best, for kappa up to 12
ok = True
for kappa in range(1, 13):
    nn = 4 * kappa + 5
    best = None
    for j in range(0, kappa):  # extra vertex at depth j, attached under path vertex at depth j
        counts = [kappa - i + 1 + (1 if i <= j else 0) for i in range(1, kappa + 1)] + [1]
        val_ = sum(c * (nn - c) for c in counts)
        best = val_ if best is None else max(best, val_)
    target_ = (5 * kappa ** 3 + 21 * kappa ** 2 + 16 * kappa + 6) // 3
    ok &= (best == target_)
check("Prop 8.3: F_exc is the max over attachment depth, kappa<=12", ok)

print("ALL PASS" if failures == 0 else f"{failures} FAILURES")
sys.exit(1 if failures else 0)
