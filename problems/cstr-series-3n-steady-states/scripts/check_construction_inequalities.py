#!/usr/bin/env python3
"""Rigorous check (ball arithmetic, python-flint/Arb) of the numerical inequalities behind Theorem 2.

Explicit family (manuscript, Theorem 2), for n >= 1 and eps > 0, m = max(1, eps):
  reduced parameters  X = 8, K = 1/(4n), lambda = 1/(200 m), z_c = 1/20, Delta = X/(lambda K) = 6400 n m,
  window points       a = 3/20, b = 3,  inlet bounds c_- = (1+K)^-(n-1), zeta_+ = max(eps, 9).
Conditions of the window lemma (Lemma 4) and why they hold uniformly in n and eps:
  (W0) zeta_+ >= eps and (z_c + X/(1+K))/(1-lambda) <= zeta_+     [uses K > 0, lambda <= 1/200]
  (W1) lambda*zeta_+ + z_c + X e^{-1/a} < a                          [uses lambda*zeta_+ <= 9/200]
  (W2) X c_- e^{-1/b} > (1+K) b  <=  X e^{-1/b} > b e^{1/4}           [uses (1+K)^n <= e^{1/4}]
  (W3) (1+K) X psi(a) < 1 and (1+K) X psi(b) < 1, psi(z) = e^{-1/z}/z^2  [uses 1+K <= 5/4]
The reductions in brackets are elementary (proved in the manuscript); this script checks the resulting
n- and eps-free numerical inequalities with Arb, checks the reductions themselves on sample values
(n = 1..2000, several eps), and prints the thesis parameters of the family for n = 3, 4 and eps = 1/20.
Each claim prints PASS/FAIL; exit code 1 on any FAIL.
Environment: python-flint 0.9 (Arb), WSL Fedora 44 ~/envs/sci. Runtime: under 5 s.
"""
import sys
from flint import arb, fmpq, ctx

ctx.prec = 200
fails = 0


def report(name, ok):
    global fails
    print(("PASS " if ok else "FAIL ") + name)
    if not ok:
        fails += 1


def A(p, q=1):
    return arb(fmpq(p, q))


X, zc, a, b = A(8), A(1, 20), A(3, 20), A(3)
psi = lambda t: (-1/t).exp()/t**2

# ---- n- and eps-free inequalities ----
lhs_W0 = (zc + X)/(1 - A(1, 200))
report("W0  (1/20 + 8)/(1 - 1/200) = %s <= 9" % lhs_W0.str(10), lhs_W0 < 9)
lhs_W1 = A(9, 200) + zc + X*(-1/a).exp()
report("W1  9/200 + 1/20 + 8 e^{-20/3} = %s < 3/20" % lhs_W1.str(10), lhs_W1 < a)
l2, r2 = X*(-1/b).exp(), b*A(1, 4).exp()
report("W2  8 e^{-1/3} = %s > 3 e^{1/4} = %s" % (l2.str(10), r2.str(10)), l2 > r2)
s_a = A(5, 4)*X*psi(a)
report("W3a (5/4) 8 psi(3/20) = %s < 1" % s_a.str(10), s_a < 1)
s_b = A(5, 4)*X*psi(b)
report("W3b (5/4) 8 psi(3) = %s < 1" % s_b.str(10), s_b < 1)
report("W3c window points satisfy 0 < a < 1/2 < b (psi increasing on (0,1/2], decreasing on [1/2,oo))",
       (a > 0) and (a < A(1, 2)) and (b > A(1, 2)))

# ---- the elementary reductions, checked on samples (they are proved in the manuscript) ----
ok = True
for n in range(1, 2001):
    K = A(1, 4*n)
    if not ((1 + K)**n <= A(1, 4).exp()):
        ok = False
        break
report("R1  (1 + 1/(4n))^n <= e^{1/4} for n = 1..2000 (proof: n log(1+x/n) <= x)", ok)
ok = True
for eps in [fmpq(1, 1000), fmpq(1, 20), fmpq(1, 2), fmpq(1), fmpq(3, 2), fmpq(9), fmpq(10), fmpq(1000)]:
    m = max(eps, fmpq(1))
    lam = 1/(200*m)                      # exact rationals
    zp = max(eps, fmpq(9))
    ok &= (lam*zp <= fmpq(9, 200)) and (zp >= eps) and (lam <= fmpq(1, 200))
report("R2  lambda*zeta_+ <= 9/200, lambda <= 1/200, zeta_+ >= eps (exact, sample eps in [1/1000, 1000])", ok)
ok = True
for n in [1, 2, 3, 4, 10, 100, 2000]:
    K = A(1, 4*n)
    cminus = (1 + K)**(-(n - 1))
    ok &= bool(X*cminus*(-1/b).exp() > (1 + K)*b)
    ok &= bool((1 + K)*X*psi(a) < 1) and bool((1 + K)*X*psi(b) < 1)
report("R3  (W2),(W3) in their original n-dependent form for n in {1,2,3,4,10,100,2000}", ok)

# ---- thesis parameters of the family ----
print("\nThesis parameters (Bisbis eqs. (5)-(8), all reactors identical, tau = 1, alpha_k = 1):")
for n in (3, 4):
    for eps in (fmpq(1, 20), fmpq(1), fmpq(5)):
        m = max(eps, fmpq(1))
        Bv = 6400*n*m/eps**2
        beta = 200*m - 1
        etav = (1 - 20*eps)/(20*eps**2)
        Da0 = (-1/arb(eps)).exp()/(4*n)
        print("  n=%d eps=%s: B=%s, Da0=e^{-1/eps}/%d=%s, beta_k=%s, eta_k=%s"
              % (n, eps, Bv, 4*n, Da0.str(8), beta, etav))
print("\nSUMMARY: %s (%d failure(s))" % ("ALL PASS" if fails == 0 else "FAILURES", fails))
sys.exit(1 if fails else 0)
