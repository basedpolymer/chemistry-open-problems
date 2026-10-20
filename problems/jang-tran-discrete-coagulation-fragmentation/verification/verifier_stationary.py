#!/usr/bin/env python3
"""Independent verifier checks of the stationary state (recursion (1.9)) -- Theorems 2 and 3.

Written by the independent verifier; does not reuse the solvers' code.

  S1  exact rational arithmetic: rho~_n > 0 and number remainder N - sum_{j<=n} rho~_j > 0 for n <= 40,
      at m = 1/2, 51/100, 3/4, 9/10, 99/100, 999/1000; and rho~ == 0 at m = 1.
  S2  Arb ball arithmetic (python-flint), n <= 1500, m = 1/2 + 1/1024, 3/4, 19/20: rho~_n > 0, remainders
      > 0, partial mass < m (rigorous FINITE statements only).
  S3  evidence only: the radius R(m) is computed independently by integrating the regularised ODE for
      r = (m - P)/(1 - z) across z = 1, then the inverse function z(P) up to the turning point P = m + 1/2;
      with Arb coefficients (m = 0.6, n <= 3000; m = 0.75, n <= 6000) checks rho~_n R^n n^{5/2} -> b/(2 sqrt pi)
      (b of Theorem 3) with an O(1/n) error, and sum n rho~_n (partial sum + asymptotic tail) = m;
      m = 0.9 in floating point (n <= 30000).

Environment: Python 3, numpy, scipy, python-flint. Runtime: about 20-60 seconds, one thread.
Deterministic. Prints PASS/FAIL per check, exits 1 on failure.
"""
import os
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
import sys
from fractions import Fraction as Q
import numpy as np
from scipy.integrate import solve_ivp
from flint import arb, ctx

FAILS = 0


def report(name, ok, info=""):
    global FAILS
    print(("PASS " if ok else "FAIL ") + name + (("  [" + info + "]") if info else ""), flush=True)
    if not ok:
        FAILS += 1


def recursion_generic(m, n_max, zero, conv_fn=None):
    """rho_n from (1.9) written as ((2m+1)n+1) rho_n = 2 * (remaining number) + sum i(n-i) rho_i rho_{n-i}."""
    N = m * (1 - m)
    rho = [zero]
    rem = N
    rems = []
    for n in range(1, n_max + 1):
        conv = zero
        for i in range(1, n):
            conv += i * (n - i) * rho[i] * rho[n - i]
        c = (2 * rem + conv) / ((2 * m + 1) * n + 1)
        rho.append(c)
        rem = rem - c
        rems.append(rem)
    return rho[1:], rems


# ---------------------------------------------------------------- S1 exact
for m in [Q(1, 2), Q(51, 100), Q(3, 4), Q(9, 10), Q(99, 100), Q(999, 1000)]:
    rho, rems = recursion_generic(m, 40, Q(0))
    # consistency with the literal form of (1.9)
    lit_ok = True
    for n in range(1, 41):
        num = 2 * m * (1 - m) + sum(i * (n - i) * rho[i - 1] * rho[n - i - 1] for i in range(1, n)) \
            - 2 * sum(rho[:n - 1])
        if num / ((2 * m + 1) * n + 1) != rho[n - 1]:
            lit_ok = False
    pm = sum((n + 1) * r for n, r in enumerate(rho))
    report("S1 m=%s: (1.9) exact, rho_n > 0 and remainders > 0 for n <= 40, partial mass < m" % m,
           lit_ok and all(r > 0 for r in rho) and all(r > 0 for r in rems) and pm < m,
           "partial mass %.6f" % float(pm))
rho1, _ = recursion_generic(Q(1), 30, Q(0))
report("S1 m=1: (1.9) gives rho == 0 (no stationary state of mass 1)", all(r == 0 for r in rho1))

# ---------------------------------------------------------------- S2 Arb
ctx.prec = 300
for (p, q) in [(513, 1024), (3, 4), (19, 20)]:
    m = arb(p) / q
    nmax = 1500
    rho = [arb(0)] * (nmax + 1)
    rem = m * (1 - m)
    pm = arb(0)
    okpos, okrem = True, True
    for n in range(1, nmax + 1):
        conv = arb(0)
        for i in range(1, n // 2 + 1):
            t = i * (n - i) * rho[i] * rho[n - i]
            conv += t if 2 * i == n else 2 * t
        c = (2 * rem + conv) / ((2 * m + 1) * n + 1)
        rho[n] = c
        rem -= c
        pm += n * c
        okpos = okpos and bool(c > 0)
        okrem = okrem and bool(rem > 0)
    report("S2 m=%d/%d (Arb, 300 bits): rho_n > 0, remainders > 0 for n <= %d, partial mass < m" % (p, q, nmax),
           okpos and okrem and bool(pm < m), "m - partial mass ~ %.3e" % float((m - pm).mid()))

# ---------------------------------------------------------------- S3 floating point
def recursion_float(m, n_max):
    rho = np.zeros(n_max + 1)
    a = np.zeros(n_max + 1)          # a_i = i rho_i
    rem = m * (1 - m)
    for n in range(1, n_max + 1):
        conv = np.dot(a[1:n], a[n - 1:0:-1]) if n > 1 else 0.0
        c = (2 * rem + conv) / ((2 * m + 1) * n + 1)
        rho[n] = c
        a[n] = n * c
        rem -= c
    return rho


def radius_from_ode(m, rho):
    """Integrate r' = (m - r + 2z(z+2)r^2)/(z(1+z)(1+2(1-z)r)) from z0 = 0.5 (series start) to the
    turning point 1 + 2(1-z) r = 0."""
    z0 = 0.5
    n = np.arange(len(rho))
    P0 = np.sum(n * rho * z0 ** n)
    r0 = (m - P0) / (1 - z0)

    def f(z, y):
        r = y[0]
        return [(m - r + 2 * z * (z + 2) * r * r) / (z * (1 + z) * (1 + 2 * (1 - z) * r))]

    # stop well before the turning point, where a - 2P = 1 + 2(1-z) r = 0.05 (this happens for z > 1)
    def ev(z, y):
        return 1 + 2 * (1 - z) * y[0] - 0.05
    ev.terminal = True
    sol = solve_ivp(f, (z0, 50.0), [r0], method='LSODA', rtol=1e-13, atol=1e-15, events=ev)
    z1 = sol.t_events[0][0]
    r1 = sol.y_events[0][0][0]
    P1 = m - (1 - z1) * r1
    a = 2 * m + 1

    # then integrate the inverse function z(P), regular up to P = a/2:  dz/dP = z(1-z^2)(a-2P)/D(z,P)
    def g(P, y):
        z = y[0]
        D = 2 * z * (m - P) * (P + 1 - m) - (1 - z) ** 2 * P
        return [z * (1 - z * z) * (a - 2 * P) / D]
    sol2 = solve_ivp(g, (P1, a / 2), [z1], method='LSODA', rtol=1e-13, atol=1e-15)
    return sol2.y[0, -1]


def recursion_arb(m_arb, nmax):
    rho = [arb(0)] * (nmax + 1)
    rem = m_arb * (1 - m_arb)
    for n in range(1, nmax + 1):
        conv = arb(0)
        for i in range(1, n // 2 + 1):
            t = i * (n - i) * rho[i] * rho[n - i]
            conv += t if 2 * i == n else 2 * t
        rho[n] = (2 * rem + conv) / ((2 * m_arb + 1) * n + 1)
        rem -= rho[n]
    return rho, rem


# Arb coefficients (no cancellation problem) for the asymptotic and mass checks at m = 0.6, 0.75
ctx.prec = 600
for (p, q, nmax) in [(3, 5, 3000), (3, 4, 6000)]:
    m = p / q
    rho_a, rem_a = recursion_arb(arb(p) / q, nmax)
    ok_rad = all(bool(r.rad() < 1e-30 * abs(r.mid())) for r in rho_a[1:])
    rho = np.array([0.0] + [float(r.mid()) for r in rho_a[1:]])
    R = radius_from_ode(m, rho[:400])
    b = np.sqrt(((R - 1) ** 2 * (m + 0.5) + 1.5 * R) / (R * R - 1))
    const = b / (2 * np.sqrt(np.pi))
    nn = np.arange(nmax + 1, dtype=float)
    scaled = rho[1:] * np.exp(nn[1:] * np.log(R)) * nn[1:] ** 2.5 / const
    e_half, e_end = scaled[nmax // 2 - 1] - 1, scaled[-1] - 1
    richardson = 2 * scaled[-1] - scaled[nmax // 2 - 1] - 1      # removes the O(1/n) term
    report("S3 m=%.2f: rho_n R^n n^(5/2) -> b/(2 sqrt pi) with error ~ c/n (R from the r-ODE, Arb coefficients)" % m,
           ok_rad and abs(e_end) < 5e-2 and 1.6 < e_half / e_end < 2.4 and abs(richardson) < 5e-3,
           "R=%.10f, ratio-1 at n=%d: %.3e, at n=%d: %.3e, Richardson: %.1e"
           % (R, nmax // 2, e_half, nmax, e_end, richardson))
    k = np.arange(nmax + 1, nmax + 400000, dtype=float)
    tail = const * np.sum(k ** -1.5 * np.exp(-k * np.log(R)))
    pm = sum((n * rho_a[n] for n in range(1, nmax + 1)), arb(0))
    mass = float(pm.mid()) + tail
    report("S3 m=%.2f: sum_{n<=%d} n rho~_n + asymptotic tail = m; number remainder > 0" % (m, nmax),
           abs(mass - m) < 1e-9 and bool(rem_a > 0),
           "mass-m=%.2e, tail used=%.2e" % (mass - m, tail))
    if p == 3 and q == 4:
        big = float(sum((n * rho_a[n] for n in range(3001, nmax + 1)), arb(0)).mid()) + tail
        report("S3 m=0.75: mass carried by sizes > 3000 is ~2e-9 (paper, Section 4)", 1e-9 < big < 4e-9,
               "%.3e" % big)

# m = 0.9: floating point only (the stationary tail extends to sizes ~ 1/(R-1) ~ 4000)
m = 0.9
nmax = 30000
rho = recursion_float(m, nmax)
R = radius_from_ode(m, rho[:400])
b = np.sqrt(((R - 1) ** 2 * (m + 0.5) + 1.5 * R) / (R * R - 1))
const = b / (2 * np.sqrt(np.pi))
nn = np.arange(nmax + 1, dtype=float)
k = np.arange(nmax + 1, nmax + 400000, dtype=float)
tail = const * np.sum(k ** -1.5 * np.exp(-k * np.log(R)))
mass = np.dot(nn, rho) + tail
report("S3 m=0.90 (floating point): sum_{n<=30000} n rho~_n + asymptotic tail = m within 1e-6",
       abs(mass - m) < 1e-6 and rho[1:].min() > 0,
       "R=%.10f, mass-m=%.2e, tail used=%.2e" % (R, mass - m, tail))

print("\nSUMMARY:", "ALL PASS" if FAILS == 0 else "%d FAIL" % FAILS)
sys.exit(1 if FAILS else 0)
