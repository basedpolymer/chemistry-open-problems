#!/usr/bin/env python3
"""Exact exploratory rational functions; not a proof of all-order positivity.

Environment: Python 3, sympy 1.14. At most one computational thread.
Typical runtime: under a minute for six coefficients. Run from this folder.
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import sympy as s

m = s.symbols("m")
rho = [s.Integer(0)]
tail = m*(1-m)
for n in range(1, 7):
    c = s.cancel((2*tail + sum(i*(n-i)*rho[i]*rho[n-i]
                             for i in range(1, n))) / ((2*m+1)*n+1))
    residual = s.cancel(((2*m+1)*n+1)*c - 2*tail
                        - sum(i*(n-i)*rho[i]*rho[n-i] for i in range(1,n)))
    if residual != 0:
        print("FAIL recursion identity:", n, flush=True)
        raise AssertionError(n)
    rho.append(c)
    tail = s.cancel(tail-c)
    print(n, s.factor(c), flush=True)
    num, den = s.fraction(s.cancel(c/(m*(1-m))))
    print("numerator coefficients:", s.Poly(num, m).all_coeffs(), flush=True)
    print("PASS recursion identity:", n, flush=True)
