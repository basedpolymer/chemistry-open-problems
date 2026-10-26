#!/usr/bin/env python3
"""Exact Hurwitz proof: no Hopf for 3 <= M <= 14 and ANY K > 0.

Python 3.12 or later / SymPy. About a minute. Jacobians are derived from the original
quadratic vector field, restricted to the mass hyperplane. Each leading
Hurwitz determinant is a polynomial in K. Positive coefficients suffice;
otherwise exact real-root isolation and a positive value at K=1 are required.
"""
from pathlib import Path
import sys
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "verification"))
from v02_dense_certify import reduced_int

def main():
    K, x = sp.symbols("K x")
    for M in range(3, 15):
        A0, A1, _, _ = reduced_int(M)
        P = (sp.Matrix(A0) + K * sp.Matrix(A1)).charpoly(x).as_expr()
        a = sp.Poly(P, x).all_coeffs()
        n = M - 1
        H = sp.Matrix(n, n, lambda i, j: a[2*j-i+1] if 0 <= 2*j-i+1 <= n else 0)
        for k in range(1, n + 1):
            d = sp.Poly(H[:k, :k].det(method="domain-ge"), K)
            positive_coefficients = all(c >= 0 for c in d.all_coeffs()) and not d.is_zero
            ok = positive_coefficients or (d.eval(1) > 0 and d.count_roots(0, sp.oo) == 0)
            if not ok:
                print("FAIL", M, k, flush=True)
                return 1
        print(f"PASS M={M}: all {n} Hurwitz determinants positive for K>0", flush=True)
    print("SUMMARY PASS: M=3..14, every K>0, exact polynomial arithmetic", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
