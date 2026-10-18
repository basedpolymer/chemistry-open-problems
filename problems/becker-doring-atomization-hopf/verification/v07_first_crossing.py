#!/usr/bin/env python3
"""v07_first_crossing.py -- exact (PARI/GP) location of ALL purely imaginary eigenvalue crossings for small M.

For the reduced Jacobian A_hat(K) = A0 + K A1 on the mass hyperplane (exact integer matrices, built from the
vector field by v02_dense_certify.reduced_int), let P(lam, K) = charpoly(A_hat(K)) in Z[K][lam].  Writing
P(i y) = Pe(y^2) + i y Po(y^2), a nonzero purely imaginary eigenvalue at K requires Res_u(Pe(u), Po(u)) (K) = 0.
Moreover P(0, K) = det A_hat(K) != 0 for K > 0 (checked: no positive real root).  Hence the number of eigenvalues
with Re > 0 is constant on each interval of (0, inf) between consecutive positive real roots of the resultant.
The script prints all positive real roots K* (rigorous real-root isolation, polrootsreal) and, for each, whether
the corresponding common root u = y^2 is positive (a genuine imaginary pair) -- by evaluating the eigenvalues of
A_hat(K*) at 100 digits -- and the number of unstable eigenvalues on each interval (eigenvalues at an interior
rational point, 100 digits; their real parts are far from 0 there).
Claims: M=14: no imaginary crossing for any K>0 and stable everywhere; M=15,16: exactly one crossing K_1(M),
stable on (0, K_1) => the Hopf point of Theorem A is the first loss of stability (and the only crossing).
Usage: python v07_first_crossing.py M [M ...]      Environment: WSL Fedora, PARI/GP 2.15+.
"""
import sys, os, subprocess, tempfile, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v02_dense_certify import reduced_int

GP = r"""
default(parisize, 400000000); default(realprecision, 100);
A0 = %s; A1 = %s; n = matsize(A0)[1];
Ah = A0 + K*A1;
P = charpoly(Ah, x);
c(k) = polcoeff(P, k, x);
Pe = sum(k=0, n\2, c(2*k)*(-1)^k*u^k);
Po = sum(k=0, (n-1)\2, c(2*k+1)*(-1)^k*u^k);
R = polresultant(Pe, Po, u);
R = R/content(R);
d0 = subst(P, x, 0);
r0 = polrootsreal(d0); np0 = #select(z->z>0, r0);
rr = polrootsreal(R); pos = select(z->z>0, rr);
print("M=", n+1, " deg_K(resultant)=", poldegree(R), " positive roots of det(A_hat)=", np0, " positive real roots of resultant: ", #pos);
cnt(Kv) = { my(E = polroots(subst(P, K, Kv))); [#select(z->real(z)>0, E), vecmin(apply(z->abs(real(z)), E))] };
genuine(Kv) = { my(E = polroots(subst(P, K, Kv)), i0 = 1); for(i=2,#E, if(abs(real(E[i])) < abs(real(E[i0])), i0 = i)); [abs(real(E[i0])), abs(imag(E[i0]))] };
for(i=1, #pos, print("  K* = ", pos[i], "   eigenvalue closest to iR at K*: [|Re|, |Im|] = ", genuine(pos[i])));
"""

def routh_count(M, Kq):
    """exact number of eigenvalues of A_hat(Kq) with Re > 0 (Routh array over Q; Kq rational)"""
    import sympy as sp
    A0, A1, _, _ = reduced_int(M)
    Kq = sp.Rational(Kq)
    Ah = sp.Matrix(A0) + Kq*sp.Matrix(A1)
    x = sp.symbols("x")
    c = sp.Poly(Ah.charpoly(x).as_expr(), x).all_coeffs()        # highest degree first
    n = len(c) - 1
    r1 = c[0::2]; r2 = c[1::2]
    col = [r1[0], r2[0]]
    rows = [r1, r2]
    while len(rows) < n + 1:
        a, b = rows[-2], rows[-1]
        if b[0] == 0:
            return None                                             # degenerate Routh array (not expected)
        new = [(b[0]*(a[i + 1] if i + 1 < len(a) else 0) - a[0]*(b[i + 1] if i + 1 < len(b) else 0))/b[0]
               for i in range(max(len(a) - 1, 1))]
        rows.append(new); col.append(new[0])
    return sum(1 for u, v in zip(col, col[1:]) if (u > 0) != (v > 0)), c[-1] != 0

def mat_gp(Mx):
    return "[" + ";".join(",".join(str(v) for v in row) for row in Mx) + "]"

def main():
    for M in [int(a) for a in sys.argv[1:]]:
        t0 = time.time()
        A0, A1, _, _ = reduced_int(M)
        src = GP % (mat_gp(A0), mat_gp(A1))
        with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False) as fh:
            fh.write(src + "\nquit;\n")
            fn = fh.name
        out = subprocess.run(["gp", "-q", "-f", fn], capture_output=True, text=True, timeout=36000)
        os.unlink(fn)
        print(out.stdout.strip())
        if out.stderr.strip():
            print("stderr:", out.stderr.strip()[:2000])
        # exact Routh-Hurwitz counts of unstable eigenvalues on each interval between the crossings
        import re
        ks = [float(v) for v in re.findall(r"K\* = ([0-9.eE+-]+)", out.stdout)]
        from fractions import Fraction
        pts = [0.0] + ks
        for i, a in enumerate(pts):
            b = pts[i + 1] if i + 1 < len(pts) else None
            m = Fraction(a + b).limit_denominator(1000)/2 if b is not None else Fraction(2*a + 1).limit_denominator(1000)
            cnt = routh_count(M, m)
            print("  interval (%s, %s): exact Routh count of eigenvalues with Re>0 at K = %s : %s"
                  % (a, b if b is not None else "inf", m, cnt[0] if cnt else "degenerate"))
        print("  (%.1fs)" % (time.time() - t0), flush=True)

if __name__ == "__main__":
    main()
