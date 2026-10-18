"""
polar_check.py -- independent (non-rigorous, floating-point) cross-check of the focal values.

Method (different from the Lyapunov-function recursion of focal_mp.py / focal_exact.py):
for the time-normalised system X' = -Y + P, Y' = X + Q, in polar coordinates
    dr/dtheta = sum_{k>=2} C_k(theta) r^k,
the solution with r(0) = r0 is r(theta) = r0 + sum_{k>=2} u_k(theta) r0^k, and the displacement map
on the ray theta = 0 is d(r0) = sum_k u_k(2 pi) r0^k.  The u_k solve a triangular ODE system which we
integrate with scipy's DOP853 (rtol = atol = 1e-13).  If V_1 = ... = V_{k-1} = 0 then
u_{2k+1}(2 pi) = pi V_k and u_j(2 pi) = 0 for j <= 2k.

Usage: python polar_check.py            (runtime ~ 1 min; environment: Python 3.12, numpy, scipy, mpmath)
"""
import sys
import os
import numpy as np
from scipy.integrate import solve_ivp
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from focal_mp import taylor_monomials, normalise, to_complex, lyapunov, bh_system  # noqa: E402


def displacement_coeffs(P, Q, M):
    """Return u_k(2 pi), k = 2..M, for X' = -Y + P, Y' = X + Q (P, Q dicts of floats)."""
    Pn = {}
    Qn = {}
    for (i, j), c in P.items():
        Pn.setdefault(i + j, []).append((i, j, float(c)))
    for (i, j), c in Q.items():
        Qn.setdefault(i + j, []).append((i, j, float(c)))

    def rhs(th, uvec):
        co, si = np.cos(th), np.sin(th)
        A = np.zeros(M + 1)   # A[n]: coefficient of r^n in r'
        B = np.zeros(M + 1)   # B[m]: coefficient of r^m in theta'
        B[0] = 1.0
        for n in range(2, M + 1):
            pn = sum(c * co**i * si**j for (i, j, c) in Pn.get(n, []))
            qn = sum(c * co**i * si**j for (i, j, c) in Qn.get(n, []))
            A[n] = co * pn + si * qn
            if n - 1 <= M:
                B[n - 1] = co * qn - si * pn
        # C = A / B as power series in r, up to r^M
        C = np.zeros(M + 1)
        for k in range(M + 1):
            s = A[k] - sum(C[j] * B[k - j] for j in range(k))
            C[k] = s / B[0]
        # r(theta) as series in r0: rs[1] = 1, rs[k] = u_k
        rs = np.zeros(M + 1)
        rs[1] = 1.0
        rs[2:] = uvec
        powk = rs.copy()
        out = np.zeros(M + 1)
        for k in range(2, M + 1):
            powk = np.convolve(powk, rs)[:M + 1]
            out += C[k] * powk
        return out[2:]

    sol = solve_ivp(rhs, (0.0, 2 * np.pi), np.zeros(M - 1), method='DOP853', rtol=1e-13, atol=1e-13)
    return sol.y[:, -1]


def normalised_PQ(a, b, d, N):
    tf, tg = bh_system(a, b, d)
    Tf = taylor_monomials(tf, N)
    Tg = taylor_monomials(tg, N)
    P, Q, om, tr, A = normalise(Tf, Tg, N)
    return P, Q


def run():
    mp.mp.dps = 40
    ok = True
    M = 9
    # point with L1 = 0 exactly (closed form of the paper), generic otherwise
    a, d = mp.mpf('1.1'), mp.mpf('2.7')
    b = (-2 + a * (1 + mp.sqrt(1 + 8 * d))) / (2 * d)
    cases = [("generic (a,b,d)=(1.2,0.3,2.5)", (mp.mpf('1.2'), mp.mpf('0.3'), mp.mpf('2.5')), 1),
             ("on L1=0: a=1.1, d=2.7", (a, b, d), 2)]
    # point with L1 = L2 = 0 (solve V2 = 0 in a along b(a,d), d = 3.1)

    def V(a, b, d, k):
        tf, tg = bh_system(a, b, d)
        N = 2 * k + 2
        P, Q, om, tr, A = normalise(taylor_monomials(tf, N), taylor_monomials(tg, N), N)
        return [v.real for v in lyapunov(to_complex(P, Q, N), k)]
    d2 = mp.mpf('3.1')
    bb = lambda a: (-2 + a * (1 + mp.sqrt(1 + 8 * d2))) / (2 * d2)
    a2 = mp.findroot(lambda a: V(a, bb(a), d2, 2)[1], mp.mpf('1.0'))
    cases.append(("on L1=L2=0: d=3.1", (a2, bb(a2), d2), 3))
    # triple zero L1=L2=L3=0
    sol = mp.findroot(lambda a, b, d: V(a, b, d, 3), (mp.mpf('1.01282'), mp.mpf('0.65463'), mp.mpf('3.28862')))
    cases.append(("triple zero L1=L2=L3=0", tuple(sol), 4))
    for name, (a, b, d), k in cases:
        Vs = V(a, b, d, k)
        P, Q = normalised_PQ(a, b, d, M)
        u = displacement_coeffs(P, Q, M)  # u[0] = u_2, ..., u[M-2] = u_M
        print("==", name, " (a,b,d) =", mp.nstr(a, 12), mp.nstr(b, 12), mp.nstr(d, 12))
        print("   u_j(2pi), j=2..%d:" % M, " ".join("%.3e" % x for x in u))
        print("   pi*V_1..V_%d:" % k, " ".join("%.10e" % float(mp.pi * v) for v in Vs))
        lead = u[2 * k + 1 - 2]
        target = float(mp.pi * Vs[k - 1])
        lower = max(abs(x) for x in u[:2 * k - 1])
        rel = abs(lead - target) / abs(target)
        good = rel < 1e-5 and lower < 1e-6 * abs(target) + 1e-11
        print("   u_%d(2pi) = %.10e  vs  pi V_%d = %.10e  rel.err %.2e ; max|u_j|, j<%d: %.2e  -> %s"
              % (2 * k + 1, lead, k, target, rel, 2 * k + 1, lower, "PASS" if good else "FAIL"))
        ok = ok and good
    return ok


if __name__ == "__main__":
    sys.exit(0 if run() else 1)
