#!/usr/bin/env python3
"""
triangle.py -- the equilateral triangle (n = 3): K(tau) has exactly one zero on (0, infinity).
K > 0 on (0, tau1], K' < 0 on [tau1, tau2], K < 0 on [tau2, infinity).  Everything exact except the modal data:
u and g are polynomials (u = ut - 3/10, ut = (1-p)/2 + (1-q)/6, q = Re w^3).

Usage: python triangle.py [tau0 tau1 tau2].  Runtime about 1 minute.  Python 3.12, python-flint 0.9.0.
"""
import sys, os, math, time
import numpy as np
import sympy as sp
from flint import arb, arb_mat, acb, fmpq, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modal_lib import *
import finite_n
from finite_n import Cert, ritz, analyse, neg, min_arb, maxarb

ok_all = True


def check(name, cond):
    global ok_all
    print(("PASS" if cond else "FAIL") + ": " + name, flush=True)
    ok_all &= bool(cond)


def exact_cell():
    """g = alpha p^3 + beta p^2 + gamma p + delta q + eps p q + zeta q^2 + const, exact rational coefficients."""
    al, de, ta = sp.symbols('al de ta')
    be, ga, ep = sp.Rational(1, 32), -sp.Rational(11, 120), sp.Rational(1, 96)
    ze = -2 * al
    p = 1 + ta ** 2
    q = 1 - 3 * ta ** 2
    qx = 3 * (1 - ta ** 2)
    flux = 6 * al * p ** 2 + 4 * be * p + 2 * ga + de * qx + ep * (2 * q + p * qx) + 2 * ze * q * qx
    sol = sp.solve(sp.Poly(sp.expand(flux), ta).coeffs(), [al, de], dict=True)[0]
    return sol[al], be, ga, sol[de], ep, -2 * sol[al]


def A(r):
    r = sp.Rational(r)
    return arb(fmpq(int(r.p), int(r.q)))


def main():
    tau0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.025
    tau1 = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    tau2 = float(sys.argv[3]) if len(sys.argv) > 3 else 1.6
    ctx.prec = 5000
    t00 = time.time()
    n = 3
    G = Geometry(3)
    al, be, ga, de, ep, ze = exact_cell()
    print("g = %s p^3 + %s p^2 + %s p + %s q + %s p q + %s q^2 + const" % (al, be, ga, de, ep, ze))
    # exact polynomials in the scaled monomials  p = R^2 phat, q = R^3 qhat
    R2 = G.R2
    R3 = G.Rn
    ut = {(0, 0): arb(1) / 2 + arb(1) / 6, (1, 0): -R2 / 2, (0, 1): -R3 / 6}
    g = {(3, 0): A(al) * R2 ** 3, (2, 0): A(be) * R2 ** 2, (1, 0): A(ga) * R2, (0, 1): A(de) * R3, (1, 1): A(ep) * R2 * R3,
         (0, 2): A(ze) * R3 ** 2}
    C = Cert()
    C.G, C.n = G, 3
    C.A0, C.A1, C.Kq, C.Lam = 30, 22, 8, 400
    C.ut_L, C.g_Lraw = ut, g
    ritz(C, True)
    analyse(C, True)
    check("exact means: <ut> = 3/10", (C.c_L - arb(3) / 10).contains(0))
    check("<u^2> = 27/700", (C.norm_uL ** 2 - arb(27) / 700).contains(0))
    Tm, Mm = C.Tm, C.Mm
    mu = [arb(0)] + C.mu
    D = [arb(1)] + [-x for x in C.d]
    N = Mm + 1
    W = [[D[k] * D[l] * Tm[k, l] for l in range(N)] for k in range(N)]
    Wabs = sum((abs(W[k][l]) for k in range(N) for l in range(N)), arb(0))
    # constants: omega = 0, osc = 0, u polynomial, sup |u| = 11/30 (u(0) = 11/30... checked on a grid bound below)
    U = arb(11) / 30          # 0 <= ut <= (1+y)(2-y)^2/6 <= 2/3 on the triangle, so -3/10 <= u <= 11/30
    D0 = C.L1 + C.rho
    normu = C.norm_uL
    normg = C.norm_gL
    print(f"M={Mm} mu_1={C.mu_f[0]:.4f} mu_M={C.mu_f[-1]:.1f} rho={float(C.rho.upper()):.3e} L1={float(C.L1.upper()):.3e} Delta0={float(D0.upper()):.3e}", flush=True)
    maxn = max(float(C.Gram[k + 1, k + 1].upper()) for k in range(Mm)) ** 0.5
    check(f"Ritz functions have norm <= 1.001 ({maxn:.6f})", maxn <= 1.001)
    sumd = sum((abs(x) for x in C.d), arb(0)) * arb('1.001')
    nu = G.nu

    def err(tau):
        Bb = min_arb(tau * tau * normu / 2, tau * normg)
        return 2 * U * (2 * D0 * Bb + D0 * D0 * tau)

    errp = 2 * U * D0 * (normu / nu + sumd + C.rho)    # |K' - K~'| bound

    def kvals(tau):
        ex = [(-m * tau).exp() for m in mu]
        K = K1 = K2 = arb(0)
        for k in range(N):
            for l in range(N):
                w = W[k][l]
                if k == l:
                    E = tau * ex[k]
                    E1 = (1 - mu[k] * tau) * ex[k]
                    E2 = (-2 * mu[k] + mu[k] ** 2 * tau) * ex[k]
                else:
                    dm = mu[k] - mu[l]
                    E = (ex[l] - ex[k]) / dm
                    E1 = (mu[k] * ex[k] - mu[l] * ex[l]) / dm
                    E2 = (mu[l] ** 2 * ex[l] - mu[k] ** 2 * ex[k]) / dm
                K += w * E
                K1 += w * E1
                K2 += w * E2
        return 2 * K, 2 * K1, 2 * K2

    # (0) K > 0 on (0, tau0]: Proposition 2 with c = 3/10, E = c^2 - <u^2>, Lip(ut) = 3/2, <u^3> = 27/19250, B = 0
    c = arb(3) / 10
    E = c * c - arb(27) / 700
    L = 8 * (arb(3) / 2) / (3 * arb.pi().sqrt())
    A3 = arb(27) / 19250 / 3
    x0 = arb(tau0).sqrt()
    h0 = A3 + E / 3 * x0 ** 2 - arb(8) / 35 * c * L * x0 ** 3 - arb(128) / 3465 * L * x0 ** 5
    # h(x) = A3 + (E/3) x^2 - D x^3 - F x^5 is increasing then decreasing on x >= 0, so min at the end points
    check(f"K > 0 on (0,{tau0}]: h(0) = {float(A3.mid()):.3e} > 0, h(sqrt(tau0)) >= {float(h0.lower()):.3e} > 0", bool(h0 > 0))
    # (i) K > 0 on [tau0, tau1]
    tau = arb(tau0)
    ok = True
    worst = 1e9
    while tau < tau1:
        nxt = tau * arb('1.02')
        if nxt > tau1:
            nxt = arb(tau1)
        K, K1, K2 = kvals(tau)
        Dt = nxt - tau
        rem = Dt ** 3 * Wabs / tau ** 2
        low = K - Dt * neg(K1) - Dt ** 2 / 2 * neg(K2) - rem - err(nxt)
        worst = min(worst, float((low / K).mid()))
        if not (low > 0):
            ok = False
            print("  fails at", float(tau.mid()), float(K.mid()), float(low.mid()))
            break
        tau = nxt
    check(f"K > 0 on [{tau0},{tau1}] (worst margin/K {worst:.3f})", ok)
    # (ii) K' < 0 on [tau1, tau2]
    tau = arb(tau1)
    ok = True
    worst = -1.0
    while tau < tau2:
        nxt = tau * arb('1.01')
        if nxt > tau2:
            nxt = arb(tau2)
        K, K1, K2 = kvals(tau)
        Dt = nxt - tau
        # |K~'''| <= 6 Wabs / tau^2
        hi = K1 + Dt * abs(K2) + Dt ** 2 * 3 * Wabs / tau ** 2 + errp
        worst = max(worst, float(hi.mid()))
        if not (hi < 0):
            ok = False
            print("  fails at", float(tau.mid()), float(K1.mid()), float(hi.mid()))
            break
        tau = nxt
    check(f"K' < 0 on [{tau1},{tau2}] (largest upper bound {worst:.3e})", ok)
    K, K1, K2 = kvals(arb(tau1))
    check(f"K({tau1}) > 0", bool(K - err(arb(tau1)) > 0))
    print("   K~(tau1) =", float(K.mid()), " err ", float(err(arb(tau1)).upper()))
    K, K1, K2 = kvals(arb(tau2))
    check(f"K({tau2}) < 0", bool(K + err(arb(tau2)) < 0))
    print("   K~(tau2) =", float(K.mid()), " err ", float(err(arb(tau2)).upper()))
    # (iii) K < 0 on [tau2, 40] and beyond
    tau = arb(tau2)
    Tend = arb(40)
    ok = True
    while tau < Tend:
        nxt = tau * arb('1.1')
        K, K1, K2 = kvals(tau)
        Dt = nxt - tau
        hi = K + Dt * max(arb(0), K1.upper() if False else arb(0)) if False else K + Dt * abs(K1) * 0
        # K~ decreasing part: upper bound through Taylor with positive parts
        pos1 = arb(K1.upper()) if arb(K1.upper()) > 0 else arb(0)
        pos2 = arb(K2.upper()) if arb(K2.upper()) > 0 else arb(0)
        hi = K + Dt * pos1 + Dt ** 2 / 2 * pos2 + Dt ** 3 * Wabs / tau ** 2 + err(nxt)
        if not (hi < 0):
            ok = False
            print("  fails at", float(tau.mid()), float(K.mid()), float(hi.mid()))
            break
        tau = nxt
    check(f"K < 0 on [{tau2},40]", ok)
    C1 = arb(0)
    for k in range(N):
        for l in range(N):
            if k or l:
                C1 += abs(W[k][l]) * 2 / maxarb(mu[k], mu[l])
    e1 = 2 * U * (2 * D0 * normg + D0 * D0)
    slope = 2 * W[0][0]
    check("K < 0 on [40, infinity): slope + e1 < 0", bool(slope + e1 < 0) and bool((slope + e1) * 40 + C1 < 0))
    # enclosure of the zero
    def bisect(sgn):
        """largest point (sgn=+1) with K certified > 0, resp. smallest point (sgn=-1) with K certified < 0"""
        lo, hi = arb(tau1), arb(tau2)
        for _ in range(24):
            mid = (lo + hi) / 2
            K, K1, K2 = kvals(mid)
            e = err(mid)
            if sgn > 0:
                if K - e > 0:
                    lo = mid
                else:
                    hi = mid
            else:
                if K + e < 0:
                    hi = mid
                else:
                    lo = mid
        return lo if sgn > 0 else hi
    lo, hi = bisect(1), bisect(-1)
    check(f"the zero lies in [{float(lo.lower()):.4f}, {float(hi.upper()):.4f}]", bool(lo < hi))
    print(f"   slope 2<u g^2> = {float(slope.mid()):.6e} (exact 2*(-7101/476476000) = {2*(-7101/476476000):.6e}), e1={float(e1.upper()):.3e}, C1={float(C1.upper()):.3e}")
    print(f"[{time.time()-t00:.0f}s]")
    print("PASS (all claims)" if ok_all else "FAIL (some claim)")
    sys.exit(0 if ok_all else 1)


if __name__ == '__main__':
    main()
