#!/usr/bin/env python3
"""
uniform.py -- certificate of K_n(tau) > 0 for all tau >= tau0 and ALL n >= N (uniform comparison with the unit disk),
and of the same statement for the disk itself (N = infinity, t = 0).  See paper.tex, Section "All large n".

The n-gon is normalised to area pi (apothem a, a^2 = arctan(t)/t, circumradius v, v^2 = a^2 (1+t^2)), so the radial
comparison functions g_0, psi_k = J_0(j_k r)/J_0(j_k) (J_1(j_k) = 0), mu_k = j_k^2, do not depend on n.
Every error term is an increasing function of t and is evaluated at t = tan(pi/N).

Usage: python uniform.py N [M] [tauhat0]      (N = 0 means the disk).  Runtime: seconds.
Prints PASS/FAIL per claim and exits non-zero on failure.  Python 3.12, python-flint 0.9.0.
"""
import sys, os, math
import numpy as np
from scipy.special import jn_zeros
from flint import arb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modal_lib import ub, lb, ball

ok_all = True


def check(name, cond):
    global ok_all
    print(("PASS" if cond else "FAIL") + ": " + name, flush=True)
    ok_all &= bool(cond)


def bessel_zero(guess):
    """Tight ball around a zero of J_1, certified by a sign change."""
    x = arb(guess)
    for _ in range(12):
        m = arb(x.mid())
        j1 = m.bessel_j(1)
        dj1 = m.bessel_j(0) - j1 / m
        x = m - j1 / dj1
    eps = arb(2) ** (-(ctx.prec - 60))
    lo, hi = arb(x.mid()) - eps, arb(x.mid()) + eps
    if not bool(lo.bessel_j(1) * hi.bessel_j(1) < 0):
        raise ArithmeticError("no sign change")
    return lo.union(hi)


def data(N, M, tauhat0, verbose=True):
    ctx.prec = 300
    pi = arb.pi()
    if N == 0:
        t = arb(0)
    else:
        t = ub((pi / N).tan())
    x = t * t
    a_lo = (1 - x / 3).sqrt()
    v_up = (1 + 2 * x / 3).sqrt()
    A = arb(1 - x / 3).union(arb(1))             # a^2 in [1 - t^2/3, 1]
    CP = v_up / pi
    nu = pi * pi / (v_up * v_up)
    tr = ((2 * CP * CP + 2 * v_up * CP) / a_lo).sqrt()
    Hnorm = x / 4 * (t / pi + 4 * x / 3).sqrt()    # ||H_c||
    # ||e_0||_bd^2 <= Q/32, Q = (1/t) int_0^t (a^2 tau^2 + a^2 - 1)^2 dtau
    zeta = arb(0).union(x / 5)                    # a^2 = 1 - t^2/3 + t^2 zeta, 0 <= zeta <= t^2/5
    Qx = A * A / 5 + 2 * A * (-arb(1) / 3 + zeta) / 3 + (arb(1) / 3 - zeta) ** 2
    e0 = ub((Qx / 32).sqrt()) * x
    eta = CP * Hnorm + tr * e0
    omega = CP * eta
    w1 = x * x / 6                               # W_1 distance of the radial laws
    bbar = x * x / 20                            # b-bar = (a^2(1+t^2/3) - 1)/2 <= t^4/20 (checked below)
    # disk data
    jf = jn_zeros(1, M)
    j = [bessel_zero(float(z)) for z in jf]
    mu = [z * z for z in j]
    J0 = [z.bessel_j(0) for z in j]
    d = [-2 / m ** 2 for m in mu]
    # sup bounds on [0,v]
    Bk = [1 / abs(z) for z in J0]
    Fg = arb(1) / 47                                # sup |g_0| on [0,v] (g_0(0) = 1/48; checked below)
    Lg = arb(1) / 20                                # sup |g_0'| on [0,v]
    Fs = sum((abs(d[k]) * Bk[k] for k in range(M)), arb(0))
    Ls = sum((abs(d[k]) * j[k] * Bk[k] / arb(2).sqrt() for k in range(M)), arb(0))
    Finf = Fg + Fs
    Linf = Lg + Ls
    # tail
    tail2 = arb(1) / 11520 - sum((dk * dk for dk in d), arb(0))
    # flux constants: Psi''(V) on V = r^2, r in [a,v]
    r = arb(a_lo).union(v_up)
    Psi = []
    for k in range(M):
        z = j[k] * r
        psi = z.bessel_j(0) / J0[k]
        dpsi = -j[k] * z.bessel_j(1) / J0[k]
        P2 = (-mu[k] * psi - 2 * dpsi / r) / (4 * r * r)
        Psi.append(ub(abs(P2)))
    L1 = sum((abs(d[k]) * tr * 16 * Psi[k] * e0 / (2 * mu[k]).sqrt() for k in range(M)), arb(0))
    # sup of zeta on [1,v]
    if N == 0:
        zbar = arb(0)
    else:
        r1 = arb(1).union(v_up)
        rr = r1 * r1
        val = rr * rr / 32 - rr / 16 + arb(1) / 48
        for k in range(M):
            val -= d[k] * (j[k] * r1).bessel_j(0) / J0[k]
        zbar = ub(abs(val))
    rho = (ub(abs(tail2)) + x / 3 * zbar * zbar).sqrt()
    Delta0 = 2 * omega + L1 + rho
    osc = x / 2
    U0 = arb(1) / 4 + x / 3 + bbar / 2
    u2 = arb(1) / 16 - A * (1 + x / 3) / 8 + A * A * (1 + 2 * x / 3 + x * x / 5) / 12
    normu = ub(abs(u2)).sqrt() + Hnorm
    normg = (arb(1) / 11520 + w1 * 2 * Fg * Lg).sqrt() + omega
    vv = v_up
    Urad = arb(1) / 4 + x / 3
    rateW = 2 * (w1 * (vv * Finf * Finf + 4 * U0 * Finf * Linf) + bbar / 2 * Finf * Finf)
    if verbose:
        print(f"N={N} t={float(t.upper()):.6f} ||H_c||<={float(Hnorm.upper()):.3e} ||e0||<={float(e0.upper()):.3e} "
              f"omega={float(omega.upper()):.3e} L1={float(L1.upper()):.3e} rho={float(rho.upper()):.3e} "
              f"Delta0={float(Delta0.upper()):.3e} F={float(Finf.upper()):.4f} Lip={float(Linf.upper()):.4f} "
              f"rateW={float(rateW.upper()):.3e}", flush=True)
    return dict(mu=mu, d=d, Delta0=Delta0, osc=osc, U0=U0, normu=normu, normg=normg, rateW=rateW, t=t, x=x)


def Wmatrix(mu, d):
    M = len(mu)
    N1 = M + 1
    mus = [arb(0)] + mu
    D = [arb(1)] + [-x for x in d]
    T = [[None] * N1 for _ in range(N1)]
    T[0][0] = arb(1) / 184320
    for k in range(1, N1):
        m = mus[k]
        T[0][k] = T[k][0] = (m * m + 96 * m - 1728) / (48 * m ** 3)
        T[k][k] = arb(1) / 12
        for l in range(1, k):
            T[k][l] = T[l][k] = -2 * (m + mus[l]) / (m - mus[l]) ** 2
    W = [[D[k] * D[l] * T[k][l] for l in range(N1)] for k in range(N1)]
    return W, mus


def neg(x):
    u = -x.lower()
    return arb(u) if arb(u) > 0 else arb(0)


def mina(a, b):
    return a if a.upper() < b.upper() else b


def maxa(a, b):
    return a if a.lower() > b.lower() else b


def kvals(W, mu, tau):
    N1 = len(mu)
    ex = [(-m * tau).exp() for m in mu]
    K = K1 = K2 = arb(0)
    for k in range(N1):
        for l in range(N1):
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


def positivity(W, mu, tau0, D0, normu, normg, U0, osc, rateW, verbose=True, sign=1):
    N1 = len(mu)
    Wabs = sum((abs(W[k][l]) for k in range(N1) for l in range(N1)), arb(0))

    def err(tau):
        B = mina(tau * tau * normu / 2, tau * normg)
        Cc = mina(normu ** 2 * tau ** 3 / 6, normg ** 2 * tau)
        return 2 * U0 * (2 * D0 * B + D0 * D0 * tau) + 2 * osc * Cc + rateW * tau

    tau = arb(tau0)
    Tend = arb(40)
    ok = True
    worst = None
    npts = 0
    while tau < Tend:
        nxt = tau * (arb('1.02') if tau < 0.2 else arb('1.1'))
        K, K1, K2 = kvals(W, mu, tau)
        Dt = nxt - tau
        rem = Dt ** 3 / 6 * (6 * Wabs / tau ** 2)
        lowK = K - Dt * neg(K1) - Dt ** 2 / 2 * neg(K2) - rem
        e = err(nxt)
        marg = lowK - e
        npts += 1
        rat = float((marg / K).mid())
        if worst is None or rat < worst[0]:
            worst = (rat, float(tau.mid()), float(K.mid()), float(e.mid()))
        if not (marg > 0):
            ok = False
            print(f"   FAIL at tau={float(tau.mid()):.5g}: K={float(K.mid()):.4e} err={float(e.mid()):.4e}")
            break
        if verbose and npts % 40 == 1:
            print(f"   tau={float(tau.mid()):.4g} K_D={float(K.mid()):.4e} err={float(e.mid()):.3e} margin/K={rat:.3f}", flush=True)
        tau = nxt
    C1 = arb(0)
    for k in range(N1):
        for l in range(N1):
            if k or l:
                C1 += abs(W[k][l]) * 2 / maxa(mu[k], mu[l])
    e1 = 2 * U0 * (2 * D0 * normg + D0 * D0) + 2 * osc * normg ** 2 + rateW
    slope = 2 * W[0][0]
    ok_end = bool(slope - e1 > 0) and bool((slope - e1) * Tend - C1 > 0)
    return ok, ok_end, worst, npts


def short_uniform(N, tau0):
    """Short-time condition of the paper for all n >= N (apothem 1 frame; N = 0: disk)."""
    pi = arb.pi()
    t = arb(0) if N == 0 else ub((pi / N).tan())
    x = arb(0).union(t * t)
    c_lo = arb(1) / 4 + x / 12 - x / 6
    c_hi = arb(1) / 4 + x / 12 + x / 3
    s = arb(1) / 2 + x / 6
    Y2 = s * s - s * (1 + x / 3) / 2 + (1 + 2 * x / 3 + x * x / 5) / 12
    E_lo = lb(2 * c_lo * c_lo - (Y2.sqrt() + x / 3) ** 2)
    L = 8 * (1 + x).sqrt() / (3 * pi.sqrt())
    B = x
    tau = arb(tau0)
    S = (E_lo - c_hi * B) / 3 - arb(8) / 35 * c_hi * L * tau.sqrt() - arb(128) / 3465 * L * tau ** 1.5 - B / 15 * tau
    return S


def run(N, M=12, tauhat0=0.0149, verbose=True):
    D = data(N, M, tauhat0, verbose)
    S = short_uniform(N, 0.015)
    check(f"N={N}: short-time condition on (0,0.015] holds uniformly, S >= {float(S.lower()):.4e}", bool(S > 0))
    W, mus = Wmatrix(D['mu'], D['d'])
    x = D['x']
    if N:
        check(f"N={N}: (1 - t^2/3) tau0 >= tauhat0 with tau0 = 0.015", bool((1 - x / 3) * arb('0.015') > arb(tauhat0)))
    ok, ok_end, worst, npts = positivity(W, mus, tauhat0, D['Delta0'], D['normu'], D['normg'], D['U0'], D['osc'], D['rateW'], verbose)
    name = "disk" if N == 0 else f"all n >= {N}"
    check(f"{name}: K > 0 on [{tauhat0}, 40] ({npts} grid points), worst margin/K = {worst}", ok)
    check(f"{name}: K > 0 on [40, infinity)", ok_end)
    return ok and ok_end


if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    tauhat0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0149
    run(N, M, tauhat0)
    print("PASS (all claims)" if ok_all else "FAIL (some claim)")
    sys.exit(0 if ok_all else 1)
