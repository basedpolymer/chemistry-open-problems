#!/usr/bin/env python3
"""
finite_n.py -- certificate of K_n(tau) > 0 for tau >= tau0, for one regular n-gon (n >= 4).
See paper.tex, Section "Certified computation".

Usage: python finite_n.py n [tau0]
Prints PASS/FAIL, exit code 1 on FAIL.  Writes ../certificates/finite_n_<n>.json.
Runtime: 1-5 minutes per n on one core.  Environment: Python 3.12, python-flint 0.9.0, numpy.
Floating point is used only to choose coefficients; all inequalities are checked in arb ball arithmetic.
"""
import sys, os, math, time, json
import numpy as np
from flint import arb, arb_mat, acb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modal_lib import *


def params(n):
    """(A0, A1, Kq, Lam, KH, KL): basis ranges, Ritz cut-off, harmonic terms; any choice is rigorous."""
    table = [(4, (38, 26, 8, 800, 40, 8)), (5, (38, 24, 8, 800, 40, 8)), (6, (38, 22, 8, 800, 40, 8)),
             (8, (40, 18, 8, 800, 40, 8)), (12, (38, 12, 10, 800, 40, 8)), (20, (38, 10, 14, 800, 36, 6)),
             (32, (38, 8, 12, 800, 32, 4)), (64, (38, 6, 8, 800, 24, 3))]
    for nmax, p in table:
        if n <= nmax:
            return p
    return (38, 4, 6, 800, 24, 3)


def float_fits(n, KH):
    """floating-point least squares for a_k (torsion) and b_j (cell problem) on the half side."""
    t = math.tan(math.pi / n)
    R = math.sqrt(1 + t * t)
    M = max(6000, 30 * n * KH)
    th = np.linspace(0, np.pi, M)
    tau = t * .5 * (1 - np.cos(th))
    wt = np.gradient(tau)
    W = np.sqrt(wt)
    lw = np.log((1 + 1j * tau) / R)
    rew = lambda m: np.real(np.exp(m * lw))
    dxw = lambda m: np.real(m * np.exp((m - 1) * lw)) / R
    A = np.array([rew(n * k) for k in range(KH + 1)]).T
    rhs = -(1 - (1 + tau ** 2)) / 2
    a, *_ = np.linalg.lstsq(A * W[:, None], rhs * W, rcond=1e-14)

    def avg_h(m):
        return (2 / t) * (1 / (m + 2)) * (rew(m) * wt).sum()
    mu = .5 - (1 + t * t / 3) / 4 + sum(a[k] * avg_h(n * k) for k in range(KH + 1))
    C = .5 - mu + a[0]
    r2 = 1 + tau ** 2
    dxg = (-C / 4) * 2 + (1 / 32) * 4 * r2
    for k in range(1, KH + 1):
        m = n * k
        c = -a[k] / (4 * (1 + m))
        dxg = dxg + c * (2 * rew(m) + r2 * dxw(m))
    B = np.array([dxw(n * j) for j in range(1, KH + 1)]).T
    b, *_ = np.linalg.lstsq(B * W[:, None], -dxg * W, rcond=1e-14)
    return a, b


class Cert:
    pass


def setup_functions(n, verbose=True):
    """Torsion and cell approximants.  High accuracy pair (ut_H, g_H) used only through 1D integrals;
    low-degree pair (ut_L, g_L) used in the two-dimensional integrals."""
    A0, A1, Kq, Lam, KH, KL = params(n)
    G = Geometry(n)
    aH_f, bH_f = float_fits(n, KH)
    aL_f, bL_f = float_fits(n, KL)
    aH = [ball(x) for x in aH_f]
    bH = [arb(0)] + [ball(x) for x in bH_f]
    aL = [ball(x) for x in aL_f] + [arb(0)] * (KH - KL)
    bL = [arb(0)] + [ball(x) for x in bL_f] + [arb(0)] * (KH - KL)
    t, R = G.t, G.R
    m1 = n * KH + 6
    taus, wts = side_nodes(G, m1)
    y = [acb(1, tau) / R for tau in taus]
    yinv = [1 / v for v in y]
    z = [v ** n for v in y]
    hval = [[arb(1)] * m1]            # h_k(1,tau_j)
    dxh = [[arb(0)] * m1]             # d/dx h_k(1,tau_j)
    zk = [acb(1)] * m1
    for k in range(1, KH + 1):
        zk = [u * v for u, v in zip(zk, z)]
        hval.append([u.real for u in zk])
        c = arb(n * k) / R
        dxh.append([(u * v).real * c for u, v in zip(zk, yinv)])
    mean_h = [sum((w * hval[k][j] for j, w in enumerate(wts)), arb(0)) * 2 / (n * k + 2) for k in range(KH + 1)]
    r2vals = [1 + tau * tau for tau in taus]
    mean_p = (1 + t * t / 3) / 2
    alpha = arb(1) / 32

    def pair(a, b):
        c_ = arb(1) / 2 - mean_p / 2 + sum((a[k] * mean_h[k] for k in range(KH + 1)), arb(0))
        beta_ = -(arb(1) / 2 - c_ + a[0]) / 4
        gamma_ = [arb(0)] + [-a[k] / (4 * (1 + n * k)) for k in range(1, KH + 1)]
        return c_, beta_, gamma_
    c_H, beta_H, gamma_H = pair(aH, bH)
    c_Lf, beta_L, gamma_L = pair(aL, bL)
    # boundary values of ut_H and flux of g_H
    ut_side = []
    flux = []
    for j in range(m1):
        p = r2vals[j]
        s = (1 - p) / 2
        f = 4 * alpha * p + 2 * beta_H
        for k in range(0, KH + 1):
            s += aH[k] * hval[k][j]
            if k >= 1:
                f += gamma_H[k] * (2 * hval[k][j] + p * dxh[k][j]) + bH[k] * dxh[k][j]
        ut_side.append(s)
        flux.append(f)
    bd_ut = (2 * sum((w * u * u for w, u in zip(wts, ut_side)), arb(0))).sqrt()     # ||ut_H||_dP
    bd_e = (2 * sum((w * u * u for w, u in zip(wts, flux)), arb(0))).sqrt()         # ||e||_dP
    normE = ((1 + R) / 2).sqrt() * bd_ut
    eta = G.CP * normE + G.tr * bd_e
    omega_H = G.CP * eta
    # || g_H - g_L ||: Gram matrix of p, p h_k, h_k through the cone formula
    comps = [(2, [r2vals[j] for j in range(m1)], beta_H - beta_L)]
    for k in range(1, KH + 1):
        comps.append((2 + n * k, [r2vals[j] * hval[k][j] for j in range(m1)], gamma_H[k] - gamma_L[k]))
        comps.append((n * k, hval[k], bH[k] - bL[k]))
    Nc = len(comps)
    V = arb_mat(Nc, m1, [e for c in comps for e in c[1]])
    Vw = arb_mat(Nc, m1, [c[1][j] * wts[j] for c in comps for j in range(m1)])
    Gc = V * Vw.transpose()
    s = arb(0)
    for i in range(Nc):
        for j in range(Nc):
            s += comps[i][2] * comps[j][2] * Gc[i, j] * 2 / (comps[i][0] + comps[j][0] + 2)
    dist_HL = ub(abs(s)).sqrt()
    omega_L = omega_H + dist_HL
    C = Cert()
    for key, val in locals().items():
        if key not in ('C', 'hval', 'dxh', 'zk', 'y', 'yinv', 'z', 'taus', 'wts', 'ut_side', 'flux', 'r2vals',
                       'comps', 'V', 'Vw', 'Gc', 'pair', 's', 'f', 'p'):
            setattr(C, key, val)
    C.a, C.b, C.gamma, C.beta = aL, bL, gamma_L, beta_L
    if verbose:
        print(f"n={n}: KH={KH} KL={KL} ||ut_H||_bd={float(bd_ut.upper()):.3e} ||e||_bd={float(bd_e.upper()):.3e} "
              f"||E||<={float(normE.upper()):.3e} eta={float(eta.upper()):.3e} omega_H={float(omega_H.upper()):.3e} "
              f"dist_HL={float(dist_HL.upper()):.3e}", flush=True)
    return C


def side_range(G, fval, dfball, pieces):
    """rigorous enclosure of the range of a function on [0,t], from values at piece midpoints and a derivative"""
    lo, hi = None, None
    for i in range(pieces):
        c = G.t * (2 * i + 1) / (2 * pieces)
        rho = G.t / (2 * pieces)
        v = fval(c)
        I = c.union(c + rho).union(c - rho)
        d = abs(dfball(I)) * rho
        l, u = (v - d).lower(), (v + d).upper()
        lo = l if lo is None or l < lo else lo
        hi = u if hi is None or u > hi else hi
    return arb(lo), arb(hi)


def wall_quantities(C, pieces=None):
    """range of ut_L on the wall, Lipschitz bound of ut_L, oscillation."""
    G = C.G
    n, KL, a = G.n, C.KL, C.a
    R = G.R
    pieces = pieces or max(300, 12 * n * max(KL, 1))

    def utL(tau):
        y = acb(1, tau) / R
        z = y ** n
        s = (-tau * tau) / 2
        zk = acb(1)
        for k in range(KL + 1):
            s += a[k] * zk.real
            zk *= z
        return s

    def dutL(tau):                  # d/dtau ut_L(1,tau) = d/dy on the side
        y = acb(1, tau) / R
        z = y ** n
        s = -tau
        zk = z
        for k in range(1, KL + 1):
            s += a[k] * ((zk / y) * acb(0, arb(n * k) / R)).real
            zk *= z
        return s

    def d2utL(tau):
        y = acb(1, tau) / R
        z = y ** n
        s = arb(-1)
        zk = z
        for k in range(1, KL + 1):
            m = n * k
            s += a[k] * ((zk / (y * y)) * (-arb(m * (m - 1)) / G.R2)).real
            zk *= z
        return s

    def dxutL(tau):                 # d/dx ut_L(1,tau)
        y = acb(1, tau) / R
        z = y ** n
        s = arb(-1)
        zk = z
        for k in range(1, KL + 1):
            s += a[k] * (zk / y).real * (arb(n * k) / R)
            zk *= z
        return s

    def ddxutL(tau):                # d/dtau of the previous
        y = acb(1, tau) / R
        z = y ** n
        s = arb(0)
        zk = z
        for k in range(1, KL + 1):
            m = n * k
            s += a[k] * ((zk / (y * y)) * acb(0, arb(m * (m - 1)) / G.R2)).real
            zk *= z
        return s

    lo, hi = side_range(G, utL, dutL, pieces)
    C.ut_lo, C.ut_hi = lo, hi
    C.delta_L = arb(max(abs(lo).upper(), abs(hi).upper()))      # sup_wall |ut_L|
    C.osc_L = hi - lo
    gx_lo, gx_hi = side_range(G, dxutL, ddxutL, pieces)
    gy_lo, gy_hi = side_range(G, dutL, d2utL, pieces)
    C.Lip_L = (arb(max(abs(gx_lo).upper(), abs(gx_hi).upper())) ** 2 +
               arb(max(abs(gy_lo).upper(), abs(gy_hi).upper())) ** 2).sqrt()
    return C


def build_polys(C):
    """ut_L, u_L, g_L as polynomials in the scaled monomials."""
    G, n, KL = C.G, C.n, C.KL
    hs = harmonic_polys(G, KL)
    a, gamma, b = C.a, C.gamma, C.b
    phat = {(1, 0): G.R2}                        # p = R^2 phat
    ut = {(0, 0): arb(1) / 2, (1, 0): -G.R2 / 2}
    for k in range(KL + 1):
        ut = padd(ut, hs[k], a[k])
    g = {(2, 0): C.alpha * G.R2 * G.R2, (1, 0): C.beta * G.R2}
    for k in range(1, KL + 1):
        g = padd(g, pmul(phat, hs[k]), gamma[k])
        g = padd(g, hs[k], b[k])
    C.ut_L, C.g_Lraw = ut, g
    C.hs = hs
    return C


def ritz(C, verbose=True):
    G, n = C.G, C.n
    A0, A1, Kq, Lam = C.A0, C.A1, C.Kq, C.Lam
    S = [(a, 0) for a in range(A0 + 1)] + [(a, k) for k in range(1, Kq + 1) for a in range(A1 + 1)]
    C.S = S
    idx = {m: i for i, m in enumerate(S)}
    aS = max(m[0] for m in S)
    kS = max(m[1] for m in S)
    au = max(m[0] for m in C.ut_L)
    ku = max(m[1] for m in C.ut_L)
    ag = max(m[0] for m in C.g_Lraw)
    kg = max(m[1] for m in C.g_Lraw)
    aU, kU = max(aS, ag), max(kS, kg)
    amax = max(2 * aU + au, 2 * (aS + n), 2 * aU) + 2
    kmax = max(2 * kU + ku, 2 * kS) + 1
    t0 = time.time()
    M = Moments(G, amax, kmax)
    C.M = M
    if verbose:
        print(f"   moments amax={amax} kmax={kmax} nodes={M.m} [{time.time()-t0:.1f}s]", flush=True)
    t0 = time.time()
    N = len(S)
    I = M.I
    iR2 = 1 / G.R2
    Gm = [[None] * N for _ in range(N)]
    St = [[None] * N for _ in range(N)]
    for i, (ai, ki) in enumerate(S):
        for j in range(i + 1):
            aj, kj = S[j]
            A, K = ai + aj, ki + kj
            g = I[A][K]
            Gm[i][j] = Gm[j][i] = g
            coef = 4 * ai * aj + 2 * n * (ai * kj + aj * ki)
            s = arb(0)
            if A >= 1 and coef != 0:
                s += I[A - 1][K] * coef
            if ki > 0 and kj > 0:
                s += I[A + n - 1][K - 2] * (n * n * ki * kj)
            St[i][j] = St[j][i] = s * iR2
    L = cholesky(Gm)
    Lm = tomat(L)
    Li = Lm.inv()
    Stm = tomat(St)
    St2 = Li * Stm * Li.transpose()
    if verbose:
        print(f"   gram+chol N={N} [{time.time()-t0:.1f}s]", flush=True)
    A = np.array([[float(St2[i, j]) for j in range(N)] for i in range(N)])
    A = (A + A.T) / 2
    mu, V = np.linalg.eigh(A)
    sel = [i for i in range(N) if 1e-3 < mu[i] <= Lam]
    C.mu_f = [float(mu[i]) for i in sel]
    Vm = arb_mat(N, len(sel), [float(V[i, j]) for i in range(N) for j in sel])
    Cm = Li.transpose() * Vm
    # exact function: midpoints of the coefficients
    Mm = len(sel)
    coeff = [[arb(Cm[i, k].mid()) for k in range(Mm)] for i in range(N)]
    # mean-zero correction
    means = [sum((coeff[i][k] * I[S[i][0]][S[i][1]] for i in range(N)), arb(0)) for k in range(Mm)]
    i0 = idx[(0, 0)]
    for k in range(Mm):
        coeff[i0][k] = coeff[i0][k] - means[k]
    C.coeff = coeff
    C.Mm = Mm
    C.mu = [ball(x) for x in C.mu_f]
    return C


def analyse(C, verbose=True):
    G, n, S, coeff, Mm = C.G, C.n, C.S, C.coeff, C.Mm
    M = C.M
    I = M.I
    mu = C.mu
    # mean-zero u_L and g_L
    C.c_L = mean_poly(M, C.ut_L)
    u_L = padd(C.ut_L, {(0, 0): -C.c_L})
    m_g = mean_poly(M, C.g_Lraw)
    g_L = padd(C.g_Lraw, {(0, 0): -m_g})
    C.u_L, C.g_L = u_L, g_L
    # unified monomial list
    U = list(S) + [m for m in g_L if m not in set(S)]
    uidx = {m: i for i, m in enumerate(U)}
    NU = len(U)
    t0 = time.time()
    # Gram on U, and M^u
    Gr = [[None] * NU for _ in range(NU)]
    Mu = [[None] * NU for _ in range(NU)]
    ulist = list(u_L.items())
    for i, (ai, ki) in enumerate(U):
        for j in range(i + 1):
            aj, kj = U[j]
            A, K = ai + aj, ki + kj
            Gr[i][j] = Gr[j][i] = I[A][K]
            s = arb(0)
            for (b, l), c in ulist:
                s += c * I[A + b][K + l]
            Mu[i][j] = Mu[j][i] = s
    Grm, Mum = tomat(Gr), tomat(Mu)
    # coefficient matrix on U: column 0 = g_L, columns 1..Mm = modes
    CU = [[arb(0)] * (Mm + 1) for _ in range(NU)]
    for m, c in g_L.items():
        CU[uidx[m]][0] = c
    for i in range(len(S)):
        for k in range(Mm):
            CU[i][k + 1] = coeff[i][k]
    CUm = tomat(CU)
    GrC = Grm * CUm
    Gram = CUm.transpose() * GrC
    Tm = CUm.transpose() * (Mum * CUm)
    if verbose:
        print(f"   T matrix NU={NU} [{time.time()-t0:.1f}s]", flush=True)
    d = [float((Gram[0, k + 1] / Gram[k + 1, k + 1]).mid()) for k in range(Mm)]
    C.d = [ball(x) for x in d]
    # rho^2 = || g_L - sum d_k psi_k ||^2
    vec = [arb(1)] + [-x for x in C.d]
    rho2 = arb(0)
    for k in range(Mm + 1):
        for l in range(Mm + 1):
            rho2 += vec[k] * vec[l] * Gram[k, l]
    C.rho = rho2.sqrt() if rho2 > 0 else arb(0, float(ub(abs(rho2)))).sqrt()
    C.Gram, C.Tm = Gram, Tm
    C.norm_uL = inner_poly(M, u_L, u_L).sqrt()
    C.norm_gL = Gram[0, 0].sqrt()
    # residuals: interior
    t0 = time.time()
    polys = []
    for k in range(Mm):
        f = {S[i]: coeff[i][k] for i in range(len(S))}
        Lf = lap_poly(G, f)
        R = padd(pscale(Lf, -1), f, -mu[k])
        polys.append(R)
    keys = sorted(set(m for R in polys for m in R))
    kid = {m: i for i, m in enumerate(keys)}
    NR = len(keys)
    GR = [[None] * NR for _ in range(NR)]
    for i, (ai, ki) in enumerate(keys):
        for j in range(i + 1):
            aj, kj = keys[j]
            GR[i][j] = GR[j][i] = I[ai + aj][ki + kj]
    means = [I[a][k] for a, k in keys]
    CR = [[arb(0)] * Mm for _ in range(NR)]
    for k, R in enumerate(polys):
        for m, c in R.items():
            CR[kid[m]][k] = c
    CRm = tomat(CR)
    GR2 = CRm.transpose() * (tomat(GR) * CRm)
    nR = []
    for k in range(Mm):
        mean = sum((CR[i][k] * means[i] for i in range(NR)), arb(0))
        v = GR2[k, k] - mean * mean
        nR.append(ub(abs(v)).sqrt())
    if verbose:
        print(f"   interior residuals NR={NR} [{time.time()-t0:.1f}s]", flush=True)
    # residuals: flux on the side
    t0 = time.time()
    deg = max(2 * a + n * k for a, k in S)
    m2 = deg + n // 2 + 4
    taus, wts = side_nodes(G, m2)
    aS = max(m[0] for m in S)
    kS = max(m[1] for m in S)
    rows = []
    for tau in taus:
        p = (1 + tau * tau) / G.R2
        w = acb(1, tau)
        wn1 = w ** (n - 1)
        q = (wn1 * w).real / G.Rn
        Q1 = wn1.real * n / G.Rn
        pw = [arb(1)]
        for a in range(aS):
            pw.append(pw[-1] * p)
        qw = [arb(1)]
        for k in range(kS):
            qw.append(qw[-1] * q)
        row = []
        for (a, k) in S:
            v = arb(0)
            if a > 0:
                v += a * pw[a - 1] * qw[k] * (2 / G.R2)
            if k > 0:
                v += k * pw[a] * qw[k - 1] * Q1
            row.append(v)
        rows.append(row)
    Fm = tomat(rows)                         # (m2, N)
    Cc = tomat(coeff)                        # (N, Mm)
    Bm = Fm * Cc                             # (m2, Mm)
    nb = []
    for k in range(Mm):
        s = arb(0)
        for j in range(m2):
            s += wts[j] * Bm[j, k] * Bm[j, k]
        nb.append(ub(2 * s).sqrt())
    if verbose:
        print(f"   flux residuals m2={m2} [{time.time()-t0:.1f}s]", flush=True)
    ell = [G.CP * nR[k] + G.tr * nb[k] for k in range(Mm)]
    C.ell, C.nR, C.nb = ell, nR, nb
    L1 = arb(0)
    for k in range(Mm):
        L1 += abs(C.d[k]) * ell[k] / (2 * mu[k]).sqrt()
    C.L1 = L1
    return C


def grid_positivity(C, tau0, verbose=True):
    """K~(tau) - err(tau) > 0 on [tau0, infinity)."""
    G, Mm = C.G, C.Mm
    T = C.Tm
    mu = [arb(0)] + C.mu
    D = [arb(1)] + [-x for x in C.d]
    N = Mm + 1
    # coefficient matrix W_kl = D_k D_l T_kl (symmetrised is automatic)
    Wm = [[D[k] * D[l] * T[k, l] for l in range(N)] for k in range(N)]
    Wabs = sum((abs(Wm[k][l]) for k in range(N) for l in range(N)), arb(0))
    # constants
    omega = C.omega_L
    Dinf = omega + C.L1
    D1 = omega + C.rho
    D0 = Dinf + D1
    normu = C.norm_uL + C.osc_L
    normg = C.norm_gL + omega
    UL = C.U_L
    osc = C.osc_L
    C.Delta0 = D0

    def err(tau):
        B = min_arb(tau * tau * normu / 2, tau * normg)
        Cc = min_arb(normu ** 2 * tau ** 3 / 6, normg ** 2 * tau)
        return 2 * UL * (2 * D0 * B + D0 * D0 * tau) + 2 * osc * Cc

    def kvals(tau):
        ex = [(-mk * tau).exp() for mk in mu]
        K = K1 = K2 = arb(0)
        for k in range(N):
            for l in range(N):
                w = Wm[k][l]
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

    # geometric grid
    tau = arb(tau0)
    worst = None
    ok = True
    table = []
    ratio = arb('1.02')
    Tend = arb(40)
    npts = 0
    while tau < Tend:
        nxt = tau * ratio
        if tau > 0.2:
            nxt = tau * arb('1.1')
        K, K1, K2 = kvals(tau)
        Dt = nxt - tau
        rem = Dt ** 3 / 6 * (6 * Wabs / tau ** 2)    # |K~'''| <= 2 Wabs * 3/tau^2
        lowK = K - Dt * neg(K1) - Dt ** 2 / 2 * neg(K2) - rem
        e = err(nxt)
        marg = lowK - e
        npts += 1
        if not (marg > 0):
            ok = False
            print(f"   FAIL at tau={float(tau.mid()):.5g}: K~={float(K.mid()):.4e} lowK={float(lowK.mid()):.4e} err={float(e.mid()):.4e}", flush=True)
            break
        ratio_m = float((marg / K).mid())
        if worst is None or ratio_m < worst[0]:
            worst = (ratio_m, float(tau.mid()), float(K.mid()), float(e.mid()))
        if npts % 40 == 1 and verbose:
            print(f"   tau={float(tau.mid()):.4g} K~={float(K.mid()):.4e} K'={float(K1.mid()):.3e} err={float(e.mid()):.3e} margin/K={ratio_m:.3f}", flush=True)
        tau = nxt
    # tau >= Tend: K~ >= 2 T_00 tau - C1, err <= e1 tau
    C1 = arb(0)
    for k in range(N):
        for l in range(N):
            if k == 0 and l == 0:
                continue
            C1 += abs(Wm[k][l]) * 2 / maxarb(mu[k], mu[l])
    e1 = 2 * UL * (2 * D0 * normg + D0 * D0) + 2 * osc * normg ** 2
    slope = 2 * Wm[0][0]
    ok_end = bool(slope - e1 > 0) and bool((slope - e1) * Tend - C1 > 0)
    C.grid = dict(ok=ok, ok_end=ok_end, worst=worst, npts=npts, slope=float(slope.mid()), e1=float(e1.upper()), C1=float(C1.upper()))
    return ok and ok_end


def neg(x):
    """upper bound of the negative part of a ball"""
    u = -x.lower()
    return arb(u) if arb(u) > 0 else arb(0)


def min_arb(a, b):
    return a if a.upper() < b.upper() else b


def maxarb(a, b):
    return a if a.lower() > b.lower() else b


def sup_u(C):
    """Upper bound of sup_P |u_L| using a grid of the half-sector and the Lipschitz bound."""
    G, n, KL = C.G, C.n, C.KL
    Ns, Nv = 120, 24
    best = None
    tf = float(G.t.mid())
    for i in range(Ns + 1):
        s = (i + 0.5) / (Ns + 1)
        for j in range(Nv + 1):
            v = (j + 0.5) / (Nv + 1)
            x = arb(s)
            y = arb(s * v) * G.t
            w = acb(x, y) / G.R
            z = w ** n
            val = (1 - x * x - y * y) / 2
            zk = acb(1)
            for k in range(KL + 1):
                val += C.a[k] * zk.real
                zk *= z
            u = val.upper()
            best = u if best is None or u > best else best
    # every point of T is within h of a grid point
    h = arb(1) / (2 * (Ns + 1)) * G.R + G.t / (2 * (Nv + 1))
    supU = arb(best) + C.Lip_L * h
    C.supU = supU
    lower = -(C.delta_L)                                  # ut_L superharmonic: min on the wall
    C.U_L = arb(max((C.c_L - lower).upper(), (supU - C.c_L).upper()))
    return C


def short_time(C, tau0):
    """Sufficient condition (paper, Prop. short) for K > 0 on (0, tau0], given <u^3> >= 0."""
    c_lo = C.c_L - C.ut_hi
    c_hi = C.c_L - C.ut_lo
    normu = C.norm_uL + C.osc_L
    E_lo = lb(c_lo * c_lo - normu * normu)
    B = 2 * C.osc_L
    L = 8 * C.Lip_L / (3 * arb.pi().sqrt())
    tau = arb(tau0)
    S = (E_lo - c_hi * B) / 3 - arb(8) / 35 * c_hi * L * tau.sqrt() - arb(128) / 3465 * L * tau ** 1.5 - B / 15 * tau
    C.short = dict(c_lo=float(c_lo.lower()), c_hi=float(c_hi.upper()), E_lo=float(E_lo.lower()), L=float(L.upper()),
                   B=float(B.upper()), S=float(S.lower()))
    return bool(S > 0)


def certify(n, tau0=0.015, verbose=True):
    t00 = time.time()
    ctx.prec = 5000 if n <= 6 else 3000
    Gauss.cache.clear()
    C = setup_functions(n, verbose)
    wall_quantities(C)
    build_polys(C)
    if verbose:
        print(f"   delta_L={float(C.delta_L.upper()):.3e} osc_L={float(C.osc_L.upper()):.3e} Lip_L={float(C.Lip_L.upper()):.4f} "
              f"[{time.time()-t00:.1f}s]", flush=True)
    ritz(C, verbose)
    analyse(C, verbose)
    sup_u(C)
    if verbose:
        print(f"   M={C.Mm} mu_1={C.mu_f[0]:.4f} mu_M={C.mu_f[-1]:.1f} rho={float(C.rho.upper()):.3e} L1={float(C.L1.upper()):.3e} "
              f"omega_L={float(C.omega_L.upper()):.3e} U_L={float(C.U_L.upper()):.4f} ||u||={float(C.norm_uL.upper()):.4f} ||g||={float(C.norm_gL.upper()):.4e}",
              flush=True)
        print(f"   ell: {[round(float(x.upper()),5) for x in C.ell[:6]]} ... {[round(float(x.upper()),3) for x in C.ell[-3:]]}", flush=True)
    ok_grid = grid_positivity(C, tau0, verbose)
    ok_short = short_time(C, tau0)
    ok = ok_grid and ok_short
    print(("PASS" if ok_grid else "FAIL") + f": n={n}: K(tau) > 0 for tau >= {tau0} (Delta0={float(C.Delta0.upper()):.3e}, "
          f"worst margin/K={C.grid['worst']}, {time.time()-t00:.0f}s)", flush=True)
    print(("PASS" if ok_short else "FAIL") + f": n={n}: short-time condition on (0,{tau0}]: S={C.short['S']:.3e}", flush=True)
    return ok, C


if __name__ == '__main__':
    n = int(sys.argv[1])
    tau0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.015
    ok, C = certify(n, tau0)
    sys.exit(0 if ok else 1)
