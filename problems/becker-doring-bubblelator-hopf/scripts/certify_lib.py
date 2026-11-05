"""Certification routines for the Hopf point of the limit model at fixed (beta, nu).

Environment: python-flint 0.9 (Arb), Python 3.  Deterministic.

Notation (see paper.tex):  G = G_{beta,nu},  H(t) = Re G(it),  H'(t) = Im G_{beta,nu+1}(it),
theta(t) = -Im G(it)/t  (the value of vartheta at which it is a characteristic root, when H(t)=0).

certify_point(beta, nu, t0float) proves, in ball arithmetic:
  (C1) there is t_lo < t_hi (a tiny interval T around t0float) with H(t_lo) > 0 > H(t_hi);
  (C2) H' < 0 on a larger interval T' containing T  -> H has exactly one zero t0 in T',
       it lies in T, and H'(t0) < 0 (transversality, simple root);
  (C3) Im G(it) < 0 on T  -> vartheta0 = theta(t0) > 0;
  (C4) scan of (0, G(0)/vartheta0_low] minus T': on every piece J, either H has no zero on J,
       or -Im G(it) < vartheta0_low * t on J.  Hence for vartheta >= vartheta0 the only purely
       imaginary characteristic roots are +-i t0 at vartheta = vartheta0 (this contains the
       nonresonance condition at every multiple k t0, k >= 2);
  (C5) Re ghat_1(1) > 0 (transversality in the form used by the implicit function theorem),
       |L(2)| > 0, and an enclosure of eta2 (direction coefficient) that excludes 0.
"""
from flint import arb, acb, ctx
import bdhopf as B


def iv(a, b):
    """Ball containing the closed interval [a, b] (a <= b floats or arbs)."""
    return arb(a).union(arb(b))


def H(beta, nu, t):
    return B.G(beta, nu, acb(0, t)).real


def Hp(beta, nu, t):
    return B.G(beta, nu + 1, acb(0, t)).imag


def find_first_zero(beta, nu, t_start=0.02, tmax=500.0, fac=1.03):
    """Floating-point search (exploration only) for the first sign change of H from + to -."""
    t = t_start
    h0 = float(H(beta, nu, arb(t)).mid())
    while t < tmax:
        t2 = t * fac
        h2 = float(H(beta, nu, arb(t2)).mid())
        if h0 > 0 and h2 < 0:
            a, b = t, t2
            for _ in range(70):
                m = 0.5 * (a + b)
                if float(H(beta, nu, arb(m)).mid()) > 0:
                    a = m
                else:
                    b = m
            return 0.5 * (a + b)
        t, h0 = t2, h2
    return None


def refine_zero(beta, nu, t0, it=60):
    a, b = t0 * (1 - 1e-6), t0 * (1 + 1e-6)
    for _ in range(it):
        m = 0.5 * (a + b)
        if float(H(beta, nu, arb(m)).mid()) > 0:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def certify_point(beta, nu, t0f, **kw):
    """beta, nu: arb balls (exact parameters inside).  t0f: float approximation of the zero."""
    r = max(t0f, 1.0) * 2.0 ** -40
    return certify_T(beta, nu, t0f - r, t0f + r, **kw)


def certify_T(beta, nu, tlo, thi, scan_depth=40, npieces=64):
    """beta, nu: arb balls; [tlo, thi]: interval that should contain the zero t0(beta, nu) for
    every parameter value in the balls.  All checks are uniform over the balls."""
    out = {}
    ok = True
    # (C1)
    Hlo, Hhi = H(beta, nu, arb(tlo)), H(beta, nu, arb(thi))
    c1 = bool(Hlo > 0) and bool(Hhi < 0)
    out['C1'] = c1
    ok &= c1
    T = iv(tlo, thi)
    t0f = 0.5 * (tlo + thi)
    # (C2) derivative negative on T' = [t0 - d, t0 + d], pieces
    d = 0.05 * t0f
    while True:
        pieces = 16
        good = True
        A, Bq = tlo - d, thi + d
        for j in range(pieces):
            a = A + (Bq - A) * j / pieces
            b = A + (Bq - A) * (j + 1) / pieces if j < pieces - 1 else Bq
            if not bool(Hp(beta, nu, iv(a, b)) < 0):
                good = False
                break
        if good or d < 1e-7 * t0f:
            break
        d /= 2
    c2 = good
    out['C2'] = c2
    out['Tprime'] = (tlo - d, thi + d)
    ok &= c2
    # (C3)
    GT = B.G(beta, nu, acb(0, T))
    c3 = bool(GT.imag < 0)
    out['C3'] = c3
    ok &= c3
    vth0 = -GT.imag / T
    vlow = vth0.lower()
    g0 = B.G0(beta, nu)
    Tmax = float((g0 / vlow).upper()) * (1 + 1e-12)
    out['Tmax'] = Tmax
    # (C4) scan
    stack = []
    # left part (0, t0-d], right part [t0+d, Tmax]
    # [0, a00]: |G(it) - G(0)| <= t G_{beta,nu+1}(0), so H > 0 there if G(0) - a00 G_{beta,nu+1}(0) > 0
    a00 = 2.0 ** -20
    c0 = bool(g0 - arb(a00) * B.G0(beta, nu + 1) > 0)
    out['C4_origin'] = c0
    ok &= c0
    nl = npieces
    a0, b0 = a00, tlo - d
    for j in range(nl):
        stack.append((a0 + (b0 - a0) * j / nl, a0 + (b0 - a0) * (j + 1) / nl if j < nl - 1 else b0, 0))
    if Tmax > thi + d:
        nr = npieces
        a0, b0 = thi + d, Tmax
        for j in range(nr):
            stack.append((a0 + (b0 - a0) * j / nr, a0 + (b0 - a0) * (j + 1) / nr if j < nr - 1 else b0, 0))
    c4 = True
    nev = 0
    while stack:
        a, b, dep = stack.pop()
        J = iv(a, b)
        GJ = B.G(beta, nu, acb(0, J))
        nev += 1
        if bool(GJ.real > 0) or bool(GJ.real < 0):
            continue
        if bool(-GJ.imag < vlow * arb(a)):
            continue
        if dep >= scan_depth:
            c4 = False
            out['C4_fail_at'] = (a, b)
            break
        m = 0.5 * (a + b)
        stack.append((a, m, dep + 1))
        stack.append((m, b, dep + 1))
    out['C4'] = c4
    out['C4_evals'] = nev
    ok &= c4
    # (C5)
    q = B.hopf_quantities(beta, nu, T)
    c5 = bool(q['ghat1_1'].real > 0) and bool(abs(q['L2']) > 0) and (bool(q['eta2'] < 0) or bool(q['eta2'] > 0))
    out['C5'] = c5
    ok &= c5
    out['q'] = q
    out['T'] = T
    out['ok'] = ok
    return out
