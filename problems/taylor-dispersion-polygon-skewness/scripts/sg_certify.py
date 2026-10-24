#!/usr/bin/env python3
"""
sg_certify.py -- computer-assisted proof that the GEOMETRIC (short-time) skewness
    S^G_n = <u^3>/<u^2>^(3/2),   u = ut - <ut>,  -Lap ut = 2 in P_n, ut = 0 on dP_n,
of the regular n-gon is POSITIVE, for each n in a given range (result.md, Section 6, Lemma 8).

Notation (apothem-1 n-gon, t = tan(pi/n), R = sec(pi/n)):
  ut_N = (1-r^2)/2 + sum_{k=0}^K a_k Re((w/S)^(nk))      (floating-point least squares, then exact),
  u_N  = ut_N - <ut_N>,  E = ut - ut_N harmonic,  E = e_b := -ut_N on dP,  f = E - <E> = u - u_N,
  phi  = u_N^2 - <u_N^2>,  Phi = Phi_p + Phi_h  with -Lap Phi_p = phi exactly (polynomial particular
         solution) and Phi_h a harmonic polynomial fitted so that rho := Phi|_dP is small.
Exact decomposition:
  int u^3 = int u_N^3 + 3 int phi E + 3 int u_N f^2 + int f^3,
  int phi E = int_dP rho d_nE - int_dP e_b d_nPhi            (Green, E harmonic, int_dP d_nE = 0),
Bounds:
  |int_dP rho d_nE| <= ||rho||_dP ||d_nE||_dP <= ||rho||_dP (t+R) ||d_s e_b||_dP     (Rellich, Lemma 7),
  ||f||_P <= ||E||_P <= eps_E := sqrt((1+R)/2) ||e_b||_dP                             (Lemma 3),
  |3 int u_N f^2 + int f^3| <= (3 Uinf + osc) eps_E^2.
Checked: L := int u_N^3 - 3 int_dP e_b d_nPhi - 3 ||rho|| (t+R) ||d_s e_b|| - (3 Uinf + osc) eps_E^2 > 0.

Usage: python3 sg_certify.py nmin nmax [K]      (default 4 30)
Runtime: ~1-30 s per n (one core).  PASS/FAIL per n, exit 1 on any FAIL.
Writes certificates/sg_certificates_<nmin>_<nmax>.json.
"""
import sys, os, json, math, time
import numpy as np
from flint import arb, arb_poly, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polylib import Terms, PF, range_on
from lt_certify import float_eval, set_prec


def default_K(n):
    if n <= 6:
        return 16
    if n <= 12:
        return 10
    return 6


def lsq_fit(n, basis_ms, target, deriv, S, M=4000):
    """Float least squares on side 0 (tau in [0,t]) of sum c_j Re((w/S)^m_j) (or its d/dx) ~ target(tau)."""
    tf = math.tan(math.pi / n)
    sgrid = (1 - np.cos(np.linspace(0, np.pi, M))) / 2
    tau = tf * sgrid
    W = np.sqrt(np.gradient(tau))
    A = np.array([float_eval({(0, m): S ** (-m)}, tau, n_deriv=deriv) for m in basis_ms]).T
    c, *_ = np.linalg.lstsq(A * W[:, None], target(tau) * W, rcond=None)
    return c


def terms_float(T):
    return {k: float(v.mid()) for k, v in T.d.items()}


def certify_sg(n, K=None, J=None, pieces=200, verbose=True):
    K = default_K(n) if K is None else K
    J = K if J is None else J
    t0 = time.time()
    set_prec(n, 6 * n * max(K, J) + 16)
    S = float(1 / math.cos(math.pi / n))
    Sa = arb(S)
    t = (arb.pi() / n).tan()
    R = 1 / (arb.pi() / n).cos()
    area = n * t
    # --- ut_N ----------------------------------------------------------------------------------------
    a_coef = lsq_fit(n, [n * k for k in range(K + 1)], lambda tau: tau ** 2 / 2, False, S)
    ut = Terms({(0, 0): arb(0.5), (1, 0): arb(-0.5)})
    for k in range(K + 1):
        ut.add(0, n * k, arb(float(a_coef[k])) / Sa ** (n * k))
    UT = ut.to_pf()
    m_u = 2 * n * UT.int_T(t) / area
    uN = ut.copy().add(0, 0, -m_u)
    U = uN.to_pf()
    # --- phi and Phi ------------------------------------------------------------------------------
    u2 = uN.mul(uN)
    m_u2 = 2 * n * u2.to_pf().int_T(t) / area
    phi = u2.copy().add(0, 0, -m_u2)
    Phip = phi.neg_lap_inverse()
    assert Phip.laplacian().scale(-1).equals(phi)
    Phip_f = terms_float(Phip)
    c_coef = lsq_fit(n, [n * j for j in range(J + 1)], lambda tau: -float_eval(Phip_f, tau), False, S)
    harm = Terms()
    for j in range(J + 1):
        harm.add(0, n * j, arb(float(c_coef[j])) / Sa ** (n * j))
    assert harm.laplacian().d == {}
    PhiT = Phip + harm
    PHI = PhiT.to_pf()
    # --- boundary quantities ------------------------------------------------------------------------
    eb = -UT.boundary()                       # boundary values of E on side 0
    rmin, rmax = range_on(UT.boundary(), arb(0), t, pieces)
    delta = arb(max(abs(rmin).upper(), abs(rmax).upper()))
    osc = arb((rmax - rmin).upper())
    eb2 = 2 * n * (eb * eb).integral()(t)
    epsE = ((1 + R) / 2 * eb2).sqrt()
    deb = eb.derivative()
    deb2 = 2 * n * (deb * deb).integral()(t)
    rho = PHI.boundary()
    rho2 = 2 * n * (rho * rho).integral()(t)
    dnPhi = PHI.dx_boundary()
    B1 = 2 * n * (eb * dnPhi).integral()(t)  # int_dP e_b d_n Phi (exact)
    # --- interior integrals -----------------------------------------------------------------------
    U2 = U * U
    I3 = 2 * n * (U2 * U).int_T(t)
    I2 = 2 * n * U2.int_T(t)
    Uinf = arb(max((m_u + delta).upper(), (R * R / 2 + delta - m_u).upper()))
    rell = 3 * rho2.sqrt() * (t + R) * deb2.sqrt()
    second = (3 * Uinf + osc) * epsE * epsE
    L = I3 - 3 * B1 - rell - second
    ok = bool(L > 0)
    # enclosure of S^G: int u^3 in [L, I3 - 3 B1 + rell + second];
    # int u^2 = int u_N^2 + 2 int u_N E + int f^2,  |2 int u_N E| <= 2 ||u_N|| eps_E
    i3e = arb(((L + (I3 - 3 * B1 + rell + second)) / 2).mid(), (((I3 - 3 * B1 + rell + second) - L) / 2).upper())
    e2 = 2 * I2.sqrt() * epsE + epsE * epsE
    i2e = arb(I2.mid(), (e2 + I2.rad()).upper())
    SG = (i3e / area) / (i2e / area) ** 1.5
    dt = time.time() - t0
    res = dict(n=n, K=K, J=J, prec=ctx.prec, delta=float(delta.upper()), L2_E=float(epsE.upper()),
               L2_rho=float(rho2.sqrt().upper()), L2_dseb=float(deb2.sqrt().upper()),
               int_uN3=I3.str(10), boundary_term=(3 * B1).str(6), rellich=float(rell.upper()),
               second=float(second.upper()), lower=float(L.lower()), SG=SG.str(6),
               a_coef=[float(x).hex() for x in a_coef], c_coef=[float(x).hex() for x in c_coef],
               S=S.hex(), PASS=ok, seconds=round(dt, 1))
    if verbose:
        print(f"n={n:3d} K={K:2d} ||E||<={res['L2_E']:.1e} ||rho||<={res['L2_rho']:.1e} int u_N^3={float(I3.mid()):.5e} "
              f"3B1={float(B1.mid())*3:.2e} rell<={res['rellich']:.1e} 2nd<={res['second']:.1e} "
              f"L>={res['lower']:.4e} S^G in {SG.str(4)} [{dt:.1f}s] {'PASS' if ok else 'FAIL'}", flush=True)
    return res


def main():
    nmin = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    Kf = int(sys.argv[3]) if len(sys.argv) > 3 else None
    print(f"sg_certify.py: certify S^G_n > 0 for {nmin} <= n <= {nmax}")
    allok = True
    allres = []
    for n in range(nmin, nmax + 1):
        r = certify_sg(n, Kf)
        allres.append(r)
        allok &= r['PASS']
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, '..', 'certificates', f'sg_certificates_{nmin}_{nmax}.json'), 'w') as f:
        json.dump(allres, f, indent=1)
    print(("PASS" if allok else "FAIL") + f": geometric skewness positive for all {nmin} <= n <= {nmax}")
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
