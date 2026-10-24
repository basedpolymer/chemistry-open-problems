#!/usr/bin/env python3
"""
lt_certify.py -- computer-assisted proof that the long-time skewness coefficient of Taylor dispersion
in the regular n-gon is POSITIVE, for each n in a given range (default 4..N0 with N0 = 24).

Claim checked for each n (apothem-1 n-gon P, see result.md, Theorem 2 and Section 4):
    I_n := int_P u g^2 dA > 0,
where  -Lap ut = 2 in P, ut = 0 on dP,  u = ut - <ut>,  -Lap g = u in P, d_n g = 0 on dP, <g> = 0.
(By Proposition 1 of result.md, <u g_2> = 2 <u g^2>, so sign(long-time skewness) = sign(I_n).)

Method (all certified quantities in arb ball arithmetic, python-flint):
  * approximations ut_N = (1-r^2)/2 + sum_k a_k Re((w/S)^(nk)) and
    g_N = particular solution of -Lap g_N = u_N := ut_N - <ut_N> + sum_j b_j Re((w/S)^(nj)) - const,
    coefficients a_k, b_j found by floating-point least squares (they are then exact binary numbers;
    how they were found is irrelevant for rigour);
  * E = ut - ut_N is harmonic with boundary values -ut_N:  |E| <= delta := max_dP |ut_N|,
    ||E||_{L2(P)}^2 <= ((1+R)/2) ||ut_N||_{L2(dP)}^2   (Lemma 3);
  * w = g - g_N:  ||grad w|| <= eta := p ||E|| + tau_tr ||e||_{L2(dP)},  ||w|| <= omega := p eta,
    p = 2R/pi (Payne-Weinberger), tau_tr = sqrt(2 p^2 + 2 R p) (trace inequality, x.n = 1 on dP),
    e = d_n g_N on dP (Lemma 4);
  * |I_n - I_N| <= Err (Lemma 5), and we check I_N - Err > 0.
Also prints an enclosure of kappa_n = 6<u g^2>/(2<u g>)^(3/2) (apothem 1) and of the
circumradius-1 value kappa_n cos(pi/n) for comparison with Table 1 of Aminian-Camassa-McLaughlin.

Usage:  python3 lt_certify.py [nmin nmax]      (default 4 24)
Runtime: about 1-3 s per n on one core (python3.13, python-flint 0.9.0, numpy 2.x).
Output: PASS/FAIL per n; exit code 1 on any FAIL. Writes certificates/lt_certificates.json.
"""
import sys, os, json, math, time
import numpy as np
from flint import arb, arb_poly, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polylib import Terms, PF, sup_abs_on, range_on, ub, lb, sqrt_ub, lo_hi

ctx.prec = 320


def set_prec(n, maxdeg):
    """Working precision: the tau-monomial representation of Re((1+i tau)^m) has cancellation of about
    m log2(1+t) bits; we add a safety margin.  (Any precision is rigorous; too low only widens balls.)"""
    ctx.prec = 160 + int(1.25 * maxdeg * math.log2(1 + math.tan(math.pi / n)))


def choose_K(n):
    # number of harmonic correction terms (heuristic choice; any choice is rigorous)
    if n <= 5:
        return 16
    if n <= 8:
        return 10
    if n <= 14:
        return 6
    return 4


def float_eval(terms_f, tau, n_deriv=False):
    """Float evaluation on side x = 1 of sum c r^(2a) Re(w^m) (or its x-derivative)."""
    w = 1 + 1j * tau
    r2 = 1 + tau ** 2
    out = np.zeros_like(tau)
    for (a, m), c in terms_f.items():
        if not n_deriv:
            out += c * np.real(r2 ** a * w ** m)
        else:
            val = 0
            if a > 0:
                val = val + 2 * a * r2 ** (a - 1) * w ** m
            if m > 0:
                val = val + m * r2 ** a * w ** (m - 1)
            out += c * np.real(val)
    return out


def certify(n, K=None, J=None, pieces=240, verbose=True, raw=False):
    K = choose_K(n) if K is None else K
    J = K if J is None else J
    t0 = time.time()
    set_prec(n, 4 * (n * max(K, J) + 4))
    h = math.pi / n
    tf = math.tan(h)
    S = float(1 / math.cos(h))          # exact binary scale factor for the basis (any value is fine)
    Sa = arb(S)
    # ---- floating-point fit of ut_N -------------------------------------------------------------
    M = 3000
    sgrid = (1 - np.cos(np.linspace(0, np.pi, M))) / 2
    tau = tf * sgrid
    wts = np.gradient(tau)
    W = np.sqrt(wts)
    basis = [{(0, n * k): S ** (-n * k)} for k in range(0, K + 1)]
    A = np.array([float_eval(b, tau) for b in basis]).T
    rhs = -(1 - (1 + tau ** 2)) / 2   # -(1 - r^2)/2 on x=1
    a_coef, *_ = np.linalg.lstsq(A * W[:, None], rhs * W, rcond=None)
    # ---- rigorous objects --------------------------------------------------------------------
    t = arb.pi() / n
    t = t.tan()
    R = 1 / (arb.pi() / n).cos()
    area = n * t                       # |P|
    ut = Terms({(0, 0): arb(0.5), (1, 0): arb(-0.5)})
    for k in range(0, K + 1):
        ut.add(0, n * k, arb(float(a_coef[k])) / Sa ** (n * k))
    UT = ut.to_pf()
    m_u = 2 * n * UT.int_T(t) / area
    uN = ut.copy().add(0, 0, -m_u)
    gp = uN.neg_lap_inverse()
    assert gp.laplacian().scale(-1).equals(uN), "particular solution check failed"
    # ---- floating-point fit of the harmonic part of g_N ---------------------------------------
    gp_f = {k: float(v.mid()) for k, v in gp.d.items()}
    rhs2 = -float_eval(gp_f, tau, n_deriv=True)
    B = np.array([float_eval({(0, n * j): S ** (-n * j)}, tau, n_deriv=True) for j in range(1, J + 1)]).T
    b_coef, *_ = np.linalg.lstsq(B * W[:, None], rhs2 * W, rcond=None)
    harm = Terms()
    for j in range(1, J + 1):
        harm.add(0, n * j, arb(float(b_coef[j - 1])) / Sa ** (n * j))
    assert harm.laplacian().d == {}, "harmonic part not harmonic"
    gN = gp + harm
    G0 = gN.to_pf()
    m_g = 2 * n * G0.int_T(t) / area
    gN.add(0, 0, -m_g)
    U = uN.to_pf()
    G = gN.to_pf()
    # ---- boundary quantities -------------------------------------------------------------------
    bu = UT.boundary()                  # ut_N(1, tau)
    rmin, rmax = range_on(bu, arb(0), t, pieces)
    delta = arb(max(abs(rmin).upper(), abs(rmax).upper()))
    osc = arb((rmax - rmin).upper())   # oscillation of E = -ut_N on dP  (>= sup |E - <E>|)
    bu2 = 2 * n * (bu * bu).integral()(t)     # ||ut_N||^2_{L2(dP)}
    e = G.dx_boundary()                 # normal derivative of g_N on side 0
    e2 = 2 * n * (e * e).integral()(t)
    Gb = G.boundary()
    # ---- interior integrals ---------------------------------------------------------------------
    G2 = G * G
    UG = U * G
    I_N = 2 * n * (U * G2).int_T(t)
    intG = 2 * n * G.int_T(t)
    intG2 = 2 * n * G2.int_T(t)
    intG4 = 2 * n * (G2 * G2).int_T(t)
    intUG = 2 * n * UG.int_T(t)
    intUG2sq = 2 * n * (UG * UG).int_T(t)
    gradG2 = intUG + 2 * n * (Gb * e).integral()(t)    # ||grad g_N||^2 = int g_N u_N + int_dP g_N e
    # ---- constants -------------------------------------------------------------------------------
    p = 2 * R / arb.pi()
    tau_tr = (2 * p * p + 2 * R * p).sqrt()
    epsE = ((1 + R) / 2 * bu2).sqrt()
    eN = e2.sqrt()
    eta = p * epsE + tau_tr * eN
    omega = p * eta
    sd_G2 = (intG4 - intG2 * intG2 / area).sqrt()
    sd_UG = (intUG2sq - intUG * intUG / area).sqrt()
    normG = intG2.sqrt()
    Uinf = arb(max((m_u + delta).upper(), (R * R / 2 + delta - m_u).upper()))
    T1 = epsE * sd_G2
    T2 = osc * (2 * normG * omega + omega * omega)
    T3 = 2 * sd_UG * omega
    T4 = Uinf * omega * omega
    Err = T1 + T2 + T3 + T4
    lower = I_N - Err
    ok = bool(lower > 0)
    # ---- enclosure of kappa_n -------------------------------------------------------------------
    I_enc = arb(I_N.mid(), (Err + I_N.rad()).upper())
    gn = gradG2.sqrt()
    grad_enc = arb(gn.mid(), (eta + gn.rad()).upper())
    ug_enc = grad_enc * grad_enc / area
    kappa = 6 * (I_enc / area) / (2 * ug_enc) ** 1.5
    kappa_circ = kappa * (arb.pi() / n).cos()
    dt = time.time() - t0
    res = dict(n=n, K=K, J=J, delta=float(delta.upper()), L2_bdry_ut=float(bu2.sqrt().upper()),
               L2_E=float(epsE.upper()), L2_e=float(eN.upper()), eta=float(eta.upper()), omega=float(omega.upper()),
               I_N=I_N.str(12), Err=float(Err.upper()), I_lower=float(lower.lower()),
               terms=[float(T1.upper()), float(T2.upper()), float(T3.upper()), float(T4.upper())],
               mean_g_check=intG.str(5), kappa_apothem1=lo_hi(kappa), kappa_circumradius1=lo_hi(kappa_circ),
               a_coef=[float(x).hex() for x in a_coef], b_coef=[float(x).hex() for x in b_coef],
               S=S.hex(), pieces=pieces, prec=ctx.prec, seconds=round(dt, 2), PASS=ok)
    if verbose:
        print(f"n={n:3d} K={K:2d} delta={res['delta']:.2e} ||E||<={res['L2_E']:.2e} ||e||<={res['L2_e']:.2e} "
              f"I_N={float(I_N.mid()):.6e} Err<={res['Err']:.2e} I>={res['I_lower']:.4e} "
              f"kappa(R=1) in [{lo_hi(kappa_circ)[0]:.4f}, {lo_hi(kappa_circ)[1]:.4f}]  [{dt:.1f}s]  {'PASS' if ok else 'FAIL'}")
    if raw:
        # raw ball quantities for alltime_windows.py (all integrals over P, not averages)
        intU2 = 2 * n * (U * U).int_T(t)
        intU4 = 2 * n * ((U * U) * (U * U)).int_T(t)
        intU3 = 2 * n * ((U * U) * U).int_T(t)
        res['raw'] = dict(area=area, R=R, t=t, m_u=m_u, delta=delta, osc=osc, epsE=epsE, eta=eta, omega=omega,
                          Uinf=Uinf, I_N=I_N, Err=Err, intG2=intG2, intUG=intUG, sd_UG=sd_UG, gradG2=gradG2,
                          intU2=intU2, intU3=intU3, intU4=intU4)
    return res


def main():
    nmin, nmax = 4, 24
    if len(sys.argv) >= 3:
        nmin, nmax = int(sys.argv[1]), int(sys.argv[2])
    print(f"lt_certify.py: certify I_n = int_P u g^2 > 0 for {nmin} <= n <= {nmax} (prec {ctx.prec} bits)")
    allres = []
    allok = True
    for n in range(nmin, nmax + 1):
        r = certify(n)
        allres.append(r)
        allok &= r['PASS']
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, '..', 'certificates', f'lt_certificates_{nmin}_{nmax}.json')
    with open(out, 'w') as f:
        json.dump(allres, f, indent=1)
    print(("PASS" if allok else "FAIL") + f": long-time coefficient positive for all {nmin} <= n <= {nmax}")
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
