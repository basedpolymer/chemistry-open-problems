#!/usr/bin/env python3
"""
sg_certify.py -- computer-assisted proof that the GEOMETRIC (short-time) skewness
    S^G_n = <u^3>/<u^2>^(3/2),   u = ut - <ut>,  -Lap ut = 2 in P_n, ut = 0 on dP_n,
of the regular n-gon is POSITIVE, for each n in a given range.

Method (result.md, Section 6): ut_N = (1-r^2)/2 + sum_{k=0}^K a_k Re((w/S)^(nk)) (floating-point
least-squares coefficients, then exact binary numbers), u_N = ut_N - <ut_N>, f = u - u_N = E - <E>,
E = ut - ut_N harmonic with boundary values -ut_N.  Then (Lemma 3)
    ||f||_{L2(P)} <= ||E||_{L2(P)} <= eps_E := sqrt((1+R)/2) ||ut_N||_{L2(dP)},   |f| <= osc_dP(ut_N),
and, since int f = 0,
    | int u^3 - int u_N^3 | <= 3 eps_E sd(u_N^2) + 3 Uinf eps_E^2 + osc eps_E^2,
with sd(phi) = ||phi - <phi>||_{L2(P)} and Uinf >= sup_P |u_N|.  We check int u_N^3 > that bound.
Also prints a rigorous enclosure of S^G_n (circumradius-independent).

Usage: python3 sg_certify.py nmin nmax [K]      (default 4 12)
Runtime: seconds to minutes per n (grows with K n).  Exit 1 on any FAIL.
"""
import sys, os, json, math, time
import numpy as np
from flint import arb, arb_poly, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from polylib import Terms, PF, range_on, ub
from lt_certify import float_eval

ctx.prec = 320


def default_K(n):
    return {4: 40, 5: 40, 6: 40, 7: 36, 8: 32}.get(n, 30)


def fit_ut(n, K, M=6000):
    h = math.pi / n
    tf = math.tan(h)
    S = float(1 / math.cos(h))
    sgrid = (1 - np.cos(np.linspace(0, np.pi, M))) / 2
    tau = tf * sgrid
    wts = np.gradient(tau)
    W = np.sqrt(wts)
    A = np.array([float_eval({(0, n * k): S ** (-n * k)}, tau) for k in range(0, K + 1)]).T
    rhs = -(1 - (1 + tau ** 2)) / 2
    a_coef, *_ = np.linalg.lstsq(A * W[:, None], rhs * W, rcond=None)
    return S, a_coef


def certify_sg(n, K=None, pieces=400, verbose=True):
    K = default_K(n) if K is None else K
    t0 = time.time()
    S, a_coef = fit_ut(n, K)
    Sa = arb(S)
    t = (arb.pi() / n).tan()
    R = 1 / (arb.pi() / n).cos()
    area = n * t
    ut = Terms({(0, 0): arb(0.5), (1, 0): arb(-0.5)})
    for k in range(0, K + 1):
        ut.add(0, n * k, arb(float(a_coef[k])) / Sa ** (n * k))
    UT = ut.to_pf()
    m_u = 2 * n * UT.int_T(t) / area
    U = UT.add_const(-m_u)
    bu = UT.boundary()
    rmin, rmax = range_on(bu, arb(0), t, pieces)
    delta = arb(max(abs(rmin).upper(), abs(rmax).upper()))
    osc = arb((rmax - rmin).upper())
    bu2 = 2 * n * (bu * bu).integral()(t)
    epsE = ((1 + R) / 2 * bu2).sqrt()
    U2 = U * U
    i2 = 2 * n * U2.int_T(t)
    i3 = 2 * n * (U2 * U).int_T(t)
    i4 = 2 * n * (U2 * U2).int_T(t)
    sd_U2 = (i4 - i2 * i2 / area).sqrt()
    Uinf = arb(max((m_u + delta).upper(), (R * R / 2 + delta - m_u).upper()))
    Err3 = 3 * epsE * sd_U2 + 3 * Uinf * epsE * epsE + osc * epsE * epsE
    ok = bool(i3 - Err3 > 0)
    # enclosure of S^G: |int u^2 - int u_N^2| <= 2 eps_E sd(u_N) + eps_E^2  (int f = 0)
    sd_U = (i2 - (2 * n * U.int_T(t)) ** 2 / area).sqrt()
    Err2 = 2 * epsE * sd_U + epsE * epsE
    i3e = arb(i3.mid(), (Err3 + i3.rad()).upper())
    i2e = arb(i2.mid(), (Err2 + i2.rad()).upper())
    SG = (i3e / area) / (i2e / area) ** 1.5
    dt = time.time() - t0
    res = dict(n=n, K=K, delta=float(delta.upper()), L2_E=float(epsE.upper()), int_uN3=i3.str(10),
               Err3=float(Err3.upper()), lower=float((i3 - Err3).lower()), SG=SG.str(6),
               a_coef=[float(x).hex() for x in a_coef], S=S.hex(), PASS=ok, seconds=round(dt, 1))
    if verbose:
        print(f"n={n:3d} K={K:3d} delta={res['delta']:.2e} ||E||<={res['L2_E']:.2e} int u_N^3={float(i3.mid()):.6e} "
              f"Err<={res['Err3']:.2e}  S^G in {SG.str(5)}  [{dt:.1f}s] {'PASS' if ok else 'FAIL'}", flush=True)
    return res


def main():
    nmin = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    Kf = int(sys.argv[3]) if len(sys.argv) > 3 else None
    print(f"sg_certify.py: certify S^G_n > 0 for {nmin} <= n <= {nmax} (prec {ctx.prec})")
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
