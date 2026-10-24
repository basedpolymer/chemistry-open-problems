#!/usr/bin/env python3
"""
alltime_windows.py -- explicit certified time windows for the SIGN of the skewness at finite times
(result.md, Section 7, Theorem D).  Partial result towards "positivity for all t > 0".

For cross-sectionally uniform initial data with a symmetric longitudinal profile (e.g. a delta release)
and ANY Peclet number, sign Sk(tau) = sign K(tau) (Proposition 3), where
    K(tau) = 2 int_0^tau int_0^sigma < (e^{(sigma-s)Lap} u) u A(s) > ds dsigma,  A(s) = int_0^s e^{r Lap} u dr,
<.> = cross-sectional average, Lap = Neumann Laplacian, tau = time in units of a^2/kappa (a = apothem).
  Short times (Lemma 11):  |K(tau) - <u^3> tau^3/3| <= C_a tau^(7/2) + C_b tau^4,
      C_a = (32/105) c0 ||grad u|| ||u^2||,  C_b = (pi/48) c0^2 ||u||_inf ||grad u||^2,
      c0 = sup_{x>0} (1-e^{-x})/sqrt(x) <= 0.64 (checked below);  norms w.r.t. the average.
  Long times (Lemma 12):   |K(tau) - <u g2> tau| <= a1/mu2 + 2 b1/mu2^2,
      a1 = ||u|| ||g2||, ||g2|| <= 2 sd(u g1)/mu2, b1 = ||u|| ||u||_inf ||g1||, <u g2> = 2 <u g1^2>,
      mu2 = first nonzero Neumann eigenvalue (square: pi^2/4 exactly; others: Payne-Weinberger pi^2/D^2).
Output: t0 (K > 0 on (0, t0]) and T1 (sign K = sign <u g2> on [T1, inf)) for n = 3 (exact data) and n = 4..8.
Runtime: ~1 min.  PASS/FAIL per claim, exit 1 on failure.
"""
import sys, os, math
import sympy as sp
from flint import arb, fmpq, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ok_all = True


def check(name, cond):
    global ok_all
    print(("PASS" if cond else "FAIL") + ": " + name, flush=True)
    ok_all &= bool(cond)


ctx.prec = 200
PI = arb.pi()
# --- c0 = sup (1-e^{-x})/sqrt(x) <= C0 ---------------------------------------------------------------
C0 = arb('0.64')
okc = True
# x in (0, 0.4]: (1-e^{-x})/sqrt(x) <= sqrt(x) <= sqrt(0.4) < 0.6383 ; x >= 4: <= 1/2
N = 3600
for i in range(N):
    a = arb('0.4') + arb('3.6') * i / N
    b = arb('0.4') + arb('3.6') * (i + 1) / N
    x = a.union(b)
    v = (1 - (-x).exp()) / x.sqrt()
    if not (v < C0):
        okc = False
check("c0 = sup_{x>0} (1-e^{-x})/sqrt(x) <= 0.64", okc and arb('0.4').sqrt() < C0)


def windows(name, avg_u3_lo, grad_u, u2norm, uinf, unorm, g1norm, sd_ug1, ug1sq_enc, mu2):
    """All inputs are arb balls (upper bounds except avg_u3_lo = lower bound of <u^3> > 0)."""
    Ca = arb(32) / 105 * C0 * grad_u * u2norm
    Cb = PI / 48 * C0 * C0 * uinf * grad_u * grad_u
    c = avg_u3_lo / 3
    x = (-Ca + (Ca * Ca + 4 * Cb * c).sqrt()) / (2 * Cb)
    t0 = arb(fmpq(int(float(x.lower()) ** 2 * 0.999 * 1e12), 10 ** 12))   # exact rational slightly below x^2
    marg = c - Ca * t0.sqrt() - Cb * t0
    ok0 = bool(marg > 0)
    g2 = 2 * sd_ug1 / mu2
    a1 = unorm * g2
    b1 = unorm * uinf * g1norm
    const = a1 / mu2 + 2 * b1 / (mu2 * mu2)
    ug2 = 2 * ug1sq_enc
    if ug2 > 0:
        sgn = 'positive'
        T1 = const / ug2
    elif ug2 < 0:
        sgn = 'negative'
        T1 = const / (-ug2)
    else:
        sgn = 'undetermined'
        T1 = arb('inf')
    T1u = float(arb(T1.upper()))
    check(f"{name}: Sk(tau) > 0 for 0 < tau <= t0 = {float(t0):.4e}", ok0)
    check(f"{name}: Sk(tau) {sgn} for tau >= T1 = {T1u:.4g}", sgn != 'undetermined')
    return float(t0), T1u, sgn


rows = []
# --- n = 3, exact polynomial data (triangle with inradius 1: average over the triangle) -------------
y, z = sp.symbols('y z', real=True)
s3 = sp.sqrt(3)
ut = (1 + y) * ((2 - y) ** 2 - 3 * z ** 2) / 6
R_ = sp.Rational
g1 = (y ** 6 / 960 - y ** 5 / 96 - y ** 4 * z ** 2 / 64 + y ** 4 / 32 + y ** 3 * z ** 2 / 48 + y ** 2 * z ** 4 / 64
      + y ** 2 * z ** 2 / 16 - R_(11, 120) * y ** 2 + y * z ** 4 / 32 - z ** 6 / 960 + z ** 4 / 32 - R_(11, 120) * z ** 2)


def avg(f):
    I = sp.integrate(sp.integrate(sp.expand(f), (z, -(2 - y) / s3, (2 - y) / s3)), (y, -1, 2))
    return sp.nsimplify(sp.simplify(I / (3 * s3)))


u = sp.expand(ut - avg(ut))
g1 = sp.expand(g1 - avg(g1))
q = {k: avg(v) for k, v in dict(u2=u ** 2, u3=u ** 3, u4=u ** 4, ut=ut, g12=g1 ** 2, ug=u * g1, ug2=u * g1 ** 2,
                                 ugsq=(u * g1) ** 2).items()}
check("triangle: <u^3> = 27/19250, <u g1^2> = -7101/476476000",
      q['u3'] == R_(27, 19250) and q['ug2'] == R_(-7101, 476476000))


def A(r):
    r = sp.Rational(r)
    return arb(fmpq(int(r.p), int(r.q)))


tri = dict(avg_u3_lo=A(q['u3']), grad_u=(2 * A(q['ut'])).sqrt(), u2norm=A(q['u4']).sqrt(), uinf=A(R_(11, 30)),
           unorm=A(q['u2']).sqrt(), g1norm=A(q['g12']).sqrt(), sd_ug1=(A(q['ugsq']) - A(q['ug']) ** 2).sqrt(),
           ug1sq_enc=A(q['ug2']), mu2=PI * PI / 12)       # Payne-Weinberger, diameter 2 sqrt(3)
rows.append(("n=3",) + windows("n=3 (equilateral triangle)", **tri))

# --- n >= 4: certified enclosures from the a-posteriori machinery ----------------------------------------
from lt_certify import certify
from sg_certify import certify_sg
for n, K in [(4, 20), (5, 16), (6, 12), (8, 10)]:
    r = certify(n, K=K, verbose=False, raw=True)['raw']
    rs = certify_sg(n, K=16, verbose=False, raw=True)['raw']
    P = r['area']
    sP = P.sqrt()
    eps = r['epsE']
    ut_avg_hi = r['m_u'] + eps / sP
    ut_avg_lo = r['m_u'] - eps / sP
    grad_u = arb((2 * ut_avg_hi).upper()).sqrt()
    uinf = arb(max(ut_avg_hi.upper(), (r['R'] * r['R'] / 2 - ut_avg_lo).upper()))
    u2norm = (r['intU4'].sqrt() + eps * (2 * r['Uinf'] + r['osc'])) / sP
    unorm = (r['intU2'].sqrt() + eps) / sP
    g1norm = (r['intG2'].sqrt() + r['omega']) / sP
    sd_ug1 = (r['sd_UG'] + r['osc'] * (r['intG2'].sqrt() + r['omega']) + r['Uinf'] * r['omega']) / sP
    ug1sq = arb((r['I_N'] - r['Err']).lower()) / P            # lower bound of <u g1^2>
    avg_u3_lo = arb(rs['L'].lower()) / P
    if n == 4:
        mu2 = PI * PI / 4                                       # square [-1,1]^2, Neumann, exact
    else:
        D = 2 * r['R'] if n % 2 == 0 else 2 * r['R'] * (PI / (2 * n)).cos()
        mu2 = PI * PI / (D * D)
    rows.append((f"n={n}",) + windows(f"n={n}", avg_u3_lo, grad_u, u2norm, uinf, unorm, g1norm, sd_ug1, ug1sq, mu2))

print("\n  n  | Sk > 0 on (0, t0] | sign on [T1, inf)")
for row in rows:
    print(f"  {row[0]:4s} | t0 = {row[1]:.3e}     | {row[3]} for tau >= {row[2]:.3g}")
print("PASS (all claims)" if ok_all else "FAIL (some claim)")
sys.exit(0 if ok_all else 1)
