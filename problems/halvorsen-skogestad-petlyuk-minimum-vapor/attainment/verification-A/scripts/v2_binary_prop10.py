#!/usr/bin/env python3
"""V2 (independent verifier): check of sections.md Proposition 10 (N = 2), used by theorem.md §7 for N = 2.

NUMERICAL/SYMBOLIC EVIDENCE ONLY; the proof is the verifier's line-by-line re-derivation in verification.md.
Written from scratch by V2 (no code of S, D, L, A or C is imported).

(1) sympy: with 1 - q eliminated by the feed equation alpha z/(alpha - theta) + (1 - z)/(1 - theta) = 1 - q,
    f(s) := alpha s (s + 1 - q - z) - (s + 1 - q)(s + 1 - z) factors as (alpha - 1)(s - H_1)(s - s_-) with
    H_1 = (1 - z)/(theta - 1), s_- = 1/(alpha - 1) - z theta/(alpha - theta); and the lower bound
    H_1 - s_- >= z alpha (theta - 1)/((alpha - theta)(alpha - 1)) (exact difference = (1-z)/(theta-1) - (1-z)/(alpha-1)).
    Also: x_T-hat < x_B-hat  <=>  f(s) > 0 (cleared denominators), and H_1 + 1 - q = alpha z/(alpha - theta).
(2) mpmath (60 digits): the shooting construction of the proof (bisection on eta = x_B) for several (alpha, z, q)
    including q < 0 and q > 1, s = H_1 + delta, equal and UNEQUAL n_T, n_B; for each solution, every stage, the feed
    stage, the reboiler and the condenser balance is recomputed from the streams, every mole fraction is checked to be
    in (0, 1) (Theorem 1 claims positivity also for N = 2), and the purities are printed; R = s.
Environment: Python 3, sympy, mpmath (WSL Fedora sci env). Runtime: about 10 s. Prints PASS/FAIL; exit 1 on FAIL.
"""
import sys
import sympy as sp
import mpmath as mp

FAILS = []


def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- (1) symbolic
a, z, th, s = sp.symbols('alpha z theta s', positive=True)
omq = a * z / (a - th) + (1 - z) / (1 - th)            # 1 - q
q = 1 - omq
H1 = (1 - z) / (th - 1)
sm = 1 / (a - 1) - z * th / (a - th)
f = a * s * (s + 1 - q - z) - (s + 1 - q) * (s + 1 - z)
check('f(s) = (alpha-1)(s-H1)(s-s_-) identically', sp.simplify(sp.expand(f - (a - 1) * (s - H1) * (s - sm))) == 0)
check('H1 + 1 - q = alpha z/(alpha-theta)', sp.simplify(H1 + omq - a * z / (a - th)) == 0)
diff = sp.simplify(H1 - sm - z * a * (th - 1) / ((a - th) * (a - 1)))
check('H1 - s_- - z alpha (theta-1)/((alpha-theta)(alpha-1)) = (1-z)/(theta-1) - (1-z)/(alpha-1)',
      sp.simplify(diff - ((1 - z) / (th - 1) - (1 - z) / (a - 1))) == 0)
# x_T-hat < x_B-hat <=> V_T/L_T < alpha s/L_B <=> V_T L_B < alpha s L_T (L_T, L_B > 0); difference = f
VT, LT, LB = s + 1 - q, s + 1 - q - z, s + 1 - z
check('alpha s L_T - V_T L_B = f(s)', sp.simplify(a * s * LT - VT * LB - f) == 0)

# ---------------------------------------------------------------- (2) numerical shooting
mp.mp.dps = 60


def theta_of(al, zz, qq):
    g = lambda t: al * zz / (al - t) + (1 - zz) / (1 - t) - (1 - qq)
    lo, hi = mp.mpf(1), mp.mpf(al)
    for _ in range(400):          # g increases from -inf to +inf on (1, alpha)
        mid = (lo + hi) / 2
        if g(mid) > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def profiles(al, zz, qq, sv, nT, nB, eta):
    e = lambda x: al * x / (1 + (al - 1) * x)
    einv = lambda y: y / (al - (al - 1) * y)
    D, B = zz, 1 - zz
    VT, LB = sv + 1 - qq, sv + 1 - zz
    LT = VT - D
    xD = 1 - (1 - zz) * eta / zz
    xb = [eta]                               # reboiler liquid, then stages 1..nB, then the feed stage
    for _ in range(nB + 1):
        xb.append((sv * e(xb[-1]) + B * eta) / LB)
    xt = [einv(xD)]                          # top stage nT liquid, going down to the feed stage
    for _ in range(nT):
        xt.append(einv((LT * xt[-1] + D * xD) / VT))
    return xb, xt, xD, (VT, LT, LB, D, B, e)


def solve(al, zz, qq, sv, nT, nB):
    eta1 = mp.mpf('1e-55')
    while True:
        xb, xt, _, _ = profiles(al, zz, qq, sv, nT, nB, eta1)
        if xb[-1] > xt[-1]:
            break
        eta1 *= 2
        if eta1 > min(1, zz / (1 - zz)):
            return None
    lo, hi = mp.mpf(0), eta1
    for _ in range(260):
        mid = (lo + hi) / 2
        xb, xt, _, _ = profiles(al, zz, qq, sv, nT, nB, mid)
        if xb[-1] > xt[-1]:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def balances(al, zz, qq, sv, nT, nB, eta):
    """Rebuild the column (light component) and return the worst component-balance residual and min/max x."""
    xb, xt, xD, (VT, LT, LB, D, B, e) = profiles(al, zz, qq, sv, nT, nB, eta)
    xf = (xb[-1] + xt[-1]) / 2
    # liquid compositions bottom -> top: reboiler, nB stages, feed stage, nT stages; condenser liquid xD
    liq = xb[:-1] + [xf] + list(reversed(xt[:-1]))
    nb = nB
    res = []
    VB = sv
    # reboiler (unit 0): in L_B x_1 ; out V_B e(x_0) + B x_0
    res.append(LB * liq[1] - VB * e(liq[0]) - B * liq[0])
    for k in range(1, len(liq)):
        below, above = liq[k - 1], (liq[k + 1] if k + 1 < len(liq) else xD)
        if k <= nb:            # bottom-section stage
            r = LB * above + VB * e(below) - LB * liq[k] - VB * e(liq[k])
        elif k == nb + 1:      # feed stage (feed light flow z), liquid down L_B, vapour up V_T
            r = LT * above + VB * e(below) + zz - LB * liq[k] - VT * e(liq[k])
        else:                  # top-section stage; the liquid above the top stage is the reflux xD
            r = LT * above + VT * e(below) - LT * liq[k] - VT * e(liq[k])
        res.append(r)
    # condenser: vapour V_T e(x_top stage) -> liquid (L_T + D) xD
    res.append(VT * e(liq[-1]) - (LT + D) * xD)
    allx = liq + [xD]
    return max(abs(r) for r in res), min(allx), max(allx), abs(xb[-1] - xt[-1])


cases = [(2, '0.5', 1, '0.02'), ('1.5', '0.3', 0, '0.05'), (3, '0.7', '0.5', '0.01'), ('1.2', '0.5', '1.3', '0.05'),
         ('2.5', '0.4', '-0.3', '0.02'), (4, '0.1', '1.6', '0.03')]
for (al, zz, qq, dl) in cases:
    al, zz, qq, dl = mp.mpf(al), mp.mpf(zz), mp.mpf(qq), mp.mpf(dl)
    thv = theta_of(al, zz, qq)
    h1 = (1 - zz) / (thv - 1)
    sv = h1 + dl
    LT = sv + 1 - qq - zz
    check('case alpha=%s z=%s q=%s: L_T > 0, V_T > 0' % (mp.nstr(al, 3), mp.nstr(zz, 3), mp.nstr(qq, 3)),
          LT > 0 and sv + 1 - qq > 0, 'H1 = %s, s = %s' % (mp.nstr(h1, 8), mp.nstr(sv, 8)))
    prev = None
    for (nT, nB) in [(20, 20), (40, 40), (80, 80), (160, 160), (40, 120), (120, 40)]:
        eta = solve(al, zz, qq, sv, nT, nB)
        if eta is None:
            check('  n_T=%d n_B=%d: bracket found' % (nT, nB), False)
            continue
        worst, xmin, xmax, gap = balances(al, zz, qq, sv, nT, nB, eta)
        imp_D = (1 - zz) * eta / zz
        ok = worst < mp.mpf('1e-40') and xmin > 0 and xmax < 1
        check('  n_T=%3d n_B=%3d: steady state, all x in (0,1)' % (nT, nB), ok,
              'impurities x_B = %s, 1-x_D = %s | worst balance %s | min x %s' %
              (mp.nstr(eta, 3), mp.nstr(imp_D, 3), mp.nstr(worst, 2), mp.nstr(xmin, 3)))
    # purity -> 0 as n grows (equal lengths)
    etas = [solve(al, zz, qq, sv, n, n) for n in (40, 80, 160)]
    check('  impurity decreases with n (40, 80, 160)', etas[0] > etas[1] > etas[2],
          ', '.join(mp.nstr(x, 3) for x in etas))

print('=' * 80)
print('ALL PASS' if not FAILS else 'FAILURES: %s' % FAILS)
sys.exit(1 if FAILS else 0)
