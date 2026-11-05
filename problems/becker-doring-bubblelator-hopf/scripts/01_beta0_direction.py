"""Case beta = 0 (size-independent removal, r = 0): certified sign of the direction coefficient.

Claims checked (each prints PASS/FAIL; exit code 1 on any FAIL):
  B0-1  For beta = 0 the closed form of eta_2 (paper, Prop. 5.2) agrees with the general
        integral-based enclosure of bdhopf.hopf_quantities at nu = 1/9, 1/2, 1, 9 (consistency).
  B0-2  eta_2(nu) < 0 for every nu > 0, i.e. for every omega = pi/(2(nu+1)) in (0, pi/2):
        - omega in [0, 0.05]          : monotone enclosures of sinc, tanc, Lg, atanc (limit nu -> oo
                                         included, eta_2 -> 1/(10 pi) - 3/20),
        - omega in [0.05, pi/2 - 0.05]: direct ball arithmetic,
        - eps = pi/2 - omega in [0, 0.05]: rewritten form, limit nu -> 0 included (eta_2 -> -1/36).
  B0-3  the limits: eta_2 -> -1/36 (nu -> 0) and eta_2 -> 1/(10 pi) - 3/20 (nu -> oo), by evaluating the
        two end-region formulas at the end points (enclosures of width < 1e-6 around the exact values).

Runtime: about 1 minute.  Environment: WSL Fedora, Python 3, python-flint 0.9.  Deterministic.
"""
import sys, os
from flint import arb, acb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bdhopf as B

ctx.prec = 128
PI = arb.pi()
FAIL = []


def report(name, ok, msg=''):
    print(('PASS' if ok else 'FAIL'), name, msg)
    if not ok:
        FAIL.append(name)


def iv(a, b):
    return arb(a).union(arb(b))


# ---- monotone helper functions on intervals [a, b] (a, b exact floats/arbs, 0 <= a <= b) ----
def sinc_iv(a, b):            # sin x / x, decreasing on [0, pi]
    hi = arb(1) if a == 0 else arb(a).sin() / arb(a)
    lo = arb(b).sin() / arb(b) if b > 0 else arb(1)
    return lo.union(hi)


def tanc_iv(a, b):            # tan x / x, increasing on [0, pi/2)
    lo = arb(1) if a == 0 else arb(a).tan() / arb(a)
    hi = arb(b).tan() / arb(b) if b > 0 else arb(1)
    return lo.union(hi)


def Lg_iv(ya, yb):            # log(1+y)/y, decreasing on [0, oo); ya, yb arb >= 0
    hi = arb(1) if float(ya.upper()) == 0 else None
    if hi is None:
        hi = (1 + ya.lower()).log() / ya.lower()
    lo = arb(1) if float(yb.upper()) == 0 else (1 + yb.upper()).log() / yb.upper()
    return lo.union(hi)


def atanc_iv(xa, xb):         # atan(x)/x, decreasing on [0, oo)
    hi = arb(1) if float(xa.upper()) == 0 else xa.lower().atan() / xa.lower()
    lo = arb(1) if float(xb.upper()) == 0 else xb.upper().atan() / xb.upper()
    return lo.union(hi)


def eta2_from(N, s, c_over, ReA, ImA, Ns, Ns2):
    """eta_2 = -(Re A / c)(1 + N s^2)/(4 N s |2i+A|^2) - 1/4 + (2 + Im A)/(4|2i+A|^2),
    with c_over = (Re A)/c supplied directly when needed."""
    D = ReA ** 2 + (2 + ImA) ** 2
    return -c_over * (1 + Ns2) / (4 * Ns * D) - arb(1) / 4 + (2 + ImA) / (4 * D)


def eta2_direct(w):
    """Ball evaluation on a ball w inside (0, pi/2)."""
    N = PI / (2 * w)
    s, c = w.sin(), w.cos()
    A = acb(c, 2 * s) ** (-N)
    return eta2_from(N, s, A.real / c, A.real, A.imag, N * s, N * s * s)


def eta2_small_omega(a, b):
    """Enclosure on omega in [a, b] subset [0, 0.05], a may be 0."""
    w = iv(a, b)
    sw = sinc_iv(a, b)                      # sin w / w
    s = iv(arb(a).sin(), arb(b).sin())      # sin w (increasing)
    c = iv(arb(b).cos(), arb(a).cos())
    tw = iv(arb(a).tan(), arb(b).tan())
    Ns = PI / 2 * sw
    Ns2 = Ns * s
    y = 3 * s * s
    ylo, yhi = 3 * arb(a).sin() ** 2, 3 * arb(b).sin() ** 2
    ReL = PI / 4 * (3 * s * sw) * Lg_iv(ylo, yhi)                 # N * 0.5 log(1+3 s^2)
    x_lo, x_hi = 2 * arb(a).tan(), 2 * arb(b).tan()
    ImL = PI * tanc_iv(a, b) * atanc_iv(x_lo, x_hi)              # N * atan(2 tan w)
    A = (-acb(ReL, ImL)).exp()
    N_dummy = None
    return eta2_from(N_dummy, s, A.real / c, A.real, A.imag, Ns, Ns2)


def eta2_near_half_pi(a, b):
    """Enclosure on eps = pi/2 - omega in [a, b] subset [0, 0.05]."""
    e = iv(a, b)
    N = 1 / (1 - 2 * e / PI)
    se = iv(arb(a).sin(), arb(b).sin())     # sin eps = c
    s = iv(arb(b).cos(), arb(a).cos())      # cos eps = s
    te_lo, te_hi = arb(a).tan(), arb(b).tan()
    ReW = -N / 2 * (1 - arb(3) / 4 * se * se).log()
    q = -N * (1 - tanc_iv(a, b) * atanc_iv(te_lo / 2, te_hi / 2) / 2)      # Im W / eps
    ImW = q * e
    # |Im W| range for sinc: |q| * eps <= |q|_max * b
    absImW_hi = float((abs(q) * arb(b)).upper())
    sincImW = sinc_iv(0.0, absImW_hi) if a == 0 else sinc_iv(0.0, absImW_hi)
    pref = arb(2) ** (-N) * ReW.exp()
    ReA_over_c = pref * sincImW * q / sinc_iv(a, b)
    ReA = pref * ImW.sin()
    ImA = -pref * ImW.cos()
    return eta2_from(N, s, ReA_over_c, ReA, ImA, N * s, N * s * s)


def finite(x):
    return x.is_finite() if hasattr(x, 'is_finite') else True


def close(x, y, tol=1e-9):
    """strict check: both finite, they overlap, and the enclosure is narrow"""
    return x.is_finite() and y.is_finite() and x.overlaps(y) and float(x.rad()) < tol


def certify_region(fn, a, b, depth=0, maxdepth=30, stats=None):
    v = fn(a, b)
    if v.is_finite() and v < 0:
        if stats is not None:
            stats['n'] += 1
            stats['lo'] = min(stats['lo'], float(v.lower()))
            stats['hi'] = max(stats['hi'], float(v.upper()))
        return True
    if depth >= maxdepth:
        print('   cannot certify on', a, b, v)
        return False
    m = 0.5 * (a + b)
    return certify_region(fn, a, m, depth + 1, maxdepth, stats) and certify_region(fn, m, b, depth + 1, maxdepth, stats)


def main():
    # B0-1 consistency with the integral-based routine
    B.TOL_BITS = 80
    for (p, qd) in [(1, 9), (1, 2), (1, 1), (9, 1)]:
        nu = arb(p) / qd
        N = nu + 1
        w = PI / (2 * N)
        t0 = w.tan()
        Q = B.hopf_quantities(arb(0), nu, t0)
        closed = eta2_direct(w)
        ok = Q['eta2'].overlaps(closed) and Q['eta0'].overlaps(w.cos() ** (N + 1) / w.sin()) \
            and Q['kappa0'].overlaps(w.cos() ** N) and abs(Q['G1'].real) < arb(10) ** -15
        report('B0-1 nu=%d/%d' % (p, qd), ok, 'eta2 = %s (closed form %s)' % (Q['eta2'].str(12), closed.str(12)))
    # B0-2 sign on (0, pi/2)
    lims = {}
    st = {'n': 0, 'lo': 1e9, 'hi': -1e9}
    okL = certify_region(eta2_small_omega, 0.0, 0.05, stats=st)
    report('B0-2a omega in [0, 0.05] (nu >= 30.4)', okL, 'accepted boxes=%d, eta2 in [%.5f, %.5f]' % (st['n'], st['lo'], st['hi']))
    lims['L'] = (st['lo'], st['hi'])
    st = {'n': 0, 'lo': 1e9, 'hi': -1e9}
    wa, wb = 0.05, float((PI / 2).mid()) - 0.05
    okM = certify_region(lambda a, b: eta2_direct(iv(a, b)), wa, wb, stats=st)
    report('B0-2b omega in [0.05, pi/2 - 0.05]', okM, 'accepted boxes=%d, eta2 in [%.5f, %.5f]' % (st['n'], st['lo'], st['hi']))
    lims['M'] = (st['lo'], st['hi'])
    st = {'n': 0, 'lo': 1e9, 'hi': -1e9}
    # eps in [0, pi/2 - wb]; pi/2 - wb >= 0.05 - 1e-16, use 0.0500001 to overlap
    okR = certify_region(eta2_near_half_pi, 0.0, 0.0500001, stats=st)
    report('B0-2c eps = pi/2 - omega in [0, 0.05] (nu <= 0.033)', okR, 'accepted boxes=%d, eta2 in [%.5f, %.5f]' % (st['n'], st['lo'], st['hi']))
    lims['R'] = (st['lo'], st['hi'])
    # overlap of the regions: pi/2 - wb must be <= 0.0500001
    report('B0-2d regions cover (0, pi/2)', bool(PI / 2 - arb(wb) < arb(0.0500001)))
    # limits
    l0 = eta2_near_half_pi(0.0, 0.0)
    linf = eta2_small_omega(0.0, 0.0)
    report('B0-3 limits', close(l0, arb(-1) / 36, 1e-6) and close(linf, 1 / (10 * PI) - arb(3) / 20, 1e-6),
           'nu->0: %s (=-1/36), nu->oo: %s (=1/(10pi)-3/20)' % (l0.str(10), linf.str(10)))
    print('summary: eta2 enclosure over all nu > 0 within [%.5f, %.5f]' % (min(v[0] for v in lims.values()), max(v[1] for v in lims.values())))
    if FAIL:
        print('FAILED:', FAIL)
        sys.exit(1)
    print('ALL PASS')


if __name__ == '__main__':
    main()
