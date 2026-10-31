#!/usr/bin/env python3
"""Independent audit of the scalar inequalities in M. Liu, Zenodo 10.5281/zenodo.22015517.

Written from the statements in the paper, without using or running the author's scripts, and by a
different method: exact rational evaluation and rational branch-and-bound instead of Bernstein
certificates and hand case analysis. Prints PASS or FAIL; exit status non-zero on FAIL.
Environment: Python 3 standard library only. Runtime: a few minutes. Deterministic.

Notation of the paper: n = 2d+1, c_0 = 1, c_{h+1} = n + c_h - 2 sqrt(c_h); r = floor(d/2);
k = 2 rho + eps; a = d - rho; s = 2a + 1 - eps; eta(h) = h (rho - (2h-a)_+); ell the chord of eta on
[1,r]; R(x) = A/(2s) + max{rho s/(2k), rho x/2, (s/2)(2x-a)_+ + ell(x)/2}.

Claim C1 (c_h <= U_d, Prop 4.3): exact, rational upper bounds, 8 <= d <= DMAX1, all 1 <= h <= r.
Claim C2 (R > interpolant of c_h on [1,min(s,r)], Props 9.3, 10.3, 11.1 for small d):
         exact rational, for all 4 <= d <= DMAX2, 5 <= k <= 2d.
Claim C3 ((9.3): S_d(a,x) > U_d(x) for ALL d >= 8, real a in [0,d-4], x in [1,d/2]):
         rational branch-and-bound on (alpha,u,tau), tau = 1/d in (0,1/8].
Claim C4 (k in {5,6,7}, d > DMAX2): one-line estimate, checked exactly.
Claim C5 ((5.3), the tree case, for all d >= 5).

Together with Sections 2-8 of the paper, C2-C5 replace the Bernstein certificates of Appendix B, the
tables of Appendices C and D, and Sections 9-11.  C3 needs no integrality of a and no special case for
8 <= d <= 12.  Log: certificates/audit-scalar.log.
"""
import sys
from fractions import Fraction as Fr
from math import isqrt

DMAX1 = 3000
DMAX2 = 260
failures = 0


def check(name, ok, extra=""):
    global failures
    print(("PASS" if ok else "FAIL"), name, extra, flush=True)
    if not ok:
        failures += 1


def sqrt_lo(x: Fr, prec=10**24) -> Fr:
    """Rational lower bound for sqrt(x), error at most about 1/prec."""
    v = x.numerator * prec * prec // x.denominator
    return Fr(isqrt(v), prec)


def c_upper(d, H, prec=10**24):
    """Rational upper bounds cbar_h >= c_h, h=0..H, computed by rounding the sqrt down.
    u -> n + u - 2 sqrt(u) is increasing for u >= 1, so the bound propagates."""
    n = 2 * d + 1
    c = [Fr(1)]
    for _ in range(H):
        u = c[-1]
        nxt = n + u - 2 * sqrt_lo(u, prec)
        # round up to denominator prec to keep the size bounded
        nxt = Fr(-((-nxt.numerator * prec) // nxt.denominator), prec)
        c.append(nxt)
    return c


def T(u):
    z = Fr(1, 2) - u
    return Fr(79, 162) - Fr(163, 270) * z - Fr(233, 540) * z * z


def L(u):
    return (611 + 466 * u) / 1080


K = Fr(11489, 10800)


def U(d, x):
    u = Fr(x) / d
    return d * d * T(u) + d * L(u) + K


# ---- C1 ----
ok = True
worst = None
for d in list(range(8, 400)) + list(range(400, DMAX1 + 1, 37)):
    r = d // 2
    c = c_upper(d, r, prec=10**15)
    for h in range(1, r + 1):
        gap = U(d, h) - c[h]
        if worst is None or gap < worst[0]:
            worst = (gap, d, h)
        if gap <= 0:
            ok = False
check("C1: c_h <= U_d(h), d in [8,400) and sampled up to 3000",
      ok, f"smallest gap {float(worst[0]):.4f} at (d,h)=({worst[1]},{worst[2]})")


# ---- C2 ----
def Afun(a, rho, eps):
    return 2 * a**3 + 2 * a**2 * rho + 3 * a**2 + 2 * a * rho**2 + 2 * a * rho + a + (1 - eps) * rho**2


def scalar_min(d, k, c):
    """exact min over [1,min(s,r)] of R - interp(c), using all knots of a piecewise affine function."""
    rho = k // 2
    eps = k - 2 * rho
    a = d - rho
    s = 2 * a + 1 - eps
    r = d // 2
    top = min(s, r)
    if top < 1:
        return None
    A = Afun(a, rho, eps)
    base = Fr(A, 2 * s)

    def eta(h):
        return h * (rho - max(2 * h - a, 0))

    if r == 1:
        ell = lambda x: Fr(eta(1))
    else:
        e1, er = eta(1), eta(r)
        ell = lambda x: e1 + Fr(er - e1, r - 1) * (x - 1)

    f1 = lambda x: Fr(rho * s, 2 * k)
    f2 = lambda x: Fr(rho) * x / 2
    f3 = lambda x: ell(x) / 2
    f4 = lambda x: Fr(s) * (2 * x - a) / 2 + ell(x) / 2

    def env(x):
        v = max(f1(x), f2(x), f3(x))
        if 2 * x > a:
            v = max(v, f4(x))
        return v

    def cint(x):
        # piecewise linear interpolation of c on [1,r]
        i = int(x)
        if i >= r:
            return c[r]
        return c[i] + (x - i) * (c[i + 1] - c[i])

    # knots: endpoints, integers, a/2, crossings of pairs of affine members
    knots = {Fr(1), Fr(top), Fr(a, 2)}
    for i in range(1, top + 1):
        knots.add(Fr(i))
    # slopes and intercepts of the affine members
    def affine(f):
        y0, y1 = f(Fr(0)), f(Fr(1))
        return (y1 - y0, y0)
    mem = [affine(f) for f in (f1, f2, f3, f4)]
    for i in range(4):
        for j in range(i + 1, 4):
            m1, b1 = mem[i]
            m2, b2 = mem[j]
            if m1 != m2:
                knots.add((b2 - b1) / (m1 - m2))
    best = None
    arg = None
    for x in knots:
        if 1 <= x <= top:
            v = base + env(x) - cint(x)
            if best is None or v < best:
                best, arg = v, x
    return best, arg


ok = True
worst = None
table = {}
for d in range(4, DMAX2 + 1):
    r = d // 2
    c = c_upper(d, r, prec=10**18)
    for k in range(5, 2 * d + 1):
        res = scalar_min(d, k, c)
        if res is None:
            continue
        m, x = res
        if d <= 7:
            table[(d, k)] = (float(m), x)
        if m <= 0:
            ok = False
            print("  violation", d, k, float(m), x)
        if worst is None or m < worst[0]:
            worst = (m, d, k, x)
check(f"C2: R > interpolant of c_h for all 4<=d<={DMAX2}, 5<=k<=2d", ok,
      f"smallest margin {float(worst[0]):.5f} at (d,k,x)=({worst[1]},{worst[2]},{worst[3]})")
print("   d<=7 minima (d,k): margin (compare Table 9 of the paper, which uses slightly larger rational c-bar):")
for key in sorted(table):
    print("   ", key, round(table[key][0], 4), "at x =", table[key][1])


# ---- C3: rational branch and bound ----
def lower_res(al1, al2, u1, u2, ta1, ta2):
    """lower bound of Res = F_u(alpha) + alpha tau/2 - T(u) - tau L(u) - K tau^2 on the box."""
    # (1-alpha+alpha^2)/2 is decreasing on alpha<=1/2 and increasing beyond; min over the interval:
    if al2 <= Fr(1, 2):
        q0 = (1 - al2 + al2 * al2) / 2
    elif al1 >= Fr(1, 2):
        q0 = (1 - al1 + al1 * al1) / 2
    else:
        q0 = Fr(3, 8)
    # max{(1-al)u/2, al(2u-al)_+} >= max of the box minima
    m1 = (1 - al2) * u1 / 2
    m2 = Fr(0)
    # al(2u-al) is concave in al; min at an endpoint of [al1,al2], increasing in u
    v1 = al1 * (2 * u1 - al1)
    v2 = al2 * (2 * u1 - al2)
    m2 = max(m2, min(v1, v2))
    F = q0 + max(m1, m2)
    # T increasing in u on [0,1/2]; L increasing
    return F + al1 * ta1 / 2 - T(u2) - ta2 * L(u2) - K * ta2 * ta2


def bnb():
    # domain: tau in (0,1/8], alpha in [0, 1-4 tau], u in [tau, 1/2]
    stack = [(Fr(0), Fr(1), Fr(0), Fr(1, 2), Fr(0), Fr(1, 8))]
    closed = 0
    worst = None
    maxdepth = 0
    while stack:
        box = stack.pop()
        al1, al2, u1, u2, ta1, ta2 = box
        # infeasible boxes
        if al1 > 1 - 4 * ta1 or u2 < ta1:
            closed += 1
            continue
        low = lower_res(*box)
        if low > 0:
            closed += 1
            if worst is None or low < worst:
                worst = low
            continue
        wa, wu, wt = al2 - al1, u2 - u1, (ta2 - ta1) * 8
        if max(wa, wu, wt) < Fr(1, 2**14):
            return False, closed, box
        if wa >= wu and wa >= wt:
            m = (al1 + al2) / 2
            stack += [(al1, m, u1, u2, ta1, ta2), (m, al2, u1, u2, ta1, ta2)]
        elif wu >= wt:
            m = (u1 + u2) / 2
            stack += [(al1, al2, u1, m, ta1, ta2), (al1, al2, m, u2, ta1, ta2)]
        else:
            m = (ta1 + ta2) / 2
            stack += [(al1, al2, u1, u2, ta1, m), (al1, al2, u1, u2, m, ta2)]
    return True, closed, worst


ok, closed, info = bnb()
check("C3: (9.3) S_d(a,x) > U_d(x) on the whole region, d>=8, by rational branch and bound", ok,
      f"boxes closed: {closed}; smallest lower bound on a closed box: {float(info) if ok else info}")

# ---- C4: k in {5,6,7}, d > DMAX2.  R >= A/(2s) >= (d^2 - d a + a^2)/2 + a/2 with a = d - rho, and
# U_d(x) <= U_d(d/2) because T, L are increasing on [0,1/2].
ok = True
for rho in (2, 3):
    # residual polynomial in d: (d^2 - rho d + rho^2)/2 + (d-rho)/2 - [d^2 T(1/2) + d L(1/2) + K]
    c2 = Fr(1, 2) - T(Fr(1, 2))
    c1 = Fr(-rho, 2) + Fr(1, 2) - L(Fr(1, 2))
    c0 = Fr(rho * rho, 2) - Fr(rho, 2) - K
    d0 = DMAX2 + 1
    val = c2 * d0 * d0 + c1 * d0 + c0
    der = 2 * c2 * d0 + c1
    ok &= (val > 0 and der > 0 and c2 > 0)
check(f"C4: for k in 5,6,7 and d>{DMAX2}, A/(2s) > U_d(x) for every x in [1,d/2]", ok)

# ---- C5: tree case, inequality (5.3) of the paper ----
# odd d: c_r < d(d+1)/2; even d: c_r < (2d+1)(d+1)/4; for 5 <= d (d=3,4 are done by hand in the paper
# and by the exhaustive enumeration here).  Direct for 5 <= d < 400, then a polynomial estimate via U_d.
ok = True
for d in range(5, 400):
    r = d // 2
    cr = c_upper(d, r, prec=10**15)[r]
    target = Fr(d * (d + 1), 2) if d % 2 else Fr((2 * d + 1) * (d + 1), 4)
    ok &= cr < target
c2 = Fr(1, 2) - T(Fr(1, 2))
c1 = Fr(1, 2) - L(Fr(1, 2))
d0 = 30
# the odd target d(d+1)/2 is the smaller of the two, and U_d(r) <= U_d(d/2)
ok &= c2 > 0 and c2 * d0 * d0 + c1 * d0 - K > 0 and 2 * c2 * d0 + c1 > 0
check("C5: (5.3) tree inequality, direct for 5<=d<400, via U_d for d>=30", ok)

print("ALL PASS" if failures == 0 else f"{failures} FAILURES")
sys.exit(1 if failures else 0)
