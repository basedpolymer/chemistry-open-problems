"""
polylib.py -- exact/ball-arithmetic calculus for D_n-invariant polynomial functions on the
regular n-gon P_n with apothem 1 (sides tangent to the unit circle, side 0 = {x = 1, |y| <= t},
t = tan(pi/n), circumradius R = sec(pi/n)).

Representation used for rigorous computations (python-flint arb balls):

  F(x, y) = sum_d F_d(x, y),  F_d homogeneous of degree d,
  stored as {d: P_d} with P_d(tau) = F_d(1, tau)  (an arb_poly in tau).

Facts used (all elementary):
  * products: (FG)_d(1,tau) = sum_{d1+d2=d} P^F_{d1}(tau) P^G_{d2}(tau);
  * the half-sector T = {0 <= x <= 1, 0 <= y <= t x} is the cone over the segment
    {x = 1, 0 <= y <= t}, so for homogeneous F_d:  int_T F_d dA = (1/(d+2)) int_0^t P_d(tau) dtau;
  * a D_n-invariant function satisfies int_P F = 2n int_T F, and its restriction to the boundary
    is determined by its restriction to {x = 1, 0 <= tau <= t};
  * Euler: x F_x + y F_y = d F for F homogeneous of degree d, so on x = 1:
    dF_d/dx (1,tau) = d P_d(tau) - tau P_d'(tau)  (this is the outward normal derivative on side 0).

Functions are built from the polar monomials r^(2a) Re(w^m), w = x + i y, which are homogeneous of
degree 2a + m with P(tau) = (1 + tau^2)^a Re((1 + i tau)^m), and
  Laplace(r^(2a) Re w^m) = 4 a (a + m) r^(2a-2) Re w^m.
The class Terms keeps such linear combinations symbolically (exact Laplacian, exact particular
solutions); Terms.to_pf() converts to the homogeneous representation.

Environment: python3.13, python-flint 0.9.0.  No floating point is used in this module.
"""
from flint import arb, arb_poly, fmpz, ctx

_binom_cache = {}


def re_pow_poly(m):
    """arb_poly for Re((1 + i tau)^m) = sum_j (-1)^j C(m, 2j) tau^(2j) (exact integer coefficients)."""
    if m in _binom_cache:
        return _binom_cache[m]
    co = [0] * (m + 1)
    from math import comb
    for j in range(0, m // 2 + 1):
        co[2 * j] = comb(m, 2 * j) * (-1) ** j
    p = arb_poly(co)
    _binom_cache[m] = p
    return p


def r2_pow_poly(a):
    """arb_poly for (1 + tau^2)^a."""
    p = arb_poly([1])
    q = arb_poly([1, 0, 1])
    for _ in range(a):
        p = p * q
    return p


class PF:
    """Homogeneous-component representation {d: arb_poly P_d}."""

    def __init__(self, comps=None):
        self.c = {}
        if comps:
            for d, p in comps.items():
                self.c[d] = p

    def copy(self):
        return PF(dict(self.c))

    def __add__(self, o):
        r = dict(self.c)
        for d, p in o.c.items():
            r[d] = r[d] + p if d in r else p
        return PF(r)

    def __neg__(self):
        return PF({d: -p for d, p in self.c.items()})

    def __sub__(self, o):
        return self + (-o)

    def scale(self, s):
        return PF({d: p * s for d, p in self.c.items()})

    def __mul__(self, o):
        if not isinstance(o, PF):
            return self.scale(o)
        r = {}
        for d1, p1 in self.c.items():
            for d2, p2 in o.c.items():
                d = d1 + d2
                pr = p1 * p2
                r[d] = r[d] + pr if d in r else pr
        return PF(r)

    def add_const(self, c):
        r = dict(self.c)
        r[0] = r[0] + arb_poly([c]) if 0 in r else arb_poly([c])
        return PF(r)

    def boundary(self):
        """F(1, tau) as an arb_poly."""
        s = arb_poly([0])
        for p in self.c.values():
            s = s + p
        return s

    def dx_boundary(self):
        """dF/dx (1, tau) = sum_d (d P_d - tau P_d') as an arb_poly (outward normal derivative on side 0)."""
        s = arb_poly([0])
        tau = arb_poly([0, 1])
        for d, p in self.c.items():
            s = s + p * d - tau * p.derivative()
        return s

    def int_T(self, t):
        """Integral over the half-sector T (vertices (0,0), (1,0), (1,t))."""
        s = arb(0)
        for d, p in self.c.items():
            s += p.integral()(t) / (d + 2)
        return s


class Terms:
    """Linear combination sum c_{a,m} r^(2a) Re(w^m) with arb coefficients (keys (a, m))."""

    def __init__(self, d=None):
        self.d = {}
        if d:
            for k, v in d.items():
                self.d[k] = arb(v)

    def copy(self):
        return Terms(dict(self.d))

    def add(self, a, m, c):
        k = (a, m)
        self.d[k] = self.d.get(k, arb(0)) + c
        return self

    def __add__(self, o):
        r = self.copy()
        for k, v in o.d.items():
            r.add(k[0], k[1], v)
        return r

    def scale(self, s):
        return Terms({k: v * s for k, v in self.d.items()})

    def __sub__(self, o):
        return self + o.scale(-1)

    def laplacian(self):
        r = Terms()
        for (a, m), c in self.d.items():
            if a > 0:
                r.add(a - 1, m, c * (4 * a * (a + m)))
        return r

    def neg_lap_inverse(self):
        """Particular solution p with -Laplace(p) = self (term-by-term)."""
        r = Terms()
        for (a, m), c in self.d.items():
            r.add(a + 1, m, -c / (4 * (a + 1) * (a + 1 + m)))
        return r

    def to_pf(self):
        comps = {}
        for (a, m), c in self.d.items():
            if c == 0:
                continue
            p = r2_pow_poly(a) * re_pow_poly(m) * c
            deg = 2 * a + m
            comps[deg] = comps[deg] + p if deg in comps else p
        return PF(comps)

    def equals(self, o):
        """Rigorous check that two Terms are identical: every coefficient difference contains 0
        and has zero radius (coefficients are exact) -- or contains zero (if balls)."""
        keys = set(self.d) | set(o.d)
        ok = True
        for k in keys:
            diff = self.d.get(k, arb(0)) - o.d.get(k, arb(0))
            if not diff.contains(0):
                ok = False
        return ok


def sup_abs_on(poly, lo, hi, pieces):
    """Rigorous upper bound for max_{lo <= tau <= hi} |poly(tau)| (lo, hi arb; uses upper(hi)).
    On each piece [c - rho, c + rho] the Taylor expansion poly(c + s) = sum b_k s^k is computed in ball
    arithmetic and |poly| <= |b_0| + sum_{k>=1} |b_k| rho^k."""
    lo_u = arb(lo.lower()) if isinstance(lo, arb) else arb(lo)
    hi_u = arb(hi.upper()) if isinstance(hi, arb) else arb(hi)
    best = arb(0)
    L = hi_u - lo_u
    for i in range(pieces):
        a = lo_u + L * i / pieces
        b = lo_u + L * (i + 1) / pieces
        c = (a + b) / 2
        rho = arb(((b - a) / 2).upper())
        sh = poly(arb_poly([c, 1]))
        co = sh.coeffs()
        s = abs(co[0]) if co else arb(0)
        rk = arb(1)
        for k in range(1, len(co)):
            rk = rk * rho
            s = s + abs(co[k]) * rk
        u = arb(s.upper())
        if u > best:
            best = u
    return best


def range_on(poly, lo, hi, pieces):
    """Rigorous enclosure [m, M] of the range of poly on [lo, hi] (as two arbs lower/upper)."""
    lo_u = arb(lo.lower()) if isinstance(lo, arb) else arb(lo)
    hi_u = arb(hi.upper()) if isinstance(hi, arb) else arb(hi)
    L = hi_u - lo_u
    mn = None
    mx = None
    for i in range(pieces):
        a = lo_u + L * i / pieces
        b = lo_u + L * (i + 1) / pieces
        c = (a + b) / 2
        rho = arb(((b - a) / 2).upper())
        sh = poly(arb_poly([c, 1]))
        co = sh.coeffs()
        rest = arb(0)
        rk = arb(1)
        for k in range(1, len(co)):
            rk = rk * rho
            rest = rest + abs(co[k]) * rk
        lo_v = (co[0] - rest).lower()
        hi_v = (co[0] + rest).upper()
        mn = lo_v if mn is None or lo_v < mn else mn
        mx = hi_v if mx is None or hi_v > mx else mx
    return arb(mn), arb(mx)


def ub(x):
    """Upper bound (as exact arb) of a ball."""
    return arb(x.upper())


def lb(x):
    return arb(x.lower())


def sqrt_ub(x):
    return arb(arb(x.upper()).sqrt().upper())
