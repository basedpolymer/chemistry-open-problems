"""
ivl.py -- minimal rigorous interval arithmetic on numpy arrays (Agent P, Claude Opus 5.5, effort xhigh).

Every elementary operation is computed in IEEE-754 binary64 round-to-nearest and the result is widened by one ulp
in each direction (np.nextafter), which encloses the exact result of the operation on the endpoints.  Only the
operations needed by certify_ss.py are provided (+, -, *, / by positive intervals, sums).
"""
import numpy as np

INF = np.inf


def dn(a):
    return np.nextafter(a, -INF)


def up(a):
    return np.nextafter(a, INF)


class I:
    __array_priority__ = 1000

    def __init__(self, lo, hi=None):
        lo = np.asarray(lo, dtype=float)
        self.lo = lo
        self.hi = lo.copy() if hi is None else np.asarray(hi, dtype=float)

    @staticmethod
    def point(v):
        """Enclosure of exact floats (no widening needed)."""
        v = np.asarray(v, dtype=float)
        return I(v, v.copy())

    def __add__(self, o):
        o = _as(o)
        return I(dn(self.lo + o.lo), up(self.hi + o.hi))

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, o):
        o = _as(o)
        return I(dn(self.lo - o.hi), up(self.hi - o.lo))

    def __rsub__(self, o):
        return _as(o) - self

    def __mul__(self, o):
        o = _as(o)
        p = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return I(dn(np.minimum.reduce(p)), up(np.maximum.reduce(p)))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = _as(o)
        if np.any(o.lo <= 0):
            raise ValueError('division by an interval not strictly positive')
        p = [self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi]
        return I(dn(np.minimum.reduce(p)), up(np.maximum.reduce(p)))

    def __getitem__(self, k):
        return I(self.lo[k], self.hi[k])

    def mag(self):
        return np.maximum(np.abs(self.lo), np.abs(self.hi))

    def mid(self):
        return 0.5 * (self.lo + self.hi)

    def rad_up(self):
        return up(0.5 * (self.hi - self.lo))

    @property
    def shape(self):
        return self.lo.shape


def _as(o):
    return o if isinstance(o, I) else I.point(o)


def isum(x, axis):
    """Rigorous sum along an axis (sequential, each partial sum widened)."""
    lo = np.moveaxis(x.lo, axis, 0); hi = np.moveaxis(x.hi, axis, 0)
    slo = lo[0].copy(); shi = hi[0].copy()
    for k in range(1, lo.shape[0]):
        slo = dn(slo + lo[k]); shi = up(shi + hi[k])
    return I(slo, shi)


def frac_enclosure(fr):
    """Interval of floats enclosing a Fraction (or int)."""
    v = float(fr)
    from fractions import Fraction
    F = Fraction(v)
    lo = v if F <= fr else dn(v)
    hi = v if F >= fr else up(v)
    return lo, hi
