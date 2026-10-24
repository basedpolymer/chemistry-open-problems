"""ad -- forward-mode automatic differentiation over Arb complex balls (acb), n real directions.

AD(v, d): value v (acb) and partial derivatives d[0..n-1] (acb) with respect to n REAL variables.
Used for the mean value form  f(X) in f(X_m) + sum_i (d_i f)(box) (X_i - X_m,i).
Special functions are given by their power series (with rigorous tail bounds) so that they are valid for balls
containing 0:  L(u) = log(1+u)/u,  M1(u) = L'(u),  Ex(b) = (e^b - 1)/b,  Ex1(b) = Ex'(b),  Ex2(b) = Ex1'(b),
N1(u) = (1 - L(u))/u,  N1'(u).
"""
from flint import arb, acb, ctx

def _tiny():
    return arb(2) ** (-ctx.prec)

def _absup(x):
    return arb(acb(x).abs_upper())

class AD:
    __slots__ = ("v", "d")
    n = 4
    def __init__(self, v, d=None):
        self.v = acb(v)
        self.d = [acb(0)] * AD.n if d is None else [acb(x) for x in d]
    @staticmethod
    def var(v, i):
        d = [acb(0)] * AD.n; d[i] = acb(1); return AD(v, d)
    @staticmethod
    def c(x):
        return x if isinstance(x, AD) else AD(x)
    def __add__(s, o):
        if isinstance(o, AD):
            return AD(s.v + o.v, [a + b for a, b in zip(s.d, o.d)])
        return AD(s.v + o, s.d)
    __radd__ = __add__
    def __neg__(s):
        return AD(-s.v, [-a for a in s.d])
    def __sub__(s, o):
        if isinstance(o, AD):
            return AD(s.v - o.v, [a - b for a, b in zip(s.d, o.d)])
        return AD(s.v - o, s.d)
    def __rsub__(s, o):
        return AD(o - s.v, [-a for a in s.d])
    def __mul__(s, o):
        if isinstance(o, AD):
            return AD(s.v * o.v, [a * o.v + s.v * b for a, b in zip(s.d, o.d)])
        return AD(s.v * o, [a * o for a in s.d])
    __rmul__ = __mul__
    def inv(s):
        iv = 1 / s.v; m = -iv * iv
        return AD(iv, [a * m for a in s.d])
    def __truediv__(s, o):
        if isinstance(o, AD):
            return s * o.inv()
        return AD(s.v / o, [a / o for a in s.d])
    def __rtruediv__(s, o):
        return s.inv() * o
    def __pow__(s, n):
        n = int(n)
        if n == 0:
            return AD(1)
        if n < 0:
            return (s ** (-n)).inv()
        p = s.v ** (n - 1); m = n * p
        return AD(p * s.v, [a * m for a in s.d])
    def exp(s):
        e = s.v.exp(); return AD(e, [a * e for a in s.d])
    def sqrt(s):
        r = s.v.sqrt(); m = 1 / (2 * r); return AD(r, [a * m for a in s.d])
    def conjugate(s):
        return AD(s.v.conjugate(), [a.conjugate() for a in s.d])
    @property
    def real(s):
        return AD(acb(s.v.real), [acb(a.real) for a in s.d])
    @property
    def imag(s):
        return AD(acb(s.v.imag), [acb(a.imag) for a in s.d])
    def abs_upper(s):
        return s.v.abs_upper()
    def abs_lower(s):
        return s.v.abs_lower()
    def union(s, o):
        o = o.v if isinstance(o, AD) else acb(o)
        return s.v.union(o)

def val(x):
    return x.v if isinstance(x, AD) else acb(x)

def exp(x):
    return x.exp() if isinstance(x, AD) else acb(x).exp()

def sqrt(x):
    return x.sqrt() if isinstance(x, AD) else acb(x).sqrt()

def conj(x):
    return x.conjugate() if isinstance(x, AD) else acb(x).conjugate()

def re(x):
    return x.real if isinstance(x, AD) else acb(acb(x).real)

def im(x):
    return x.imag if isinstance(x, AD) else acb(acb(x).imag)

# ------------------------------------------------------------------------------------------------ series functions
def _series_eval(coef, x, tailbound):
    """sum_k coef(k) x^k with explicit tail bound tailbound(n, r) for the terms k >= n; returns (acb)"""
    x = acb(x); r = _absup(x)
    n = 2
    while not (tailbound(n, r) < _tiny()) and n < 5000:
        n += 1
    s = acb(0); p = acb(1)
    for k in range(n):
        s += coef(k) * p
        p *= x
    t = tailbound(n, r)
    return s + acb(arb(0, t), arb(0, t))

def _geom_tail(cmax):
    """tail bound for |coef(k)| <= cmax(k) nonincreasing, |x| = r < 1:  sum_{k>=n} cmax(n) r^k = cmax(n) r^n/(1-r)"""
    def tb(n, r):
        if not (r < 1):
            return arb("inf")
        return cmax(n) * r ** n / (1 - r)
    return tb

# L(u) = sum (-1)^k u^k/(k+1)
def _L_val(u):
    return _series_eval(lambda k: arb((-1) ** k) / (k + 1), u, _geom_tail(lambda n: arb(1) / (n + 1)))
# M1(u) = L'(u) = sum_{k>=1} (-1)^k k u^(k-1)/(k+1) = sum_{j>=0} (-1)^(j+1) (j+1) u^j/(j+2)
def _M1_val(u):
    return _series_eval(lambda j: arb((-1) ** (j + 1)) * (j + 1) / (j + 2), u, _geom_tail(lambda n: arb(1)))
# N1(u) = (1 - L)/u = sum_{j>=0} (-1)^j u^j/(j+2)
def _N1_val(u):
    return _series_eval(lambda j: arb((-1) ** j) / (j + 2), u, _geom_tail(lambda n: arb(1) / (n + 2)))
# N1'(u) = sum_{j>=1} (-1)^j j u^(j-1)/(j+2) = sum_{i>=0} (-1)^(i+1) (i+1) u^i/(i+3)
def _N1p_val(u):
    return _series_eval(lambda i: arb((-1) ** (i + 1)) * (i + 1) / (i + 3), u, _geom_tail(lambda n: arb(1)))

def _fac(n):
    return arb.fac_ui(n)

# Ex(b) = sum b^k/(k+1)!,  Ex1 = Ex' = sum (k+1) b^k/(k+2)!,  Ex2 = Ex1' = sum (k+1)(k+2) b^k/(k+3)!
def _exp_tail(shift, poly):
    """terms poly(k) b^k/(k+shift)! with poly(k) <= (k+shift)^deg: tail <= sum_{k>=n} r^k (k+shift)^deg/(k+shift)!
       <= r^n e^r * poly(n) / (n+shift)! * 2  (valid once n + shift > 2 r + deg + 2)"""
    def tb(n, r):
        if not (arb(n + shift) > 2 * r + 6):
            return arb("inf")
        return 2 * r ** n * r.exp() * poly(n) / _fac(n + shift)
    return tb
def _Ex_val(b):
    return _series_eval(lambda k: 1 / _fac(k + 1), b, _exp_tail(1, lambda n: arb(1)))
def _Ex1_val(b):
    return _series_eval(lambda k: arb(k + 1) / _fac(k + 2), b, _exp_tail(2, lambda n: arb(n + 1)))
def _Ex2_val(b):
    return _series_eval(lambda k: arb((k + 1) * (k + 2)) / _fac(k + 3), b, _exp_tail(3, lambda n: arb((n + 1) * (n + 2))))

def _wrap(fval, fder):
    def f(x):
        if isinstance(x, AD):
            v = fval(x.v); dv = fder(x.v)
            return AD(v, [a * dv for a in x.d])
        return fval(acb(x))
    return f

Lfun = _wrap(_L_val, _M1_val)
N1fun = _wrap(_N1_val, _N1p_val)
Exfun = _wrap(_Ex_val, _Ex1_val)
Ex1fun = _wrap(_Ex1_val, _Ex2_val)
