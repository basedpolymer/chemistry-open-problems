"""l1reg -- first Lyapunov coefficient of the Becker-Doring/atomization system in an eps-regular closed form.

Uses the transfer-function form (checked against result 8's closed form, see check_l1reg.py):
  with q normalized by q_1 = 1, p the adjoint (p B = lambda p, <p, q> = 1), S the quadratic map (f = Bv + v_1 S v),
  r(mu) = (mu - B)^{-1} S q  (for mu = 0: the solution in the mass hyperplane H), sigma(y) = p S y,
      l1 = (1/omega) Re[ sigma(r(0)) + sigma(conj r(0)) + 2 Re(r(0)_1) sigma(q) + sigma(r(2 i omega))
                         + r(2 i omega)_1 sigma(conj q) ].
Every vector is a finite sum of 'atoms' (exponential sequences in l = 1..M):
  Slow(c, bt) : c beta^l,            beta = 1 + eps^2 bt      (b := M log beta = bt L(eps^2 bt))
  BM(c, gt)   : c gamma^(l-M),       gamma = A (1 + eps^2 gt) (boundary layer at l = M)
  Dd(c, g1, g2): c (gamma1^(l-M) - gamma2^(l-M))/(gamma1 - gamma2)   (two nearly equal boundary layers)
  B1(c, pt)   : c psi'^l,            psi' = (1 + eps^2 pt)/A  (boundary layer at l = 1; only in p S)
Coefficients are 'scaled numbers' SN(v, p) = v eps^p with v regular at eps = 0 and p an integer; sums over l of
products of atoms are evaluated by closed forms whose regular parts are explicit SNs and whose exponentially small
parts (powers of T = A^-M) are enclosed in disks.  The final quantity
      Lhat := eps^2 l1(q_1 = 1) = l1(q_1 = 1)/M
is obtained as a regular expression, so eps may be a ball containing 0 (eps = 0 is the limit M = infinity).
"""
from flint import arb, acb, ctx
import allm_core as C

_E = None          # current eps (acb)
_EPSBALL = None    # current eps (arb)
_KAP = None

class SN:
    __slots__ = ("v", "p")
    def __init__(s, v, p):
        s.v = acb(v); s.p = int(p)
    def at(s, p):
        """value at power p <= s.p:  v eps^(s.p - p)"""
        assert p <= s.p, "irregular: requested power %d > available %d" % (p, s.p)
        return s.v * _E ** (s.p - p)
    def __add__(s, o):
        if isinstance(o, int) and o == 0:
            return s
        p = min(s.p, o.p); return SN(s.at(p) + o.at(p), p)
    __radd__ = __add__
    def __neg__(s):
        return SN(-s.v, s.p)
    def __sub__(s, o):
        return s + (-o)
    def __mul__(s, o):
        if isinstance(o, SN):
            return SN(s.v * o.v, s.p + o.p)
        return SN(s.v * acb(o), s.p)
    __rmul__ = __mul__
    def __truediv__(s, o):
        if isinstance(o, SN):
            return SN(s.v / o.v, s.p - o.p)
        return SN(s.v / acb(o), s.p)
    def conj(s):
        return SN(s.v.conjugate(), s.p)
    def re(s):
        return SN(acb(s.v.real), s.p)

def snsum(lst):
    out = 0
    for x in lst:
        out = out + x
    return out

def expsmall(fac, j, p):
    """T eps^-j * (complex number of modulus <= fac), represented as an SN at power p (T = A^-M)"""
    k = j + p
    fac = C.arb_up(fac)
    if k >= 0:
        t = C.T_over_eps_bound(_EPSBALL, _KAP, k)
    else:
        t = C.T_over_eps_bound(_EPSBALL, _KAP, 0) * C.arb_up(_EPSBALL) ** (-k)
    return SN(C.disk(t * fac), p)

# ------------------------------------------------------------------------------------------------ atoms
class Slow:
    __slots__ = ("c", "bt")
    def __init__(s, c, bt): s.c = c; s.bt = acb(bt)
    def conj(s): return Slow(s.c.conj(), s.bt.conjugate())
    def scale(s, f): return Slow(s.c * f, s.bt)
class BM:
    __slots__ = ("c", "gt")
    def __init__(s, c, gt): s.c = c; s.gt = acb(gt)
    def conj(s): return BM(s.c.conj(), s.gt.conjugate())
    def scale(s, f): return BM(s.c * f, s.gt)
class Dd:
    __slots__ = ("c", "g1", "g2")
    def __init__(s, c, g1, g2): s.c = c; s.g1 = acb(g1); s.g2 = acb(g2)
    def conj(s): return Dd(s.c.conj(), s.g1.conjugate(), s.g2.conjugate())
    def scale(s, f): return Dd(s.c * f, s.g1, s.g2)
class B1:
    __slots__ = ("c", "pt")
    def __init__(s, c, pt): s.c = c; s.pt = acb(pt)
    def conj(s): return B1(s.c.conj(), s.pt.conjugate())
    def scale(s, f): return B1(s.c * f, s.pt)

def seq_conj(seq): return [a.conj() for a in seq]
def seq_scale(seq, f): return [a.scale(f) for a in seq]

class Ctx:
    """regular data at (eps, kappa); sets the module eps"""
    def __init__(self, eps, kap):
        global _E, _EPSBALL, _KAP
        self.eps = arb(eps); self.kap = arb(kap)
        _E = acb(self.eps); _EPSBALL = self.eps; _KAP = self.kap
        self.e = acb(self.eps); self.K = self.kap * self.eps; self.A = 1 + self.K

def _beta(cx, bt): return 1 + cx.e ** 2 * bt
def _b(cx, bt): return bt * C.Lfun(cx.e ** 2 * bt)
def _gam(cx, gt): return cx.A * (1 + cx.e ** 2 * gt)
def _lg(cx, gt): return gt * C.Lfun(cx.e ** 2 * gt)          # M log(gamma/A); |gamma^-M| = T exp(-Re lg)
def _dg(cx, gt): return cx.kap + cx.A * cx.e * gt                # (gamma - 1)/eps

# ------------------------------------------------------------------------------------------------ sums
def psum(cx, X, Y):
    """sum_{k=1}^{N} X_k Y_k  (N = M - 1), X in {Slow, B1}, Y in {Slow, BM, Dd}"""
    e = cx.e; kap = cx.kap; A = cx.A
    cc = X.c * Y.c
    if isinstance(X, Slow):
        bX = _b(cx, X.bt); beta = _beta(cx, X.bt)
        if isinstance(Y, Slow):
            bt = X.bt + Y.bt + e * e * X.bt * Y.bt
            b = bX + _b(cx, Y.bt)
            val = C.Exfun(b) * C.Lfun(e * e * bt) - e * e
            return cc * SN(val, -2)
        if isinstance(Y, BM):
            d = kap + A * e * (X.bt + Y.gt + e * e * X.bt * Y.gt)
            reg = SN(bX.exp() / d, -1)
            gam = _gam(cx, Y.gt)
            fac = acb(beta).abs_upper() * acb(gam).abs_upper() * (-(_lg(cx, Y.gt).real)).exp() / acb(d).abs_lower()
            return cc * (reg + expsmall(fac, 1, 7))
        if isinstance(Y, Dd):
            d1 = kap + A * e * (X.bt + Y.g1 + e * e * X.bt * Y.g1)
            d2 = kap + A * e * (X.bt + Y.g2 + e * e * X.bt * Y.g2)
            reg = SN(-beta * bX.exp() / (d1 * d2), -2)
            gh = Y.g1.union(Y.g2)
            dh = kap + A * e * (X.bt + gh + e * e * X.bt * gh)
            gam = _gam(cx, gh)
            inner = (-(_lg(cx, gh).real)).exp() * (1 / acb(dh).abs_lower()
                     + C.arb_up(cx.eps) * acb(beta).abs_upper() * acb(gam).abs_upper() / acb(dh).abs_lower() ** 2)
            fac = acb(beta).abs_upper() * inner
            return cc * (reg + expsmall(fac, 3, 6))
    if isinstance(X, B1):
        if isinstance(Y, Slow):
            beta = _beta(cx, Y.bt); psi = (1 + e * e * X.pt) / A
            dp = kap - e * (X.pt + Y.bt + e * e * X.pt * Y.bt)
            reg = SN(psi * beta * A / dp, -1)
            fac = A * ((_lg(cx, X.pt).real + _b(cx, Y.bt).real).exp()) / acb(dp).abs_lower()
            return cc * (reg + expsmall(fac, 1, 7))
        if isinstance(Y, BM):
            c2 = X.pt + Y.gt + e * e * X.pt * Y.gt
            mlog = (c2 * C.Lfun(e * e * c2)).real
            fac = (C.arb_up(mlog).max(arb(0)) - _lg(cx, Y.gt).real).exp()
            return cc * expsmall(fac, 2, 6)
        if isinstance(Y, Dd):
            gh = Y.g1.union(Y.g2)
            c2 = X.pt + gh + e * e * X.pt * gh
            mlog = (c2 * C.Lfun(e * e * c2)).real
            gam = _gam(cx, gh)
            fac = (C.arb_up(mlog).max(arb(0)) - _lg(cx, gh).real).exp() / acb(gam).abs_lower()
            return cc * expsmall(fac, 4, 6)
    raise TypeError((type(X), type(Y)))

def pair(cx, Xs, Ys):
    return snsum([psum(cx, X, Y) for X in Xs for Y in Ys])

def mass(cx, Y):
    """sum_{l=1}^{M} l Y_l"""
    e = cx.e
    if isinstance(Y, Slow):
        beta = _beta(cx, Y.bt); b = _b(cx, Y.bt); u = e * e * Y.bt
        L = C.Lfun(u)
        return Y.c * SN(beta * (L * L * Ex1fun(b) + e * e * b.exp() * N1fun(u)), -4)
    if isinstance(Y, BM):
        gam = _gam(cx, Y.gt); dg = _dg(cx, Y.gt)
        reg = SN(gam / dg - e * gam / dg ** 2, -3)
        fac = acb(gam).abs_upper() * (-(_lg(cx, Y.gt).real)).exp() / acb(dg).abs_lower() ** 2
        return Y.c * (reg + expsmall(fac, 2, 5))
    if isinstance(Y, Dd):
        d1 = _dg(cx, Y.g1); d2 = _dg(cx, Y.g2)
        reg = SN(-(1 - e * e) / (d1 * d2) + e * (d1 + d2) / (d1 * d1 * d2 * d2), -4)
        gh = Y.g1.union(Y.g2); dh = _dg(cx, gh); gam = _gam(cx, gh)
        fac = (-(_lg(cx, gh).real)).exp() * (1 / acb(dh).abs_lower() ** 2
                                              + 2 * C.arb_up(cx.eps) * acb(gam).abs_upper() / acb(dh).abs_lower() ** 3)
        return Y.c * (reg + expsmall(fac, 4, 4))
    raise TypeError(type(Y))

def val(cx, Y, where):
    """Y_l at l = 1, N, M"""
    e = cx.e
    if isinstance(Y, Slow):
        beta = _beta(cx, Y.bt); eb = _b(cx, Y.bt).exp()
        v = {"1": beta, "N": eb / beta, "M": eb}[where]
        return Y.c * SN(v, 0)
    if isinstance(Y, BM):
        gam = _gam(cx, Y.gt)
        if where == "M": return Y.c * SN(1, 0)
        if where == "N": return Y.c * SN(1 / gam, 0)
        fac = acb(gam).abs_upper() * (-(_lg(cx, Y.gt).real)).exp()
        return Y.c * expsmall(fac, 0, 8)
    if isinstance(Y, Dd):
        if where == "M": return Y.c * SN(0, 0)
        if where == "N": return Y.c * SN(-1 / (_gam(cx, Y.g1) * _gam(cx, Y.g2)), 0)
        gh = Y.g1.union(Y.g2)
        fac = (-(_lg(cx, gh).real)).exp()
        return Y.c * expsmall(fac, 2, 8)
    raise TypeError(type(Y))

def seqval(cx, seq, where): return snsum([val(cx, Y, where) for Y in seq])
def seqmass(cx, seq): return snsum([mass(cx, Y) for Y in seq])

# ------------------------------------------------------------------------------------------------ entire helpers
def Ex1fun(b):
    """(1 - e^b + b e^b)/b^2 = sum_{n>=0} (n+1) b^n/(n+2)!  (entire); |remainder| <= r^n e^r (n+1)/(n+2)! style bound"""
    b = acb(b); r = C.acb_absup(b)
    if not b.contains(0) and r > arb("1e-3"):
        return (1 - b.exp() + b * b.exp()) / (b * b)
    n = 1; tiny = C.tiny(); f = arb(6); rn = r; er = r.exp()       # f = (n+2)!
    while not (rn * er * (n + 2) / f < tiny) and n < 4000:
        n += 1; f *= (n + 2); rn *= r
    s = acb(0); p = acb(1); ff = arb(2)
    for k in range(n):
        s += (k + 1) * p / ff
        ff *= (k + 3); p *= b
    # tail sum_{k>=n} (k+1) r^k/(k+2)! <= r^n e^r (n+2)/(n+2)!  (since (k+1)/(k+2)! <= (n+2)/(n+2)! * 1/(k-n)!)
    bound = rn * er * (n + 2) / f
    return s + acb(arb(0, bound), arb(0, bound))

def N1fun(u):
    """(1 - L(u))/u = sum_{k>=1} (-1)^(k+1) u^(k-1)/(k+1)"""
    u = acb(u); r = C.acb_absup(u)
    if r < arb("0.5"):
        n = 1; tiny = C.tiny(); rn = arb(1)
        while not (rn / (1 - r) < tiny) and n < 4000:
            n += 1; rn *= r
        s = acb(0); p = acb(1)
        for k in range(1, n + 1):
            s += (-1) ** (k + 1) * p / (k + 1)
            p *= u
        bound = rn / (1 - r)
        return s + acb(arb(0, bound), arb(0, bound))
    return (1 - C.Lfun(u)) / u

# ------------------------------------------------------------------------------------------------ the blocks
def blocks(eps, kap, z):
    """all building blocks at the Hopf point (kap, z) (balls) and eps (ball, eps >= 0)."""
    cx = Ctx(eps, kap)
    e = cx.e; A = cx.A; kap = cx.kap; z = acb(z)
    phi = 1 + e * e * z
    g = kap + A * e * z
    Lh = g * z / phi                       # Lambda-hat: lambda = eps^3 Lh
    E = (z * C.Lfun(e * e * z)).exp()      # phi^M
    Tq = C.disk(C.T_over_eps_bound(cx.eps, kap, 1))   # T/eps
    # eigenvector q (q_1 = 1)
    C1 = (g - A * A * phi * Tq / E) / (g * E / phi + e * e * z * A * Tq / E)
    C2h = (A * E + e * e * z) / (g * E + e * e * z * A * phi * Tq / E)
    q = [Slow(SN(C1 * E, 0), -z / phi), BM(SN(C2h, -1), z)]
    # adjoint (unnormalized)
    th = Lh / (kap + e * e * Lh)
    cA = (kap + e) * th + kap
    psi = 1 / (A * phi)
    a1 = cA * psi + e * e * Tq * (A * z - e * Lh) / E
    a2 = -(cA * phi - E * (g / phi + e * e * Lh))
    S_ = a1 * phi + a2 * psi
    pc = -(1 - e * e * th) * S_
    p = [Slow(SN(a1, 0), z), B1(SN(a2, 0), -z / phi), Slow(SN(pc, 0), 0)]
    # normalization <p, q> = sum_{k=1}^{N} p_k q_k + p_M q_M
    pM = SN(a1 * E + pc, 0) + SN(a2 * Tq / E, 1)
    qM = SN(C1, 0) + SN(C2h, -1)
    Np = pair(cx, p, q) + pM * qM
    # p S (k = 1..N), normalized
    pS = [Slow(SN(a1 * z, 2) / Np, z), B1(SN(-a2 * g / (A * phi), 1) / Np, -z / phi), Slow(SN(-th * S_, 2) / Np, 0)]
    out = dict(cx=cx, Lh=Lh, q=q, Np=Np, pS=pS, E=E, C1=C1, C2h=C2h, a1=a1, a2=a2)
    out["sig_q"] = pair(cx, pS, q)
    out["sig_qb"] = pair(cx, pS, seq_conj(q))
    # resolvents
    def resolvent(m, zeta):
        """r = (mu - B)^{-1} S q, mu = eps^3 m, s1 = 1 + eps^2 zeta (slow root), s2 = A/s1"""
        gt2 = -zeta / (1 + e * e * zeta)
        Xi = (kap - e * zeta) / (1 + e * e * zeta) + A * e * z - e * e * m
        Rs = Slow(SN(z * C1 * E / (m - Lh), -1), -z / phi)
        Dp = Dd(SN(C2h * g / Xi, -1), z, gt2)
        base = [Rs, Dp]
        s1a = Slow(SN(1, 0), zeta)
        s2a = BM(SN(1, 0), gt2)
        mu = SN(m, 3)
        def rowM(seq):
            return seqval(cx, seq, "1") * SN(A, 0) + seqval(cx, seq, "N") * SN(A, 0) - seqval(cx, seq, "M") * (SN(A, 0) + mu)
        qN = seqval(cx, q, "N")
        a11 = rowM([s1a]).at(0)
        # row M on s2^(l-M):  A s2^(1-M) + A/s2 - (A + mu) = (s1 - A - mu) + A s2^(1-M) = -eps (kap - eps zeta + eps^2 m) + exp-small
        g2v = _gam(cx, gt2)
        fac12 = acb(A).abs_upper() * acb(g2v).abs_upper() * (-(_lg(cx, gt2).real)).exp()
        a12 = (SN(-(kap - e * zeta + e * e * m), 1) + expsmall(fac12, 0, 9)).at(1)
        rhs1 = (-qN - rowM(base)).at(-1)
        a21 = mass(cx, s1a).at(-4); a22 = mass(cx, s2a).at(-3)
        rhs2 = (-seqmass(cx, base)).at(-5)
        det = a11 * a22 - a12 * a21
        h1 = (rhs1 * a22 - a12 * rhs2) / det
        H2 = (a11 * rhs2 - a21 * rhs1) / det
        r = base + [Slow(SN(h1, -1), zeta), BM(SN(H2, -2), gt2)]
        return r, det
    r0, det0 = resolvent(acb(0), acb(0))
    m2 = 2 * Lh
    cm = kap + e * e * m2
    zeta2 = -2 * m2 / (cm + (cm * cm + 4 * e * m2).sqrt())
    r2, det2 = resolvent(m2, zeta2)
    out.update(r0=r0, r2=r2, det0=det0, det2=det2, zeta2=zeta2)
    out["r0_1"] = seqval(cx, r0, "1")
    out["r2_1"] = seqval(cx, r2, "1")
    out["T1"] = pair(cx, pS, r0)
    out["T2"] = pair(cx, pS, seq_conj(r0))
    out["T3"] = out["r0_1"].re() * 2 * out["sig_q"]
    out["T4"] = pair(cx, pS, r2)
    out["T5"] = out["r2_1"] * out["sig_qb"]
    br = out["T1"] + out["T2"] + out["T3"] + out["T4"] + out["T5"]
    out["bracket"] = br
    if br.p >= 1:
        bhat = br.at(1)
        out["Lhat"] = bhat.real / Lh.imag      # = eps^2 l1(q1=1) = l1(q1=1)/M
    else:
        out["Lhat"] = None
    return out
