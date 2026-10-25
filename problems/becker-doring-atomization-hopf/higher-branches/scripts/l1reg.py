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
import ad
from ad import AD

_E = None          # current eps (acb)
_EPSBALL = None    # current eps (arb)
_KAP = None

class SN:
    __slots__ = ("v", "p")
    def __init__(s, v, p):
        s.v = v if isinstance(v, AD) else acb(v); s.p = int(p)
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
        return SN(s.v * (o if isinstance(o, AD) else acb(o)), s.p)
    __rmul__ = __mul__
    def __truediv__(s, o):
        if isinstance(o, SN):
            return SN(s.v / o.v, s.p - o.p)
        return SN(s.v / (o if isinstance(o, AD) else acb(o)), s.p)
    def conj(s):
        return SN(s.v.conjugate(), s.p)
    def re(s):
        return SN(ad.re(s.v), s.p)

def snsum(lst):
    out = 0
    for x in lst:
        out = out + x
    return out

_DISK_MODE = None      # None: compute; "record": compute and store radii; "replay": reuse stored radii in order
_DISK_LOG = []

def expsmall(fac, j, p):
    """T eps^-j * (complex number of modulus <= fac), represented as an SN at power p (T = A^-M).
    In "replay" mode the radius computed in the preceding "record" pass (box evaluation) is reused, so that the
    midpoint evaluation of the mean value form uses disks valid on the whole box (result.md, Lemma 5.2)."""
    global _DISK_LOG
    if _DISK_MODE == "replay":
        jj, pp, rad = _DISK_LOG.pop(0)
        assert (jj, pp) == (j, p), "replay mismatch"
        return SN(C.disk(rad), p)
    k = j + p
    fac = C.arb_up(fac)
    if k >= 0:
        t = C.T_over_eps_bound(_EPSBALL, _KAP, k)
    else:
        t = C.T_over_eps_bound(_EPSBALL, _KAP, 0) * C.arb_up(_EPSBALL) ** (-k)
    rad = C.arb_up(t * fac)
    if _DISK_MODE == "record":
        _DISK_LOG.append((j, p, rad))
    return SN(C.disk(rad), p)

# ------------------------------------------------------------------------------------------------ atoms
def _cv(x):
    return x if isinstance(x, AD) else acb(x)

def _up(x):
    """arb upper bound of |x| (x acb or AD; uses the value only)"""
    return arb(ad.val(x).abs_upper())

def _lo(x):
    return arb(ad.val(x).abs_lower())

def _hull(a, b):
    return ad.val(a).union(ad.val(b))

def _expre_up(x):
    """upper bound of exp(Re x)"""
    return arb(arb(ad.val(x).real.upper()).exp().upper())

class Slow:
    __slots__ = ("c", "bt")
    def __init__(s, c, bt): s.c = c; s.bt = _cv(bt)
    def conj(s): return Slow(s.c.conj(), ad.conj(s.bt))
    def scale(s, f): return Slow(s.c * f, s.bt)
class BM:
    __slots__ = ("c", "gt")
    def __init__(s, c, gt): s.c = c; s.gt = _cv(gt)
    def conj(s): return BM(s.c.conj(), ad.conj(s.gt))
    def scale(s, f): return BM(s.c * f, s.gt)
class Dd:
    __slots__ = ("c", "g1", "g2")
    def __init__(s, c, g1, g2): s.c = c; s.g1 = _cv(g1); s.g2 = _cv(g2)
    def conj(s): return Dd(s.c.conj(), ad.conj(s.g1), ad.conj(s.g2))
    def scale(s, f): return Dd(s.c * f, s.g1, s.g2)
class B1:
    __slots__ = ("c", "pt")
    def __init__(s, c, pt): s.c = c; s.pt = _cv(pt)
    def conj(s): return B1(s.c.conj(), ad.conj(s.pt))
    def scale(s, f): return B1(s.c * f, s.pt)

def seq_conj(seq): return [a.conj() for a in seq]
def seq_scale(seq, f): return [a.scale(f) for a in seq]

class Ctx:
    """regular data at the evaluation point (eps, kappa) -- acb/arb balls or AD -- and the parameter BOXES
    (eps_box, kap_box) used for all bounds of exponentially small terms (T = A^-M)."""
    def __init__(self, eps, kap, eps_box=None, kap_box=None):
        global _E, _EPSBALL, _KAP
        self.eps_box = arb(eps_box if eps_box is not None else ad.val(eps).real)
        self.kap_box = arb(kap_box if kap_box is not None else ad.val(kap).real)
        self.e = eps if isinstance(eps, AD) else acb(eps)
        self.kap = kap if isinstance(kap, AD) else acb(kap)
        self.eps = self.eps_box
        _E = self.e; _EPSBALL = self.eps_box; _KAP = self.kap_box
        self.K = self.kap * self.e; self.A = 1 + self.K

def _beta(cx, bt): return 1 + cx.e ** 2 * bt
def _b(cx, bt): return bt * ad.Lfun(cx.e ** 2 * bt)
def _gam(cx, gt): return cx.A * (1 + cx.e ** 2 * gt)
def _lg(cx, gt): return gt * ad.Lfun(cx.e ** 2 * gt)          # M log(gamma/A); |gamma^-M| = T exp(-Re lg)
def _dg(cx, gt): return cx.kap + cx.A * cx.e * gt                # (gamma - 1)/eps

# ------------------------------------------------------------------------------------------------ sums
# Exponentially small parts: |gamma^-M| = T exp(-Re lg(gt)), etc.  All bounds use VALUES only (constants for AD).

def psum(cx, X, Y):
    """sum_{k=1}^{N} X_k Y_k  (N = M - 1), X in {Slow, B1}, Y in {Slow, BM, Dd}"""
    e = cx.e; kap = cx.kap; A = cx.A
    epsup = C.arb_up(cx.eps_box)
    cc = X.c * Y.c
    if isinstance(X, Slow):
        bX = _b(cx, X.bt); beta = _beta(cx, X.bt)
        if isinstance(Y, Slow):
            bt = X.bt + Y.bt + e * e * X.bt * Y.bt
            b = bX + _b(cx, Y.bt)
            v = ad.Exfun(b) * ad.Lfun(e * e * bt) - e * e
            return cc * SN(v, -2)
        if isinstance(Y, BM):
            d = kap + A * e * (X.bt + Y.gt + e * e * X.bt * Y.gt)
            reg = SN(ad.exp(bX) / d, -1)
            gam = _gam(cx, Y.gt)
            fac = _up(beta) * _up(gam) * _expre_up(-_lg(cx, Y.gt)) / _lo(d)
            return cc * (reg + expsmall(fac, 1, 7))
        if isinstance(Y, Dd):
            d1 = kap + A * e * (X.bt + Y.g1 + e * e * X.bt * Y.g1)
            d2 = kap + A * e * (X.bt + Y.g2 + e * e * X.bt * Y.g2)
            reg = SN(-beta * ad.exp(bX) / (d1 * d2), -2)
            gh = _hull(Y.g1, Y.g2)
            ev = ad.val(e); btv = ad.val(X.bt); kv = ad.val(kap); Av = ad.val(A)
            dh = kv + Av * ev * (btv + gh + ev * ev * btv * gh)
            gam = Av * (1 + ev * ev * gh)
            lgh = gh * ad.Lfun(ev * ev * gh)
            inner = _expre_up(-lgh) * (1 / _lo(dh) + epsup * _up(beta) * _up(gam) / _lo(dh) ** 2)
            fac = _up(beta) * inner
            return cc * (reg + expsmall(fac, 3, 6))
    if isinstance(X, B1):
        if isinstance(Y, Slow):
            beta = _beta(cx, Y.bt); psi = (1 + e * e * X.pt) / A
            dp = kap - e * (X.pt + Y.bt + e * e * X.pt * Y.bt)
            reg = SN(psi * beta * A / dp, -1)
            fac = _up(A) * _expre_up(_lg(cx, X.pt) + _b(cx, Y.bt)) / _lo(dp)
            return cc * (reg + expsmall(fac, 1, 7))
        if isinstance(Y, BM):
            c2 = ad.val(X.pt + Y.gt + e * e * X.pt * Y.gt); ev = ad.val(e)
            mlog = c2 * ad.Lfun(ev * ev * c2)
            fac = _expre_up(mlog).max(arb(1)) * _expre_up(-_lg(cx, Y.gt))
            return cc * expsmall(fac, 2, 6)
        if isinstance(Y, Dd):
            gh = _hull(Y.g1, Y.g2); ev = ad.val(e); ptv = ad.val(X.pt); Av = ad.val(A)
            c2 = ptv + gh + ev * ev * ptv * gh
            mlog = c2 * ad.Lfun(ev * ev * c2)
            gam = Av * (1 + ev * ev * gh)
            lgh = gh * ad.Lfun(ev * ev * gh)
            fac = _expre_up(mlog).max(arb(1)) * _expre_up(-lgh) / _lo(gam)
            return cc * expsmall(fac, 4, 6)
    raise TypeError((type(X), type(Y)))

def pair(cx, Xs, Ys):
    return snsum([psum(cx, X, Y) for X in Xs for Y in Ys])

def mass(cx, Y):
    """sum_{l=1}^{M} l Y_l"""
    e = cx.e
    epsup = C.arb_up(cx.eps_box)
    if isinstance(Y, Slow):
        beta = _beta(cx, Y.bt); b = _b(cx, Y.bt); u = e * e * Y.bt
        L = ad.Lfun(u)
        return Y.c * SN(beta * (L * L * ad.Ex1fun(b) + e * e * ad.exp(b) * ad.N1fun(u)), -4)
    if isinstance(Y, BM):
        gam = _gam(cx, Y.gt); dg = _dg(cx, Y.gt)
        reg = SN(gam / dg - e * gam / dg ** 2, -3)
        fac = _up(gam) * _expre_up(-_lg(cx, Y.gt)) / _lo(dg) ** 2
        return Y.c * (reg + expsmall(fac, 2, 5))
    if isinstance(Y, Dd):
        d1 = _dg(cx, Y.g1); d2 = _dg(cx, Y.g2)
        reg = SN(-(1 - e * e) / (d1 * d2) + e * (d1 + d2) / (d1 * d1 * d2 * d2), -4)
        gh = _hull(Y.g1, Y.g2); ev = ad.val(e); kv = ad.val(cx.kap); Av = ad.val(cx.A)
        dh = kv + Av * ev * gh; gam = Av * (1 + ev * ev * gh); lgh = gh * ad.Lfun(ev * ev * gh)
        fac = _expre_up(-lgh) * (1 / _lo(dh) ** 2 + 2 * epsup * _up(gam) / _lo(dh) ** 3)
        return Y.c * (reg + expsmall(fac, 4, 4))
    raise TypeError(type(Y))

def val(cx, Y, where):
    """Y_l at l = 1, N, M"""
    e = cx.e
    if isinstance(Y, Slow):
        beta = _beta(cx, Y.bt); eb = ad.exp(_b(cx, Y.bt))
        v = {"1": beta, "N": eb / beta, "M": eb}[where]
        return Y.c * SN(v, 0)
    if isinstance(Y, BM):
        gam = _gam(cx, Y.gt)
        if where == "M": return Y.c * SN(1, 0)
        if where == "N": return Y.c * SN(1 / gam, 0)
        fac = _up(gam) * _expre_up(-_lg(cx, Y.gt))
        return Y.c * expsmall(fac, 0, 8)
    if isinstance(Y, Dd):
        if where == "M": return Y.c * SN(0, 0)
        if where == "N": return Y.c * SN(-1 / (_gam(cx, Y.g1) * _gam(cx, Y.g2)), 0)
        gh = _hull(Y.g1, Y.g2); ev = ad.val(e)
        fac = _expre_up(-(gh * ad.Lfun(ev * ev * gh)))
        return Y.c * expsmall(fac, 2, 8)
    raise TypeError(type(Y))

def seqval(cx, seq, where): return snsum([val(cx, Y, where) for Y in seq])
def seqmass(cx, seq): return snsum([mass(cx, Y) for Y in seq])

# ------------------------------------------------------------------------------------------------ the blocks
def blocks(eps, kap, z, eps_box=None, kap_box=None):
    """All building blocks at the Hopf point (kap, z) and eps (balls or AD).  eps_box, kap_box: parameter boxes
    used for the bounds of the exponentially small terms (default: eps, kap themselves)."""
    cx = Ctx(eps, kap, eps_box, kap_box)
    e = cx.e; A = cx.A; kap = cx.kap
    z = z if isinstance(z, AD) else acb(z)
    phi = 1 + e * e * z
    g = kap + A * e * z
    Lh = g * z / phi                       # Lambda-hat: lambda = eps^3 Lh
    E = ad.exp(z * ad.Lfun(e * e * z))     # phi^M
    Tq = C.disk(C.T_over_eps_bound(cx.eps_box, cx.kap_box, 1))   # T/eps (over the boxes)
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
    def resolvent(m, zeta):
        """r = (mu - B)^{-1} S q, mu = eps^3 m, s1 = 1 + eps^2 zeta (slow root), s2 = A/s1;
        rows 2..N by construction, row M and the mass condition sum l r_l = 0 (equivalent to row 1) solved."""
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
        # row M on s2^(l-M):  A s2^(1-M) + A/s2 - (A + mu) = (s1 - A - mu) + A s2^(1-M)
        #                                               = -eps (kap - eps zeta + eps^2 m) + exp-small
        g2v = _gam(cx, gt2)
        fac12 = _up(A) * _up(g2v) * _expre_up(-_lg(cx, gt2))
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
    zeta2 = -2 * m2 / (cm + ad.sqrt(cm * cm + 4 * e * m2))
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
        out["Lhat"] = ad.re(bhat) / ad.im(Lh)      # = eps^2 l1(q1=1) = l1(q1=1)/M  (real, as acb/AD)
    else:
        out["Lhat"] = None
    return out

def Lhat_box(eps_box, kap_box, x_box, y_box):
    """Rigorous enclosure of Lhat(eps, kappa, x + i y) for all (eps, kappa, x, y) in the box, by the mean value
    form: value at the midpoint (exponentially small terms as disks over the whole box) + gradient over the box
    (forward AD in the 4 real directions, the exponentially small disks being held fixed) times the half-widths.
    Returns an arb."""
    global _DISK_MODE, _DISK_LOG
    eps_box = arb(eps_box); kap_box = arb(kap_box); x_box = arb(x_box); y_box = arb(y_box)
    em, km, xm, ym = [arb(b.mid()) for b in (eps_box, kap_box, x_box, y_box)]
    AD.n = 4
    kv = AD.var(kap_box, 0); xv = AD.var(x_box, 1); yv = AD.var(y_box, 2); ev = AD.var(eps_box, 3)
    zv = xv + yv * acb(0, 1)
    # box evaluation (values and gradient over the whole box); record the exponentially small disks
    _DISK_MODE = "record"; _DISK_LOG = []
    try:
        G = blocks(ev, kv, zv, eps_box, kap_box)["Lhat"]
    finally:
        _DISK_MODE = None
    if G is None:
        return None
    # midpoint evaluation with the SAME disks (valid on the whole box)
    _DISK_MODE = "replay"
    try:
        F = blocks(em, km, acb(xm, ym), eps_box, kap_box)["Lhat"]
        assert not _DISK_LOG, "replay incomplete"
    finally:
        _DISK_MODE = None
    if F is None:
        return None
    val = acb(F)
    for i, (b, m) in enumerate(zip((kap_box, x_box, y_box, eps_box), (km, xm, ym, em))):
        val += G.d[i] * (b - m)
    return val.real
