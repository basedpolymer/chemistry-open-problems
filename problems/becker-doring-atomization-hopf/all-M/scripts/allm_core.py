"""allm_core -- eps-regular formulation of the Becker-Doring/atomization Hopf problem (all-M extension of result 8).

Model: Pego & Velazquez, Nonlinearity 33 (2020) 1812-1846, eqs (1)-(4).  Matrix B (27), characteristic function
F (71)/(76), scaled variables kappa = K sqrt(M), z = M (phi - 1), eps = M^(-1/2) (their Sec. 6-9).

Everything here is written as an explicit function of (eps, kappa, z) that is REGULAR at eps = 0, so that eps can be
a ball (interval) containing 0.  This file is self-contained (it does not import result 8's bdlib).

Notation (all exact identities, M = eps^-2, K = kappa eps, A = 1 + K, phi = 1 + eps^2 z):
    A phi - 1   = eps g,        g = kappa + A eps z
    A phi^2 - 1 = eps h,        h = g phi + eps z
    lambda      = (A - 1/phi)(phi - 1) = eps^3 g z / phi,     Lambda := M lambda / K = z + eps z^2/(kappa phi)
    E           = phi^M = exp(z L(eps^2 z)),  L(u) = log(1+u)/u      (dE/dz = E/phi)
    T           = A^-M  = exp(-(kappa/eps) L(kappa eps))   (exponentially small)
    Q^eps(z)    = K^-3 F(1 + z/M) = E P2t - P1t + (T/eps) W / kappa^3
       P2t = (phi g^3 - eps A^2 z^2 + A^2 z^2 g)/kappa^3
       P1t = (kappa (g^2 phi^2 - eps^2 z^2) + (kappa + eps) z g h)/kappa^3
       W   = -A^2 (g^2 phi^2 - eps^2 z^2) + (A eps^4 z^3 + A^2 phi g^2 (z + phi))/E
    At eps = 0:  Q^0(z) = e^z (1 + z^2/kappa^2) - (1 + z)   (paper (87)).
"""
from flint import arb, acb, arb_mat, ctx

# ------------------------------------------------------------------------------------------------ helpers
def arb_up(x):
    """upper bound (arb, exact) of a real ball"""
    return arb(arb(x).upper())

def acb_absup(x):
    return arb(acb(x).abs_upper())

def tiny():
    return arb(2) ** (-ctx.prec)

def Lfun(u):
    """L(u) = log(1+u)/u, analytic for |u| < 1; valid for balls containing 0."""
    u = acb(u)
    r = acb_absup(u)
    if r < arb("0.5"):
        n = 1
        eps = tiny()
        rn = r
        while not (rn / ((n + 1) * (1 - r)) < eps) and n < 4000:
            n += 1
            rn = rn * r
        s = acb(0)
        p = acb(1)
        for k in range(n):
            s += p / (k + 1)
            p = -p * u
        bound = rn / ((n + 1) * (1 - r))
        return s + acb(arb(0, bound), arb(0, bound))
    return acb.log1p(u) / u

def Exfun(b):
    """Ex(b) = (e^b - 1)/b (entire); valid for balls containing 0."""
    b = acb(b)
    r = acb_absup(b)
    if not b.contains(0) and r > arb("1e-3"):
        return acb.expm1(b) / b
    n = 1
    eps = tiny()
    fact = arb(2)        # (n+1)!
    rn = r
    er = r.exp()
    while not (rn * er / fact < eps) and n < 4000:
        n += 1
        fact *= (n + 1)
        rn *= r
    s = acb(0)
    p = acb(1)
    f = arb(1)
    for k in range(n):
        f *= (k + 1)          # (k+1)!
        s += p / f
        p *= b
    bound = rn * er / fact
    return s + acb(arb(0, bound), arb(0, bound))

def eps_interval(lo, hi, k=None):
    """ball containing the real interval [lo, hi] (0 <= lo < hi): endpoints rounded outwards to multiples of 2^-k,
    so that midpoint and radius are exact at the working precision (>= k + 10 bits) and the ball is exactly
    [lo', hi'] with lo' >= 0."""
    if k is None:
        k = ctx.prec - 20
    assert ctx.prec >= k + 10
    s = arb(2) ** k
    a = (arb(arb(lo).lower()) * s).floor().unique_fmpz()
    b = (arb(arb(hi).upper()) * s).ceil().unique_fmpz()
    if a is None or b is None:
        raise ValueError("endpoint too imprecise")
    if a < 0:
        a = 0
    m = arb(a + b) / (2 * s); r = arb(b - a) / (2 * s)
    assert m.rad() == 0 and r.rad() == 0
    # note: arb radii have a 30-bit mantissa, so the ball may extend below lo' by a relative 2^-30;
    # all eps-formulas used with such balls are analytic in eps there, and the only non-analytic part
    # (T = A^-M) is bounded using the upper endpoint only, so the enclosures remain valid on [lo', hi'].
    return arb(m, r)

def disk(r):
    """complex ball centred at 0 with radius r (arb upper bound)"""
    r = arb_up(r)
    return acb(arb(0, r), arb(0, r))

def T_over_eps_bound(eps, kap, k):
    """rigorous upper bound of sup T(eps', kap')/eps'^k over eps' in the ball eps (eps' > 0), kap' in kap,
    T = A^-M = (1 + kap eps)^(-1/eps^2).  Uses log(1+x) >= x - x^2/2 (x >= 0):  T <= exp(-kap/eps + kap^2/2),
    and eps -> exp(-c/eps)/eps^k is increasing on (0, c/k)."""
    eps = arb(eps); kap = arb(kap)
    e_hi = arb_up(eps)
    k_lo = arb(kap.lower()); k_hi = arb(kap.upper())
    assert e_hi < k_lo / max(k, 1), "eps too large for the monotonicity argument"
    return arb_up((-(k_lo / e_hi) + k_hi * k_hi / 2).exp() / e_hi ** k)

# ------------------------------------------------------------------------------------------------ dual numbers
class DN:
    """value + holomorphic derivative in z + derivatives in kappa and eps (all acb balls)"""
    __slots__ = ("v", "dz", "dk", "de")
    def __init__(self, v, dz=0, dk=0, de=0):
        self.v = acb(v); self.dz = acb(dz); self.dk = acb(dk); self.de = acb(de)
    @staticmethod
    def c(x):
        return x if isinstance(x, DN) else DN(x)
    def __add__(s, o):
        o = DN.c(o); return DN(s.v + o.v, s.dz + o.dz, s.dk + o.dk, s.de + o.de)
    __radd__ = __add__
    def __neg__(s):
        return DN(-s.v, -s.dz, -s.dk, -s.de)
    def __sub__(s, o):
        return s + (-DN.c(o))
    def __rsub__(s, o):
        return DN.c(o) - s
    def __mul__(s, o):
        o = DN.c(o)
        return DN(s.v * o.v, s.dz * o.v + s.v * o.dz, s.dk * o.v + s.v * o.dk, s.de * o.v + s.v * o.de)
    __rmul__ = __mul__
    def inv(s):
        iv = 1 / s.v; iv2 = -iv * iv
        return DN(iv, s.dz * iv2, s.dk * iv2, s.de * iv2)
    def __truediv__(s, o):
        return s * DN.c(o).inv()
    def __rtruediv__(s, o):
        return DN.c(o) * s.inv()
    def __pow__(s, n):
        n = int(n); p = s.v ** (n - 1); m = n * p
        return DN(p * s.v, m * s.dz, m * s.dk, m * s.de)

def M1fun(u):
    """M1(u) = (1/(1+u) - L(u))/u = sum_{k>=1} (-1)^k k/(k+1) u^(k-1); valid near 0 (|u| < 1/2)."""
    u = acb(u)
    r = acb_absup(u)
    if r < arb("0.5"):
        n = 1; eps = tiny(); rn = arb(1)
        while not (rn / (1 - r) < eps) and n < 4000:
            n += 1; rn = rn * r
        s = acb(0); p = acb(1)
        for k in range(1, n + 1):
            s += (-1) ** k * p * k / (k + 1)
            p = p * u
        bound = rn / (1 - r)
        return s + acb(arb(0, bound), arb(0, bound))
    return (1 / (1 + u) - Lfun(u)) / u

def Efun(z, eps):
    """E = phi^M = exp(z L(eps^2 z)) as a dual number: dE/dz = E/phi, dE/deps = E * 2 eps z^2 M1(eps^2 z)."""
    z = acb(z); e = acb(arb(eps)); u = e * e * z
    Ev = (z * Lfun(u)).exp()
    return DN(Ev, Ev / (1 + u), 0, Ev * 2 * e * z * z * M1fun(u))

# ------------------------------------------------------------------------------------------------ Q^eps and Lambda
def QL(kap, z, eps, split=False):
    """Dual numbers (Q^eps, Lambda) at (kappa, z) for eps a real ball (may contain 0; eps >= 0 is assumed).
    With split=True returns (Qreg, Xdisk, Lambda): Q^eps = Qreg + X with X the exponentially small part,
    enclosed (value and derivatives) in disks."""
    eps = arb(eps)
    kd = DN(kap, 0, 1, 0)
    zd = DN(z, 1, 0, 0)
    e = DN(acb(eps), 0, 0, 1)
    A = 1 + kd * e
    phi = 1 + e * e * zd
    g = kd + A * e * zd
    h = g * phi + e * zd
    E = Efun(z, eps)
    k3 = kd ** 3
    P2t = (phi * g ** 3 - e * A * A * zd * zd + A * A * zd * zd * g) / k3
    P1t = (kd * (g * g * phi * phi - e * e * zd * zd) + (kd + e) * zd * g * h) / k3
    Qreg = E * P2t - P1t
    if arb_up(eps) > 0:
        W = -A * A * (g * g * phi * phi - e * e * zd * zd) + (A * e ** 4 * zd ** 3 + A * A * phi * g * g * (zd + phi)) / E
        # (T/eps): value <= t1; d/dkappa = -T/(eps^2 A), |.| <= t2; d/deps: |.| <= (3 kappa + eps) T/eps^3
        t1 = T_over_eps_bound(eps, kap, 1)
        t2 = T_over_eps_bound(eps, kap, 2)
        t3 = T_over_eps_bound(eps, kap, 3) * (3 * arb_up(kap) + arb_up(eps))
        Td = DN(disk(t1), 0, disk(t2), disk(t3))
        X = Td * W / k3
    else:
        X = DN(0)
    Lam = zd + e * zd * zd / (kd * phi)
    if split:
        return Qreg, X, Lam
    return Qreg + X, Lam

def Q_centered(zbox, kapbox, epsbox):
    """Enclosure of Q^eps(z; kappa) over a box (z, kappa, eps) by the mean value form around the midpoint
    (exponentially small part added as a disk)."""
    zbox = acb(zbox); kapbox = arb(kapbox); epsbox = arb(epsbox)
    zm = acb(zbox.mid()); km = arb(kapbox.mid()); em = arb(epsbox.mid())
    if arb(epsbox.lower()) < -arb(2) ** -20 * arb_up(epsbox):
        raise ValueError("eps must be >= 0")
    Qm, Xm, _ = QL(km, zm, em, split=True)
    Qb, Xb, _ = QL(kapbox, zbox, epsbox, split=True)
    val = Qm.v + Qb.dz * (zbox - zm) + Qb.dk * (kapbox - km) + Qb.de * (epsbox - em) + Xb.v
    return val

def Q0(z, kap):
    """limit function Q(z; kappa) = e^z (1 + z^2/kappa^2) - (1 + z)"""
    z = acb(z); kap = acb(kap)
    return z.exp() * (1 + z * z / (kap * kap)) - (1 + z)

# ------------------------------------------------------------------------------------------------ Hopf system
def fvec_jac(X, eps):
    """f = (Re Q, Im Q, Re Lambda) and its Jacobian in (kappa, x, y) on the box X, for all eps' in the ball eps.
    f is enclosed by the mean value form in eps around mid(eps) (less overestimation for wide eps)."""
    kap, x, y = X
    Q, Lam = QL(kap, acb(x, y), eps)
    eps = arb(eps)
    if eps.rad() > 0:
        em = arb(eps.mid())
        Qm, Lm = QL(kap, acb(x, y), em)
        Qv = Qm.v + Q.de * (eps - em); Lv = Lm.v + Lam.de * (eps - em)
        Qv = acb(Qv.real.intersection(Q.v.real), Qv.imag.intersection(Q.v.imag)) if hasattr(arb, "intersection") else Qv
    else:
        Qv, Lv = Q.v, Lam.v
    f = [Qv.real, Qv.imag, Lv.real]
    J = [[Q.dk.real, Q.dz.real, -Q.dz.imag],
         [Q.dk.imag, Q.dz.imag, Q.dz.real],
         [Lam.dk.real, Lam.dz.real, -Lam.dz.imag]]
    return f, J

def newton(X0, eps, it=100):
    """floating (midpoint) Newton iteration at eps = mid(eps)"""
    em = arb(arb(eps).mid())
    X = [arb(arb(x).mid()) for x in X0]
    for _ in range(it):
        f, J = fvec_jac(X, em)
        Jm = arb_mat([[c.mid() for c in row] for row in J])
        dX = Jm.solve(arb_mat([[c.mid()] for c in f]))
        X = [arb((X[i] - dX[i, 0]).mid()) for i in range(3)]
        if max(abs(float(dX[i, 0].mid())) for i in range(3)) < 2.0 ** (-ctx.prec + 40) * (1 + abs(float(X[0].mid()))):
            break
    return X

def krawczyk(Xm, rad, eps):
    """Krawczyk test on the box Xm +- rad, uniformly for all eps' in the ball eps.
    Returns the Krawczyk image (subset of the box) or None."""
    box = [arb(Xm[i], rad[i]) for i in range(3)]
    fm, Jm = fvec_jac(Xm, eps)
    _, JX = fvec_jac(box, eps)
    Jmid = fvec_jac(Xm, arb(arb(eps).mid()))[1]
    Y = arb_mat([[c.mid() for c in row] for row in Jmid]).inv()
    Y = arb_mat([[Y[i, j].mid() for j in range(3)] for i in range(3)])
    I3 = arb_mat([[1 if i == j else 0 for j in range(3)] for i in range(3)])
    dev = arb_mat([[box[i] - Xm[i]] for i in range(3)])
    Kx = arb_mat([[Xm[i]] for i in range(3)]) - Y * arb_mat([[c] for c in fm]) + (I3 - Y * arb_mat(JX)) * dev
    if all(box[i].contains_interior(Kx[i, 0]) for i in range(3)):
        return [Kx[i, 0] for i in range(3)], box
    return None, box

def certify_hopf(X0, eps, rels=(1e-30, 1e-24, 1e-18, 1e-12, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3)):
    """Newton at mid(eps), then Krawczyk on growing boxes, uniformly in eps.  Returns (enclosure, box) or (None,None)."""
    X = newton(X0, eps)
    w = arb_up(arb(eps).rad())
    for r in rels:
        rad = [arb(abs(float(x.mid())) * r + r) for x in X]
        enc, box = krawczyk(X, rad, eps)
        if enc is not None:
            return enc, box
    return None, None

def hopf_checks(box, eps):
    """On the whole Krawczyk box (all eps' in eps):
       simple : 0 not in dQ/dz  (simple root of F -> simple eigenvalue, paper Lemma 5.4 + f(phi) != 0 check)
       transv : lower bound of d Re(Lambda)/d kappa along the root curve; Re dlambda/dkappa = eps^3 (Re Lambda +
                kappa d Re Lambda/dkappa) = eps^3 kappa dReLambda/dkappa at Re Lambda = 0
       fphi   : |phi f(phi)| = |A E + eps^2 z| > 0  (f(phi) = A phi^(M-1) + 1 - 1/phi, f(phi) phi = A E + phi - 1)
       Lam_im : enclosure of Im Lambda (omega = kappa eps^3 Im Lambda)"""
    kap, x, y = box
    z = acb(x, y)
    Q, Lam = QL(kap, z, eps)
    simple = not Q.dz.contains(0)
    dzdk = -Q.dk / Q.dz
    dLam = Lam.dz * dzdk + Lam.dk
    tr = dLam.real
    e = acb(arb(eps))
    phi = 1 + e * e * z
    Ev = (z * Lfun(e * e * z)).exp()
    fphi = (1 + kap * arb(eps)) * Ev + e * e * z
    return simple, tr, fphi, Lam.v.imag

# ------------------------------------------------------------------------------------------------ H2: localization
# Notation: w = phi - 1, a = |w|, u = Re w.  For an eigenvalue lambda = (A - 1/phi)(phi - 1) with Re lambda >= 0
# and |phi| >= 1 one has u >= 0 (Re phi >= 1), and Re lambda >= 0  <=>  a^2 (1 - A u) <= u (K + 2 A u),
# hence u >= u_min(a) (root of 2A u^2 + (K + A a^2) u - a^2 = 0), u_min increasing in a, decreasing in K.
# Gershgorin (K < 1): Re lambda >= 0 => |lambda| <= 2A + 2 => a <= A_MAX = 2 + 2 sqrt 2.

def u_min(a, K):
    """lower bound of u = Re(phi - 1) on {|phi - 1| = a, Re lambda >= 0} (valid for the K-ball: uses K_hi)"""
    a = arb(a); Kh = arb_up(K); Ah = 1 + Kh
    b = Kh + Ah * a * a
    disc = b * b + 8 * Ah * a * a
    return arb(((-b + disc.sqrt()) / (4 * Ah)).lower())

def loc_check_I(eps, kap, R0, a1):
    """Sum of the ratio bounds rho_a..rho_f (see result.md, Lemma L1) on a in [R0/M, a1]; must be < 1.
    eps, kap: balls (eps >= 0); returns an upper bound (arb)."""
    eps = arb(eps); kap = arb(kap); R0 = arb(R0); a1 = arb(a1)
    e = arb(0).union(arb_up(eps)) if arb(eps.lower()) <= 0 else eps
    K = kap * e; A = 1 + K
    rho_a = (1 + a1) * (kap / (A * R0) + e) ** 2
    rho_b = e / kap
    rho_c = (1 + a1) ** 2 * (kap ** 2 / (A * A * R0 * R0) + kap * e / (A * R0))
    rho_d = e * e / (A * A)
    rho_e = (kap ** 2 + kap * e) / (A * A * R0) + (kap * e + e * e) * (2 + a1) / A
    Ka = K + A * a1
    Nf0 = A * A * (Ka ** 2 * (1 + a1) ** 2 + a1 ** 2) + A * a1 ** 3 + A * A * (1 + a1) * Ka ** 2 * (1 + a1)
    Nf1 = A * A * (1 + a1) * Ka ** 2 * a1
    if arb_up(eps) > 0:
        rho_f = (Nf0 * T_over_eps_bound(eps, kap, 3) + Nf1 * T_over_eps_bound(eps, kap, 5)) / (A * A * R0 * R0 * kap)
    else:
        rho_f = arb(0)
    s = rho_a + rho_b + rho_c + rho_d + rho_e + rho_f
    return arb_up(s)

def loc_check_II(eps, kap, a_lo, a_hi):
    """On a in [a_lo, a_hi] (a_lo > 0): returns a lower bound of
         |phi|^M * LB(|P2|)/M - UB(|P1|)/M - UB(exp-small)/M
    over all u in [u_min(a_lo), a_hi], all eps' in eps (eps' > 0, M = 1/eps'^2) and kappa' in kap; a positive value
    proves that F has no zero phi with |phi - 1| in [a_lo, a_hi], |phi| >= 1, Re lambda >= 0.
    Every quantity is monotone in (K, u, a, 1/M) >= 0, so lower/upper bounds are evaluated at the endpoints:
       |K + A w|^2 = K^2 + 2 K A u + A^2 a^2,   |phi|^2 = 1 + 2u + a^2,   |A phi^2 - 1| <= K + 2Aa + Aa^2,
       |P2| >= M A^2 a^2 |K+Aw| - |phi||K+Aw|^3 - A^2 a^2,
       |P1| <= K(|K+Aw|^2 |phi|^2 + a^2) + (KM+1) a |K+Aw| |A phi^2 - 1|,
       A^-M(|R1| + |R2|) <= T [A^2(|K+Aw|^2|phi|^2 + a^2) + A a^3 + A^2 |phi| |K+Aw|^2 (M a + |phi|)],
       |phi|^M >= (1 + 2 u_min(a_lo) + a_lo^2)^(M_lo/2)."""
    eps = arb(eps); kap = arb(kap)
    eh = arb_up(eps); el = arb(eps.lower())
    if el < 0:
        el = arb(0)
    kl = arb(kap.lower()); kh = arb_up(kap)
    Kl = kl * el; Kh = arb_up(kh * eh)
    Al = 1 + Kl; Ah = 1 + Kh
    al = arb(a_lo); ah = arb_up(a_hi)
    ul = u_min(al, Kh); uh = ah
    invMl = el * el; invMh = arb_up(eh * eh)
    Mlo = 1 / invMh
    kw_l = (Kl * Kl + 2 * Kl * Al * ul + Al * Al * al * al).sqrt()
    kw_h = arb_up((Kh * Kh + 2 * Kh * Ah * uh + Ah * Ah * ah * ah).sqrt())
    ph2_l = 1 + 2 * ul + al * al
    ph2_h = arb_up(1 + 2 * uh + ah * ah)
    ph_h = arb_up(ph2_h.sqrt())
    LB2 = arb((Al * Al * al * al * kw_l - (ph_h * kw_h ** 3 + Ah * Ah * ah * ah) * invMh).lower())
    UB1 = arb_up(Kh * (kw_h * kw_h * ph2_h + ah * ah) * invMh + (Kh + invMh) * ah * kw_h * (Kh + 2 * Ah * ah + Ah * ah * ah))
    if eh > 0:
        T = T_over_eps_bound(eps, kap, 0)
        EXP = arb_up(T * ((Ah * Ah * (kw_h * kw_h * ph2_h + ah * ah) + Ah * ah ** 3 + Ah * Ah * ph2_h * kw_h * kw_h) * invMh
                          + Ah * Ah * ph_h * kw_h * kw_h * ah))
    else:
        EXP = arb(0)
    if not (LB2 > 0):
        return arb(-1)
    grow = arb((ph2_l ** (Mlo / 2)).lower())
    return arb((grow * LB2 - UB1 - EXP).lower())

def loc_check_II_range(eps, kap, a1, amax=None, n0=64, depth=12):
    """verify loc_check_II > 0 on [a1, amax] by adaptive subdivision; returns (ok, number of intervals)"""
    if amax is None:
        amax = 2 + 2 * arb(2).sqrt()
    amax = arb_up(amax)
    stack = []
    a1 = arb(a1)
    for i in range(n0):
        lo = a1 + (amax - a1) * i / n0; hi = a1 + (amax - a1) * (i + 1) / n0
        stack.append((arb(lo.lower()), arb_up(hi), 0))
    cnt = 0
    while stack:
        lo, hi, d = stack.pop()
        cnt += 1
        if loc_check_II(eps, kap, lo, hi) > 0:
            continue
        if d >= depth:
            return False, cnt
        mid = arb(((lo + hi) / 2).mid())
        stack.append((lo, mid, d + 1)); stack.append((mid, hi, d + 1))
    return True, cnt

def region_check(eps, kap, R0, sigma):
    """upper bound of Re Lambda on {|z| <= R0, Re z <= -sigma}; must be < 0"""
    eps = arb(eps); kap = arb(kap); R0 = arb(R0); sigma = arb(sigma)
    e = arb(0).union(arb_up(eps))
    kl = arb(kap.lower())
    q = e * e * R0
    assert arb_up(q) < 1
    corr = e ** 3 * R0 ** 3 / (kl * (1 - q))
    c1 = -sigma + e * sigma ** 2 / kl
    c2 = -R0 + e * R0 ** 2 / kl
    return arb_up(arb_up(c1).max(arb_up(c2)) + corr)

# ------------------------------------------------------------------------------------------------ argument principle
import math

def _seg_box(z0, z1):
    """complex ball containing the segment [z0, z1]"""
    return acb(z0).union(acb(z1))

def winding_count(path_pts, kap, eps, max_depth=30, Qfun=None):
    """Number of zeros (with multiplicity) of Q^eps(.; kappa) inside the closed polygon path_pts (list of complex
    points, counterclockwise, closed: last == first), valid for all kappa in kap, eps in eps.
    Each edge is subdivided until the mean-value enclosure of Q over the edge (and the parameter boxes) excludes 0;
    then the argument increment along the edge is the principal difference of the arguments of Q at the
    endpoints (at the midpoint parameters).  Returns (count, number of segments) or (None, info) on failure."""
    if Qfun is None:
        Qfun = Q_centered
    km = arb(arb(kap).mid()); em = arb(arb(eps).mid())
    total = 0.0
    nseg = 0
    def val(zz):
        return QL(km, acb(zz), em)[0].v
    for i in range(len(path_pts) - 1):
        stack = [(acb(path_pts[i]), acb(path_pts[i + 1]), 0)]
        while stack:
            z0, z1, d = stack.pop()
            B = Qfun(_seg_box(z0, z1), kap, eps)
            if B.contains(0) or not B.is_finite():
                if d >= max_depth:
                    return None, ("zero not excluded near", complex(z0.mid()), complex(z1.mid()))
                zm = acb(((z0 + z1) / 2).mid())
                stack.append((zm, z1, d + 1)); stack.append((z0, zm, d + 1))   # z0..zm processed first
                continue
            v0 = val(z0); v1 = val(z1)
            a0 = math.atan2(float(v0.imag.mid()), float(v0.real.mid()))
            a1 = math.atan2(float(v1.imag.mid()), float(v1.real.mid()))
            da = a1 - a0
            while da > math.pi: da -= 2 * math.pi
            while da < -math.pi: da += 2 * math.pi
            # rigour of the principal difference: B excludes 0, so |true increment| < pi; the float error is
            # bounded by the relative radii of v0, v1 (checked small) plus 1e-12
            for v in (v0, v1):
                rel = float(acb_absup(acb(arb(0, v.real.rad()), arb(0, v.imag.rad())))) / max(float(v.abs_lower()), 1e-300)
                if rel > 1e-6:
                    return None, ("endpoint value not precise", complex(z0.mid()))
            total += da
            nseg += 1
    w = total / (2 * math.pi)
    n = round(w)
    if abs(w - n) > 1e-3:
        return None, ("non-integer winding", w)
    return n, nseg

def contour_D(R0, sigma, n_arc=96, n_seg=24):
    """Counterclockwise convex polygon P containing D = {|z| <= R0, Re z >= -sigma}: vertices on the circle of
    radius R1 = R0/cos(Delta/2) (Delta = angular step), so that every chord is at distance >= R0 from 0, at angles
    in [-th1, th1] with R1 cos th1 = -sigma, closed by the vertical chord Re z = -sigma.  Vertices are rounded
    outwards slightly (radius * (1 + 1e-12)) to absorb float rounding of the vertex coordinates."""
    R0 = float(R0); sigma = float(sigma)
    # fixed point: Delta depends on th1 which depends on R1
    R1 = R0
    for _ in range(50):
        th1 = math.pi - math.acos(sigma / R1)
        Delta = 2 * th1 / n_arc
        R1 = R0 / math.cos(Delta / 2) * (1 + 1e-12)
    th1 = math.pi - math.acos(sigma / R1)
    Delta = 2 * th1 / n_arc
    assert R1 * math.cos(Delta / 2) >= R0 * (1 + 1e-13)
    pts = []
    for k in range(n_arc + 1):
        th = -th1 + Delta * k
        pts.append(complex(R1 * math.cos(th), R1 * math.sin(th)))
    xs = pts[-1].real
    Y1 = pts[-1].imag
    for k in range(1, n_seg):
        pts.append(complex(xs, Y1 - 2 * Y1 * k / n_seg))
    pts.append(pts[0])
    return pts, R1
