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
    assert e_hi * k_hi < 1
    return arb_up((-(k_lo / e_hi) + k_hi * k_hi / 2).exp() / e_hi ** k)

# ------------------------------------------------------------------------------------------------ dual numbers
class DN:
    """value + holomorphic derivative in z + derivative in kappa (all acb balls)"""
    __slots__ = ("v", "dz", "dk")
    def __init__(self, v, dz=0, dk=0):
        self.v = acb(v); self.dz = acb(dz); self.dk = acb(dk)
    @staticmethod
    def c(x):
        return x if isinstance(x, DN) else DN(x)
    def __add__(s, o):
        o = DN.c(o); return DN(s.v + o.v, s.dz + o.dz, s.dk + o.dk)
    __radd__ = __add__
    def __neg__(s):
        return DN(-s.v, -s.dz, -s.dk)
    def __sub__(s, o):
        return s + (-DN.c(o))
    def __rsub__(s, o):
        return DN.c(o) - s
    def __mul__(s, o):
        o = DN.c(o); return DN(s.v * o.v, s.dz * o.v + s.v * o.dz, s.dk * o.v + s.v * o.dk)
    __rmul__ = __mul__
    def inv(s):
        iv = 1 / s.v; return DN(iv, -s.dz * iv * iv, -s.dk * iv * iv)
    def __truediv__(s, o):
        return s * DN.c(o).inv()
    def __rtruediv__(s, o):
        return DN.c(o) * s.inv()
    def __pow__(s, n):
        n = int(n); p = s.v ** (n - 1)
        return DN(p * s.v, n * p * s.dz, n * p * s.dk)

# ------------------------------------------------------------------------------------------------ Q^eps and Lambda
def QL(kap, z, eps):
    """Dual numbers (Q^eps, Lambda) at (kappa, z) for eps a real ball (may contain 0, then eps >= 0 is assumed)."""
    eps = arb(eps)
    kd = DN(kap, 0, 1)
    zd = DN(z, 1, 0)
    e = DN(acb(eps))
    A = 1 + kd * e
    phi = 1 + e * e * zd
    g = kd + A * e * zd
    h = g * phi + e * zd
    Ev = (acb(z) * Lfun(acb(eps) ** 2 * acb(z))).exp()
    E = DN(Ev, Ev / phi.v, 0)
    k3 = kd ** 3
    P2t = (phi * g ** 3 - e * A * A * zd * zd + A * A * zd * zd * g) / k3
    P1t = (kd * (g * g * phi * phi - e * e * zd * zd) + (kd + e) * zd * g * h) / k3
    W = -A * A * (g * g * phi * phi - e * e * zd * zd) + (A * e ** 4 * zd ** 3 + A * A * phi * g * g * (zd + phi)) / E
    # exponentially small part (T/eps) W / kappa^3 :  T/eps enclosed in a disk, d/dkappa (T/eps) = -T/(eps^2 A)
    if arb_up(eps) > 0:
        t1 = T_over_eps_bound(eps, kap, 1)
        t2 = T_over_eps_bound(eps, kap, 2)
        Td = DN(disk(t1), 0, disk(t2))
        X = Td * W / k3
    else:
        X = DN(0)
    Q = E * P2t - P1t + X
    Lam = zd + e * zd * zd / (kd * phi)
    return Q, Lam

def Q0(z, kap):
    """limit function Q(z; kappa) = e^z (1 + z^2/kappa^2) - (1 + z)"""
    z = acb(z); kap = acb(kap)
    return z.exp() * (1 + z * z / (kap * kap)) - (1 + z)

# ------------------------------------------------------------------------------------------------ Hopf system
def fvec_jac(X, eps):
    kap, x, y = X
    Q, Lam = QL(kap, acb(x, y), eps)
    f = [Q.v.real, Q.v.imag, Lam.v.real]
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
