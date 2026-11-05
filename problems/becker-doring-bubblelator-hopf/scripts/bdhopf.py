"""Shared library: rigorous (Arb ball arithmetic) evaluation of
    G_{beta,nu}(lam) = int_0^oo exp(-z^(beta+1)/(beta+1) - lam z) z^nu dz ,   Re lam >= 0, Im lam >= 0,
and of the Hopf quantities of the limit model of Niethammer-Pego-Schlichting-Velazquez.

Environment: python-flint 0.9 (Arb), Python 3.  Deterministic.

Method.  For lam in the closed first quadrant the integrand is holomorphic in the sector
-pi < arg z < pi and, on every ray z = r e^{-i theta} with 0 <= theta <= phi < pi/(2(beta+1)),
    |exp(-z^(beta+1)/(beta+1) - lam z) z^nu| <= r^nu exp(-r^(beta+1) cos((beta+1)phi)/(beta+1)),
because Re(lam z) = r (Re lam cos theta + Im lam sin theta) >= 0.  By Cauchy's theorem the
integral may therefore be taken along the ray arg z = -phi (the arcs vanish).  On the ray we use
    [0, 2^-K]           : |integral| <= 2^{-K(nu+1)}/(nu+1)   (bound |integrand| <= r^nu),
    [2^-k-1, 2^-k], ... : Arb's rigorous Gauss-Legendre integrator acb.integral,
    [R, oo)             : |integral| <= c^{-s} Gamma(s, c R^{beta+1})/(beta+1),
                          c = cos((beta+1)phi)/(beta+1), s = (nu+1)/(beta+1).
All parameters (beta, nu, lam) may be balls; the result then encloses the integral for every
parameter value in the balls.
"""
from flint import arb, acb, ctx

TOL_BITS = 90


def _arb(x):
    return x if isinstance(x, arb) else arb(x)


def G(beta, nu, lam, tol_bits=None):
    """Ball enclosure of G_{beta,nu}(lam); beta, nu arb (beta >= 0, nu >= 0), lam acb in the
    closed first quadrant."""
    if tol_bits is None:
        tol_bits = TOL_BITS
    beta = _arb(beta)
    nu = _arb(nu)
    lam = acb(lam)
    assert float(lam.real.lower()) >= 0 and float(lam.imag.lower()) >= 0, "lam must lie in the closed first quadrant"
    b1 = beta + 1
    # exact dyadic angle with (beta+1) phi = 0.75 (approximately) < pi/2
    phi = arb(float(0.75 / float(b1.mid())))
    cphi, sphi = phi.cos(), phi.sin()
    e = acb(cphi, -sphi)                       # e^{-i phi}
    ang = b1 * phi
    eb = acb(ang.cos(), -ang.sin())            # e^{-i (beta+1) phi}
    c = ang.cos() / b1
    assert c > 0
    tol = arb(2) ** (-tol_bits)

    def f(r, a):
        return (-(r.pow(b1, analytic=a) * eb) / b1 - lam * r * e).exp() * r.pow(nu, analytic=a)

    # cut-off R: c R^{beta+1} >= tol_bits*ln2 + 30
    Rf = float(((arb(tol_bits * 0.7 + 30) / c) ** (1 / b1)).upper()) + 1.0
    R = arb(Rf)
    # near zero: dyadic segments down to 2^-K with 2^{-K(nu+1)}/(nu+1) <= tol
    nulow = float(nu.lower())
    K = int(tol_bits / (nulow + 1.0)) + 2
    tot = acb(0)
    for k in range(1, K):
        a_, b_ = arb(2) ** (-k - 1), arb(2) ** (-k)
        tot += acb.integral(f, a_, b_, rel_tol=tol, abs_tol=tol)
    # main part
    pts = [arb(0.5)]
    x = 1.0
    while x < Rf:
        pts.append(arb(x))
        x *= 2
    pts.append(R)
    for p, q in zip(pts[:-1], pts[1:]):
        tot += acb.integral(f, p, q, rel_tol=tol, abs_tol=tol)
    # error terms
    eps = arb(2) ** (-K)
    head = eps ** (nu + 1) / (nu + 1)
    s = (nu + 1) / b1
    tail = c ** (-s) * (c * R ** b1).gamma_upper(s) / b1
    err = (head + tail).upper()
    tot = tot * (e ** (nu + 1))
    tot += acb(arb(0, err), arb(0, err))
    return tot


def G0(beta, nu):
    """G_{beta,nu}(0) = (beta+1)^((nu-beta)/(beta+1)) Gamma((nu+1)/(beta+1)) (exact formula)."""
    beta = _arb(beta)
    nu = _arb(nu)
    b1 = beta + 1
    return b1 ** ((nu - beta) / b1) * ((nu + 1) / b1).gamma()


def hopf_quantities(beta, nu, t):
    """All Hopf quantities at the (ball) frequency t (rescaled), assuming Re G(it) = 0 at the
    true t0 inside the ball.  Returns a dict of balls."""
    beta = _arb(beta)
    nu = _arb(nu)
    b1 = beta + 1
    I = acb(0, 1)
    G1 = G(beta, nu, acb(0, t))
    G2 = G(beta, nu, acb(0, 2 * t))
    Gb = G(beta, nu + beta + 1, acb(0, t))
    g0 = G0(beta, nu)
    vth0 = -G1.imag / t
    mu0 = 1 / (b1 * t ** b1)
    delta0 = t ** (nu + 2) * vth0
    tn1 = t ** (nu + 1)
    gh01 = tn1 * G1
    gh02 = tn1 * G2
    gh11 = -(t ** (nu + beta + 2)) * Gb
    L2 = 2 * I * delta0 + gh02
    absL2sq = L2.real ** 2 + L2.imag ** 2
    mu2 = -delta0 ** 2 * gh02.real / (4 * absL2sq * gh11.real)
    delta2 = -mu2 * gh11.imag - delta0 / 4 + delta0 ** 2 * (2 * delta0 + gh02.imag) / (4 * absL2sq)
    eta2 = (nu + 2) * mu2 / mu0 + b1 * delta2 / delta0
    kappa2 = (nu + 1) / b1 * mu2 / mu0 + delta2 / delta0
    eta0 = (vth0 / g0) ** b1
    kappa0 = t * vth0 / g0
    return dict(G1=G1, G2=G2, Gb=Gb, G0=g0, vartheta0=vth0, mu0=mu0, delta0=delta0,
                ghat0_1=gh01, ghat0_2=gh02, ghat1_1=gh11, L2=L2, mu2=mu2, delta2=delta2,
                eta0=eta0, kappa0=kappa0, eta2=eta2, kappa2=kappa2)


def H3_line(t, deriv=0, tol_bits=None):
    """Case beta = 3, nu = 0 only.  Enclosure of H(t) = Re G_{3,0}(it) (deriv=0) or of
    H'(t) (deriv=1) for a real ball t > 0, through the representation
        H(t) = 1/2 int_R exp(-y^4/4 - i t y) dy ,   H'(t) = 1/2 int_R (-i y) exp(-y^4/4 - i t y) dy,
    shifted (Cauchy) to the line y = x - i c, c = tc^(1/3)/2 (tc = midpoint of t, exact dyadic),
    which passes through the saddle points: the integrand there has the size of H itself, so the
    ball enclosure has a small *relative* error even for t-balls of moderate width.
    Tail |x| >= X with X >= sqrt(12) c:  Re (x - i c)^4 >= x^4/2, so
        |tail| <= e^{-t c} (1 + c/X) 2 e^{-X^4/8} / X^(3-deriv)   (both tails, factor 1/2 included).
    """
    if tol_bits is None:
        tol_bits = TOL_BITS
    t = _arb(t)
    tc = float(t.mid())
    c = arb(float(tc ** (1.0 / 3.0) / 2.0))           # exact dyadic shift
    X = max(float((arb(12).sqrt() * c).upper()) + 1.0, (8 * (tol_bits * 0.7 + 10 + float(tc ** (4.0 / 3.0)))) ** 0.25 + 1.0)
    Xa = arb(X)
    tol = arb(2) ** (-tol_bits) * (-(t.mid() * c * 0.375)).exp()
    ic = acb(0, c)

    def f(x, a):
        y = x - ic
        v = (-(y ** 4) / 4 - acb(0, 1) * t * y).exp()
        if deriv == 1:
            v = v * (-acb(0, 1) * y)
        return v

    pts = [-X] + [-X + 2 * X * k / 16 for k in range(1, 16)] + [X]
    tot = acb(0)
    for p, q in zip(pts[:-1], pts[1:]):
        tot += acb.integral(f, arb(p), arb(q), rel_tol=tol, abs_tol=tol)
    tot = tot / 2
    tail = (-(t * c)).exp() * (1 + c / Xa) * 2 * (-(Xa ** 4) / 8).exp() / Xa ** (3 - deriv)
    err = tail.upper()
    val = tot.real if deriv == 0 else tot.real   # true value is real
    return val + arb(0, err)

