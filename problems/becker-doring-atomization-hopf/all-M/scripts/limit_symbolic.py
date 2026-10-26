#!/usr/bin/env python3
"""limit_symbolic.py -- the limit M -> infinity of the first Lyapunov coefficient, for EVERY branch j, in closed form.

Environment: Python 3.12 or later, sympy, python-flint 0.9, mpmath.  One thread.  Runtime: about 1-3 min.

The all-M regular formula (all-M result.md, Section 5; scripts/l1reg.py) gives Lhat = l1(q_1 = 1)/M as an expression
that is continuous in eps = M^-1/2 up to eps = 0.  This script evaluates that SAME algorithm (l1reg.blocks) at eps = 0
symbolically: every scaled number v eps^p keeps only its value at eps = 0 (a term of higher power vanishes; the
exponentially small terms vanish), with A = 1, phi = 1, g = kappa, E = phi^M = exp(z), z = i t at the Hopf point,
and exp(n i t) = w^n, w = e^{it}.
Claims (PASS/FAIL):
  S1  (numerical consistency diagnostic, not an algebraic proof) at j = 1, ..., 8 the symbolic
      expression agrees within 1e-25 with the ball evaluation of l1reg.blocks at eps = 0.
  S2  (identity) with the Hopf relations R = 1 + kappa^2, t^2 = R^2 - 1, w = -(1 + i t)/R (PV (89)-(91)), written in
      the rational parametrization R = (4 tau^4 + 1)/(4 tau^2), t = (4 tau^4 - 1)/(4 tau^2),
      kappa = (2 tau^2 - 1)/(2 tau)  (tau > 2^-1/2), the expression equals result 8's closed form
      L = -R (R-1)(7R^2 - 34R + 23) / [t (R+1)(5R-3)((R-1)^3 + R + 1)]   identically in tau.
  S3  (no singularity) the quantities the algorithm divides by -- the normalization <p, q>, the 2x2 determinants of
      the two resolvents, Im Lambda-hat = kappa t, and m - Lambda-hat -- are nonzero at eps = 0 for every
      tau >= 2 (all Hopf points have t >= t_1 > 4.49, i.e. tau > 2.12): each is a rational function of tau whose
      numerator has no real root >= 2 (exact real-root isolation).
  S4  L < 0 for tau >= 2 (7R^2 - 34R + 23 > 0 for R > 4.05; all other factors positive), so for every j:
      lim_{M -> infinity} l1(q_1 = 1)/M = L(t_j) < 0.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp

fails = 0
def report(name, ok, extra=""):
    global fails
    fails += (not ok)
    print("%-4s %s %s" % (name, "PASS" if ok else "FAIL", extra), flush=True)

t, k = sp.symbols("t kappa", positive=True)
w = sp.symbols("w")
I = sp.I
BIG = 10 ** 6

def conj(e):
    return sp.sympify(e).xreplace({I: -I, w: 1 / w})

def expo(n):            # exp(n i t)
    return w ** n

def Ex(n):              # (e^b - 1)/b, b = n i t
    return sp.Integer(1) if n == 0 else (w ** n - 1) / (n * I * t)

def Ex1(n):             # (1 - e^b + b e^b)/b^2
    if n == 0:
        return sp.Rational(1, 2)
    b = n * I * t
    return (1 - w ** n + b * w ** n) / b ** 2

class SN:
    """value of v eps^p at eps = 0, keeping the power p (v is the eps = 0 value of the regular coefficient)"""
    __slots__ = ("v", "p")
    def __init__(s, v, p):
        s.v = sp.sympify(v); s.p = int(p)
    def at(s, p):
        assert p <= s.p, "irregular"
        return s.v if s.p == p else sp.Integer(0)
    def __add__(s, o):
        if isinstance(o, int) and o == 0:
            return s
        if s.p < o.p: return SN(s.v, s.p)
        if o.p < s.p: return SN(o.v, o.p)
        return SN(s.v + o.v, s.p)
    __radd__ = __add__
    def __neg__(s): return SN(-s.v, s.p)
    def __sub__(s, o): return s + (-o)
    def __mul__(s, o):
        if isinstance(o, SN): return SN(s.v * o.v, s.p + o.p)
        return SN(s.v * o, s.p)
    __rmul__ = __mul__
    def __truediv__(s, o):
        if isinstance(o, SN): return SN(s.v / o.v, s.p - o.p)
        return SN(s.v / o, s.p)
    def conj(s): return SN(conj(s.v), s.p)
    def re(s): return SN((s.v + conj(s.v)) / 2, s.p)

ZERO = SN(0, BIG)

def snsum(lst):
    out = 0
    for x in lst:
        out = out + x
    return out

# atoms: (kind, c, n...) with bt = n i t
class Atom:
    def __init__(s, kind, c, n1, n2=None):
        s.kind = kind; s.c = c; s.n1 = n1; s.n2 = n2
    def conj(s):
        return Atom(s.kind, s.c.conj(), -s.n1, None if s.n2 is None else -s.n2)

def seq_conj(seq): return [a.conj() for a in seq]

def psum(X, Y):
    cc = X.c * Y.c
    if X.kind == "Slow":
        if Y.kind == "Slow": return cc * SN(Ex(X.n1 + Y.n1), -2)
        if Y.kind == "BM":   return cc * SN(expo(X.n1) / k, -1)
        if Y.kind == "Dd":   return cc * SN(-expo(X.n1) / k ** 2, -2)
    if X.kind == "B1":
        if Y.kind == "Slow": return cc * SN(1 / k, -1)
        return ZERO                                   # exponentially small
    raise TypeError

def pair(Xs, Ys): return snsum([psum(X, Y) for X in Xs for Y in Ys])

def mass(Y):
    if Y.kind == "Slow": return Y.c * SN(Ex1(Y.n1), -4)
    if Y.kind == "BM":   return Y.c * SN(1 / k, -3)
    if Y.kind == "Dd":   return Y.c * SN(-1 / k ** 2, -4)

def val(Y, where):
    if Y.kind == "Slow":
        return Y.c * SN(1 if where == "1" else expo(Y.n1), 0)
    if Y.kind == "BM":
        return Y.c * SN(1, 0) if where in ("M", "N") else ZERO
    if Y.kind == "Dd":
        if where == "M": return Y.c * SN(0, 0)
        if where == "N": return Y.c * SN(-1, 0)
        return ZERO

def seqval(seq, where): return snsum([val(Y, where) for Y in seq])
def seqmass(seq): return snsum([mass(Y) for Y in seq])

def blocks0():
    """l1reg.blocks at eps = 0 (symbolic); returns dict with Lhat and the denominators"""
    z = I * t
    E = w
    g = k
    Lh = k * z
    C1 = 1 / E
    C2h = 1 / k
    q = [Atom("Slow", SN(C1 * E, 0), -1), Atom("BM", SN(C2h, -1), 1)]
    th = z
    cA = k * th + k
    a1 = cA
    a2 = -(cA - E * g)
    S_ = a1 + a2
    pc = -S_
    p = [Atom("Slow", SN(a1, 0), 1), Atom("B1", SN(a2, 0), -1), Atom("Slow", SN(pc, 0), 0)]
    pM = SN(a1 * E + pc, 0)
    qM = SN(C1, 0) + SN(C2h, -1)
    Np = pair(p, q) + pM * qM
    pS = [Atom("Slow", SN(a1 * z, 2) / Np, 1), Atom("B1", SN(-a2 * g, 1) / Np, -1), Atom("Slow", SN(-th * S_, 2) / Np, 0)]
    out = dict(Np=Np)
    sig_q = pair(pS, q); sig_qb = pair(pS, seq_conj(q))
    def resolvent(m, nz):
        Rs = Atom("Slow", SN(z * C1 * E / (m - Lh), -1), -1)
        Dp = Atom("Dd", SN(C2h * g / k, -1), 1, -nz)
        base = [Rs, Dp]
        s1a = Atom("Slow", SN(1, 0), nz)
        s2a = Atom("BM", SN(1, 0), -nz)
        mu = SN(m, 3)
        def rowM(seq):
            return seqval(seq, "1") * SN(1, 0) + seqval(seq, "N") * SN(1, 0) - seqval(seq, "M") * (SN(1, 0) + mu)
        qN = seqval(q, "N")
        a11 = rowM([s1a]).at(0)
        a12 = SN(-k, 1).at(1)
        rhs1 = (-qN - rowM(base)).at(-1)
        a21 = mass(s1a).at(-4); a22 = mass(s2a).at(-3)
        rhs2 = (-seqmass(base)).at(-5)
        det = a11 * a22 - a12 * a21
        h1 = (rhs1 * a22 - a12 * rhs2) / det
        H2 = (a11 * rhs2 - a21 * rhs1) / det
        return base + [Atom("Slow", SN(h1, -1), nz), Atom("BM", SN(H2, -2), -nz)], det, m - Lh
    r0, det0, mL0 = resolvent(0, 0)
    m2 = 2 * Lh                                       # zeta_2 = -2 m2/(c + sqrt(c^2 + 4 eps m2)) = -2 z at eps = 0
    r2, det2, mL2 = resolvent(m2, -2)
    r0_1 = seqval(r0, "1"); r2_1 = seqval(r2, "1")
    T1 = pair(pS, r0); T2 = pair(pS, seq_conj(r0)); T3 = r0_1.re() * 2 * sig_q
    T4 = pair(pS, r2); T5 = r2_1 * sig_qb
    br = T1 + T2 + T3 + T4 + T5
    bhat = br.at(1)
    out.update(Lhat=((bhat + conj(bhat)) / 2) / (k * t), det0=det0, det2=det2, mL0=mL0, mL2=mL2, br_p=br.p)
    return out

def main():
    import mpmath as mp
    from flint import arb, acb, ctx
    import allm_core as C
    import l1reg as Rg
    out = blocks0()
    Lsym = out["Lhat"]
    print("bracket power at eps = 0: %d (must be <= 1; the value at power 1 is used)" % out["br_p"], flush=True)
    # S1: numerical consistency with l1reg at eps = 0
    ctx.prec = 256; mp.mp.dps = 60
    ok1 = True
    for j in range(1, 9):
        c = (2 * j - mp.mpf(1) / 2) * mp.pi
        tj = mp.findroot(lambda u: mp.tan(u) - u, c - 1 / c)
        kj = mp.sqrt(mp.sqrt(1 + tj * tj) - 1)
        X0 = [arb(mp.nstr(kj, 70)), arb(0), arb(mp.nstr(tj, 70))]
        enc, box = C.certify_hopf(X0, arb(0))
        L0 = acb(Rg.blocks(arb(0), enc[0], acb(enc[1], enc[2]))["Lhat"]).real
        num = sp.N(Lsym.subs({t: sp.Float(mp.nstr(tj, 55), 55), k: sp.Float(mp.nstr(kj, 55), 55),
                              w: sp.exp(I * sp.Float(mp.nstr(tj, 55), 55))}), 45)
        nre, nim = sp.re(num), sp.im(num)
        okj = abs(nim) < sp.Float("1e-35") and bool(abs(L0 - arb(str(sp.N(nre, 45)))) < arb("1e-25"))
        ok1 &= okj
        print("     j=%d: symbolic %s, l1reg at eps=0 %s %s" % (j, sp.N(nre, 28), L0.str(25), "ok" if okj else "MISMATCH"), flush=True)
    report("S1", ok1, "the symbolic eps = 0 port agrees with l1reg.blocks(eps = 0) for j = 1..8")
    # S2: identity in the rational parametrization
    tau = sp.symbols("tau", positive=True)
    Rv = (4 * tau ** 4 + 1) / (4 * tau ** 2); tv = (4 * tau ** 4 - 1) / (4 * tau ** 2); kv = (2 * tau ** 2 - 1) / (2 * tau)
    wv = -(1 + I * tv) / Rv
    sub = {t: tv, k: kv, w: wv}
    Lt = sp.cancel(sp.together(Lsym.subs(sub)))
    Lr8 = -Rv * (Rv - 1) * (7 * Rv ** 2 - 34 * Rv + 23) / (tv * (Rv + 1) * (5 * Rv - 3) * ((Rv - 1) ** 3 + Rv + 1))
    diff = sp.cancel(sp.together(Lt - Lr8))
    report("S2", diff == 0, "Lhat(eps = 0) = result 8's L(R, t) identically on the Hopf curve (rational in tau)")
    # S3: denominators nonzero for tau >= 2
    ok3 = True
    for name in ("Np", "det0", "det2", "mL0", "mL2"):
        e = out[name]
        e = e.v if isinstance(e, SN) else e
        e = sp.cancel(sp.together(sp.sympify(e).subs(sub)))
        num, den = sp.fraction(e)
        for part in (sp.re(sp.expand(num)), sp.im(sp.expand(num))):
            pass
        # the numerator is a polynomial in tau with complex coefficients: it vanishes only where both its real and
        # imaginary parts vanish; it suffices that one of them has no real root >= 2
        nr = sp.Poly(sp.expand(sp.re(sp.expand(num))), tau); ni = sp.Poly(sp.expand(sp.im(sp.expand(num))), tau)
        def no_root_ge2(P):
            if P.is_zero:
                return False
            return all(not (iv[1] >= 2) for iv, m in P.intervals()) if P.degree() > 0 else True
        okn = no_root_ge2(nr) or no_root_ge2(ni)
        dr = sp.Poly(sp.expand(sp.re(sp.expand(den))), tau); di = sp.Poly(sp.expand(sp.im(sp.expand(den))), tau)
        okd = no_root_ge2(dr) or no_root_ge2(di)
        ok3 &= okn and okd
        print("     %s: numerator and denominator nonvanishing for tau >= 2: %s, %s" % (name, okn, okd), flush=True)
    report("S3", ok3, "the divisions of the algorithm are by nonzero quantities at eps = 0 for every Hopf point")
    # S4: sign
    Rs_ = sp.symbols("R", positive=True)
    q7 = sp.Poly(7 * Rs_ ** 2 - 34 * Rs_ + 23, Rs_)
    roots = [float(r) for r in sp.real_roots(q7)]
    # tau >= 2 gives R >= (4*16+1)/16 = 4.0625 > max root (17 + 8 sqrt 2)/7 = 4.0447
    ok4 = all(r < sp.Rational(65, 16) for r in sp.real_roots(q7)) and all(
        sp.Poly(f, tau).intervals() == [] or all(iv[1] < 2 for iv, m in sp.Poly(f, tau).intervals())
        for f in [sp.numer(sp.together(Rv - 1)), sp.numer(sp.together(Rv + 1)), sp.numer(sp.together(5 * Rv - 3)),
                  sp.numer(sp.together((Rv - 1) ** 3 + Rv + 1)), sp.numer(sp.together(tv))])
    report("S4", ok4, "L < 0 for all tau >= 2 (roots of 7R^2 - 34R + 23: %s; R(tau = 2) = 65/16)" % roots)
    print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
    sys.exit(0 if fails == 0 else 1)

if __name__ == "__main__":
    main()
