"""bdlib -- Becker-Doring equations with linear atomization (Pego & Velazquez, Nonlinearity 33 (2020)).

Library used by all scripts of this folder.

Model (paper eqs. (1)-(4)), sizes l = 1..M, N = M-1, J_l = n_l n_1 - n_{l+1}:
    dn_l/dt = J_{l-1} - J_l           (2 <= l <= N)
    dn_M/dt = J_N - K n_M
    dn_1/dt = -J_1 - sum_{l=1}^N J_l + M K n_M
Constant equilibrium n = A = 1+K.  With n = A + v the vector field is exactly
    dv/dt = B v + v_1 S v                  (purely quadratic nonlinearity)
with B the paper's matrix (27) and S the linear map
    (S v)_1 = -2 v_1 - sum_{l=2}^N v_l,  (S v)_l = v_{l-1} - v_l (2<=l<=N),  (S v)_M = v_N.
Kuznetsov's bilinear form is therefore Bil(x, y) = x_1 S y + y_1 S x and the trilinear form is 0.

Two independent implementations of the first Lyapunov coefficient l1 are provided:
  * dense floating point (numpy/scipy), for small M;
  * closed form in ball arithmetic (python-flint acb), cost O(1) in M: every vector that occurs is an
    exponential sum  v_l = sum_e c_e beta_e^l  plus point corrections at l in {1, 2, N, M}.
Hopf points are enclosed with the Krawczyk operator in the scaled unknowns
    kappa = K sqrt(M),  z = M (phi - 1),   lambda = (A - 1/phi)(phi - 1)   (paper (31), (36), Sec. 6),
for the system  K^-3 F(phi) = 0 (complex),  Re(M lambda / K) = 0, with F the paper function (71)/(76).
"""
import numpy as np
from flint import acb, arb, acb_mat, arb_mat, ctx

# ----------------------------------------------------------------------------------------------
# Float model
# ----------------------------------------------------------------------------------------------

def field(n, K):
    M = len(n); N = M - 1
    J = n[:N]*n[0] - n[1:]
    f = np.zeros(M)
    f[1:N] = J[0:N-1] - J[1:N]
    f[M-1] = J[N-1] - K*n[M-1]
    f[0] = -J[0] - J.sum() + M*K*n[M-1]
    return f

def Bmat(K, M):
    A = 1 + K; N = M - 1
    B = np.zeros((M, M))
    B[0, 0] = -A*(M+2); B[0, 1] = 1 - K; B[0, 2:N] = -K; B[0, M-1] = M*K + 1
    for l in range(1, N):
        B[l, l-1] += A; B[l, l] += -A - 1; B[l, l+1] += 1
    B[M-1, 0] += A; B[M-1, N-1] += A; B[M-1, M-1] += -A
    return B

def Smat(M):
    N = M - 1
    S = np.zeros((M, M))
    S[0, 0] = -2; S[0, 1:N] = -1
    for l in range(1, N):
        S[l, l-1] = 1; S[l, l] = -1
    S[M-1, N-1] = 1
    return S

def nullvec(K, M):
    """right null vector of B (paper (29)), scaled to avoid overflow"""
    A = 1 + K; l = np.arange(1, M+1)
    return 1 + (A**(l - M + 1.0) - A**(2.0 - M))/K

def jac_full(n, K):
    """Jacobian of field at n (exact formula B + v_1 S + (S v) e_1^T)"""
    M = len(n); A = 1 + K; v = n - A
    S = Smat(M)
    J = Bmat(K, M) + v[0]*S
    J[:, 0] += S @ v
    return J

def lyap_dense(K, M, lam_target):
    """Kuznetsov first Lyapunov coefficient, dense float; q normalized by q^H q = 1, <p,q> = 1."""
    import scipy.linalg as sla
    B = Bmat(K, M); S = Smat(M)
    w, VL, VR = sla.eig(B, left=True, right=True)
    i = np.argmin(np.abs(w - lam_target))
    lam = w[i]; q = VR[:, i]; q = q/np.sqrt(np.vdot(q, q).real)
    p = VL[:, i]; p = p/np.conj(np.vdot(p, q))
    om = lam.imag
    Bil = lambda x, y: x[0]*(S@y) + y[0]*(S@x)
    ell = np.arange(1, M+1); vbar = nullvec(K, M)
    Mb = np.zeros((M+1, M+1)); Mb[:M, :M] = B; Mb[:M, M] = vbar; Mb[M, :M] = ell
    w11 = np.linalg.solve(Mb.astype(complex), np.concatenate([Bil(q, np.conj(q)), [0]]))[:M]
    w20 = np.linalg.solve(2j*om*np.eye(M) - B, Bil(q, q))
    c = -2*np.vdot(p, Bil(q, w11)) + np.vdot(p, Bil(np.conj(q), w20))
    return lam, c.real/(2*om), q

def unstable_pairs(K, M):
    w = np.linalg.eigvals(Bmat(K, M))
    return int(np.sum((w.real > 1e-12*max(1.0, K)) & (w.imag > 1e-12)))

def dense_crossings(M, Kmin=1e-3, Kmax=1e5, n=900):
    """float scan of K for changes in the number of unstable complex pairs"""
    Ks = np.exp(np.linspace(np.log(Kmin), np.log(Kmax), n))
    res = []; prev = unstable_pairs(Ks[0], M)
    for a, b in zip(Ks[:-1], Ks[1:]):
        c = unstable_pairs(b, M)
        if c != prev:
            res.append((a, b, prev, c))
        prev = c
    return res

def dense_hopf(M, a, b, nprev):
    """refine the crossing in [a,b] where the (nprev+1)-th pair becomes unstable (brentq on Re)"""
    import scipy.optimize as so
    def lam_new(K):
        w = np.linalg.eigvals(Bmat(K, M))
        w = w[w.imag > 1e-12]
        w = w[np.argsort(-w.real)]
        return w[nprev]
    K = so.brentq(lambda K: lam_new(K).real, a, b, xtol=1e-15, rtol=1e-15)
    lam = lam_new(K)
    A = 1 + K
    r = np.roots([A, -(A + 1 + lam), 1]); phi = r[np.argmax(np.abs(r))]
    return K, lam, phi

# ----------------------------------------------------------------------------------------------
# Closed form in ball arithmetic
# ----------------------------------------------------------------------------------------------

ZERO = (0, 0, 0, 0, 0)

def eadd(e, f):
    return tuple(a + b for a, b in zip(e, f))

class Problem:
    """Generators of the bases: (phi, conj(phi), A, s, conj(s)); a base is a monomial with integer exponents."""
    def __init__(self, M, K, phi):
        self.M = int(M); self.N = self.M - 1
        self.K = acb(K); self.A = 1 + self.K
        self.phi = acb(phi)
        self.gens = [self.phi, self.phi.conjugate(), self.A, None, None]
        self._base = {}; self._pow = {}
    def set_s(self, s):
        self.gens[3] = acb(s); self.gens[4] = acb(s).conjugate()
        self._base = {}; self._pow = {}
    def base(self, e):
        if e not in self._base:
            b = acb(1)
            for g, k in zip(self.gens, e):
                if k != 0:
                    b = b*(g**k)
            self._base[e] = b
        return self._base[e]
    def bpow(self, e, l):
        key = (e, l)
        if key not in self._pow:
            b = acb(1)
            for g, k in zip(self.gens, e):
                if k != 0:
                    b = b*(g**(k*l))
            self._pow[key] = b
        return self._pow[key]
    def geom(self, e, a, b):
        """sum_{l=a}^{b} beta_e^l  (exact count when the base is identically 1)"""
        if b < a:
            return acb(0)
        if e == ZERO:
            return acb(b - a + 1)
        be = self.base(e)
        return (self.bpow(e, b + 1) - self.bpow(e, a))/(be - 1)
    def lgeom(self, e):
        """sum_{l=1}^{M} l beta_e^l"""
        M = self.M
        if e == ZERO:
            return acb(M*(M + 1)//2)
        be = self.base(e)
        return be*(1 - (M + 1)*self.bpow(e, M) + M*self.bpow(e, M + 1))/(1 - be)**2

class Seq:
    """v_l = sum_e c_e beta_e^l + pts[l], l = 1..M"""
    def __init__(self, P, terms=None, pts=None):
        self.P = P
        self.terms = {k: acb(v) for k, v in (terms or {}).items()}
        self.pts = {k: acb(v) for k, v in (pts or {}).items()}
    def copy(self):
        return Seq(self.P, self.terms, self.pts)
    def exp_val(self, l):
        s = acb(0)
        for e, c in self.terms.items():
            s += c*self.P.bpow(e, l)
        return s
    def val(self, l):
        return self.exp_val(l) + self.pts.get(l, acb(0))
    def rsum(self, a, b):
        s = acb(0)
        for e, c in self.terms.items():
            s += c*self.P.geom(e, a, b)
        for k, v in self.pts.items():
            if a <= k <= b:
                s += v
        return s
    def mass(self):
        s = acb(0)
        for e, c in self.terms.items():
            s += c*self.P.lgeom(e)
        for k, v in self.pts.items():
            s += k*v
        return s
    def __add__(self, o):
        r = self.copy()
        for e, c in o.terms.items():
            r.terms[e] = r.terms.get(e, acb(0)) + c
        for k, v in o.pts.items():
            r.pts[k] = r.pts.get(k, acb(0)) + v
        return r
    def scale(self, c):
        c = acb(c)
        return Seq(self.P, {e: c*v for e, v in self.terms.items()}, {k: c*v for k, v in self.pts.items()})
    def __sub__(self, o):
        return self + o.scale(-1)
    def conj(self):
        t = {(e[1], e[0], e[2], e[4], e[3]): v.conjugate() for e, v in self.terms.items()}
        return Seq(self.P, t, {k: v.conjugate() for k, v in self.pts.items()})

def _apply(w, mult, rowfun, rows):
    P = w.P
    r = Seq(P, {e: c*mult(P.base(e)) for e, c in w.terms.items()})
    for row in sorted(rows):
        if 1 <= row <= P.M:
            r.pts[row] = rowfun(w, row) - r.exp_val(row)
    return r

def B_row(w, r):
    P = w.P; M, N, A, K = P.M, P.N, P.A, P.K
    if r == 1:
        return -A*(M + 2)*w.val(1) + (1 - K)*w.val(2) - K*w.rsum(3, N) + (M*K + 1)*w.val(M)
    if r == M:
        return A*w.val(1) + A*w.val(N) - A*w.val(M)
    return A*w.val(r - 1) - (A + 1)*w.val(r) + w.val(r + 1)

def applyB(w):
    P = w.P; A = P.A
    rows = {1, P.M}
    for k in w.pts:
        rows |= {k - 1, k, k + 1}
    return _apply(w, lambda b: A/b - (A + 1) + b, B_row, rows)

def S_row(w, r):
    P = w.P; M, N = P.M, P.N
    if r == 1:
        return -2*w.val(1) - w.rsum(2, N)
    if r == M:
        return w.val(N)
    return w.val(r - 1) - w.val(r)

def applyS(w):
    P = w.P
    rows = {1, P.M}
    for k in w.pts:
        rows |= {k, k + 1}
    return _apply(w, lambda b: 1/b - 1, S_row, rows)

def Bil(x, y):
    return applyS(y).scale(x.val(1)) + applyS(x).scale(y.val(1))

def inner(u, v):
    """sum_l u_l v_l (bilinear)"""
    P = u.P; M = P.M
    s = acb(0)
    for e, c in u.terms.items():
        for f, d in v.terms.items():
            s += c*d*P.geom(eadd(e, f), 1, M)
    for k in set(u.pts) | set(v.pts):
        s += u.val(k)*v.val(k) - u.exp_val(k)*v.exp_val(k)
    return s

def solve(r, mu, hbases, singular=False):
    """Exact solution of (B - mu) w = r by the ansatz  particular + h1 s1^l + h2 s2^l + e1 delta_1 + eM delta_M.
    The particular part solves the bulk rows term by term; the 4 unknowns are fixed by rows 1, 2, N, M
    (row 1 replaced by the mass condition sum l w_l = 0 when mu = 0, i.e. w is taken in the hyperplane H)."""
    P = r.P; A = P.A; M, N = P.M, P.N
    for k in r.pts:
        assert k in (1, 2, N, M), k
    part = {}
    for e, c in r.terms.items():
        b = P.base(e)
        part[e] = c/(A/b - (A + 1) + b - mu)
    w0 = Seq(P, part)
    basis = [Seq(P, {hbases[0]: 1}), Seq(P, {hbases[1]: 1}), Seq(P, None, {1: 1}), Seq(P, None, {M: 1})]
    def res(w, inhom):
        R = applyB(w) - w.scale(mu)
        return R - r if inhom else R
    rows = [1, 2, N, M]
    R0 = res(w0, True)
    Rb = [res(b, False) for b in basis]
    G = [[Rb[j].val(row) for j in range(4)] for row in rows]
    rhs = [-R0.val(row) for row in rows]
    if singular:
        G[0] = [b.mass() for b in basis]
        rhs[0] = -w0.mass()
    X = acb_mat(G).solve(acb_mat([[x] for x in rhs]))
    w = w0
    for j in range(4):
        w = w + basis[j].scale(X[j, 0])
    return w

def eigvec_right(P):
    """q_l = c1 phi^(M-l) + c2 psi^(M-l), psi = 1/(A phi), (c1, c2) = (f(psi), -f(phi))  (paper Sec. 5, 10.3)"""
    M, A, phi = P.M, P.A, P.phi
    psi = 1/(A*phi)
    c1 = A*psi**(M - 1) + 1 - A*phi
    c2 = -(A*phi**(M - 1) + 1 - 1/phi)
    return Seq(P, {(-1, 0, 0, 0, 0): c1*phi**M, (1, 0, 1, 0, 0): c2*psi**M})

def eigvec_left(P, lam):
    """row vector u with u B = lam u:  u_k = a1 phi^k + a2 psi^k + uc,  uc = -K u_1/lam"""
    M, N, A, K = P.M, P.N, P.A, P.K
    phi = P.phi; psi = 1/(A*phi)
    def useq(a1, a2):
        uc = -K*(a1*phi + a2*psi)/(lam + K)
        return Seq(P, {(1, 0, 0, 0, 0): a1, (-1, 0, -1, 0, 0): a2, ZERO: uc})
    def eqs(u):
        e1 = -A*(M + 2)*u.val(1) + A*u.val(2) + A*u.val(M) - lam*u.val(1)
        eM = (M*K + 1)*u.val(1) + u.val(N) - A*u.val(M) - lam*u.val(M)
        return e1, eM
    m11, m21 = eqs(useq(acb(1), acb(0)))
    m12, m22 = eqs(useq(acb(0), acb(1)))
    if abs(m11.mid()) + abs(m12.mid()) >= abs(m21.mid()) + abs(m22.mid()):
        a1, a2 = m12, -m11
    else:
        a1, a2 = m22, -m21
    return useq(a1, a2)

def lyapunov_cf(M, K, phi, normalize="unit", details=False):
    """First Lyapunov coefficient (Kuznetsov (3.20), trilinear term 0) at an exact Hopf point (K, phi).
    normalize: 'unit' -> q^H q = 1 ; 'q1' -> q_1 = 1.  Ball arithmetic: the result encloses the exact value
    whenever the input balls contain the exact Hopf point."""
    P = Problem(M, K, phi)
    A = P.A
    lam = (A - 1/P.phi)*(P.phi - 1)
    om = lam.imag
    q = eigvec_right(P)
    if normalize == "unit":
        q = q.scale(1/inner(q.conj(), q).real.sqrt())
    elif normalize == "q1":
        q = q.scale(1/q.val(1))
    u = eigvec_left(P, lam)
    u = u.scale(1/inner(u, q))
    qb = q.conj()
    w11 = solve(Bil(q, qb), acb(0), [ZERO, (0, 0, 1, 0, 0)], singular=True)
    mu = 2*acb(0, 1)*om
    disc = ((A + 1 + mu)**2 - 4*A).sqrt()
    P.set_s(((A + 1 + mu) + disc)/2)
    w20 = solve(Bil(q, q).scale(-1), mu, [(0, 0, 0, 1, 0), (0, 0, 1, -1, 0)])
    t1 = inner(u, Bil(q, w11))
    t2 = inner(u, Bil(qb, w20))
    l1 = (-2*t1 + t2).real/(2*om)
    if details:
        return dict(l1=l1, lam=lam, q=q, u=u, w11=w11, w20=w20, P=P)
    return l1

def dense_vec(s):
    return [s.val(l) for l in range(1, s.P.M + 1)]

# ----------------------------------------------------------------------------------------------
# Hopf points: scaled system and Krawczyk operator
# ----------------------------------------------------------------------------------------------

class D:
    """forward-mode dual numbers over acb: partials w.r.t. z (holomorphic) and kappa"""
    __slots__ = ("v", "dz", "dk")
    def __init__(self, v, dz=0, dk=0):
        self.v = acb(v); self.dz = acb(dz); self.dk = acb(dk)
    def __add__(s, o):
        o = o if isinstance(o, D) else D(o)
        return D(s.v + o.v, s.dz + o.dz, s.dk + o.dk)
    __radd__ = __add__
    def __neg__(s):
        return D(-s.v, -s.dz, -s.dk)
    def __sub__(s, o):
        return s + (-(o if isinstance(o, D) else D(o)))
    def __rsub__(s, o):
        return D(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, D) else D(o)
        return D(s.v*o.v, s.dz*o.v + s.v*o.dz, s.dk*o.v + s.v*o.dk)
    __rmul__ = __mul__
    def inv(s):
        iv = 1/s.v
        return D(iv, -s.dz*iv*iv, -s.dk*iv*iv)
    def __truediv__(s, o):
        o = o if isinstance(o, D) else D(o)
        return s*o.inv()
    def __rtruediv__(s, o):
        return D(o)*s.inv()
    def __pow__(s, n):
        n = int(n)
        p = s.v**(n - 1)
        return D(p*s.v, n*p*s.dz, n*p*s.dk)

def hopf_system(kap, z, M):
    """returns Dual numbers (Qeps = K^-3 F(phi), Lam = M lambda / K) at kappa = kap, z"""
    M = int(M)
    kd = D(kap, 0, 1); zd = D(z, 1, 0)
    K = kd*(1/arb(M).sqrt())
    A = 1 + K
    phi = 1 + zd*(arb(1)/M)
    pM = phi**M
    Ap = A*phi
    ApM = Ap**(-M)
    f1 = (phi - 1) + A*pM
    f2 = (1 - Ap) + ApM*A*A*phi
    G1 = phi*(1 - Ap)*(K*phi + (K*M + 1)*(phi - 1))
    G2 = (1 - Ap)*(A*M*(phi - 1) + K*phi + A*phi*(phi - 1))
    H1 = (phi - 1)*(K + (K*M + 1)*(1 - Ap))
    H2 = (phi - 1)*((A + M*A*A*phi)*(1 - Ap) + K*A*phi)
    F = f1*(H1 - ApM*H2) - (G1 - pM*G2)*f2
    Qe = F/(K*K*K)
    lam = (A - 1/phi)*(phi - 1)
    Lam = lam*M/K
    return Qe, Lam

def fvec_jac(X, M):
    kap, x, y = X
    Qe, Lam = hopf_system(kap, acb(x, y), M)
    f = [Qe.v.real, Qe.v.imag, Lam.v.real]
    J = [[Qe.dk.real, Qe.dz.real, -Qe.dz.imag],
         [Qe.dk.imag, Qe.dz.imag, Qe.dz.real],
         [Lam.dk.real, Lam.dz.real, -Lam.dz.imag]]
    return f, J

def newton(X0, M, it=80):
    X = [arb(arb(x).mid()) for x in X0]
    for _ in range(it):
        f, J = fvec_jac(X, M)
        Jm = arb_mat([[c.mid() for c in row] for row in J])
        dX = Jm.solve(arb_mat([[c.mid()] for c in f]))
        X = [arb((X[i] - dX[i, 0]).mid()) for i in range(3)]
        if max(abs(float(dX[i, 0].mid())) for i in range(3)) < 2.0**(-ctx.prec + 30)*(1 + abs(float(X[0].mid()))):
            break
    return X

def krawczyk(Xm, M, rad):
    box = [arb(Xm[i], rad[i]) for i in range(3)]
    fm, Jm = fvec_jac(Xm, M)
    _, JX = fvec_jac(box, M)
    Y = arb_mat([[c.mid() for c in row] for row in Jm]).inv()
    Y = arb_mat([[Y[i, j].mid() for j in range(3)] for i in range(3)])
    I3 = arb_mat([[1 if i == j else 0 for j in range(3)] for i in range(3)])
    dev = arb_mat([[box[i] - Xm[i]] for i in range(3)])
    Kx = arb_mat([[Xm[i]] for i in range(3)]) - Y*arb_mat([[c] for c in fm]) + (I3 - Y*arb_mat(JX))*dev
    if all(box[i].contains_interior(Kx[i, 0]) for i in range(3)):
        return [Kx[i, 0] for i in range(3)], box
    return None, box

def certify_hopf(X0, M, rel=1e-25):
    """Newton + Krawczyk. Returns (enclosure, box) where the box contains a unique zero of the 3x3 system
    and the enclosure (subset of the box) contains it; or (None, None)."""
    X = newton(X0, M)
    for r in [rel, rel*1e5, rel*1e10, rel*1e15]:
        rad = [arb(abs(float(x.mid()))*r + r) for x in X]
        enc, box = krawczyk(X, M, rad)
        if enc is not None:
            return enc, box
    return None, None

def hopf_checks(enc, M):
    """Nondegeneracy at the enclosed Hopf point:
       simple: 0 not in dQeps/dz (=> simple root of F => simple eigenvalue, paper Lemma 5.4 / Sec. 10.3)
       transv: lower bound of Re(Lam + kappa dLam/dkappa) along the root curve (Re dlambda/dkappa = M^-3/2 * this)
       omega : enclosure of the frequency Im lambda."""
    kap, x, y = enc
    Qe, Lam = hopf_system(kap, acb(x, y), M)
    simple = not Qe.dz.contains(0)
    dzdk = -Qe.dk/Qe.dz
    dLam = Lam.dz*dzdk + Lam.dk
    tr = (Lam.v + kap*dLam).real
    K = kap/arb(M).sqrt()
    phi = 1 + acb(x, y)/M
    lam = (1 + K - 1/phi)*(phi - 1)
    return simple, tr, lam.imag

def certify_point(X0, M, prec=192):
    """Full certificate for one (M, branch): returns dict or None"""
    old = ctx.prec
    ctx.prec = prec
    try:
        enc, box = certify_hopf(X0, M)
        if enc is None:
            return None
        kap, x, y = enc
        K = kap/arb(M).sqrt()
        phi = 1 + acb(x, y)/M
        l1 = lyapunov_cf(M, K, phi, normalize="unit")
        l1q = lyapunov_cf(M, K, phi, normalize="q1")
        simple, tr, om = hopf_checks(enc, M)
        return dict(M=M, kappa=kap, x=x, y=y, l1=l1, l1q1=l1q, simple=simple, transv=tr, omega=om, box=box)
    finally:
        ctx.prec = old
