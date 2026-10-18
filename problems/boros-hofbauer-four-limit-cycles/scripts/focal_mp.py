"""
focal_mp.py -- floating-point (mpmath, arbitrary precision) focal values for planar systems
with generalized-monomial right-hand sides and an equilibrium at (1,1).

This module is used ONLY for exploration and cross-validation (it is not part of the certificate).

Method (complex Lyapunov-function method):
  1. Taylor-expand x^al y^be = (1+u)^al (1+v)^be at (1,1).
  2. Linear change X = u, Y = -(A11 u + A12 v)/omega and time tau = omega t give
        X' = -Y + P(X,Y),  Y' = X + Q(X,Y)            (requires tr A = 0, det A = omega^2 > 0).
  3. With z = X + iY:  z' = i z + R(z, zb),  R = P + iQ.
  4. Find F = z zb + sum_{p+q>=3} f_pq z^p zb^q with  dF/dtau = sum_k V_k (z zb)^{k+1}.
     For p != q: f_pq = i S_pq/(p-q); for p = q = k+1: V_k = S_pq, f_pp = 0, where S_pq is the
     coefficient of z^p zb^q in sum_{lower degree} f_rs (r z^{r-1} zb^s R + s z^r zb^{s-1} conj(R)).
  If V_1 = ... = V_{k-1} = 0, the displacement map on a ray is  d(r) = pi V_k r^{2k+1} + O(r^{2k+2})
  (time-normalised coordinates).
"""
import mpmath as mp


def binom(al, i):
    r = mp.mpf(1)
    for k in range(i):
        r = r * (al - k) / (k + 1)
    return r


def pmul(A, B, N):
    C = {}
    for (i1, j1), c1 in A.items():
        for (i2, j2), c2 in B.items():
            if i1 + j1 + i2 + j2 <= N:
                k = (i1 + i2, j1 + j2)
                C[k] = C.get(k, 0) + c1 * c2
    return C


def taylor_monomials(terms, N):
    """terms: list of (coef, al, be) meaning coef * x^al * y^be. Returns dict (i,j)->coef of u^i v^j."""
    T = {}
    for (cf, al, be) in terms:
        for i in range(N + 1):
            bi = binom(al, i)
            for j in range(N + 1 - i):
                T[(i, j)] = T.get((i, j), 0) + cf * bi * binom(be, j)
    return T


def normalise(Tf, Tg, N):
    """Return (P, Q, omega, trace, A) for the time-normalised system."""
    A11, A12 = Tf.get((1, 0), 0), Tf.get((0, 1), 0)
    A21, A22 = Tg.get((1, 0), 0), Tg.get((0, 1), 0)
    tr = A11 + A22
    det = A11 * A22 - A12 * A21
    om = mp.sqrt(det)
    # u = X, v = al*X + be*Y
    al = -A11 / A12
    be = -om / A12
    V = {(1, 0): al, (0, 1): be}
    Upow = [{(0, 0): mp.mpf(1)}]
    Vpow = [{(0, 0): mp.mpf(1)}]
    for k in range(1, N + 1):
        Upow.append({(k, 0): mp.mpf(1)})
        Vpow.append(pmul(Vpow[-1], V, N))

    def subst(T):
        S = {}
        for (i, j), c in T.items():
            if i + j < 2 or i + j > N or c == 0:
                continue
            M = pmul(Upow[i], Vpow[j], N)
            for k, v in M.items():
                S[k] = S.get(k, 0) + c * v
        return S
    N1, N2 = subst(Tf), subst(Tg)
    P = {k: v / om for k, v in N1.items()}
    Q = {}
    for k in set(N1) | set(N2):
        Q[k] = -(A11 * N1.get(k, 0) + A12 * N2.get(k, 0)) / om**2
    return P, Q, om, tr, (A11, A12, A21, A22)


def to_complex(P, Q, N):
    """R(z,zb) = P + iQ with X = (z+zb)/2, Y = (z-zb)/(2i)."""
    I = mp.mpc(0, 1)
    zx = {(1, 0): mp.mpf(1) / 2, (0, 1): mp.mpf(1) / 2}
    zy = {(1, 0): 1 / (2 * I), (0, 1): -1 / (2 * I)}
    Xp = [{(0, 0): mp.mpc(1)}]
    Yp = [{(0, 0): mp.mpc(1)}]
    for k in range(1, N + 1):
        Xp.append(pmul(Xp[-1], zx, N))
        Yp.append(pmul(Yp[-1], zy, N))
    R = {}
    for k in set(P) | set(Q):
        c = P.get(k, 0) + I * Q.get(k, 0)
        if c == 0:
            continue
        M = pmul(Xp[k[0]], Yp[k[1]], N)
        for kk, v in M.items():
            R[kk] = R.get(kk, 0) + c * v
    return R


def lyapunov(R, kmax):
    """Return [V_1..V_kmax] (complex numbers whose imaginary parts must vanish)."""
    N = 2 * kmax + 2
    I = mp.mpc(0, 1)
    Rb = {(q, p): mp.conj(c) for (p, q), c in R.items()}
    F = {(1, 1): mp.mpc(1)}
    Vs = []
    for n in range(3, N + 1):
        S = {}
        for (r, s), fc in F.items():
            m = r + s
            deg = n - m + 1  # degree of R needed
            for (p, q), rc in R.items():
                if p + q != deg:
                    continue
                if r > 0:
                    k = (r - 1 + p, s + q)
                    S[k] = S.get(k, 0) + fc * r * rc
            for (p, q), rc in Rb.items():
                if p + q != deg:
                    continue
                if s > 0:
                    k = (r + p, s - 1 + q)
                    S[k] = S.get(k, 0) + fc * s * rc
        for p in range(n + 1):
            q = n - p
            sv = S.get((p, q), 0)
            if p != q:
                F[(p, q)] = I * sv / (p - q)
            else:
                Vs.append(sv)
    return Vs


def focal_values(terms_f, terms_g, kmax, dps=60):
    with mp.workdps(dps):
        N = 2 * kmax + 2
        Tf = taylor_monomials(terms_f, N)
        Tg = taylor_monomials(terms_g, N)
        P, Q, om, tr, A = normalise(Tf, Tg, N)
        R = to_complex(P, Q, N)
        Vs = lyapunov(R, kmax)
        return Vs, om, tr, A


def bh_system(a, b, d, K=None):
    """Boros-Hofbauer Sec. 5.1: x' = 1/y - x^a y^b, y' = K(-(d-1) - 1/y + d x^a y^b)."""
    if K is None:
        K = a / (1 + b * d)
    tf = [(mp.mpf(1), 0, -1), (mp.mpf(-1), a, b)]
    tg = [(-K * (d - 1), 0, 0), (-K, 0, -1), (K * d, a, b)]
    return tf, tg


if __name__ == "__main__":
    mp.mp.dps = 50
    a, b, d = mp.mpf('1.01282'), mp.mpf('0.65463'), mp.mpf('3.28862')
    tf, tg = bh_system(a, b, d)
    Vs, om, tr, A = focal_values(tf, tg, 4)
    print("trace", tr, "omega", om)
    for k, v in enumerate(Vs, 1):
        print("V%d" % k, mp.nstr(v, 15))
