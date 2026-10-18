#!/usr/bin/env python3
"""v04_sparse_arb.py -- verifier's high-precision O(M) cross-check of l1 at large M (Arb midpoints, 320 bits).

Same algorithm as v03_sparse_float.py (reduced coordinates on the mass hyperplane, tridiagonal + rank-one
Jacobian, Thomas algorithm + Sherman-Morrison, bordered Newton for the eigenpair, secant in K on Re lam = 0,
Kuznetsov's formula), but every operation is done in Arb at 320 bits (default; env V04_PREC).  A first run at 160 bits agreed with the
certificate to 1e-17 at 9 of 12 points but only to ~1e-9 at three j>=2 points with M >= 77777: the midpoint
iteration loses digits there (near-singular shifted solves); at 320 bits all 12 points agree to ~1e-16.  The starting point is the float result
of v03 (itself started from a perturbed kappa), so the Hopf point is re-located here, not copied.
This is a high-precision NUMERICAL cross-check (midpoint iteration; no enclosure of the Hopf point is claimed):
it tests the solver's closed form at large M to ~20 digits instead of the ~5 digits of double precision.
Usage: python v04_sparse_arb.py M:j [M:j ...]
Environment: WSL Fedora, python-flint 0.9, numpy/scipy.  Runtime: ~1-3 min per point at M = 1e5.
"""
import sys, os, csv, gzip, math, time
from flint import arb, acb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v03_sparse_float as v3

ctx.prec = int(os.environ.get("V04_PREC", "320"))

def thomas(n, lo, up, dg_bulk, dg_last, sigma, rhss):
    """solve tridiagonal (lower = lo for rows>=1, upper = up for rows<=n-2, diag = dg_bulk (rows<n-1), dg_last
    (row n-1)) minus sigma*I, for several right-hand sides. Returns list of solutions."""
    nr = len(rhss)
    cp = [None]*n
    dps = [[None]*n for _ in range(nr)]
    d0 = dg_bulk - sigma
    den = d0
    cp[0] = up/den
    for k in range(nr):
        dps[k][0] = rhss[k][0]/den
    for i in range(1, n):
        dgi = d0 if i < n - 1 else dg_last - sigma
        den = dgi - lo*cp[i - 1]
        inv = 1/den
        if i < n - 1:
            cp[i] = up*inv
        for k in range(nr):
            dps[k][i] = (rhss[k][i] - lo*dps[k][i - 1])*inv
    xs = []
    for k in range(nr):
        x = [None]*n
        x[n - 1] = dps[k][n - 1]
        for i in range(n - 2, -1, -1):
            x[i] = dps[k][i] - cp[i]*x[i + 1]
        xs.append(x)
    return xs

def dot(a, b):
    s = acb(0)
    for x, y in zip(a, b):
        s += x*y
    return s

class RedArb:
    def __init__(self, M, K):
        self.M = M; self.n = M - 1; self.K = K; self.A = 1 + K
        self.lvec = [arb(l) for l in range(2, M + 1)]
    def us(self):
        u = [acb(0)]*self.n; u[0] = acb(1); u[self.n - 1] = acb(1); return u
    def vs(self):
        return [-self.A*l for l in self.lvec]
    def matvec(self, y):
        n, A = self.n, self.A
        r = [None]*n
        for i in range(n):
            d = -(A + 1) if i < n - 1 else -A
            s = d*y[i]
            if i >= 1:
                s += A*y[i - 1]
            if i <= n - 2:
                s += y[i + 1]
            r[i] = s
        c = dot(self.vs(), y)
        r[0] += c; r[n - 1] += c
        return r
    def solve(self, sigma, rhss, transpose=False):
        n, A = self.n, self.A
        if not transpose:
            lo, up = A, arb(1); u_, v_ = self.us(), self.vs()
        else:
            lo, up = arb(1), A; u_, v_ = self.vs(), self.us()
        sols = thomas(n, lo, up, -(A + 1), -A, sigma, list(rhss) + [u_])
        z2 = sols[-1]
        den = mid(1 + dot(v_, z2))   # midpoint: near an eigenvalue this ball may contain 0 (numerical check only)
        out = []
        for z1 in sols[:-1]:
            c = dot(v_, z1)/den
            out.append([a - b*c for a, b in zip(z1, z2)])
        return out

def mid(z):
    return acb(z.real.mid(), z.imag.mid())

def eigpair(R, lam, q, k, it=6):
    for _ in range(it):
        Aq = R.matvec(q)
        r = [-(a - lam*b) for a, b in zip(Aq, q)]
        z1, z2 = R.solve(lam, [r, q])
        dlam = -z1[k]/z2[k]
        q = [mid(a + b + dlam*c) for a, b, c in zip(q, z1, z2)]
        lam = mid(lam + dlam)
        if abs(dlam).mid() < arb(2)**(-110)*abs(lam).mid():   # next step would be below the working precision
            break
    return lam, q

def E(y, M):
    s = acb(0)
    for l, x in enumerate(y):
        s -= (l + 2)*x
    return [s] + list(y)

def bil_full(X, Y):
    M = len(X); N = M - 1
    D = [X[l]*Y[0] + X[0]*Y[l] for l in range(N)]
    f = [None]*M
    sD = acb(0)
    for d in D:
        sD += d
    f[0] = -D[0] - sD
    for l in range(1, N):
        f[l] = D[l - 1] - D[l]
    f[M - 1] = D[N - 1]
    return f

def bil_red(x, y, M):
    return bil_full(E(x, M), E(y, M))[1:]

def run(M, j, cert):
    t0 = time.time()
    r = cert[M]
    Kf, lamf, qf, Rf = v3.hopf(M, float(r["kappa_mid"]), float(r["omega_lo"]))
    k = max(range(M - 1), key=lambda i: abs(qf[i]))
    qf = qf/qf[k]
    q = [acb(complex(z).real, complex(z).imag) for z in qf]
    lam = acb(lamf.real, lamf.imag)
    # secant in K on Re lam(K) = 0
    Ks = []; vals = []
    K = arb(Kf)
    for K_ in (K, K*(1 + arb(2)**-40)):
        R = RedArb(M, K_); lam, q = eigpair(R, lam, q, k)
        Ks.append(K_); vals.append(lam.real)
    for _ in range(8):
        K_ = (Ks[-1] - vals[-1]*(Ks[-1] - Ks[-2])/(vals[-1] - vals[-2])).mid()
        K_ = arb(K_)
        R = RedArb(M, K_); lam, q = eigpair(R, lam, q, k)
        Ks.append(K_); vals.append(lam.real)
        if abs(Ks[-1] - Ks[-2]).mid() < arb(2)**(-145)*K_.mid():
            break
    K = Ks[-1]; R = RedArb(M, K)
    om = lam.imag
    n = M - 1
    ones = [acb(1)]*n
    shift = lam*(1 + arb(2)**-100)                 # slightly off the eigenvalue: inverse iteration stays finite
    u, = R.solve(shift, [ones], transpose=True)
    u = [mid(x) for x in u]
    u, = R.solve(shift, [u], transpose=True)
    u = [mid(x) for x in u]
    uq = dot(u, q)
    u = [x/uq for x in u]
    qb = [x.conjugate() for x in q]
    w11, = R.solve(acb(0), [bil_red(q, qb, M)])
    w20, = R.solve(2*acb(0, 1)*om, [[-x for x in bil_red(q, q, M)]])
    c = -2*dot(u, bil_red(q, w11, M)) + dot(u, bil_red(qb, w20, M))
    l1 = c.real/(2*om)
    Eq = E(q, M)
    nrm = arb(0)
    for x in Eq:
        nrm += abs(x)**2
    l1u = l1/nrm; l1q = l1/abs(Eq[0])**2
    kap = K*arb(M).sqrt()
    their = float(r["l1_unit_lo"]); their_q = float(r["l1_q1_lo"])
    relu = abs(float((l1u - their)/their)); relq = abs(float((l1q - their_q)/their_q))
    relk = abs(float((kap - arb(r["kappa_mid"]))/kap))
    print("M=%d j=%d kappa=%s (solver %s, rel %.1e) Re lam=%.1e  l1_unit=%s (solver %.16e, rel %.1e)  "
          "l1_q1 rel %.1e  %.0fs" % (M, j, kap.mid().str(25, radius=False), r["kappa_mid"][:27], relk,
                                     float(abs(lam.real).mid()), l1u.mid().str(18, radius=False), their, relu, relq,
                                     time.time() - t0), flush=True)
    return relu < 1e-14 and relq < 1e-14 and relk < 1e-20 and ((l1u > 0) == ((M, j) == (15, 1)))

def main():
    pts = [tuple(int(s) for s in x.split(":")) for x in sys.argv[1:]]
    d = os.path.join(HERE, "..", "certificates")
    certs = {}
    bad = 0
    for (M, j) in pts:
        if j not in certs:
            certs[j] = v3.load(d, j)
        bad += not run(M, j, certs[j])
    print("SUMMARY: %d points, %s (float64 in the csv limits the comparison to ~1e-16)"
          % (len(pts), "PASS" if bad == 0 else "FAIL (%d)" % bad))
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
