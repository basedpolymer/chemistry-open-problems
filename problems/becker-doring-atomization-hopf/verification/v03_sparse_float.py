#!/usr/bin/env python3
"""v03_sparse_float.py -- verifier's independent FLOATING-POINT cross-check of l1 for large M (O(M) per point).

Not a proof: double precision, used to cross-check the solver's closed form at large M by a different algorithm.
Method (written by the verifier, no use of bdlib.py):
  * reduced coordinates on the mass hyperplane, y = (v_2..v_M), v_1 = -sum_{l>=2} l v_l;
  * the reduced Jacobian is  A_hat = T - A (e_2 + e_M) l^T  with T tridiagonal (rows l = 2..M of the paper's
    matrix with the v_1 column removed): every solve is a banded solve plus a Sherman-Morrison correction;
  * Hopf point: Newton on the bordered eigen-system (A_hat - lam) q = 0, q_k = 1, and secant iteration in K on
    Re lam(K) = 0, started from a deliberately perturbed kappa (1e-4 relative) so that the Hopf point is
    re-located, not copied;
  * bilinear form from the second derivative of the fluxes, D^2 J_l [X, Y] = X_l Y_1 + X_1 Y_l, assembled with the
    equations (2)-(4) of the paper (checked against polarisation of the field for small M);
  * Kuznetsov's formula with u q = 1, normalisations 'unit' (|E q| = 1) and q_1 = 1.
Usage: python v03_sparse_float.py OUT.csv SEED NPERBRANCH [extra M:j ...]
Environment: WSL Fedora, numpy/scipy. Runtime: ~0.1-1 s per point.
"""
import sys, os, csv, gzip, math, random, time
import numpy as np
from scipy.linalg import solve_banded

def field(n, K):
    M = len(n); N = M - 1
    J = n[:N]*n[0] - n[1:]
    f = np.empty(M, dtype=np.result_type(n, float))
    f[0] = -J[0] - J.sum() + M*K*n[M - 1]
    f[1:N] = J[:N - 1] - J[1:]
    f[M - 1] = J[N - 1] - K*n[M - 1]
    return f

def bil_full(X, Y):
    """D^2 f [X, Y] from D^2 J_l[X,Y] = X_l Y_1 + X_1 Y_l"""
    M = len(X); N = M - 1
    D = X[:N]*Y[0] + X[0]*Y[:N]
    f = np.empty(M, dtype=complex)
    f[0] = -D[0] - D.sum()
    f[1:N] = D[:N - 1] - D[1:]
    f[M - 1] = D[N - 1]
    return f

class Red:
    """reduced Jacobian A_hat(K) = T(K) + us vs^T, us = e_2 + e_M (reduced indices 0, n-1), vs = -A l"""
    def __init__(self, M, K):
        self.M = M; self.K = K; A = 1 + K; self.A = A
        n = M - 1; self.n = n
        lower = np.full(n, A); diag = np.full(n, -(A + 1.0)); upper = np.full(n, 1.0)
        diag[n - 1] = -A          # row M: A v_N - A v_M (+ A v_1, moved to the rank-one part)
        self.lower, self.diag, self.upper = lower, diag, upper
        self.l = np.arange(2, M + 1, dtype=float)
        self.us = np.zeros(n); self.us[0] = 1; self.us[n - 1] = 1
        self.vs = -A*self.l
    def matvec(self, y):
        n = self.n
        r = self.diag*y
        r[1:] += self.lower[1:]*y[:-1]
        r[:-1] += self.upper[:-1]*y[1:]
        return r + self.us*(self.vs@y)
    def rmatvec(self, y):          # A_hat^T y
        n = self.n
        r = self.diag*y
        r[:-1] += self.lower[1:]*y[1:]
        r[1:] += self.upper[:-1]*y[:-1]
        return r + self.vs*(self.us@y)
    def solve(self, sigma, rhs, transpose=False):
        """(A_hat - sigma) x = rhs  (or transposed)"""
        n = self.n
        ab = np.zeros((3, n), dtype=complex)
        if not transpose:
            ab[0, 1:] = self.upper[:-1]; ab[1] = self.diag - sigma; ab[2, :-1] = self.lower[1:]
            u_, v_ = self.us, self.vs
        else:
            ab[0, 1:] = self.lower[1:]; ab[1] = self.diag - sigma; ab[2, :-1] = self.upper[:-1]
            u_, v_ = self.vs, self.us
        z1 = solve_banded((1, 1), ab, rhs.astype(complex))
        z2 = solve_banded((1, 1), ab, u_.astype(complex))
        return z1 - z2*(v_@z1)/(1 + v_@z2)
    def dense(self):
        n = self.n
        D = np.diag(self.diag) + np.diag(self.lower[1:], -1) + np.diag(self.upper[:-1], 1)
        return D + np.outer(self.us, self.vs)

def E(y, M):
    return np.concatenate([[-(np.arange(2, M + 1)@y)], y])

def bil_red(x, y, M):
    return bil_full(E(x, M), E(y, M))[1:]

def eigpair(R, lam0, q0=None, it=40):
    n = R.n
    q = R.solve(lam0, np.ones(n)) if q0 is None else q0
    for _ in range(3):
        q = R.solve(lam0, q); q = q/np.abs(q).max()
    k = int(np.argmax(np.abs(q))); q = q/q[k]
    lam = lam0
    for _ in range(it):
        # bordered Newton: (A - lam) dq - q dlam = -(A - lam) q ; dq_k = 0
        r = R.matvec(q) - lam*q
        z1 = R.solve(lam, -r); z2 = R.solve(lam, q)
        # dq = z1 + dlam z2 ; enforce dq_k = 0
        dlam = -z1[k]/z2[k]
        dq = z1 + dlam*z2
        q = q + dq; lam = lam + dlam
        if abs(dlam) < 1e-15*abs(lam):
            break
    return lam, q

def hopf(M, kap0, om0):
    """secant in K on Re lam(K) = 0"""
    K0 = kap0/math.sqrt(M)
    Ks = [K0*(1 + 1e-4), K0*(1 + 2e-4)]
    vals = []; q = None; lam = 1j*om0
    for K in Ks:
        R = Red(M, K); lam, q = eigpair(R, lam, None); vals.append(lam.real)
    for _ in range(60):
        K = Ks[-1] - vals[-1]*(Ks[-1] - Ks[-2])/(vals[-1] - vals[-2])
        R = Red(M, K); lam, q = eigpair(R, lam, q)
        Ks.append(K); vals.append(lam.real)
        if abs(Ks[-1] - Ks[-2]) < 1e-15*K:
            break
    return K, lam, q, R

def lyap(M, K, lam, q, R):
    om = lam.imag
    n = R.n
    u = R.solve(lam, np.ones(n), transpose=True)
    u = R.solve(lam, u, transpose=True)
    u = u/(u@q)
    qb = q.conj()
    w11 = R.solve(0.0, bil_red(q, qb, M))
    w20 = R.solve(2j*om, -bil_red(q, q, M))           # (A - 2i om) w20 = -B(q,q)
    c = -2*(u@bil_red(q, w11, M)) + u@bil_red(qb, w20, M)
    l1 = c.real/(2*om)
    Eq = E(q, M)
    return l1, l1/np.vdot(Eq, Eq).real, l1/abs(Eq[0])**2, om

def selftest():
    M, K = 9, 0.7
    rng = np.random.default_rng(0)
    X = rng.normal(size=M) + 1j*rng.normal(size=M); Y = rng.normal(size=M) + 1j*rng.normal(size=M)
    A = 1 + K
    pol = field(A + X + Y, K) - field(A + X, K) - field(A + Y, K) + field(A*np.ones(M), K)
    e1 = np.abs(pol - bil_full(X, Y)).max()
    R = Red(M, K)
    # dense reduced Jacobian from the field by exact central differences, compared with T + us vs^T
    B = np.array([(field(A + np.eye(M)[m], K) - field(A - np.eye(M)[m], K))/2 for m in range(M)]).T
    Ah = B[1:, 1:] - np.outer(B[1:, 0], np.arange(2, M + 1))
    e2 = np.abs(Ah - R.dense()).max()
    return e1, e2

def load(dirn, j):
    out = {}
    for fn in ("sweep_j%d.csv.gz" % j, "sweep_j%d_b.csv.gz" % j):
        for r in csv.DictReader(gzip.open(os.path.join(dirn, fn), "rt")):
            if r["status"] in ("PASS", "FAIL"):
                out[int(r["M"])] = r
    return out

def main():
    out, seed, nper = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    extra = [tuple(int(s) for s in x.split(":")) for x in sys.argv[4:]]
    e1, e2 = selftest()
    print("selftest: polarisation vs D^2J assembly %.1e ; reduced Jacobian %.1e" % (e1, e2))
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates")
    rnd = random.Random(seed)
    pts = []
    for j, m0 in ((1, 15), (2, 44), (3, 87)):
        pts += [(M, j) for M in sorted(rnd.sample(range(m0, 100001), nper))]
    pts += extra
    certs = {j: load(d, j) for j in (1, 2, 3)}
    nbad = 0; t0 = time.time(); worst = 0
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["M", "j", "kappa_mine", "kappa_solver", "l1_unit_mine", "l1_unit_solver", "rel_diff",
                     "l1_q1_mine", "l1_q1_solver", "sign_ok"])
        for (M, j) in pts:
            r = certs[j][M]
            kap_s = float(r["kappa_mid"]); l1s = float(r["l1_unit_lo"]); q1s = float(r["l1_q1_lo"])
            K, lam, q, R = hopf(M, kap_s, float(r["omega_lo"]))
            l1, l1u, l1q1, om = lyap(M, K, lam, q, R)
            kap = K*math.sqrt(M)
            rel = abs(l1u - l1s)/abs(l1s)
            worst = max(worst, rel)
            sign_ok = (l1u > 0) == ((M, j) == (15, 1)) and om > 0
            # double precision: Re lambda is O(M^-3/2) and the inner products cancel to relative order M^1.5, so
            # only ~5 digits of l1 survive at M = 1e5; the test is the sign plus a loose agreement
            nbad += (not sign_ok) or rel > 1e-3 or abs(kap - kap_s) > 1e-7*kap_s
            wr.writerow([M, j, "%.15g" % kap, "%.15g" % kap_s, "%.12e" % l1u, "%.12e" % l1s, "%.2e" % rel,
                         "%.12e" % l1q1, "%.12e" % q1s, sign_ok])
    print("%d points (%d random per branch, seed %d, + %d extra): max rel. diff l1 (unit) = %.2e ; %s ; %.0fs"
          % (len(pts), nper, seed, len(extra), worst, "PASS" if nbad == 0 else "FAIL (%d)" % nbad, time.time() - t0))
    sys.exit(1 if nbad else 0)

if __name__ == "__main__":
    main()
