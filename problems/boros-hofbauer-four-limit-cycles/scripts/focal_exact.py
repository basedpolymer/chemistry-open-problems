"""
focal_exact.py -- EXACT focal values (Lyapunov quantities) V~_1..V~_kmax of the Boros-Hofbauer system

    x' = 1/y - x^a y^b,   y' = K(-(d-1) - 1/y + d x^a y^b),   K = a/(1+bd)   (tr J(1,1) = 0),

as polynomials with rational coefficients in the variables (a, c, w), where

    c = 1 + b,   w = omega/a,   omega^2 = det J(1,1) = a^2 (d-1)/(1+bd),
    equivalently  d = (1+w^2)/(1-(c-1)w^2),  K = a(1-(c-1)w^2)/c.

Environment: Python 3.12, python-flint 0.9.0 (fmpq_mpoly). Runtime: ~1-3 min for kmax = 4.
Usage: python focal_exact.py [kmax] [outdir]      (writes outdir/Vtilde_k.txt)

Derivation (see notes.md, Sec. 3).
  u = x-1, v = y-1, and the linear change + time rescaling
      u = c X,   v = a (w Y - X),   tau = omega t
  turn the system into X' = -Y + P, Y' = X + Q with
      P = N_f/(a c w),   Q = (c N_I - N_f)/(a c),
  where I = 1/(1+v), E = (1+u)^a (1+v)^(c-1), f = I - E, and N_I, N_E, N_f = N_I - N_E denote the
  parts of degree >= 2 in (X, Y).  All of N_I, N_f are polynomials in (a, c, w) (coefficientwise).
  With D = a c w we set R~ = D (P + iQ) = N_f + i w (c N_I - N_f).
  Complex Lyapunov-function recursion (z = X + iY, z' = iz + R, R = R~/D):
      F = z zb + sum f_pq z^p zb^q,  F~_n := D^(n-2) F_n,
      S~_n = sum_{m=2}^{n-1} D^(n-1-m) [ F~_m (x) R~_{n-m+1} ]    (bracket = coefficient extraction of
             sum f_rs (r z^(r-1) zb^s R + s z^r zb^(s-1) conj(R)) ),
      f~_pq = i S~_pq/(p-q) (p != q),   V~_k = S~_{k+1,k+1},  f~_{k+1,k+1} = 0.
  Then dF/dtau = sum_k V_k (z zb)^(k+1) with V_k = V~_k / D^(2k) (coordinates X, Y), and the focal values
  in the normalisation of Boros-Hofbauer (u-coordinates, as checked numerically) are
      L_k = pi * V~_k / (a^(2k) c^(4k) w^(2k)).
"""
import sys
import os
import time
import json
from flint import fmpq_mpoly_ctx, fmpq, Ordering

CTX = fmpq_mpoly_ctx.get(('a', 'c', 'w'), Ordering.degrevlex)
A_, C_, W_ = CTX.gens()
ZERO = CTX.from_dict({})
ONE = CTX.from_dict({(0, 0, 0): 1})


class CP:
    """Complex number re + i*im with re, im in Q[a,c,w]."""
    __slots__ = ('re', 'im')

    def __init__(self, re, im=None):
        self.re = re
        self.im = ZERO if im is None else im

    def __add__(self, o):
        return CP(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return CP(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        if isinstance(o, CP):
            return CP(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
        return CP(self.re * o, self.im * o)  # o: fmpq or polynomial (real)

    def conj(self):
        return CP(self.re, -self.im)

    def times_i(self):
        return CP(-self.im, self.re)

    def is_zero(self):
        return self.re == 0 and self.im == 0


def series_mul(Aser, Bser, N):
    """Truncated product of bivariate series given as dicts {(i,j): poly}."""
    out = {}
    for (i1, j1), c1 in Aser.items():
        for (i2, j2), c2 in Bser.items():
            if i1 + j1 + i2 + j2 <= N:
                k = (i1 + i2, j1 + j2)
                out[k] = out.get(k, ZERO) + c1 * c2
    return out


def binom_poly(alpha, i):
    """binomial(alpha, i) for a polynomial alpha."""
    r = ONE
    for k in range(i):
        r = r * (alpha - k) * fmpq(1, k + 1)
    return r


def build_Rtilde(N):
    """R~ = N_f + i w (c N_I - N_f) as dict {(p,q): CP} in the variables z, zb, degrees 2..N."""
    a, c, w = A_, C_, W_
    U = {(1, 0): c}                    # u = c X
    V = {(1, 0): -a, (0, 1): a * w}    # v = a (w Y - X)
    Upow = [{(0, 0): ONE}]
    Vpow = [{(0, 0): ONE}]
    for k in range(1, N + 1):
        Upow.append(series_mul(Upow[-1], U, N))
        Vpow.append(series_mul(Vpow[-1], V, N))
    NI = {}
    NE = {}
    for j in range(2, N + 1):
        sgn = 1 if j % 2 == 0 else -1
        for k, val in Vpow[j].items():
            NI[k] = NI.get(k, ZERO) + val * sgn
    for i in range(0, N + 1):
        bi = binom_poly(a, i)
        for j in range(0, N + 1 - i):
            if i + j < 2:
                continue
            bj = binom_poly(c - 1, j)
            M = series_mul(Upow[i], Vpow[j], N)
            for k, val in M.items():
                NE[k] = NE.get(k, ZERO) + bi * bj * val
    keys = set(NI) | set(NE)
    Nf = {k: NI.get(k, ZERO) - NE.get(k, ZERO) for k in keys}
    Pt = Nf                                                        # real part, in X,Y
    Qt = {k: w * (c * NI.get(k, ZERO) - Nf[k]) for k in keys}       # imaginary part, in X,Y
    # X = (z+zb)/2, Y = -i (z - zb)/2 : expand X^i Y^j with Gaussian-rational coefficients
    zx = {(1, 0): CP(ONE * fmpq(1, 2)), (0, 1): CP(ONE * fmpq(1, 2))}
    zy = {(1, 0): CP(ZERO, ONE * fmpq(-1, 2)), (0, 1): CP(ZERO, ONE * fmpq(1, 2))}

    def cmul(Aser, Bser):
        out = {}
        for (i1, j1), c1 in Aser.items():
            for (i2, j2), c2 in Bser.items():
                if i1 + j1 + i2 + j2 <= N:
                    k = (i1 + i2, j1 + j2)
                    out[k] = out[k] + c1 * c2 if k in out else c1 * c2
        return out
    Xp = [{(0, 0): CP(ONE)}]
    Yp = [{(0, 0): CP(ONE)}]
    for k in range(1, N + 1):
        Xp.append(cmul(Xp[-1], zx))
        Yp.append(cmul(Yp[-1], zy))
    R = {}
    for (i, j) in keys:
        coef = CP(Pt[(i, j)], Qt[(i, j)])
        if coef.is_zero():
            continue
        M = cmul(Xp[i], Yp[j])
        for k, val in M.items():
            # val is a Gaussian rational (constant poly); coef * val
            term = coef * val
            R[k] = R[k] + term if k in R else term
    return R


def lyapunov_exact(kmax, verbose=True):
    N = 2 * kmax + 2
    t0 = time.time()
    R = build_Rtilde(N - 1)
    if verbose:
        print("R~ built (degree <= %d) in %.1fs" % (N - 1, time.time() - t0), flush=True)
    Rb = {(q, p): val.conj() for (p, q), val in R.items()}
    D = A_ * C_ * W_
    Dpow = [ONE]
    for k in range(1, N + 1):
        Dpow.append(Dpow[-1] * D)
    F = {(1, 1): CP(ONE)}
    Fdeg = {2: [(1, 1)]}
    Vt = []
    for n in range(3, N + 1):
        S = {}
        for m in range(2, n):
            deg = n - m + 1
            acc = {}
            for (r, s) in Fdeg.get(m, []):
                fc = F[(r, s)]
                if fc.is_zero():
                    continue
                for (p, q), rc in R.items():
                    if p + q != deg:
                        continue
                    if r > 0:
                        k = (r - 1 + p, s + q)
                        t = fc * rc * fmpq(r)
                        acc[k] = acc[k] + t if k in acc else t
                for (p, q), rc in Rb.items():
                    if p + q != deg:
                        continue
                    if s > 0:
                        k = (r + p, s - 1 + q)
                        t = fc * rc * fmpq(s)
                        acc[k] = acc[k] + t if k in acc else t
            e = n - 1 - m
            for k, val in acc.items():
                val = val * Dpow[e] if e > 0 else val
                S[k] = S[k] + val if k in S else val
        Fdeg[n] = []
        for p in range(n + 1):
            q = n - p
            sv = S.get((p, q), CP(ZERO))
            if p != q:
                F[(p, q)] = sv.times_i() * fmpq(1, p - q)
                Fdeg[n].append((p, q))
            else:
                if sv.im != 0:
                    raise RuntimeError("non-real focal value at degree %d" % n)
                Vt.append(sv.re)
        # consistency: S_qp = conj(S_pq)
        for p in range(n + 1):
            q = n - p
            if (p, q) in S and (q, p) in S:
                d1 = S[(p, q)] - S[(q, p)].conj()
                if not d1.is_zero():
                    raise RuntimeError("symmetry check failed")
        if verbose:
            print("degree %d done (%.1fs)%s" % (n, time.time() - t0,
                  ("  -> V~_%d has %d terms" % (len(Vt), len(Vt[-1]))) if n % 2 == 0 else ""), flush=True)
    return Vt


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                                  "..", "certificates")
    Vt = lyapunov_exact(kmax)
    os.makedirs(outdir, exist_ok=True)
    for k, V in enumerate(Vt, 1):
        save_poly(V, os.path.join(outdir, "Vtilde_%d.json" % k))
        print("V~_%d: %d terms, total degree %d, degrees (a,c,w) = %s" % (k, len(V), V.total_degree(),
                                                                          V.degrees()))


def save_poly(V, path):
    """Save a polynomial in Q[a,c,w] as JSON: {"vars": [...], "terms": [[ea, ec, ew, "p/q"], ...]}."""
    terms = [[int(e[0]), int(e[1]), int(e[2]), str(cf)] for e, cf in V.to_dict().items()]
    terms.sort()
    with open(path, "w") as fh:
        json.dump({"vars": ["a", "c", "w"], "terms": terms}, fh)


def load_poly(path):
    with open(path) as fh:
        data = json.load(fh)
    d = {}
    for t in data["terms"]:
        num, _, den = t[3].partition('/')
        d[(t[0], t[1], t[2])] = fmpq(int(num), int(den) if den else 1)
    return CTX.from_dict(d)


if __name__ == "__main__":
    main()
