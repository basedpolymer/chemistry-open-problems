#!/usr/bin/env python3
"""petnet.py -- Agent C (computation, red team), attainment wave.  Library, imported by the check scripts.

Generalized Petlyuk arrangement (FTC) for N components at FIXED flows, in junction / net-flow variables, in high
precision (mpmath).  Numerical evidence only, never a proof.

Model (class A of total/foundations/class.md): constant relative volatilities alpha_1 > ... > alpha_N, CMO,
equilibrium stages y = alpha x / S(x), S(x) = sum alpha_i x_i.  Topology (BRIEF-common section 2), 0-based indices:
  columns (i, m), 0 <= i < m <= N-1; unit (i, m) = junction (feed) stage of column (i, m); unit (p, p) = product p
  (p = 0: total condenser; p = N-1: partial reboiler; else a side-draw equilibrium stage, liquid product, no duty).
  Section T(i, m): bottom unit (i, m), top unit (i, m-1).   Section B(i, m): bottom unit (i+1, m), top unit (i, m).
  Process feed (z, q) enters unit (0, N-1).

Exact section relation (no approximation): in a section with flows V, L, upward net flows w (sum w = V - L) and
n stages, the liquid x_{k+1} entering stage k+1 from above satisfies x_{k+1} = M x_k / (L S(x_k)),
M = M(w) = V diag(alpha) - w alpha^T, where x_0 = liquid of the bottom unit (the vapor entering the section is in
equilibrium with it) and x_{n+1} = liquid of the top unit (for the condenser: the reflux = distillate).  Hence
x_Y ~ M^{n+1} x_X projectively.  Eigenvalues of M are V phi, phi the roots of
P_w(phi) = V prod_j (1 - phi/alpha_j) - sum_i w_i prod_{j != i} (1 - phi/alpha_j)  ( = prod(1 - phi/alpha) Psi(phi) );
polynomial left eigenvectors lhat(phi)_i = prod_{j != i} (1 - phi/alpha_j) (smooth through zero trace flows).
With the pinch mode k* of the section (descending order of phi) the section relation is EXACTLY equivalent (pivots
lhat_* . x != 0) to the N - 1 chart equations
   r above k*:  E_r = (l_r.x_X)(l_*.x_Y) - (phi_*/phi_r)^(n+1) (l_r.x_Y)(l_*.x_X) = 0
   r below k*:  E_r = (l_r.x_Y)(l_*.x_X) - (phi_r/phi_*)^(n+1) (l_r.x_X)(l_*.x_Y) = 0
and at n = infinity (t = 0) they become the pinch conditions of S's Lemma 11 (closed at the bottom above the
pinch, closed at the top below it).

Unknowns: x_u (N entries) for every unit, w_s (N entries) for every section.  Equations: normalization of every x_u;
sum(w_s) = V_s - L_s; component balances 0..N-2 of every unit; N - 1 chart equations per section.  Square.

Environment: Python 3.12, mpmath 1.3, numpy 2.x (WSL Fedora sci env).
"""
import math
import mpmath as mp

mpf = mp.mpf


# ----------------------------------------------------------------------------------------------- basic algebra
def bracket_root(f, lo, hi, increasing=None, tol=None, maxit=4000):
    """Root of a continuous function with a sign change on (lo, hi), bisection then secant polish."""
    flo, fhi = f(lo), f(hi)
    if flo == 0:
        return lo
    if fhi == 0:
        return hi
    if (flo > 0) == (fhi > 0):
        raise ValueError("no sign change")
    a, b, fa, fb = lo, hi, flo, fhi
    eps = mpf(2) ** (-mp.mp.prec + 4)
    for _ in range(maxit):
        m = (a + b) / 2
        if m <= a or m >= b:
            break
        fm = f(m)
        if fm == 0:
            return m
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b, fb = m, fm
        if abs(b - a) <= eps * max(abs(a), abs(b), mpf(1)):
            break
    return (a + b) / 2


def feed_fun(al, c, U):
    """phi -> sum_p al_p c_p/(al_p - phi) - U."""
    return lambda t: mp.fsum(a * x / (a - t) for a, x in zip(al, c) if x != 0) - U


def interval_roots(al, c, U, idx):
    """Roots of sum al c/(al - phi) = U in (al[r+1], al[r]) for r in idx (c > 0 on r, r+1)."""
    out = {}
    f = feed_fun(al, c, U)
    for r in idx:
        lo, hi = al[r + 1], al[r]
        w = (hi - lo) * mpf(2) ** (-(mp.mp.prec - 12))
        out[r] = bracket_root(f, lo + w, hi - w)
    return out


def lemma1_vector(al_loc, chi):
    """Astra's Lemma 1: unique normalized X on local components (volatilities al_loc, m entries) with
    G_X(chi_j) = 0 for the m - 1 interlacing points chi."""
    m = len(al_loc)
    om = []
    for i in range(m):
        num = mp.fprod([al_loc[i] - c for c in chi])
        den = al_loc[i] * mp.fprod([al_loc[i] - al_loc[k] for k in range(m) if k != i])
        om.append(num / den)
    s = mp.fsum(om)
    return [o / s for o in om]


def Gx(al, x, t):
    return mp.fsum(a * xi / (a - t) for a, xi in zip(al, x))


def poly_P(al, V, w, phi):
    """P_w(phi) = V prod(1 - phi/al) - sum_i w_i prod_{j != i}(1 - phi/al_j)."""
    fac = [1 - phi / a for a in al]
    N = len(al)
    tot = V * mp.fprod(fac)
    for i in range(N):
        if w[i] != 0:
            tot -= w[i] * mp.fprod([fac[j] for j in range(N) if j != i])
    return tot


def poly_P_d(al, V, w, phi):
    """d/dphi of P_w."""
    N = len(al)
    fac = [1 - phi / a for a in al]
    der = [-1 / a for a in al]
    # derivative of prod over set S = sum_k der_k prod_{j in S, j != k} fac_j
    def dprod(S):
        return mp.fsum(der[k] * mp.fprod([fac[j] for j in S if j != k]) for k in S)
    allS = list(range(N))
    tot = V * dprod(allS)
    for i in range(N):
        if w[i] != 0:
            tot -= w[i] * dprod([j for j in allS if j != i])
    return tot


def section_roots(al, V, w, guess=None):
    """All N roots of P_w, real and simple near the configurations used here, sorted descending.
    With a guess (previous roots), Newton-polish each; else polyroots."""
    N = len(al)
    if guess is not None:
        out = []
        for g in guess:
            phi = mpf(g)
            ok = False
            for _ in range(200):
                p = poly_P(al, V, w, phi)
                dp = poly_P_d(al, V, w, phi)
                if dp == 0:
                    break
                st = p / dp
                phi -= st
                if abs(st) <= abs(phi) * mpf(2) ** (-mp.mp.prec + 6):
                    ok = True
                    break
            if not ok:
                return section_roots(al, V, w, None)
            out.append(phi)
        # sanity: distinct
        so = sorted(out, reverse=True)
        for a, b in zip(so, so[1:]):
            if abs(a - b) < mpf(10) ** (-(mp.mp.dps // 3)) * max(abs(a), 1):
                return section_roots(al, V, w, None)
        return so
    # coefficients of P_w in phi via polynomial multiplication (exact in mp)
    def pmul(p, q):
        r = [mpf(0)] * (len(p) + len(q) - 1)
        for i, a in enumerate(p):
            for j, b in enumerate(q):
                r[i + j] += a * b
        return r
    def prod_poly(S):
        p = [mpf(1)]
        for j in S:
            p = pmul(p, [mpf(1), -1 / al[j]])   # 1 - phi/al_j, ascending coefficients
        return p
    tot = [V * c for c in prod_poly(range(N))]
    for i in range(N):
        if w[i] != 0:
            q = prod_poly([j for j in range(N) if j != i])
            for k, c in enumerate(q):
                tot[k] -= w[i] * c
    coeffs = list(reversed(tot))   # descending for polyroots
    while coeffs and coeffs[0] == 0:
        coeffs.pop(0)
    rts = mp.polyroots(coeffs, maxsteps=400, extraprec=4 * mp.mp.prec)
    rr = []
    for r in rts:
        if abs(mp.im(r)) > mpf(10) ** (-(mp.mp.dps // 2)) * max(1, abs(r)):
            raise ValueError("complex root of a section polynomial: %s" % r)
        rr.append(mp.re(r))
    return sorted(rr, reverse=True)


def lhat(al, phi):
    N = len(al)
    fac = [1 - phi / a for a in al]
    return [mp.fprod([fac[j] for j in range(N) if j != i]) for i in range(N)]


def dot(a, b):
    return mp.fsum(x * y for x, y in zip(a, b))


# ----------------------------------------------------------------------------------------------- the problem
class Problem:
    def __init__(self, alpha, z, q):
        self.N = len(alpha)
        self.al = [mpf(a) for a in alpha]
        zz = [mpf(x) for x in z]
        s = mp.fsum(zz)
        self.z = [x / s for x in zz]
        self.q = mpf(q)
        N = self.N
        self.theta = [interval_roots(self.al, self.z, 1 - self.q, [j])[j] for j in range(N - 1)]
        self.H = [mp.fsum(self.al[p] * self.z[p] / (self.theta[j] - self.al[p]) for p in range(j + 1, N))
                  for j in range(N - 1)]
        self.Hmax = max(self.H)
        # topology
        self.cols = [(i, m) for i in range(N) for m in range(i + 1, N)]
        self.units = [(i, m) for i in range(N) for m in range(i, N)]
        self.secs = [('T', i, m) for (i, m) in self.cols] + [('B', i, m) for (i, m) in self.cols]

    @staticmethod
    def sec_ends(s):
        typ, i, m = s
        if typ == 'T':
            return (i, m), (i, m - 1)      # bottom unit, top unit
        return (i + 1, m), (i, m)


# ----------------------------------------------------------------------------------------------- design
def preferred_split(pb, i, m, f, U, chi):
    """Theorem achieve: preferred split of column (i, m) for feed f (support i..m), vapor U, with interval roots
    chi[r] (r = i..m-1).  Returns d (support i..m-1), b, VT, VB."""
    al, N = pb.al, pb.N
    coef = {}
    for p in range(i, m):
        num = mp.fprod([al[p] - chi[r] for r in range(i, m)])
        den = mp.fprod([al[p] - al[pp] for pp in range(i, m) if pp != p])
        coef[p] = num / den / al[p]
    VT = f[i] / coef[i]
    d = [mpf(0)] * N
    for p in range(i, m):
        d[p] = VT * coef[p]
    d[i] = f[i]                      # exact: the lightest component goes up entirely
    b = [f[p] - d[p] for p in range(N)]
    b[i] = mpf(0)
    for p in range(N):
        if p < i or p > m:
            b[p] = mpf(0)
    VB = VT - U
    return d, b, VT, VB


def astra_h(pb, i, m, chi):
    """Astra's redistribution h (eq. (4)) for column (i, m) with middle roots chi[i+1..m-2]."""
    al, N = pb.al, pb.N
    h = [mpf(0)] * N
    mid = [chi[r] for r in range(i + 1, m - 1)]
    beta = mp.fprod([al[p] for p in range(i + 1, m)]) / mp.fprod(mid) if mid else al[i + 1]
    for p in range(i + 1, m):
        num = (al[p] - beta) * mp.fprod([al[p] - c for c in mid])
        den = al[p] * mp.fprod([al[p] - al[k] for k in range(i + 1, m) if k != p])
        h[p] = num / den
    return h, beta


class Design:
    """Operating specification F(delta, eta) and its limit net flows.

    delta  : R = max_j H_j + delta (reboiler duty, the only heat input).
    eta_nb : dict column -> Astra slack eta_[i,m] (non-binary columns); default eta_nb_default.
    eta_bin: dict binary column (p, p+1) -> (kind, eta) with kind 'B' (top total f_p - eta) or 'T' (f_p + eta),
             or kind 'N' (unforced: top total = f_p, the floating-front case).
    """

    def __init__(self, pb, delta, eta_nb=None, eta_bin=None, eta_nb_default=None, bin_default=None):
        self.pb = pb
        N, al = pb.N, pb.al
        self.delta = mpf(delta)
        self.R = pb.Hmax + self.delta
        cols = sorted(pb.cols, key=lambda c: -(c[1] - c[0]))
        self.col = {}
        eta_nb = eta_nb or {}
        eta_bin = eta_bin or {}
        for (i, m) in cols:
            f = [mpf(0)] * N
            U = mpf(0)
            if (i, m) == (0, N - 1):
                f = list(pb.z)
                U = 1 - pb.q
            if m + 1 <= N - 1:
                c = self.col[(i, m + 1)]
                f = [a + b for a, b in zip(f, c['d'])]
                U += c['VT']
            if i - 1 >= 0:
                c = self.col[(i - 1, m)]
                f = [a + b for a, b in zip(f, c['b'])]
                U -= c['VB']
            rec = dict(f=f, U=U)
            if m - i >= 2:
                chi = interval_roots(al, f, U, list(range(i, m)))
                d0, b0, VT0, VB0 = preferred_split(pb, i, m, f, U, chi)
                eta = mpf(eta_nb.get((i, m), eta_nb_default if eta_nb_default is not None else 0))
                h, beta = astra_h(pb, i, m, chi)
                d = [x + eta * y for x, y in zip(d0, h)]
                b = [x - eta * y for x, y in zip(b0, h)]
                rec.update(kind='nb', chi=chi, d0=d0, b0=b0, VT0=VT0, VB0=VB0, eta=eta, h=h, beta=beta,
                           d=d, b=b, VT=VT0 + eta, VB=VB0 + eta)
            else:
                p = i
                chi = interval_roots(al, f, U, [p])
                kind, eta = eta_bin.get((i, m), bin_default if bin_default is not None else ('N', 0))
                eta = mpf(eta)
                d = [mpf(0)] * N
                b = [mpf(0)] * N
                if kind == 'B':
                    d[p] = f[p] - eta
                    b[p] = eta
                    b[p + 1] = f[p + 1]
                elif kind == 'T':
                    d[p] = f[p]
                    d[p + 1] = eta
                    b[p + 1] = f[p + 1] - eta
                else:
                    d[p] = f[p]
                    b[p + 1] = f[p + 1]
                yT = al[p] * f[p] / (al[p] - chi[p])
                yB = al[p + 1] * f[p + 1] / (chi[p] - al[p + 1])
                rec.update(kind=kind, eta=eta, chi=chi, d=d, b=b, yT=yT, yB=yB)
            self.col[(i, m)] = rec
        # vapor of the binary columns from the bottom: V_B(N-2, N-1) = R; side draws: V_T(p,p+1) = V_B(p-1,p)
        Vb = self.R
        for p in range(N - 2, -1, -1):
            c = self.col[(p, p + 1)]
            c['VB'] = Vb
            c['VT'] = Vb + c['U']
            Vb = c['VT']
            c['excess'] = c['VT'] - c['yT']
        # liquids, totals, products
        for (i, m), c in self.col.items():
            c['Dtot'] = mp.fsum(c['d'])
            c['Btot'] = mp.fsum(c['b'])
            c['LT'] = c['VT'] - c['Dtot']
            c['LB'] = c['VB'] + c['Btot']
        self.P = []
        self.xprod = []
        for p in range(N):
            v = [mpf(0)] * N
            if p >= 1:
                v = [a + b for a, b in zip(v, self.col[(p - 1, p)]['b'])]
            if p <= N - 2:
                v = [a + b for a, b in zip(v, self.col[(p, p + 1)]['d'])]
            P = mp.fsum(v)
            self.P.append(P)
            self.xprod.append([x / P for x in v])
        # section flows and limit net flows
        self.V, self.L, self.wlim = {}, {}, {}
        for (i, m), c in self.col.items():
            self.V[('T', i, m)], self.L[('T', i, m)] = c['VT'], c['LT']
            self.V[('B', i, m)], self.L[('B', i, m)] = c['VB'], c['LB']
            self.wlim[('T', i, m)] = list(c['d'])
            self.wlim[('B', i, m)] = [-x for x in c['b']]

    def flow_checks(self):
        """Total liquid and vapor balances at every unit; positivity of flows.  Returns (max abs residual,
        min flow)."""
        pb = self.pb
        N = pb.N
        res = []
        for u in pb.units:
            Lin = Lout = Vin = Vout = mpf(0)
            for s in pb.secs:
                bu, tu = pb.sec_ends(s)
                if bu == u:   # section above u: u sends vapor V up, receives liquid L
                    Vout += self.V[s]
                    Lin += self.L[s]
                if tu == u:   # section below u: u receives vapor, sends liquid
                    Vin += self.V[s]
                    Lout += self.L[s]
            Fl = Fv = mpf(0)
            if u == (0, N - 1):
                Fl, Fv = pb.q, 1 - pb.q
            if u[0] == u[1]:
                p = u[0]
                if p == 0:          # total condenser: no vapor out, product liquid
                    res.append(Vin - (Lout + self.P[0]))
                    continue
                if p == N - 1:      # reboiler: duty R
                    res.append(Lin - (Vout + self.P[p]) + (self.R - 0) * 0)
                    res.append(Vout - self.R)
                    continue
                res.append(Lin + Vin - Lout - Vout - self.P[p])
                res.append(Vout - Vin)          # no duty at side draws
                continue
            res.append(Lin + Vin + Fl + Fv - Lout - Vout)
            res.append(Vout - Vin - Fv)         # adiabatic junction
        flows = list(self.V.values()) + list(self.L.values()) + list(self.P)
        return max(abs(r) for r in res), min(flows)


# ----------------------------------------------------------------------------------------------- limit config
def pinch_index(pb, s, kind):
    """Pinch mode index (0-based, descending phi) of section s at the limit configuration."""
    typ, i, m = s
    if kind == 'nb':
        return (m - 1) if typ == 'T' else (i + 1)
    if kind == 'B':           # binary (p, p+1) B-pinned: T pinch at its major root (index p), B at phi_1 (index p)
        return i
    if kind == 'T':           # T-pinned: T pinch at phi_2 just below alpha_{p+1} (index p+1), B at its root (p+1)
        return i + 1
    raise ValueError("unforced binary column: the pinch is ambiguous (floating front)")


def limit_config(des):
    """Independent infinite-stage limit configuration: junction compositions (Lemma 1 vectors, pinch eigenvectors),
    product compositions, section net flows, section roots and pinch indices."""
    pb = des.pb
    N, al = pb.N, pb.al
    X = {}
    info = {}
    roots = {}
    for s in pb.secs:
        roots[s] = section_roots(al, des.V[s], des.wlim[s])
    for (i, m), c in des.col.items():
        if c['kind'] == 'nb':
            rT = roots[('T', i, m)]
            rB = roots[('B', i, m)]
            # r = root of Psi_T in (al_{i+1}, al_i); s = root of Psi_B in (al_m, al_{m-1}); middle roots
            r = [x for x in rT if al[i + 1] < x < al[i]]
            sb = [x for x in rB if al[m] < x < al[m - 1]]
            assert len(r) == 1 and len(sb) == 1, (i, m, rT, rB)
            r, sb = r[0], sb[0]
            mids_T = [[x for x in rT if al[j + 1] < x < al[j]] for j in range(i + 1, m - 1)]
            mids_B = [[x for x in rB if al[j + 1] < x < al[j]] for j in range(i + 1, m - 1)]
            midgap = max([abs(a[0] - b[0]) for a, b in zip(mids_T, mids_B)] + [mpf(0)])
            chi = [r] + [a[0] for a in mids_T] + [sb]
            loc = lemma1_vector(al[i:m + 1], chi)
            x = [mpf(0)] * N
            for k, v in enumerate(loc):
                x[i + k] = v
            X[(i, m)] = x
            info[(i, m)] = dict(r=r, s=sb, midgap=midgap)
        else:
            p = i
            if c['kind'] == 'B':
                sec = ('B', i, m)
                k = pinch_index(pb, sec, 'B')
            elif c['kind'] == 'T':
                sec = ('T', i, m)
                k = pinch_index(pb, sec, 'T')
            else:
                X[(i, m)] = None
                continue
            phi = roots[sec][k]
            w = des.wlim[sec]
            e = [mpf(0)] * N
            for j in (p, p + 1):
                e[j] = w[j] / (al[j] - phi)
            se = mp.fsum(e)
            X[(i, m)] = [v / se for v in e]
            info[(i, m)] = dict(phi=phi, sec=sec)
    for p in range(N):
        X[(p, p)] = list(des.xprod[p])
    kstar = {}
    for s in pb.secs:
        typ, i, m = s
        kstar[s] = pinch_index(pb, s, des.col[(i, m)]['kind'])
    return dict(X=X, info=info, roots=roots, kstar=kstar, w=des.wlim)


# ----------------------------------------------------------------------------------------------- the system
class System:
    """Finite-stage (n per section, or dict) or limit (n = None) system in junction/net-flow variables."""

    def __init__(self, des, lim, nst=None, kstar=None, drop=None):
        self.des = des
        self.pb = des.pb
        self.lim = lim
        pb = self.pb
        N = pb.N
        self.N = N
        self.uidx = {u: k for k, u in enumerate(pb.units)}
        self.sidx = {s: k for k, s in enumerate(pb.secs)}
        self.nU, self.nS = len(pb.units), len(pb.secs)
        self.nvar = N * (self.nU + self.nS)
        self.nst = nst
        self.kstar = dict(lim['kstar'])
        if kstar:
            self.kstar.update(kstar)
        # dropped (redundant) component balance per unit: a major component (product p: p; junction (i,m): i)
        self.drop = {u: u[0] for u in pb.units}
        if drop:
            self.drop.update(drop)
        self.root_guess = {s: list(lim['roots'][s]) for s in pb.secs}
        # constant row scalings of the chart equations (from the limit configuration)
        self.scale = {}
        for s in pb.secs:
            rts = lim['roots'][s]
            ls = [lhat(pb.al, ph) for ph in rts]
            nrm = [max(abs(v) for v in l) for l in ls]
            k = self.kstar[s]
            self.scale[s] = [1 / (nrm[r] * nrm[k]) for r in range(N)]

    def n_of(self, s):
        if self.nst is None:
            return None
        if isinstance(self.nst, dict):
            return self.nst[s]
        return self.nst

    # variables
    def pack(self, X, W):
        v = []
        for u in self.pb.units:
            v += list(X[u])
        for s in self.pb.secs:
            v += list(W[s])
        return mp.matrix(v)

    def unpack(self, v):
        N = self.N
        X, W = {}, {}
        for k, u in enumerate(self.pb.units):
            X[u] = [v[N * k + j] for j in range(N)]
        off = N * self.nU
        for k, s in enumerate(self.pb.secs):
            W[s] = [v[off + N * k + j] for j in range(N)]
        return X, W

    def start(self):
        X = {u: list(self.lim['X'][u]) for u in self.pb.units}
        W = {s: list(self.lim['w'][s]) for s in self.pb.secs}
        return self.pack(X, W)

    def section_data(self, s, W, update_guess=False):
        pb = self.pb
        rts = section_roots(pb.al, self.des.V[s], W[s], self.root_guess[s])
        if update_guess:
            self.root_guess[s] = rts
        return rts

    def residual(self, v, update_guess=False, tvals=None, allrows=False):
        pb, N, des = self.pb, self.N, self.des
        X, W = self.unpack(v)
        F = []
        for u in pb.units:
            F.append(mp.fsum(X[u]) - 1)
        for s in pb.secs:
            F.append(mp.fsum(W[s]) - (des.V[s] - des.L[s]))
        # unit balances, components 0..N-2
        for u in pb.units:
            i, m = u
            bal = [mpf(0)] * N
            if i < m:
                for j in range(N):
                    bal[j] = W[('T', i, m)][j] - W[('B', i, m)][j]
                    if m + 1 <= N - 1:
                        bal[j] -= W[('T', i, m + 1)][j]
                    if i - 1 >= 0:
                        bal[j] += W[('B', i - 1, m)][j]
                    if u == (0, N - 1):
                        bal[j] -= pb.z[j]
            else:
                p = i
                for j in range(N):
                    if p + 1 <= N - 1:
                        bal[j] += W[('T', p, p + 1)][j]
                    if p - 1 >= 0:
                        bal[j] -= W[('B', p - 1, p)][j]
                    bal[j] -= des.P[p] * X[u][j]
            if allrows:
                F += bal
            else:
                F += [bal[j] for j in range(N) if j != self.drop[u]]
        # chart equations
        for s in pb.secs:
            bu, tu = pb.sec_ends(s)
            xX, xY = X[bu], X[tu]
            rts = self.section_data(s, W, update_guess)
            k = self.kstar[s]
            ls = [lhat(pb.al, ph) for ph in rts]
            aX = [dot(l, xX) for l in ls]
            aY = [dot(l, xY) for l in ls]
            n = self.n_of(s)
            for r in range(N):
                if r == k:
                    continue
                if n is None:
                    t = mpf(0)
                else:
                    t = (rts[k] / rts[r]) ** (n + 1) if r < k else (rts[r] / rts[k]) ** (n + 1)
                if tvals is not None:
                    tvals[(s, r)] = t
                if r < k:
                    E = aX[r] * aY[k] - t * aY[r] * aX[k]
                else:
                    E = aY[r] * aX[k] - t * aX[r] * aY[k]
                F.append(E * self.scale[s][r])
        return mp.matrix(F)

    def jacobian(self, v, h=None, allrows=False, fd=False):
        if not fd:
            return self.jacobian_an(v, allrows=allrows)
        if h is None:
            h = mpf(10) ** (-(mp.mp.dps // 3))
        n = len(v)
        J = None
        for j in range(n):
            vp = v.copy()
            vm = v.copy()
            vp[j] += h
            vm[j] -= h
            Fp = self.residual(vp, allrows=allrows)
            Fm = self.residual(vm, allrows=allrows)
            if J is None:
                J = mp.matrix(len(Fp), n)
            for i in range(len(Fp)):
                J[i, j] = (Fp[i] - Fm[i]) / (2 * h)
        return J

    def jacobian_an(self, v, allrows=False):
        """Analytic Jacobian of residual() (same row order)."""
        pb, N, des = self.pb, self.N, self.des
        al = pb.al
        X, W = self.unpack(v)
        nU, nS = self.nU, self.nS
        nvar = N * (nU + nS)
        cx = lambda u, j: N * self.uidx[u] + j
        cw = lambda s, j: N * nU + N * self.sidx[s] + j
        rows = []
        for u in pb.units:
            rows.append({cx(u, j): mpf(1) for j in range(N)})
        for s in pb.secs:
            rows.append({cw(s, j): mpf(1) for j in range(N)})
        for u in pb.units:
            i, m = u
            comps = range(N) if allrows else [j for j in range(N) if j != self.drop[u]]
            for j in comps:
                r = {}
                if i < m:
                    r[cw(('T', i, m), j)] = mpf(1)
                    r[cw(('B', i, m), j)] = mpf(-1)
                    if m + 1 <= N - 1:
                        r[cw(('T', i, m + 1), j)] = mpf(-1)
                    if i - 1 >= 0:
                        r[cw(('B', i - 1, m), j)] = mpf(1)
                else:
                    p = i
                    if p + 1 <= N - 1:
                        r[cw(('T', p, p + 1), j)] = mpf(1)
                    if p - 1 >= 0:
                        r[cw(('B', p - 1, p), j)] = mpf(-1)
                    r[cx(u, j)] = -des.P[p]
                rows.append(r)
        for s in pb.secs:
            bu, tu = pb.sec_ends(s)
            xX, xY = X[bu], X[tu]
            V = des.V[s]
            w = W[s]
            rts = self.section_data(s, W)
            k = self.kstar[s]
            n = self.n_of(s)
            # root derivatives: dphi/dw_i = -(dP/dw_i)/(dP/dphi), dP/dw_i = -prod_{j != i}(1 - phi/al_j)
            dphi = []
            ls, dls = [], []
            for ph in rts:
                fac = [1 - ph / a for a in al]
                dP = poly_P_d(al, V, w, ph)
                dphi.append([mp.fprod([fac[jj] for jj in range(N) if jj != ii]) / dP for ii in range(N)])
                ls.append([mp.fprod([fac[jj] for jj in range(N) if jj != ii]) for ii in range(N)])
                # d lhat_i / d phi = sum_{j != i} (-1/al_j) prod_{l != i, j} fac_l
                dls.append([mp.fsum((-1 / al[jj]) * mp.fprod([fac[ll] for ll in range(N) if ll != ii and ll != jj])
                                    for jj in range(N) if jj != ii) for ii in range(N)])
            A = [dot(l, xX) for l in ls]
            B = [dot(l, xY) for l in ls]
            dA = [dot(dl, xX) for dl in dls]     # d(l_r . x_X)/d phi_r
            dB = [dot(dl, xY) for dl in dls]
            for r in range(N):
                if r == k:
                    continue
                sc = self.scale[s][r]
                row = {}
                def add(c, val):
                    row[c] = row.get(c, mpf(0)) + sc * val
                if n is None:
                    t = mpf(0)
                    dlogt = [mpf(0)] * N
                elif r < k:
                    t = (rts[k] / rts[r]) ** (n + 1)
                    dlogt = [(n + 1) * (dphi[k][ii] / rts[k] - dphi[r][ii] / rts[r]) for ii in range(N)]
                else:
                    t = (rts[r] / rts[k]) ** (n + 1)
                    dlogt = [(n + 1) * (dphi[r][ii] / rts[r] - dphi[k][ii] / rts[k]) for ii in range(N)]
                if r < k:      # E = A_r B_k - t B_r A_k
                    for jj in range(N):
                        add(cx(bu, jj), ls[r][jj] * B[k] - t * B[r] * ls[k][jj])
                        add(cx(tu, jj), A[r] * ls[k][jj] - t * A[k] * ls[r][jj])
                    for ii in range(N):
                        dAr, dBr = dA[r] * dphi[r][ii], dB[r] * dphi[r][ii]
                        dAk, dBk = dA[k] * dphi[k][ii], dB[k] * dphi[k][ii]
                        val = dAr * B[k] + A[r] * dBk - t * dlogt[ii] * B[r] * A[k] - t * (dBr * A[k] + B[r] * dAk)
                        add(cw(s, ii), val)
                else:          # E = B_r A_k - t A_r B_k
                    for jj in range(N):
                        add(cx(bu, jj), B[r] * ls[k][jj] - t * B[k] * ls[r][jj])
                        add(cx(tu, jj), A[k] * ls[r][jj] - t * A[r] * ls[k][jj])
                    for ii in range(N):
                        dAr, dBr = dA[r] * dphi[r][ii], dB[r] * dphi[r][ii]
                        dAk, dBk = dA[k] * dphi[k][ii], dB[k] * dphi[k][ii]
                        val = dBr * A[k] + B[r] * dAk - t * dlogt[ii] * A[r] * B[k] - t * (dAr * B[k] + A[r] * dBk)
                        add(cw(s, ii), val)
                rows.append(row)
        J = mp.matrix(len(rows), nvar)
        for a, row in enumerate(rows):
            for c, val in row.items():
                J[a, c] = val
        return J

    def labels(self, allrows=False):
        """Row and column labels."""
        pb, N = self.pb, self.N
        cols = [('x', u, j) for u in pb.units for j in range(N)] + [('w', s, j) for s in pb.secs for j in range(N)]
        rows = [('norm', u) for u in pb.units] + [('sum', s) for s in pb.secs]
        for u in pb.units:
            rows += [('bal', u, j) for j in range(N) if allrows or j != self.drop[u]]
        for s in pb.secs:
            rows += [('chart', s, r) for r in range(N) if r != self.kstar[s]]
        return rows, cols

    def newton(self, v0, tol=None, maxit=60, verbose=False, jac_every=1):
        if tol is None:
            tol = mpf(10) ** (-(mp.mp.dps - 8))
        v = v0.copy()
        F = self.residual(v, update_guess=True)
        nF = mp.norm(F, mp.inf)
        J = None
        hist = []
        for it in range(maxit):
            if nF < tol:
                break
            if J is None or it % jac_every == 0:
                J = self.jacobian(v)
            dv = solve_linear(J, -F)
            lam = mpf(1)
            ok = False
            while lam > mpf(10) ** -12:
                vn = v + lam * dv
                try:
                    Fn = self.residual(vn)
                except Exception:
                    lam /= 2
                    continue
                nFn = mp.norm(Fn, mp.inf)
                if nFn < nF * (1 - lam / 10) or nFn < tol:
                    ok = True
                    break
                lam /= 2
            if not ok:
                hist.append(('stall', float(nF)))
                break
            v = vn
            F = self.residual(v, update_guess=True)
            nF = mp.norm(F, mp.inf)
            hist.append((float(lam), mp.nstr(nF, 5)))
            if verbose:
                print("    newton it %d lam %.3g |F| %s" % (it, float(lam), mp.nstr(nF, 5)), flush=True)
        return v, nF, hist


def solve_linear(A, b):
    """Solve A x = b; mpmath LU for small systems, python-flint arb_mat for larger ones."""
    n = A.rows
    if n <= 60:
        return mp.lu_solve(A, b)
    import flint
    flint.ctx.prec = mp.mp.prec + 20
    Af = flint.arb_mat(n, n, [flint.arb(mp.nstr(A[i, j], mp.mp.dps + 5)) for i in range(n) for j in range(n)])
    bf = flint.arb_mat(n, 1, [flint.arb(mp.nstr(b[i], mp.mp.dps + 5)) for i in range(n)])
    x = Af.solve(bf)
    return mp.matrix([mpf(x[i, 0].mid().str(mp.mp.dps + 5, radius=False)) for i in range(n)])


# ----------------------------------------------------------------------------------------------- profiles
def right_vec(al, w, phi, V=None):
    """Right eigenvector of M(w) for the root phi.  Near a trace mode (phi close to alpha_k) uses the exact regularized
    form e_i = w_i/(alpha_i - phi) (i != k), e_k = Psi^(k)(phi)/alpha_k (L's Lemma 2.2), which avoids 0/0."""
    N = len(al)
    k = min(range(N), key=lambda i: abs(al[i] - phi))
    if V is not None and abs(al[k] - phi) < mpf(10) ** (-(mp.mp.dps // 4)) * al[k]:
        Psik = V - mp.fsum(al[i] * w[i] / (al[i] - phi) for i in range(N) if i != k)
        e = [w[i] / (al[i] - phi) if i != k else Psik / al[k] for i in range(N)]
    else:
        e = [wi / (a - phi) for a, wi in zip(al, w)]
    mx = max(abs(x) for x in e)
    return [x / mx for x in e]


def reconstruct(sysm, v):
    """Stage profiles of every section from the chart solution (stable two-sided eigen expansion).
    Returns dict s -> list of x_k, k = 0..n+1 (x_0 = bottom unit liquid, x_{n+1} = top unit liquid)."""
    pb, N = sysm.pb, sysm.N
    X, W = sysm.unpack(v)
    prof = {}
    for s in pb.secs:
        bu, tu = pb.sec_ends(s)
        n = sysm.n_of(s)
        rts = sysm.section_data(s, W)
        k = sysm.kstar[s]
        ls = [lhat(pb.al, ph) for ph in rts]
        es = [right_vec(pb.al, W[s], ph, sysm.des.V[s]) for ph in rts]
        le = [dot(l, e) for l, e in zip(ls, es)]
        a = [dot(l, X[bu]) / c for l, c in zip(ls, le)]
        bcoef = [dot(l, X[tu]) / c for l, c in zip(ls, le)]
        out = []
        for st in range(n + 2):
            x = [mpf(0)] * N
            for r in range(N):
                if r >= k:
                    c = (a[r] / a[k]) * (rts[r] / rts[k]) ** st
                else:
                    c = (bcoef[r] / bcoef[k]) * (rts[k] / rts[r]) ** (n + 1 - st)
                for j in range(N):
                    x[j] += c * es[r][j]
            sm = mp.fsum(x)
            out.append([xx / sm for xx in x])
        # exact end values
        out[0] = list(X[bu])
        out[-1] = list(X[tu])
        prof[s] = out
    return prof


def stage_check(sysm, v, prof):
    """Every cross-section balance V y_k - L x_{k+1} = w (relative residual per component), every stage composition
    and S(x_k) > 0, and every unit balance written with the actual stream compositions.  Returns dict."""
    pb, N, des = sysm.pb, sysm.N, sysm.des
    al = pb.al
    X, W = sysm.unpack(v)
    worst_rel = mpf(0)
    minx = mpf(1)
    minS = mpf(10) ** 9
    negs = []
    for s in pb.secs:
        V, L = des.V[s], des.L[s]
        xs = prof[s]
        for st in range(len(xs) - 1):
            x, xn = xs[st], xs[st + 1]
            S = dot(al, x)
            minS = min(minS, S)
            y = [a * xi / S for a, xi in zip(al, x)]
            for j in range(N):
                lhs = V * y[j] - L * xn[j]
                den = V * abs(y[j]) + L * abs(xn[j]) + abs(W[s][j])
                rel = abs(lhs - W[s][j]) / den if den > 0 else mpf(0)
                worst_rel = max(worst_rel, rel)
        for st in range(1, len(xs) - 1):
            for j in range(N):
                if xs[st][j] < minx:
                    minx = xs[st][j]
                if xs[st][j] <= 0:
                    negs.append((s, st, j, xs[st][j]))
    # unit balances with actual streams
    unit_rel = mpf(0)
    for u in pb.units:
        i, m = u
        inflow = [mpf(0)] * N
        outflow = [mpf(0)] * N
        xu = X[u]
        Su = dot(al, xu)
        yu = [a * xi / Su for a, xi in zip(al, xu)]
        for s in pb.secs:
            bu, tu = pb.sec_ends(s)
            V, L = des.V[s], des.L[s]
            xs = prof[s]
            n = len(xs) - 2
            if bu == u:   # section above u: liquid from its stage 1 in, vapor of u out
                if n >= 1:
                    xl = xs[1]
                else:
                    xl = xs[1]
                for j in range(N):
                    inflow[j] += L * xl[j]
                    outflow[j] += V * yu[j]
            if tu == u:   # section below u: vapor from its stage n in, liquid of u out
                xv = xs[n]
                Sv = dot(al, xv)
                for j in range(N):
                    inflow[j] += V * al[j] * xv[j] / Sv
                    if u == (0, 0):
                        pass
                    outflow[j] += L * xu[j]
        if u == (0, N - 1):
            for j in range(N):
                inflow[j] += pb.z[j]
        if i == m:
            p = i
            for j in range(N):
                if p == 0:
                    # total condenser: vapor in = (reflux + distillate) * x_D; reflux counted above
                    outflow[j] += des.P[p] * xu[j]
                else:
                    outflow[j] += des.P[p] * xu[j]
        if u == (0, 0):
            # the condenser has no vapor outflow; the balance above already uses L x_D (reflux) + D x_D
            pass
        for j in range(N):
            den = abs(inflow[j]) + abs(outflow[j])
            if den > 0:
                unit_rel = max(unit_rel, abs(inflow[j] - outflow[j]) / den)
    # unit compositions positive
    for u in pb.units:
        for j in range(N):
            minx = min(minx, X[u][j])
            if X[u][j] <= 0:
                negs.append((u, None, j, X[u][j]))
        minS = min(minS, dot(al, X[u]))
    return dict(worst_rel=worst_rel, unit_rel=unit_rel, minx=minx, minS=minS, negs=negs)


def impurities(sysm, v):
    X, W = sysm.unpack(v)
    out = []
    for p in range(sysm.N):
        out.append(mp.fsum(X[(p, p)][j] for j in range(sysm.N) if j != p))
    return out


# ----------------------------------------------------------------------------------------------- analysis tools
def to_numpy(J):
    import numpy as np
    return np.array([[float(J[i, j]) for j in range(J.cols)] for i in range(J.rows)])


def svals(J):
    """Singular values (descending).  mpmath SVD (exact to working precision) for small matrices, else numpy."""
    import numpy as np
    if J.rows <= 80 and J.rows == J.cols:
        S = mp.svd_r(J, compute_uv=False)
        return sorted([abs(S[i]) for i in range(len(S))], reverse=True)
    return [mpf(x) for x in np.linalg.svd(to_numpy(J), compute_uv=False)]


def classify(sysm, vlim, J=None, allrows=False, tol=None):
    """Major / trace classification at the limit point: a column is trace iff its limit value is 0; a row is trace
    iff all its derivatives with respect to major columns vanish (|.| < tol)."""
    if tol is None:
        tol = mpf(10) ** (-(mp.mp.dps // 4))
    rows, cols = sysm.labels(allrows=allrows)
    if J is None:
        J = sysm.jacobian(vlim, allrows=allrows)
    mc = [j for j in range(len(cols)) if abs(vlim[j]) > tol]
    tc = [j for j in range(len(cols)) if abs(vlim[j]) <= tol]
    tr, mr = [], []
    for i in range(len(rows)):
        if max([abs(J[i, j]) for j in mc] + [mpf(0)]) < tol:
            tr.append(i)
        else:
            mr.append(i)
    return dict(J=J, rows=rows, cols=cols, mc=mc, tc=tc, mr=mr, tr=tr)


def submatrix(J, R, C):
    M = mp.matrix(len(R), len(C))
    for a, i in enumerate(R):
        for b, j in enumerate(C):
            M[a, b] = J[i, j]
    return M


def trace_components(cl, sysm, tol=None):
    """Split the trace block by component: for each component k, the trace columns of k and the trace rows whose
    nonzero entries lie on trace columns of k only.  Rows coupling several components are reported."""
    if tol is None:
        tol = mpf(10) ** (-(mp.mp.dps // 4))
    J, cols, rows = cl['J'], cl['cols'], cl['rows']
    comp_of_col = {j: cols[j][2] for j in cl['tc']}
    blocks = {}
    mixed = []
    for i in cl['tr']:
        ks = set(comp_of_col[j] for j in cl['tc'] if abs(J[i, j]) > tol)
        if len(ks) == 1:
            blocks.setdefault(ks.pop(), {'rows': [], 'cols': []})['rows'].append(i)
        else:
            mixed.append((rows[i], sorted(ks)))
    for j in cl['tc']:
        blocks.setdefault(comp_of_col[j], {'rows': [], 'cols': []})['cols'].append(j)
    return blocks, mixed


def kappa_check(sysm, v):
    """(H7)-type signs.  For every section s and every component k with zero LIMIT net flow in s: the end J at which
    the mode of k is closed (top end if the mode is below the pinch, bottom end if above), G^(k)_{x_J}(alpha_k),
    Psi^(k)_s(alpha_k), kappa = |G| / Psi, and whether the sign rule holds (Psi > 0, G > 0 at a top end, G < 0 at
    a bottom end), so that the limit closed-end relation x_{J,k} = -w_k G/Psi reads x_{J,k} = kappa q with q the net
    flow of k from J into s.  Also reports whether k is a trace at the OTHER end (if not, the finite-n value of
    x_{J,k} carries a source term of order t and the ratio x/(kappa q) need not tend to 1)."""
    pb, N, des = sysm.pb, sysm.N, sysm.des
    al = pb.al
    X, W = sysm.unpack(v)
    Xl = sysm.lim['X']
    out = []
    for s in pb.secs:
        bu, tu = pb.sec_ends(s)
        rts = sysm.section_data(s, W)
        k_star = sysm.kstar[s]
        wl = sysm.lim['w'][s]
        for k in range(N):
            if wl[k] != 0:
                continue
            idx = min(range(N), key=lambda r: abs(rts[r] - al[k]))
            if idx == k_star:
                out.append(dict(sec=s, k=k, end=None, note='trace mode is the pinch'))
                continue
            end = tu if idx > k_star else bu
            other = bu if end == tu else tu
            xJ = X[end]
            G = mp.fsum(al[i] * xJ[i] / (al[i] - al[k]) for i in range(N) if i != k)
            Psi = des.V[s] - mp.fsum(al[i] * W[s][i] / (al[i] - al[k]) for i in range(N) if i != k)
            q = -W[s][k] if end == tu else W[s][k]
            kap = abs(G) / Psi
            out.append(dict(sec=s, k=k, end=end, x=xJ[k], q=q, G=G, Psi=Psi, kappa=kap,
                            other_trace=(Xl[other][k] == 0),
                            ratio=(xJ[k] / (kap * q)) if q != 0 else None,
                            sign_ok=bool(Psi > 0 and ((end == tu and G > 0) or (end == bu and G < 0)))))
    return out


def right_vec_lim(al, V, w, phi, tiny=None):
    """Right eigenvector of M(w) for the mode phi, robust at a trace mode (w_k = 0, phi = alpha_k):
    e_k = 1, e_i = w_i alpha_k / ((alpha_i - alpha_k) Psi^(k)(alpha_k)) (i != k)."""
    N = len(al)
    if tiny is None:
        tiny = mpf(10) ** (-(mp.mp.dps // 2))
    for k in range(N):
        if abs(w[k]) < tiny and abs(phi - al[k]) < tiny * max(1, abs(al[k])):
            Psik = V - mp.fsum(al[i] * w[i] / (al[i] - al[k]) for i in range(N) if i != k)
            e = [w[i] * al[k] / ((al[i] - al[k]) * Psik) if i != k else mpf(1) for i in range(N)]
            return e
    return [wi / (a - phi) for a, wi in zip(al, w)]


def hyp_checks(des, lim, K=600):
    """Numerical tests of (H5)-(H8)/(A2)-(A5) at the limit configuration.
    Returns (ok, messages, data).  (A5)/(H6): limit half-profiles xi^bot_k, xi^top_k (A's formula, modes on the
    closed side dropped) componentwise >= 0 for k = 0..K and S > 0; (H8)/(A3): pivots pi_bot, pi_top > 0 with
    1.e_* = 1 = l_*.e_*; (H7)(iii)-(v)."""
    pb = des.pb
    N, al = pb.N, pb.al
    msgs = []
    ok = True
    data = {}
    tol = mpf(10) ** (-(mp.mp.dps - 15))
    Xl, Wl, R, KS = lim['X'], lim['w'], lim['roots'], lim['kstar']
    minhp = mpf(1)
    minpiv = mpf(10) ** 6
    for s in pb.secs:
        bu, tu = pb.sec_ends(s)
        V = des.V[s]
        w = Wl[s]
        rts = R[s]
        k = KS[s]
        es = [right_vec_lim(al, V, w, ph) for ph in rts]
        ls = [lhat(al, ph) for ph in rts]
        # normalize: 1.e_* = 1, l_r.e_r = 1
        se = mp.fsum(es[k])
        if se == 0:
            ok = False
            msgs.append(f"{s}: 1.e_* = 0")
            continue
        es = [[x / mp.fsum(e) if mp.fsum(e) != 0 else x for x in e] if r == k else e for r, e in enumerate(es)]
        ls = [[x / dot(l, e) for x in l] for l, e in zip(ls, es)]
        cb = [dot(l, Xl[bu]) for l in ls]
        ct = [dot(l, Xl[tu]) for l in ls]
        pib, pit = cb[k], ct[k]
        minpiv = min(minpiv, pib, pit)
        if not (pib > 0 and pit > 0):
            ok = False
            msgs.append(f"(H8)/(A3) FAIL {s}: pivots pi_bot = {mp.nstr(pib, 5)}, pi_top = {mp.nstr(pit, 5)}")
        # closed-side coefficients must vanish (A4)
        cl_b = max([abs(cb[r]) for r in range(k)] + [mpf(0)])
        cl_t = max([abs(ct[r]) for r in range(k + 1, N)] + [mpf(0)])
        if cl_b > tol or cl_t > tol:
            ok = False
            msgs.append(f"(A4) FAIL {s}: closed coefficients {mp.nstr(cl_b, 3)}, {mp.nstr(cl_t, 3)}")
        mn = mpf(1)
        for kk in range(K + 1):
            xb = [mpf(0)] * N
            xt = [mpf(0)] * N
            for r in range(N):
                if r >= k:
                    c = cb[r] / pib * (rts[r] / rts[k]) ** kk
                    for j in range(N):
                        xb[j] += c * es[r][j]
                if r <= k:
                    c = ct[r] / pit * (rts[k] / rts[r]) ** kk
                    for j in range(N):
                        xt[j] += c * es[r][j]
            for vec in (xb, xt):
                sm = mp.fsum(vec)
                if sm <= 0:
                    mn = min(mn, mpf(-1))
                    continue
                vv = [x / sm for x in vec]
                mn = min(mn, min(vv))
                if dot(al, vv) <= 0:
                    mn = min(mn, mpf(-1))
        minhp = min(minhp, mn)
        if mn < -tol:
            ok = False
            msgs.append(f"(H6)/(A5) FAIL {s}: limit half-profile has a negative entry {mp.nstr(mn, 5)}")
        data[str(s)] = dict(pi_bot=float(pib), pi_top=float(pit), min_halfprofile=float(mn))
    msgs.append(f"(H6)/(A5) min over all sections and k <= {K} of the half-profile entries: {mp.nstr(minhp, 4)} (>= 0 required; exact zeros are trace entries)")
    msgs.append(f"(H8)/(A3) smallest pivot (normalized 1.e_* = 1 = l_*.e_*): {mp.nstr(minpiv, 5)}")
    # (H7) (iii)-(v)
    for s in pb.secs:
        bu, tu = pb.sec_ends(s)
        w = Wl[s]
        for j in range(N):
            if w[j] != 0 and (Xl[bu][j] == 0 or Xl[tu][j] == 0):
                ok = False
                msgs.append(f"(H7)(v) FAIL {s}: component {j} has nonzero net flow but is absent at an end")
        k = KS[s]
        for j in range(N):
            if w[j] == 0:
                idx = min(range(N), key=lambda r: abs(R[s][r] - al[j]))
                if idx == k:
                    ok = False
                    msgs.append(f"(H5) FAIL {s}: trace mode {j} is the pinch")
                    continue
                end = tu if idx > k else bu
                if Xl[end][j] != 0:
                    ok = False
                    msgs.append(f"(H7)(iii) FAIL {s}: trace {j} closed at {end} where it is not a trace")
    for u in pb.units:
        if u[0] == u[1]:
            continue
        for j in range(N):
            if Xl[u][j] != 0:
                continue
            closed = False
            for s in pb.secs:
                bu, tu = pb.sec_ends(s)
                if u not in (bu, tu) or Wl[s][j] != 0:
                    continue
                idx = min(range(N), key=lambda r: abs(R[s][r] - al[j]))
                end = tu if idx > KS[s] else bu
                if end == u:
                    closed = True
            if not closed:
                ok = False
                msgs.append(f"(H7)(iv) FAIL: trace {j} at junction {u} is not closed by any section")
    return ok, msgs, dict(sections=data, min_halfprofile=float(minhp), min_pivot=float(minpiv))


def solve_audited(build, nst, dps, start=None, max_dps=4000, newton_kw=None):
    """Newton for G_n with precision audits: the working precision is raised until (i) every nonzero variable and
    (ii) every nonzero interior stage mole fraction of the reconstructed profiles exceeds 10^-(dps-30) in magnitude
    (so that their values and signs are resolved).  build() -> (pb, des, lim) at the current precision.
    Returns dict(sn, v, nF, conv, prof, chk, phys, dps, pb, des, lim)."""
    newton_kw = newton_kw or {}
    v_prev = start
    for _ in range(6):
        mp.mp.dps = dps
        pb, des, lim = build()
        sn = System(des, lim, nst)
        v0 = sn.start() if v_prev is None else mp.matrix([mpf(str(x)) for x in v_prev])
        v, nF, hist = sn.newton(v0, **newton_kw)
        conv = nF < mpf(10) ** (-(mp.mp.dps - 12))
        if not conv:
            break
        mn = min([abs(x) for x in v if x != 0] + [mpf(1)])
        prof = reconstruct(sn, v)
        smin = min([abs(x) for pr in prof.values() for xk in pr for x in xk if x != 0] + [mn])
        need = int(-mp.log10(smin)) + 80
        if need <= dps or dps >= max_dps:
            break
        dps = min(need, max_dps)
        v_prev = [str(x) for x in v]
    prof = reconstruct(sn, v)
    chk = stage_check(sn, v, prof)
    phys = bool(conv and chk['worst_rel'] < mpf(10) ** -25 and chk['unit_rel'] < mpf(10) ** -25
                and not chk['negs'] and chk['minS'] > 0)
    return dict(sn=sn, v=v, nF=nF, hist=hist, conv=conv, prof=prof, chk=chk, phys=phys, dps=mp.mp.dps,
                pb=pb, des=des, lim=lim)


def ratio_form_rows(sysm, v, J):
    """Row-rescale a Jacobian of the product form E_r (rows as in residual()) to L's ratio form
    E_r / (pi_bot pi_top * scale) at a zero (pivots computed with the same left eigenvectors)."""
    pb, N = sysm.pb, sysm.N
    X, W = sysm.unpack(v)
    rows, cols = sysm.labels()
    J2 = J.copy()
    for a, lab in enumerate(rows):
        if lab[0] != 'chart':
            continue
        s, r = lab[1], lab[2]
        bu, tu = pb.sec_ends(s)
        rts = sysm.section_data(s, W)
        k = sysm.kstar[s]
        lk = lhat(pb.al, rts[k])
        fac = 1 / (dot(lk, X[bu]) * dot(lk, X[tu]) * sysm.scale[s][r])
        for c in range(J.cols):
            J2[a, c] = J[a, c] * fac
    return J2


def trace_mmatrix(cl, sysm, tol=None):
    """For each component k: the trace block of DG_inf restricted to k (x-traces at units, w-traces at sections),
    the w-traces eliminated through the closed-end (chart) rows: A_k = B_x - B_w C_w^{-1} C_x (node balances in the
    trace contents).  Tests the nonsingular M-matrix property of +-A_k (Z sign pattern after making the diagonal
    positive, and entrywise nonnegative inverse).  Returns dict k -> info."""
    if tol is None:
        tol = mpf(10) ** (-(mp.mp.dps // 3))
    J, rows, cols = cl['J'], cl['rows'], cl['cols']
    out = {}
    tc = cl['tc']
    tr = cl['tr']
    comps = sorted(set(cols[j][2] for j in tc))
    for k in comps:
        xc = [j for j in tc if cols[j][0] == 'x' and cols[j][2] == k]
        wc = [j for j in tc if cols[j][0] == 'w' and cols[j][2] == k]
        rk = [i for i in tr if any(abs(J[i, j]) > tol for j in xc + wc)]
        cr = [i for i in rk if rows[i][0] == 'chart']
        br = [i for i in rk if rows[i][0] == 'bal']
        info = dict(nx=len(xc), nw=len(wc), nchart=len(cr), nbal=len(br))
        if len(cr) != len(wc) or len(br) != len(xc):
            info['ok'] = False
            info['note'] = 'counts do not allow the elimination'
            out[k] = info
            continue
        Cw = submatrix(J, cr, wc)
        Cx = submatrix(J, cr, xc)
        Bw = submatrix(J, br, wc)
        Bx = submatrix(J, br, xc)
        # orient every node balance as (outflow - inflow): residual() writes junction rows that way and product
        # rows as (inflow - outflow)
        for a_, i in enumerate(br):
            u = rows[i][1]
            if u[0] == u[1]:
                for j in range(Bw.cols):
                    Bw[a_, j] = -Bw[a_, j]
                for j in range(Bx.cols):
                    Bx[a_, j] = -Bx[a_, j]
        A = Bx - Bw * mp.inverse(Cw) * Cx
        n = A.rows
        sgn = [1 if A[i, i] > 0 else -1 for i in range(n)]
        diag_same = len(set(sgn)) == 1
        s0 = sgn[0]
        Ap = A * s0
        offmax = max([Ap[i, j] for i in range(n) for j in range(n) if i != j] + [mpf(-1)])
        Ainv = mp.inverse(Ap)
        invmin = min(Ainv[i, j] for i in range(n) for j in range(n))
        colsum = min(mp.fsum(Ap[i, j] for i in range(n)) for j in range(n))
        rowsum = min(mp.fsum(Ap[i, j] for j in range(n)) for i in range(n))
        info.update(ok=bool(diag_same and offmax <= tol and invmin >= -tol), diag_min=float(min(Ap[i, i] for i in range(n))),
                    offdiag_max=float(offmax), inverse_min=float(invmin), min_col_sum=float(colsum),
                    min_row_sum=float(rowsum), size=n)
        out[k] = info
    return out
