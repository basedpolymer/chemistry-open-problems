#!/usr/bin/env python3
"""limitsys.py -- the exact finite-stage system G_n and the limit system G_inf of the generalized Petlyuk
arrangement, in junction / net-flow variables (Agent L, attainment wave; see ../limit.md).

Module (imported by the check scripts).  Environment: Python 3.12, mpmath 1.3 (WSL Fedora sci env).

Contents
  design(alpha, z, q, eta_col, eta_bin, delta, pattern)  -> Design
      An operating specification F(delta, eta) of the arrangement and its infinite-stage limit configuration,
      built by Agent L's OWN simple recipe (HS preferred split + Astra slack in every column with >= 3
      components, relative to its actual feed; forced binary columns, pattern 'B' or 'T' per binary column;
      R = smallest reboiler duty making every binary column strictly feasible, plus delta).  It is only a test
      bed for Theorem L; the reference design is Agent D's (attainment/design/).
  System(design, n)   n = None means the limit system G_inf.
      .v_inf             the limit configuration (vector of unknowns)
      .G(v)              the residual vector (N-1 rows per section)
      .jac(v)            central finite-difference Jacobian (mpmath)
      .rows, .vars       classification (major / trace, component) used by the block-structure checks
      .unpack(v)         net flows d, b, junction compositions X, product compositions, section data

Indexing: components 0..N-1 (0 = lightest = most volatile).  Column (i, m), 0 <= i < m <= N-1.  Section
('T', i, m) runs from junction J(i,m) (bottom end) up to J(i, m-1) (product i if m-1 == i); section ('B', i, m)
runs from J(i,m) (top end) down to J(i+1, m) (product m if i+1 == m).  Product 0 = total condenser, product N-1 =
partial reboiler, products 1..N-2 = liquid side draws.
"""
import itertools
import mpmath as mp

mp.mp.dps = 60


# ----------------------------------------------------------------------------------------------- utilities
def bisect(f, lo, hi, it=400):
    flo = f(lo)
    for _ in range(it):
        mid = (lo + hi) / 2
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def feed_roots(al, f, U, i, m):
    """Roots chi_j in (al[j+1], al[j]), j = i..m-1, of sum_p al_p f_p/(al_p - chi) = U (f supported on [i,m])."""
    F = lambda t: sum(al[p] * f[p] / (al[p] - t) for p in range(i, m + 1)) - U
    out = []
    for j in range(i, m):
        eps = (al[j] - al[j + 1]) * mp.mpf(10) ** (-mp.mp.dps + 5)
        out.append(bisect(F, al[j + 1] + eps, al[j] - eps))
    return out


def psi(al, V, w, t, skip=None):
    return V - sum(al[i] * w[i] / (al[i] - t) for i in range(len(al)) if i != skip and w[i] != 0)


def lemma1_vector(al, chi, i, m):
    """Astra's Lemma 1: the normalized X on face [i,m] with G_X(chi_j) = 0 (chi interlacing)."""
    N = len(al)
    om = [mp.mpf(0)] * N
    for p in range(i, m + 1):
        num = mp.fprod([al[p] - c for c in chi])
        den = al[p] * mp.fprod([al[p] - al[k] for k in range(i, m + 1) if k != p])
        om[p] = num / den
    s = sum(om)
    return [o / s for o in om]


# ----------------------------------------------------------------------------------------------- design
class Design:
    pass


def design(alpha, z, q, eta_col=mp.mpf('1e-3'), eta_bin=mp.mpf('1e-2'), delta=mp.mpf('2e-2'), pattern=None,
           R=None):
    """Operating specification and limit configuration (Agent L's test recipe; see module docstring).

    eta_col: Astra slack of every column with >= 3 components (scalar or dict col -> value);
    eta_bin: forcing of every binary column (scalar or dict p -> value);
    pattern: dict p -> 'B' (component p forced down B[p,p+1] into product p+1) or 'T' (p+1 forced up T[p,p+1]
             into product p); default all 'B'.
    """
    al = [mp.mpf(a) for a in alpha]
    z = [mp.mpf(x) for x in z]
    q = mp.mpf(q)
    N = len(al)
    if pattern is None:
        pattern = {p: 'B' for p in range(N - 1)}
    ec = (lambda c: eta_col[c]) if isinstance(eta_col, dict) else (lambda c: mp.mpf(eta_col))
    eb = (lambda p: eta_bin[p]) if isinstance(eta_bin, dict) else (lambda p: mp.mpf(eta_bin))
    th = feed_roots(al, z, 1 - q, 0, N - 1)
    H = [sum(al[p] * z[p] / (th[r] - al[p]) for p in range(r + 1, N)) for r in range(N - 1)]
    D = Design()
    D.al, D.z, D.q, D.N, D.theta, D.H, D.pattern = al, z, q, N, th, H, dict(pattern)
    cols = sorted([(i, m) for i in range(N) for m in range(i + 1, N)], key=lambda c: (-(c[1] - c[0]), c[0]))
    D.cols = cols
    d, b, VT, VB, f, U, chi = {}, {}, {}, {}, {}, {}, {}
    zero = lambda: [mp.mpf(0)] * N
    for (i, m) in cols:
        if (i, m) == (0, N - 1):
            ff, UU = list(z), 1 - q
        else:
            ff, UU = zero(), mp.mpf(0)
            if m + 1 <= N - 1:
                ff = [ff[k] + d[(i, m + 1)][k] for k in range(N)]
                UU += VT[(i, m + 1)]
            if i - 1 >= 0:
                ff = [ff[k] + b[(i - 1, m)][k] for k in range(N)]
                UU -= VB[(i - 1, m)]
        f[(i, m)], U[(i, m)] = ff, UU
        if m - i >= 2:
            cr = feed_roots(al, ff, UU, i, m)
            chi[(i, m)] = cr

            def coef(p):
                num = mp.fprod([al[p] - c for c in cr])
                den = mp.fprod([al[p] - al[k] for k in range(i, m) if k != p])
                return num / den / al[p]
            vt0 = ff[i] / coef(i)
            d0 = zero()
            for p in range(i, m):
                d0[p] = vt0 * coef(p)
            # Astra's slack deformation (attainment-astra eq. (3)-(6)), local indices 1..mm
            a = al[i:m + 1]
            mm = len(a)
            thl = cr  # local theta_1..theta_{mm-1}
            beta = mp.fprod(a[1:mm - 1]) / mp.fprod(thl[1:mm - 2]) if mm > 3 else None
            h = zero()
            if mm > 3:
                for l in range(1, mm - 1):
                    num = (a[l] - beta) * mp.fprod([a[l] - thl[j] for j in range(1, mm - 2)])
                    den = a[l] * mp.fprod([a[l] - a[k] for k in range(1, mm - 1) if k != l])
                    h[i + l] = num / den
            e = ec((i, m))
            dd = [d0[k] + e * h[k] for k in range(N)]
            bb = [ff[k] - dd[k] for k in range(N)]
            bb[i] = mp.mpf(0)  # exact zero (rounding)
            d[(i, m)], b[(i, m)] = dd, bb
            VT[(i, m)] = vt0 + e
            VB[(i, m)] = vt0 - UU + e
    # binary columns: net flows
    for p in range(N - 1):
        c = (p, p + 1)
        ff = f[c]
        e = eb(p)
        dd, bb = zero(), zero()
        if pattern[p] == 'B':
            dd[p] = ff[p] - e
            bb[p] = e
            bb[p + 1] = ff[p + 1]
        else:
            dd[p] = ff[p]
            dd[p + 1] = e
            bb[p + 1] = ff[p + 1] - e
        d[c], b[c] = dd, bb
    # vapor chain: V_T[p,p+1] = R + sum_{p' >= p} U[p',p'+1]
    cT = {p: sum(U[(pp, pp + 1)] for pp in range(p, N - 1)) for p in range(N - 1)}
    cB = {p: cT[p] - U[(p, p + 1)] for p in range(N - 1)}

    def binary_margin(p, RR):
        """B-pinned: phi_T - phi''; T-pinned: phi''_T - phi_B (both increasing in R)."""
        c = (p, p + 1)
        vt, vb = RR + cT[p], RR + cB[p]
        if vt <= 0 or vb <= 0:
            return mp.mpf(-1)
        if pattern[p] == 'B':
            Dt = d[c][p]
            phT = al[p] * (1 - Dt / vt)
            g = lambda t: vb + al[p] * b[c][p] / (al[p] - t) + al[p + 1] * b[c][p + 1] / (al[p + 1] - t)
            eps = (al[p] - al[p + 1]) * mp.mpf(10) ** (-mp.mp.dps + 5)
            ph2 = bisect(g, al[p + 1] + eps, al[p] - eps)
            return phT - ph2
        else:
            Bt = b[c][p + 1]
            phB = al[p + 1] * (1 + Bt / vb)
            g = lambda t: vt - al[p] * d[c][p] / (al[p] - t) - al[p + 1] * d[c][p + 1] / (al[p + 1] - t)
            eps = (al[p] - al[p + 1]) * mp.mpf(10) ** (-mp.mp.dps + 5)
            ph2 = bisect(g, al[p + 1] + eps, al[p] - eps)
            return ph2 - phB
    if R is None:
        Rmin = []
        for p in range(N - 1):
            lo = max(-cT[p], -cB[p]) + mp.mpf('1e-30')
            hi = lo + 1
            while binary_margin(p, hi) <= 0:
                hi = lo + 2 * (hi - lo)
            Rmin.append(bisect(lambda RR: binary_margin(p, RR), lo, hi, it=200))
        D.Rmin_binary = Rmin
        D.Rmin = max(Rmin)
    else:                      # R prescribed (e.g. unforced test): no margin computation
        D.Rmin_binary, D.Rmin = None, None
    if R is None:
        R = D.Rmin + mp.mpf(delta)
    D.R = mp.mpf(R)
    for p in range(N - 1):
        VT[(p, p + 1)] = D.R + cT[p]
        VB[(p, p + 1)] = D.R + cB[p]
    # sections
    D.d, D.b, D.f, D.U, D.VT, D.VB, D.chi = d, b, f, U, VT, VB, chi
    D.LT = {c: VT[c] - sum(d[c]) for c in cols}
    D.LB = {c: VB[c] + sum(b[c]) for c in cols}
    D.P = []
    for p in range(N):
        tot = mp.mpf(0)
        if p - 1 >= 0:
            tot += sum(b[(p - 1, p)])
        if p + 1 <= N - 1:
            tot += sum(d[(p, p + 1)])
        D.P.append(tot)
    D.eta_col = {c: ec(c) for c in cols if c[1] - c[0] >= 2}
    D.eta_bin = {p: eb(p) for p in range(N - 1)}
    return D


# ----------------------------------------------------------------------------------------------- system
class System:
    """Exact n-stage system (n = None: limit system).  See ../limit.md Section 2-3."""

    def __init__(self, D, n=None, nstages=None):
        self.D = D
        self.n = n
        self.nst = nstages  # optional dict section -> n
        N, al = D.N, D.al
        self.N = N
        # ---- variable layout per column
        self.vars = []          # (col, 'd'|'b'|'X', k, kind, comp)  kind: 'maj' or 'tr'
        for c in D.cols:
            i, m = c
            if m - i >= 2:
                for k in range(N):
                    if k < i or k > m:
                        self.vars.append((c, 'd', k, 'tr', k))
                self.vars.append((c, 'd', m, 'tr', m))
                self.vars.append((c, 'b', i, 'tr', i))
                for k in range(i + 1, m - 1):
                    self.vars.append((c, 'd', k, 'maj', k))
            else:
                p = i
                if D.pattern[p] == 'B':
                    for k in range(N):
                        if k != p:
                            self.vars.append((c, 'd', k, 'tr', k))
                else:
                    for k in range(N):
                        if k not in (p, p + 1):
                            self.vars.append((c, 'd', k, 'tr', k))
                    self.vars.append((c, 'b', p, 'tr', p))
            for k in range(N):
                if k == i:
                    continue
                self.vars.append((c, 'X', k, 'tr' if (k < i or k > m) else 'maj', k))
        self.nv = len(self.vars)
        # ---- sections, modes, pinch, closed ends (from the limit configuration)
        lim = self._limit_quantities()
        self.lim = lim
        self.sections = []
        for c in D.cols:
            for kind in ('T', 'B'):
                self.sections.append((kind,) + c)
        self.sec_info = {}
        self.rows = []
        for s in self.sections:
            w = lim['w'][s]
            V = lim['V'][s]
            majors = [k for k in range(N) if w[k] != 0]
            traces = [k for k in range(N) if w[k] == 0]
            roots = self._major_roots(V, w, majors, s)
            modes = [('maj', r, None) for r in roots] + [('tr', al[k], k) for k in traces]
            modes.sort(key=lambda t: -t[1])
            vals = [t[1] for t in modes]
            assert all(vals[j] > vals[j + 1] for j in range(len(vals) - 1)), (s, vals)
            majroots = [t for t in modes if t[0] == 'maj']
            pin = majroots[-1] if s[0] == 'T' else majroots[0]
            ipin = modes.index(pin)
            info = dict(modes=modes, ipin=ipin, majors=majors, traces=traces)
            self.sec_info[s] = info
            for j, md in enumerate(modes):
                if j == ipin:
                    continue
                end = 'bot' if j < ipin else 'top'
                self.rows.append((s, j, md[0], md[2], end))
        assert len(self.rows) == self.nv, (len(self.rows), self.nv)
        self.v_inf = self._v_inf()

    # -- limit configuration
    def _limit_quantities(self):
        D, N = self.D, self.N
        w, V, L = {}, {}, {}
        for c in D.cols:
            w[('T',) + c] = list(D.d[c])
            w[('B',) + c] = [-x for x in D.b[c]]
            V[('T',) + c], L[('T',) + c] = D.VT[c], D.LT[c]
            V[('B',) + c], L[('B',) + c] = D.VB[c], D.LB[c]
        return dict(w=w, V=V, L=L)

    def _major_roots(self, V, w, majors, s):
        al = self.D.al
        poles = sorted([al[k] for k in majors], reverse=True)
        F = lambda t: psi(al, V, w, t)
        roots = []
        tiny = lambda a, b_: (a - b_) * mp.mpf(10) ** (-mp.mp.dps + 8)
        for j in range(len(poles) - 1):
            hi, lo = poles[j], poles[j + 1]
            roots.append(bisect(F, lo + tiny(hi, lo), hi - tiny(hi, lo)))
        if s[0] == 'T':  # extra root below the lowest pole
            lo, hi = mp.mpf(0), poles[-1]
            roots.append(bisect(F, lo, hi - tiny(hi, 0)))
        else:            # extra root above the highest pole
            lo = poles[0]
            hi = lo * 2 + 1
            while F(hi) < 0:
                hi = 2 * hi
            roots.append(bisect(F, lo + tiny(hi, lo), hi))
        return roots

    def _v_inf(self):
        D, N = self.D, self.N
        X = {}
        for c in D.cols:
            i, m = c
            if m - i >= 2:
                r = self.sec_info[('T',) + c]
                bsec = self.sec_info[('B',) + c]
                rT = [t[1] for t in r['modes'] if t[0] == 'maj']
                rB = [t[1] for t in bsec['modes'] if t[0] == 'maj']
                chi = [rT[0]] + list(D.chi[c][1:-1]) + [rB[-1]]
                X[c] = lemma1_vector(D.al, chi, i, m)
            else:
                p = i
                if D.pattern[p] == 'B':
                    info = self.sec_info[('B',) + c]
                    phi = info['modes'][info['ipin']][1]
                    w = self.lim['w'][('B',) + c]
                else:
                    info = self.sec_info[('T',) + c]
                    phi = info['modes'][info['ipin']][1]
                    w = self.lim['w'][('T',) + c]
                e = [w[k] / (D.al[k] - phi) if w[k] != 0 else mp.mpf(0) for k in range(N)]
                s = sum(e)
                X[c] = [x / s for x in e]
        self.X_inf = X
        v = []
        for (c, typ, k, kind, comp) in self.vars:
            if typ == 'd':
                v.append(D.d[c][k])
            elif typ == 'b':
                v.append(D.b[c][k])
            else:
                v.append(X[c][k])
        return v

    # -- unknowns -> derived quantities
    def unpack(self, v):
        D, N = self.D, self.N
        val = {}
        for x, (c, typ, k, kind, comp) in zip(v, self.vars):
            val[(c, typ, k)] = x
        d, b, f, X = {}, {}, {}, {}
        for c in D.cols:          # decreasing length: parents first
            i, m = c
            if c == (0, N - 1):
                ff = list(D.z)
            else:
                ff = [mp.mpf(0)] * N
                if m + 1 <= N - 1:
                    ff = [ff[k] + d[(i, m + 1)][k] for k in range(N)]
                if i - 1 >= 0:
                    ff = [ff[k] + b[(i - 1, m)][k] for k in range(N)]
            f[c] = ff
            tot = D.VT[c] - D.LT[c]
            dd = [None] * N
            for k in range(N):
                if (c, 'd', k) in val:
                    dd[k] = val[(c, 'd', k)]
            if m - i >= 2:
                dd[i] = ff[i] - val[(c, 'b', i)]
                dd[m - 1] = tot - sum(dd[k] for k in range(N) if k != m - 1)
            else:
                p = i
                if D.pattern[p] == 'B':
                    dd[p] = tot - sum(dd[k] for k in range(N) if k != p)
                else:
                    dd[p] = ff[p] - val[(c, 'b', p)]
                    dd[p + 1] = tot - sum(dd[k] for k in range(N) if k != p + 1)
            d[c] = dd
            b[c] = [ff[k] - dd[k] for k in range(N)]
            xx = [val.get((c, 'X', k)) for k in range(N)]
            xx[i] = 1 - sum(xx[k] for k in range(N) if k != i)
            X[c] = xx
        prod = []
        for p in range(N):
            tot = [mp.mpf(0)] * N
            if p - 1 >= 0:
                tot = [tot[k] + b[(p - 1, p)][k] for k in range(N)]
            if p + 1 <= N - 1:
                tot = [tot[k] + d[(p, p + 1)][k] for k in range(N)]
            prod.append([t / D.P[p] for t in tot])
        return d, b, f, X, prod

    def ends(self, s, X, prod):
        kind, i, m = s
        if kind == 'T':
            xb = X[(i, m)]
            xt = X[(i, m - 1)] if m - 1 > i else prod[i]
        else:
            xt = X[(i, m)]
            xb = X[(i + 1, m)] if i + 1 < m else prod[m]
        return xb, xt

    # -- eigen-data at current net flows
    def eigen(self, s, w):
        al, N = self.D.al, self.N
        V = self.lim['V'][s]
        info = self.sec_info[s]
        out = []
        for (typ, phi0, k) in info['modes']:
            if typ == 'maj':
                phi = phi0
                for _ in range(200):
                    F = psi(al, V, w, phi)
                    dF = -sum(al[i] * w[i] / (al[i] - phi) ** 2 for i in range(N))
                    step = F / dF
                    phi -= step
                    if abs(step) < mp.mpf(10) ** (-mp.mp.dps + 3) * abs(phi):
                        break
                ell = [al[i] / (al[i] - phi) for i in range(N)]
                out.append((phi, ell))
            else:
                # (al_k - phi) Psi^(k)(phi) - al_k w_k = 0, Newton from al_k
                phi = al[k]
                for _ in range(200):
                    Pk = psi(al, V, w, phi, skip=k)
                    dPk = -sum(al[i] * w[i] / (al[i] - phi) ** 2 for i in range(N) if i != k)
                    F = (al[k] - phi) * Pk - al[k] * w[k]
                    dF = -Pk + (al[k] - phi) * dPk
                    step = F / dF
                    phi -= step
                    if abs(step) < mp.mpf(10) ** (-mp.mp.dps + 3) * abs(phi) + mp.mpf(10) ** (-2 * mp.mp.dps):
                        break
                Pk = psi(al, V, w, phi, skip=k)
                c = w[k] / Pk          # = (al_k - phi_k)/al_k exactly
                ell = [c * al[i] / (al[i] - phi) if i != k else mp.mpf(1) for i in range(N)]
                out.append((phi, ell))
        return out

    def G(self, v):
        d, b, f, X, prod = self.unpack(v)
        res = []
        for s in self.sections:
            kind, i, m = s
            c = (i, m)
            w = d[c] if kind == 'T' else [-x for x in b[c]]
            xb, xt = self.ends(s, X, prod)
            eig = self.eigen(s, w)
            info = self.sec_info[s]
            jp = info['ipin']
            lp = eig[jp][1]
            dot = lambda l, x: sum(a * y for a, y in zip(l, x))
            pb, pt = dot(lp, xb), dot(lp, xt)
            ns = self.n if self.nst is None else self.nst[s]
            for j, (phi, ell) in enumerate(eig):
                if j == jp:
                    continue
                if j < jp:   # above the pinch: closed at the bottom
                    val = dot(ell, xb) / pb
                    if ns is not None:
                        val -= (eig[jp][0] / phi) ** (ns + 1) * dot(ell, xt) / pt
                else:        # below: closed at the top
                    val = dot(ell, xt) / pt
                    if ns is not None:
                        val -= (phi / eig[jp][0]) ** (ns + 1) * dot(ell, xb) / pb
                res.append(val)
        return res

    def jac(self, v, h=None):
        if h is None:
            h = mp.mpf(10) ** (-(mp.mp.dps // 3))
        n = len(v)
        J = mp.matrix(n, n)
        for j in range(n):
            vp = list(v)
            vm = list(v)
            vp[j] += h
            vm[j] -= h
            gp, gm = self.G(vp), self.G(vm)
            for i in range(n):
                J[i, j] = (gp[i] - gm[i]) / (2 * h)
        return J

    def row_class(self):
        """For each row: ('maj', None) or ('tr', k)."""
        return [(r[2], r[3]) for r in self.rows]

    def var_class(self):
        return [(kind, comp if kind == 'tr' else None) for (c, typ, k, kind, comp) in self.vars]


def newton(sys_, v0, tol=None, maxit=40, verbose=False):
    if tol is None:
        tol = mp.mpf(10) ** (-mp.mp.dps + 10)
    """Damped Newton (backtracking on the max-norm of the residual).  Returns (v, |G(v)|, iterations); the
    iteration count is -1 if the linear solve failed or no decrease could be found."""
    v = list(v0)
    g = sys_.G(v)
    nr = max(abs(x) for x in g)
    for it in range(maxit):
        if verbose:
            print('   newton', it, mp.nstr(nr, 5))
        if nr < tol:
            return v, nr, it
        try:
            J = sys_.jac(v)
            dv = mp.lu_solve(J, mp.matrix(g))
        except (ZeroDivisionError, ValueError):
            return v, nr, -1
        lam = mp.mpf(1)
        while lam > mp.mpf('1e-6'):
            try:
                vt = [v[i] - lam * dv[i] for i in range(len(v))]
                gt = sys_.G(vt)
                nt = max(abs(x) for x in gt)
            except (ZeroDivisionError, ValueError):
                nt = mp.inf
            if nt < nr:
                break
            lam /= 2
        else:
            return v, nr, -1
        v, g, nr = vt, gt, nt
    return v, nr, maxit
