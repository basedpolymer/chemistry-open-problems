#!/usr/bin/env python3
"""v1_independent.py -- V1 (independent verifier, fresh context): independent re-implementation and check of
Theorem D (design.md) and Theorem L / Theorem 1.4 / Corollary 4.7 (limit.md), from the statements only.

Nothing is imported from the solvers' scripts.  Differences from the solvers' implementations, on purpose:
  * the preferred split and Astra's direction h are obtained by LINEAR SOLVES of their defining conditions
    (Psi_T(theta^c_r) = 0 for all feed roots with d_i = f_i, d_m = 0;  sum h = 0 and K(theta^c_r) = 0 at the middle
    roots), not from the closed product formulas;
  * the spectra and left/right eigenvectors of the section matrices M = V diag(alpha) - w alpha^T are computed by a
    general eigen-solver (mpmath.eig), not by Underwood's formulas; left eigenvectors are normalized by a FIXED entry
    (the entry of largest modulus at the limit point; for a trace mode alpha_k this is entry k, i.e. limit.md's
    regularized normalization), which is an analytic normalization near the limit point;
  * the limit system G_inf and the exact system G_n are built from these numerical eigen-data, and the Jacobian is
    computed by central finite differences in high precision;
  * finite-stage zeros are found by the chord iteration v <- v - DG_inf(v_inf)^{-1} G_n(v) (the contraction of
    Corollary 4.7), the stage profiles are reconstructed by the recursion, and EVERY component balance of EVERY unit
    (stages, junctions, side draws, condenser, reboiler), the total-condenser composition rule, and positivity are
    checked directly on the physical equations, independently of the reformulation.
Checks (PASS/FAIL each):
  D1 flows and specification (H1), total/CMO balances with zero duty at junctions and side draws, R = H* + delta;
  D2 limit net-flow supports (H2), middle feed roots of non-binary columns = global roots (Lemma 5.3);
  D3 Astra's conditions (H3): interlacing, Psi_T(alpha_m) > 0, Psi_B(alpha_i) > 0, cone weights A,B,C,D > 0;
  D4 binary margins (H4): kappa_p >= e_p + delta/2, forcing increases the margin, forced roots located (L3, L4);
  D5 kappa_p(eta_nb; delta) - kappa_p(eta_nb; 0) = delta and kappa_p(0; delta) = e_p + delta (Prop. 5.4 (b), (c));
  D6 flux bookkeeping (Lemma 5.2 (e)) at random t;
  D7 spectra (H5): distinct positive modes, pinch rule (top: smallest major root, bottom: largest), pivots > 0 at both
     ends in the sum-normalized eigenbasis (H8), trace signs and kappa > 0 (H7);
  D8 half-profiles (H6): M^k x_bot >= 0 and M^-k x_top >= 0 for k <= K;
  D9 product impurities = (eta_{p-1}[pi_{p-1}=B] + eta_p[pi_p=T]) / P_p;
  L1 G_inf(v_inf) = 0;
  L2 zero blocks of Theorem L in the coordinates of limit.md Table 2 (trace rows vs major coordinates; trace rows of
     k vs trace coordinates of k' != k); trace rows equal Lemma 4.3's prediction;
  L3 major rows of a column depend only on the coordinates of the column and its ancestors; diagonal blocks of the
     major block nonsingular;
  L4 every trace block T_k nonsingular; A_k (Prop. 4.5) upper triangular with positive diagonal in the h-order;
  L5 sigma_min(DG_inf(v_inf)) > 0;
  F1 (optional) finite n: chord iteration converges, every unit balance holds, all mole fractions > 0, x_top reached,
     |v_n - v_inf| decreases with n.
Usage:  python v1_independent.py CASE [dps] [nlist] [extra digits for F1, default 30]
  CASE = "n3:PAT" | "n4:PAT" | "rand:N:seed[:pattern][:frac]" | "near:N:seed[:pattern][:frac]"
         (frac = eta_p / eta_hat_p, default 0.25; "near": one consecutive volatility ratio 1.05)
  nlist = comma list of section lengths for F1 (default none), e.g. 60,40,30
Environment: Python 3.12, mpmath 1.3, one thread.  Runtime: N = 3 seconds; N = 4 ~ 1 min; N = 5 several minutes.
Verifier: Claude Opus 5.5 (effort xhigh), fresh context.
"""
import sys
import time
import random
import mpmath as mp

FAIL = []


def check(name, ok, info=''):
    print(f'  [{"PASS" if ok else "FAIL"}] {name} {info}')
    if not ok:
        FAIL.append(name)
    sys.stdout.flush()


# ------------------------------------------------------------------------------------------------------------
# generic helpers
def bisect(fun, lo, hi, increasing):
    a, b = mp.mpf(lo), mp.mpf(hi)
    for _ in range(int(mp.mp.prec) + 20):
        mid = (a + b) / 2
        v = fun(mid)
        if v == 0:
            return mid
        if (v > 0) == increasing:
            b = mid
        else:
            a = mid
    return (a + b) / 2


def vec(N):
    return [mp.mpf(0)] * N


def add(u, v):
    return [a + b for a, b in zip(u, v)]


def sub(u, v):
    return [a - b for a, b in zip(u, v)]


def scal(c, u):
    return [c * a for a in u]


def dot(u, v):
    return mp.fsum(a * b for a, b in zip(u, v))


class Data:
    def __init__(self, alpha, z, q):
        self.N = len(alpha)
        self.a = [mp.mpf(x) for x in alpha]
        self.z = [mp.mpf(x) for x in z]
        s = mp.fsum(self.z)
        self.z = [x / s for x in self.z]
        self.q = mp.mpf(q)
        a, N = self.a, self.N
        F = lambda t: mp.fsum(a[p] * self.z[p] / (a[p] - t) for p in range(N)) - (1 - self.q)   # increasing
        self.theta = [bisect(F, a[j + 1], a[j], True) for j in range(N - 1)]
        self.H = [mp.fsum(a[p] * self.z[p] / (self.theta[j] - a[p]) for p in range(j + 1, N)) for j in range(N - 1)]
        self.Hs = max(self.H)
        self.e = [self.Hs - h for h in self.H]


def Psi(a, V, w, t, skip=None):
    return V - mp.fsum(a[i] * w[i] / (a[i] - t) for i in range(len(a)) if w[i] != 0 and i != skip)


def Gx(a, x, t, skip=None):
    return mp.fsum(a[i] * x[i] / (a[i] - t) for i in range(len(a)) if x[i] != 0 and i != skip)


# ------------------------------------------------------------------------------------------------------------
# design F(delta, eta, pi): own construction (linear solves of the defining conditions)
def build(D, delta, eta_nb, frac, pat, R_override=None, eta_bin=None):
    N, a = D.N, D.a
    cols = {}
    for Lc in range(N - 1, 1, -1):
        for i in range(1, N - Lc + 1):
            m = i + Lc
            c = {'kind': 'nb', 'i': i, 'm': m}
            if (i, m) == (1, N):
                f, U = list(D.z), 1 - D.q
            else:
                f, U = vec(N), mp.mpf(0)
                if m + 1 <= N:
                    f = add(f, cols[(i, m + 1)]['d'])
                    U += cols[(i, m + 1)]['VT']
                if i - 1 >= 1:
                    f = add(f, cols[(i - 1, m)]['b'])
                    U -= cols[(i - 1, m)]['VB']
            I0, M0 = i - 1, m - 1                       # 0-based light and heavy keys
            Phi = lambda t: U - mp.fsum(a[p] * f[p] / (a[p] - t) for p in range(I0, M0 + 1))   # decreasing
            th = [bisect(Phi, a[r + 1], a[r], False) for r in range(I0, M0)]               # th[r - I0]
            # preferred split: unknowns VT, d_{i+1..m-1}; d_i = f_i, d_m = 0; Psi_T(th_r) = 0 for all feed roots
            idx = list(range(I0 + 1, M0))
            nu = 1 + len(idx)
            A = mp.matrix(nu, nu)
            rhs = mp.matrix(nu, 1)
            for row, t in enumerate(th):
                A[row, 0] = 1
                for col, p in enumerate(idx):
                    A[row, col + 1] = -a[p] / (a[p] - t)
                rhs[row] = a[I0] * f[I0] / (a[I0] - t)
            sol = mp.lu_solve(A, rhs)
            VT0 = sol[0]
            d0 = vec(N)
            d0[I0] = f[I0]
            for col, p in enumerate(idx):
                d0[p] = sol[col + 1]
            # Astra's direction: h on [i+1, m-1], sum h = 0, K(th_r) = 0 at the middle feed roots
            h = vec(N)
            if len(idx) >= 2:
                nh = len(idx)
                A = mp.matrix(nh, nh)
                rhs = mp.matrix(nh, 1)
                for col in range(nh):
                    A[0, col] = 1
                for row, t in enumerate(th[1:-1]):
                    for col, p in enumerate(idx):
                        A[row + 1, col] = a[p] / (a[p] - t)
                    rhs[row + 1] = 1
                sol = mp.lu_solve(A, rhs)
                for col, p in enumerate(idx):
                    h[p] = sol[col]
            eta = mp.mpf(eta_nb[(i, m)] if isinstance(eta_nb, dict) else eta_nb)
            d = add(d0, scal(eta, h))
            b = sub(f, d)
            VT = VT0 + eta
            VB = VT - U
            c.update(f=f, U=U, th=th, d0=d0, h=h, VT0=VT0, eta=eta, d=d, b=b, VT=VT, VB=VB,
                     LT=VT - mp.fsum(d), LB=VB + mp.fsum(b), K=lambda t, h=h: 1 - mp.fsum(
                         a[p] * h[p] / (a[p] - t) for p in range(N) if h[p] != 0))
            cols[(i, m)] = c
    # binary columns: feeds, U, vapors from the reboiler up
    for p in range(1, N):
        c = {'kind': 'bin', 'i': p, 'm': p + 1, 'p': p, 'pi': pat[p - 1]}
        f, U = vec(N), mp.mpf(0)
        if p + 2 <= N:
            f = add(f, cols[(p, p + 2)]['d'])
            U += cols[(p, p + 2)]['VT']
        if p - 1 >= 1:
            f = add(f, cols[(p - 1, p + 1)]['b'])
            U -= cols[(p - 1, p + 1)]['VB']
        c.update(f=f, U=U)
        cols[(p, p + 1)] = c
    R = D.Hs + delta if R_override is None else R_override
    VB = R
    for p in range(N - 1, 0, -1):
        c = cols[(p, p + 1)]
        c['VB'] = VB
        c['VT'] = VB + c['U']
        VB = c['VT']
    for p in range(1, N):
        c = cols[(p, p + 1)]
        f, P0 = c['f'], p - 1
        Phi = lambda t: c['U'] - mp.fsum(a[k] * f[k] / (a[k] - t) for k in (P0, P0 + 1))
        th1 = bisect(Phi, a[P0 + 1], a[P0], False)
        c['th'] = [th1]
        c['kappa'] = c['VT'] - a[P0] * f[P0] / (a[P0] - th1)
        if c['pi'] == 'B':
            if p == 1:
                hat = f[P0]
            else:
                PsiB0 = c['VB'] + a[P0 + 1] * f[P0 + 1] / (a[P0 + 1] - a[P0 - 1])
                hat = min(f[P0], (a[P0 - 1] - a[P0]) * PsiB0 / a[P0])
        else:
            if p + 1 == N:
                hat = min(f[P0 + 1], c['VT'] - f[P0])
            else:
                PsiT0 = c['VT'] - a[P0] * f[P0] / (a[P0] - a[P0 + 2])
                hat = min(f[P0 + 1], (a[P0 + 1] - a[P0 + 2]) * PsiT0 / a[P0 + 1])
        c['hat'] = hat
        eta = frac * hat if eta_bin is None else mp.mpf(eta_bin[p])
        c['eta'] = eta
        d, b = vec(N), vec(N)
        if c['pi'] == 'B':
            d[P0] = f[P0] - eta
            b[P0], b[P0 + 1] = eta, f[P0 + 1]
        else:
            d[P0], d[P0 + 1] = f[P0], eta
            b[P0 + 1] = f[P0 + 1] - eta
        c.update(d=d, b=b, LT=c['VT'] - mp.fsum(d), LB=c['VB'] + mp.fsum(b))
    for c in cols.values():
        c['D'] = c['VT'] - c['LT']
        c['B'] = c['LB'] - c['VB']
    P = {}
    for p in range(1, N + 1):
        P[p] = (cols[(p - 1, p)]['B'] if p > 1 else 0) + (cols[(p, p + 1)]['D'] if p < N else 0)
    return {'D': D, 'cols': cols, 'P': P, 'R': R, 'delta': delta, 'pat': pat}


# ------------------------------------------------------------------------------------------------------------
# arrangement topology
def sections(N):
    out = []
    for i in range(1, N + 1):
        for m in range(i + 1, N + 1):
            out += [('T', i, m), ('B', i, m)]
    return out


def ends(s):          # (bottom node, top node), node = (a, b), a <= b
    k, i, m = s
    return ((i, m), (i, m - 1)) if k == 'T' else ((i + 1, m), (i, m))


def h_(J):
    return J[0] + J[1]


# ------------------------------------------------------------------------------------------------------------
# limit configuration (own computation)
def eigdata(a, V, w):
    N = len(a)
    M = mp.matrix(N, N)
    for r in range(N):
        for c in range(N):
            M[r, c] = (V * a[r] if r == c else 0) - w[r] * a[c]
    E, EL, ER = mp.eig(M, left=True, right=True)
    order = sorted(range(N), key=lambda j: -mp.re(E[j]))
    lam = [mp.re(E[j]) for j in order]
    imag = max(abs(mp.im(E[j])) for j in order)
    left = [[mp.re(EL[j, c]) for c in range(N)] for j in order]
    right = [[mp.re(ER[c, j]) for c in range(N)] for j in order]
    return M, lam, left, right, imag


def limit_config(S):
    D, cols, P = S['D'], S['cols'], S['P']
    N, a = D.N, D.a
    X = {}
    for (i, m), c in cols.items():
        I0, M0 = i - 1, m - 1
        if c['kind'] == 'nb':
            wT, wB = c['d'], scal(-1, c['b'])
            # roots r_c in (a_{i+1}, a_i) of Psi_T (decreasing) and s_c in (a_m, a_{m-1}) of Psi_B (increasing)
            rc = bisect(lambda t: Psi(a, c['VT'], wT, t), a[I0 + 1], a[I0], False)
            sc = bisect(lambda t: Psi(a, c['VB'], wB, t), a[M0], a[M0 - 1], True)
            uc = bisect(lambda t: Psi(a, c['VT'], wT, t), mp.mpf(0), a[M0 - 1], False)
            vc = bisect(lambda t: Psi(a, c['VB'], wB, t), a[I0 + 1],
                        a[I0 + 1] + mp.fsum(a[k] * abs(wB[k]) for k in range(N)) / c['VB'] + a[0] + 1, True)
            c.update(rc=rc, sc=sc, uc=uc, vc=vc)
            chis = [rc] + [D.theta[r] for r in range(I0 + 1, M0 - 1)] + [sc]
            sup = list(range(I0, M0 + 1))
            n_ = len(sup)
            A = mp.matrix(n_, n_)
            rhs = mp.matrix(n_, 1)
            for col in range(n_):
                A[0, col] = 1
            rhs[0] = 1
            for row, chi in enumerate(chis):
                for col, k in enumerate(sup):
                    A[row + 1, col] = a[k] / (a[k] - chi)
            sol = mp.lu_solve(A, rhs)
            x = vec(N)
            for col, k in enumerate(sup):
                x[k] = sol[col]
            X[(i, m)] = x
        else:
            sec = ('B', i, m) if c['pi'] == 'B' else ('T', i, m)
            V = c['VB'] if c['pi'] == 'B' else c['VT']
            w = scal(-1, c['b']) if c['pi'] == 'B' else c['d']
            _, lam, left, right, _ = eigdata(a, V, w)
            # pinch: B-pinned -> largest major root of Psi_B (> alpha_p); T-pinned -> smallest major root of Psi_T
            majors = [r for r in range(N) if not any(w[k] == 0 and abs(lam[r] / V - a[k]) < mp.mpf(10) ** (-mp.mp.dps // 2)
                                                     for k in range(N))]
            r = majors[0] if c['pi'] == 'B' else majors[-1]
            e = right[r]
            s = mp.fsum(e)
            xx = [t / s for t in e]
            # the eigen-solver returns ~1e-(dps) entries outside {p, p+1}; check, then set them to exact zeros
            off = max(abs(xx[k]) for k in range(N) if k not in (I0, I0 + 1))
            assert off < mp.mpf(10) ** (-(mp.mp.dps - 15)), off
            X[(i, m)] = [xx[k] if k in (I0, I0 + 1) else mp.mpf(0) for k in range(N)]
            c['forced_root'] = lam[r] / V
    x = {}
    for p in range(1, N + 1):
        n = vec(N)
        if p > 1:
            n = add(n, cols[(p - 1, p)]['b'])
        if p < N:
            n = add(n, cols[(p, p + 1)]['d'])
        x[p] = [t / P[p] for t in n]
    S['X'], S['x'] = X, x
    return S


def node_comp(S, J, X=None, x=None):
    X = S['X'] if X is None else X
    x = S['x'] if x is None else x
    return x[J[0]] if J[0] == J[1] else X[J]


def section_flows(S, s, cols=None):
    cols = S['cols'] if cols is None else cols
    k, i, m = s
    c = S['cols'][(i, m)]
    if k == 'T':
        return c['VT'], c['LT'], cols[(i, m)]['d']
    return c['VB'], c['LB'], scal(-1, cols[(i, m)]['b'])


# ------------------------------------------------------------------------------------------------------------
# coordinates of limit.md Table 2
def column_order(N):
    return sorted([(i, m) for i in range(1, N + 1) for m in range(i + 1, N + 1)], key=lambda c: (-(c[1] - c[0]), c))


def coord_list(S):
    N = S['D'].N
    out = []           # (col, var, comp0, class)   class = 'maj' or k (1-based trace component)
    for (i, m) in column_order(N):
        c = S['cols'][(i, m)]
        rng = range(i, m + 1)
        if c['kind'] == 'nb':
            for k in range(1, N + 1):
                if k not in rng:
                    out.append(((i, m), 'd', k - 1, k))
            out.append(((i, m), 'd', m - 1, m))
            out.append(((i, m), 'b', i - 1, i))
            for k in range(1, N + 1):
                if k not in rng:
                    out.append(((i, m), 'X', k - 1, k))
            for p in range(i + 1, m - 1):
                out.append(((i, m), 'd', p - 1, 'maj'))
            for p in range(i + 1, m + 1):
                out.append(((i, m), 'X', p - 1, 'maj'))
        else:
            p = i
            if c['pi'] == 'B':
                for k in range(1, N + 1):
                    if k != p:
                        out.append(((i, m), 'd', k - 1, k))
            else:
                for k in range(1, N + 1):
                    if k not in (p, p + 1):
                        out.append(((i, m), 'd', k - 1, k))
                out.append(((i, m), 'b', p - 1, p))
            for k in range(1, N + 1):
                if k not in (p, p + 1):
                    out.append(((i, m), 'X', k - 1, k))
            out.append(((i, m), 'X', p, 'maj'))
    return out


def state_from_coords(S, coords, vals):
    """(d_c, X_c, b_c, f_c) for all columns from coordinate values, parents first."""
    N = S['D'].N
    by = {}
    for (cd, v) in zip(coords, vals):
        by.setdefault(cd[0], []).append((cd[1], cd[2], v))
    st = {}
    for (i, m) in column_order(N):
        c = S['cols'][(i, m)]
        if (i, m) == (1, N):
            f = list(S['D'].z)
        else:
            f = vec(N)
            if m + 1 <= N:
                f = add(f, st[(i, m + 1)]['d'])
            if i - 1 >= 1:
                f = add(f, st[(i - 1, m)]['b'])
        d, X = [None] * N, [None] * N
        bco = {}
        for (var, k0, v) in by[(i, m)]:
            if var == 'd':
                d[k0] = v
            elif var == 'X':
                X[k0] = v
            else:
                bco[k0] = v
        for k0, v in bco.items():
            d[k0] = f[k0] - v
        dep = [k for k in range(N) if d[k] is None]
        assert len(dep) == 1, (i, m, dep)
        d[dep[0]] = c['D'] - mp.fsum(t for t in d if t is not None)
        depX = [k for k in range(N) if X[k] is None]
        assert len(depX) == 1
        X[depX[0]] = 1 - mp.fsum(t for t in X if t is not None)
        st[(i, m)] = {'f': f, 'd': d, 'b': sub(f, d), 'X': X}
    return st


def coords_at_limit(S, coords):
    vals = []
    for (col, var, k0, cl) in coords:
        c = S['cols'][col]
        vals.append(c['d'][k0] if var == 'd' else (c['b'][k0] if var == 'b' else S['X'][col][k0]))
    return vals


# ------------------------------------------------------------------------------------------------------------
# the systems G_n / G_inf from numerical eigen-data
class System:
    def __init__(self, S):
        self.S = S
        D = S['D']
        self.N, self.a = D.N, D.a
        self.coords = coord_list(S)
        self.v0 = coords_at_limit(S, self.coords)
        st = state_from_coords(S, self.coords, self.v0)
        # sanity: coordinates reproduce the limit configuration
        err = max(abs(st[c]['d'][k] - S['cols'][c]['d'][k]) + abs(st[c]['X'][k] - S['X'][c][k])
                  for c in st for k in range(self.N))
        self.recon_err = err
        self.secs = sections(self.N)
        self.meta = []          # per section: (r*, normalization indices, mode classes)
        tol = mp.mpf(10) ** (-mp.mp.dps // 2)
        for s in self.secs:
            V, L, w = section_flows(S, s)
            _, lam, left, right, imag = eigdata(self.a, V, w)
            cls = []
            for r in range(self.N):
                tr = [k for k in range(self.N) if w[k] == 0 and abs(lam[r] / V - self.a[k]) < tol]
                cls.append(('trace', tr[0]) if tr else ('major', None))
            majors = [r for r in range(self.N) if cls[r][0] == 'major']
            rstar = majors[-1] if s[0] == 'T' else majors[0]
            jn = []
            for r in range(self.N):
                if cls[r][0] == 'trace':
                    jn.append(cls[r][1])
                else:
                    jn.append(max(range(self.N), key=lambda j: abs(left[r][j])))
            self.meta.append({'rstar': rstar, 'jn': jn, 'cls': cls, 'imag': imag})
        # row metadata
        self.rows = []
        for si, s in enumerate(self.secs):
            bot, top = ends(s)
            mt = self.meta[si]
            for r in range(self.N):
                if r == mt['rstar']:
                    continue
                closed = bot if r < mt['rstar'] else top
                other = top if r < mt['rstar'] else bot
                cl = mt['cls'][r]
                self.rows.append({'sec': s, 'si': si, 'r': r, 'closed': closed, 'open': other,
                                  'class': 'maj' if cl[0] == 'major' else cl[1] + 1,
                                  'col': (s[1], s[2])})

    def eval(self, vals, n=None, want=False):
        S, N, a = self.S, self.N, self.a
        st = state_from_coords(S, self.coords, vals)
        X = {c: st[c]['X'] for c in st}
        x = {}
        for p in range(1, N + 1):
            nn = vec(N)
            if p > 1:
                nn = add(nn, st[(p - 1, p)]['b'])
            if p < N:
                nn = add(nn, st[(p, p + 1)]['d'])
            x[p] = [t / S['P'][p] for t in nn]
        comp = lambda J: x[J[0]] if J[0] == J[1] else X[J]
        out = []
        cache = {}
        for si, s in enumerate(self.secs):
            V, L, _ = section_flows(S, s)
            w = st[(s[1], s[2])]['d'] if s[0] == 'T' else scal(-1, st[(s[1], s[2])]['b'])
            M, lam, left, right, imag = eigdata(a, V, w)
            mt = self.meta[si]
            left = [[t / left[r][mt['jn'][r]] for t in left[r]] for r in range(N)]
            bot, top = ends(s)
            xb, xt = comp(bot), comp(top)
            rs = mt['rstar']
            pb, pt = dot(left[rs], xb), dot(left[rs], xt)
            cache[s] = (M, lam, left, right, xb, xt, w, V, L)
            for r in range(N):
                if r == rs:
                    continue
                if r < rs:
                    g = dot(left[r], xb) / pb
                    if n is not None:
                        g -= (lam[rs] / lam[r]) ** (n[s] + 1) * dot(left[r], xt) / pt
                else:
                    g = dot(left[r], xt) / pt
                    if n is not None:
                        g -= (lam[r] / lam[rs]) ** (n[s] + 1) * dot(left[r], xb) / pb
                out.append(g)
        if want:
            return out, st, X, x, cache
        return out

    def jacobian(self, vals, n=None, hstep=None):
        hstep = mp.mpf(10) ** (-(mp.mp.dps // 3)) if hstep is None else hstep
        dim = len(vals)
        J = mp.matrix(len(self.rows), dim)
        for j in range(dim):
            vp = list(vals)
            vm = list(vals)
            vp[j] += hstep
            vm[j] -= hstep
            gp, gm = self.eval(vp, n), self.eval(vm, n)
            for r in range(len(gp)):
                J[r, j] = (gp[r] - gm[r]) / (2 * hstep)
        return J


# ------------------------------------------------------------------------------------------------------------
def run_design_checks(S, eta_nb):
    D, cols, P = S['D'], S['cols'], S['P']
    N, a, delta = D.N, D.a, S['delta']
    tiny = mp.mpf(10) ** (-(mp.mp.dps - 12))
    # D1 flows and balances
    ok = all(c['VT'] > 0 and c['VB'] > 0 and c['LT'] > 0 and c['LB'] > 0 and c['D'] > 0 and c['B'] > 0
             for c in cols.values()) and all(P[p] > 0 for p in P)
    ok &= all(all(c['f'][k] > 0 for k in range(c['i'] - 1, c['m'])) and
              all(c['f'][k] == 0 for k in range(N) if not (c['i'] - 1 <= k <= c['m'] - 1)) for c in cols.values())
    # unit balances (total flow, vapor/CMO) of the specification
    bal = []
    for (i, m), c in cols.items():
        vin = c['VB'] + (cols[(i, m + 1)]['VT'] if m < N else 0) + ((1 - D.q) if (i, m) == (1, N) else 0)
        vout = c['VT'] + (cols[(i - 1, m)]['VB'] if i > 1 else 0)
        lin = c['LT'] + (cols[(i - 1, m)]['LB'] if i > 1 else 0) + (D.q if (i, m) == (1, N) else 0)
        lout = c['LB'] + (cols[(i, m + 1)]['LT'] if m < N else 0)
        bal += [vout - vin, lout - lin]
    for p in range(2, N):
        bal += [cols[(p - 1, p)]['VB'] - cols[(p, p + 1)]['VT'],
                cols[(p - 1, p)]['LB'] - cols[(p, p + 1)]['LT'] - P[p]]
    bal += [cols[(1, 2)]['VT'] - cols[(1, 2)]['LT'] - P[1], cols[(N - 1, N)]['LB'] - cols[(N - 1, N)]['VB'] - P[N]]
    ok &= max(abs(t) for t in bal) < tiny
    ok &= abs(cols[(N - 1, N)]['VB'] - (D.Hs + delta)) < tiny and abs(mp.fsum(P.values()) - 1) < tiny
    check('D1 flows > 0, specification balances (zero duty at junctions/side draws), R = H*+delta', ok,
          f'(max balance residual {mp.nstr(max(abs(t) for t in bal), 3)})')
    # D2 supports, middle roots
    ok = True
    for (i, m), c in cols.items():
        I0, M0 = i - 1, m - 1
        if c['kind'] == 'nb':
            ok &= all(c['d'][k] > 0 for k in range(I0, M0)) and all(c['d'][k] == 0 for k in range(N) if not I0 <= k < M0)
            ok &= all(c['b'][k] > 0 for k in range(I0 + 1, M0 + 1)) and all(
                abs(c['b'][k]) < tiny for k in range(N) if not I0 < k <= M0)
        else:
            p0 = I0
            if c['pi'] == 'B':
                ok &= c['d'][p0] > 0 and all(c['d'][k] == 0 for k in range(N) if k != p0)
                ok &= c['b'][p0] == c['eta'] > 0 and c['b'][p0 + 1] > 0
            else:
                ok &= c['b'][p0 + 1] > 0 and all(abs(c['b'][k]) < tiny for k in range(N) if k != p0 + 1)
                ok &= c['d'][p0 + 1] == c['eta'] > 0 and c['d'][p0] > 0
    dm = mp.mpf(0)
    for (i, m), c in cols.items():
        if c['kind'] == 'nb':
            for r in range(i, m - 2):
                dm = max(dm, abs(c['th'][r - i + 1] - D.theta[r]))
    ok &= dm < tiny
    check('D2 limit net-flow supports (H2); middle feed roots global (Lemma 5.3)', ok, f'(max |theta^c - theta| {mp.nstr(dm, 3)})')
    # D3 Astra's conditions
    ok = True
    for (i, m), c in cols.items():
        if c['kind'] != 'nb':
            continue
        I0, M0 = i - 1, m - 1
        ok &= c['eta'] > 0
        ok &= a[I0] > c['rc'] > c['th'][0] > c['vc'] > a[I0 + 1]
        ok &= a[M0 - 1] > c['uc'] > c['th'][-1] > c['sc'] > a[M0]
        wT, wB = c['d'], scal(-1, c['b'])
        for t in c['th'][1:-1]:
            ok &= abs(Psi(a, c['VT'], wT, t)) < tiny and abs(Psi(a, c['VB'], wB, t)) < tiny
        ok &= Psi(a, c['VT'], wT, a[M0]) > 0 and Psi(a, c['VB'], wB, a[I0]) > 0
        ok &= all(S['X'][(i, m)][k] > 0 for k in range(I0, M0 + 1))
    check('D3 interlacing r>th_i>v, u>th_{m-1}>s; common middle roots; Psi_T(a_m), Psi_B(a_i) > 0; X_c > 0', ok)
    # D4 binary margins
    ok = True
    for p in range(1, N):
        c = cols[(p, p + 1)]
        P0 = p - 1
        ok &= c['kappa'] >= D.e[P0] + delta / 2 - tiny and 0 < c['eta'] < c['hat']
        wT, wB = c['d'], scal(-1, c['b'])
        if c['pi'] == 'B':
            phiT = a[P0] * c['LT'] / c['VT']
            phi1 = c['forced_root']
            phi2 = bisect(lambda t: Psi(a, c['VB'], wB, t), a[P0 + 1], a[P0], True)
            ok &= Psi(a, c['VB'], wB, phiT) > c['kappa'] and phiT > phi2 and phiT > a[P0 + 1]
            ok &= phi1 > a[P0] and (p == 1 or phi1 < a[P0 - 1])
            c['phi_np'], c['phi_2'] = phiT, phi2
        else:
            phiB = a[P0 + 1] * c['LB'] / c['VB']
            phil = c['forced_root']
            phim = bisect(lambda t: Psi(a, c['VT'], wT, t), a[P0 + 1], a[P0], False)
            ok &= Psi(a, c['VT'], wT, phiB) > c['kappa'] and phiB < phim and phiB < a[P0]
            ok &= phil < a[P0 + 1] and (p + 1 == N or phil > a[P0 + 2]) and phil > 0
            c['phi_np'], c['phi_2'] = phiB, phim
    check('D4 kappa_p >= e_p + delta/2, 0 < eta_p < eta_hat_p, forcing increases the margin, (L3), (L4)', ok,
          f'(min kappa/delta {mp.nstr(min(cols[(p, p + 1)]["kappa"] / delta for p in range(1, N)), 4)})')
    # D5 kappa identities
    S0 = build(D, mp.mpf(0), eta_nb, mp.mpf('0.25'), S['pat'])
    Sz = build(D, delta, mp.mpf(0), mp.mpf('0.25'), S['pat'])
    e1 = max(abs(cols[(p, p + 1)]['kappa'] - S0['cols'][(p, p + 1)]['kappa'] - delta) for p in range(1, N))
    e2 = max(abs(Sz['cols'][(p, p + 1)]['kappa'] - D.e[p - 1] - delta) for p in range(1, N))
    check('D5 kappa_p(eta;delta) - kappa_p(eta;0) = delta and kappa_p(0;delta) = e_p + delta', max(e1, e2) < tiny,
          f'({mp.nstr(e1, 3)}, {mp.nstr(e2, 3)})')
    # D6 flux bookkeeping at random t
    rnd = random.Random(7)
    worst = mp.mpf(0)
    for _ in range(3):
        t = mp.mpf(rnd.uniform(0.05, 1.5 * float(a[0])))
        if min(abs(t - x) for x in a) < 1e-3:
            continue
        PsiS = {}
        for s in sections(N):
            V, L, w = section_flows(S, s)
            PsiS[s] = Psi(a, V, w, t)
        nodes = [(i, m) for i in range(1, N + 1) for m in range(i, N + 1)]
        for J in nodes:
            lhs = mp.fsum(PsiS[s] for s in PsiS if ends(s)[0] == J) - mp.fsum(PsiS[s] for s in PsiS if ends(s)[1] == J)
            Q = 0
            if J[0] == J[1]:
                p = J[0]
                n = scal(S['P'][p], S['x'][p])
                lhs += -mp.fsum(a[k] * n[k] / (a[k] - t) for k in range(N))
                Q = -cols[(1, 2)]['VT'] if p == 1 else (S['R'] if p == N else 0)
            rhs = Q + ((1 - D.q - mp.fsum(a[k] * D.z[k] / (a[k] - t) for k in range(N))) if J == (1, N) else 0)
            worst = max(worst, abs(lhs - rhs))
    check('D6 flux bookkeeping identity at every unit (Lemma 5.2 (e))', worst < mp.mpf(10) ** (-(mp.mp.dps - 15)),
          f'({mp.nstr(worst, 3)})')


def run_limit_checks(S, Sys, K=400):
    D, cols = S['D'], S['cols']
    N, a = D.N, D.a
    tiny = mp.mpf(10) ** (-(mp.mp.dps - 15))
    # D7 spectra, pinch rule, pivots, trace signs, kappa
    ok, okpiv, okH7 = True, True, True
    minrat, rho = mp.mpf(1), mp.mpf(0)
    for si, s in enumerate(Sys.secs):
        V, L, w = section_flows(S, s)
        M, lam, left, right, imag = eigdata(a, V, w)
        mt = Sys.meta[si]
        ok &= imag < tiny and all(t > 0 for t in lam) and all(lam[r] > lam[r + 1] for r in range(N - 1))
        rho = max(rho, max(lam[r + 1] / lam[r] for r in range(N - 1)))
        # pinch rule and traces outside the major spectrum
        maj = [r for r in range(N) if mt['cls'][r][0] == 'major']
        ok &= maj == list(range(maj[0], maj[-1] + 1))
        bot, top = ends(s)
        xb, xt = node_comp(S, bot), node_comp(S, top)
        # pivots in the sum-normalized right eigenbasis
        Rm = mp.matrix(N, N)
        for r in range(N):
            sr = mp.fsum(right[r])
            for k in range(N):
                Rm[k, r] = right[r][k] / sr
        cb = mp.lu_solve(Rm, mp.matrix(xb))
        ct = mp.lu_solve(Rm, mp.matrix(xt))
        rs = mt['rstar']
        okpiv &= cb[rs] > 0 and ct[rs] > 0
        minrat = min(minrat, cb[rs], ct[rs])
        # limit relations: closed coefficients vanish
        for r in range(N):
            if r < rs:
                ok &= abs(cb[r]) < tiny
            elif r > rs:
                ok &= abs(ct[r]) < tiny
        # H7 trace signs
        for r in range(N):
            if mt['cls'][r][0] != 'trace':
                continue
            k = mt['cls'][r][1]
            J, xJ = (bot, xb) if r < rs else (top, xt)
            Pk = Psi(a, V, w, a[k], skip=k)
            Gk = Gx(a, xJ, a[k], skip=k)
            supp = [j for j in range(N) if xJ[j] != 0]
            light = r < rs
            okH7 &= Pk > 0 and xJ[k] == 0 and (all(j > k for j in supp) if light else all(j < k for j in supp))
            okH7 &= (Gk < 0) if light else (Gk > 0)
    check('D7 spectra simple, positive, pinch = smallest major (T) / largest major (B), traces outside, limit relations',
          ok, f'(max rho {mp.nstr(rho, 6)})')
    check('D7 pivots c_r*(x_bot), c_r*(x_top) > 0 (H8)', okpiv, f'(min {mp.nstr(minrat, 4)})')
    check('D7 trace signs: Psi^(k)(a_k) > 0, x_Jk = 0, k on the correct side of supp x_J (H7)', okH7)
    # D8 half profiles, computed through the chart vectors (H6): M^k x_bot / (lambda_r*^k) = sum over the OPEN modes
    # (closed coefficients are exactly 0 by the limit relations; iterating M^-1 in floating point would amplify their
    # rounding noise by (lambda_r*/lambda_r)^k, a numerical artifact).  Also: the closed coefficients are < tiny.
    ok = True
    worst = mp.mpf(1)
    for si, s in enumerate(Sys.secs):
        V, L, w = section_flows(S, s)
        M, lam, left, right, _ = eigdata(a, V, w)
        rs = Sys.meta[si]['rstar']
        Rm = mp.matrix(N, N)
        eh = []
        for r in range(N):
            sr = mp.fsum(right[r])
            eh.append([t / sr for t in right[r]])
            for k in range(N):
                Rm[k, r] = eh[r][k]
        bot, top = ends(s)
        cb = mp.lu_solve(Rm, mp.matrix(node_comp(S, bot)))
        ct = mp.lu_solve(Rm, mp.matrix(node_comp(S, top)))
        ok &= all(abs(cb[r]) < tiny for r in range(rs)) and all(abs(ct[r]) < tiny for r in range(rs + 1, N))
        for k in list(range(0, 60)) + list(range(60, K, 10)):
            vb = [mp.fsum(cb[r] * (lam[r] / lam[rs]) ** k * eh[r][j] for r in range(rs, N)) for j in range(N)]
            vt = [mp.fsum(ct[r] * (lam[rs] / lam[r]) ** k * eh[r][j] for r in range(0, rs + 1)) for j in range(N)]
            mb, mt_ = min(vb), min(vt)
            worst = min(worst, mb / mp.fsum(vb), mt_ / mp.fsum(vt))
            if mb < -tiny or mt_ < -tiny or mp.fsum(vb) <= 0 or mp.fsum(vt) <= 0:
                ok = False
    check(f'D8 half-profiles (chart vectors) >= 0 with positive sum for k <= {K} (H6)', ok,
          f'(min normalized entry {mp.nstr(worst, 3)})')
    # D9 purities
    ok = True
    pat = S['pat']
    for p in range(1, N + 1):
        imp = (cols[(p - 1, p)]['eta'] if p > 1 and pat[p - 2] == 'B' else 0) + \
              (cols[(p, p + 1)]['eta'] if p < N and pat[p - 1] == 'T' else 0)
        ok &= abs((1 - S['x'][p][p - 1]) - imp / S['P'][p]) < tiny
    check('D9 product impurities (eta_{p-1}[B] + eta_p[T]) / P_p', ok)
    return rho


def theorem_L_checks(S, Sys):
    N = Sys.N
    a = Sys.a
    tiny = mp.mpf(10) ** (-(mp.mp.dps - 15))
    t0 = time.time()
    G0 = Sys.eval(Sys.v0)
    check('L0 coordinates of Table 2 reproduce the limit configuration', Sys.recon_err < tiny, f'({mp.nstr(Sys.recon_err, 3)})')
    check('L1 G_inf(v_inf) = 0', max(abs(g) for g in G0) < tiny, f'(max {mp.nstr(max(abs(g) for g in G0), 3)})')
    J = Sys.jacobian(Sys.v0)
    print(f'  Jacobian {J.rows}x{J.cols} in {time.time() - t0:.1f} s')
    fd = mp.mpf(10) ** (-(2 * (mp.mp.dps // 3)) + 8)       # finite-difference noise level
    rows, coords = Sys.rows, Sys.coords
    # L2 zero blocks
    z1 = z2 = mp.mpf(0)
    for ri, r in enumerate(rows):
        if r['class'] == 'maj':
            continue
        for cj, c in enumerate(coords):
            if c[3] == 'maj':
                z1 = max(z1, abs(J[ri, cj]))
            elif c[3] != r['class']:
                z2 = max(z2, abs(J[ri, cj]))
    scale = max(abs(J[i, j]) for i in range(J.rows) for j in range(J.cols))
    check('L2 d(trace rows)/d(major coords) = 0 and d(trace rows of k)/d(trace coords of k\' != k) = 0',
          max(z1, z2) < fd * scale, f'({mp.nstr(z1, 3)}, {mp.nstr(z2, 3)}; scale {mp.nstr(scale, 3)})')
    # Lemma 4.3 prediction: (dx_{J,k} + g dw_{sigma,k}) / pivot
    def lin_part(fun):
        base = fun(Sys.v0)
        out = []
        for j in range(len(Sys.v0)):
            vp = list(Sys.v0)
            vp[j] += 1
            out.append(fun(vp) - base)
        return out
    worst = mp.mpf(0)
    for ri, r in enumerate(rows):
        if r['class'] == 'maj':
            continue
        k0 = r['class'] - 1
        s = r['sec']
        J_ = r['closed']
        V, L, w = section_flows(S, s)
        xJ = node_comp(S, J_)
        g = Gx(a, xJ, a[k0], skip=k0) / Psi(a, V, w, a[k0], skip=k0)
        def xk(v, J_=J_, k0=k0):
            st = state_from_coords(S, coords, v)
            if J_[0] == J_[1]:
                p = J_[0]
                n = (st[(p - 1, p)]['b'][k0] if p > 1 else 0) + (st[(p, p + 1)]['d'][k0] if p < N else 0)
                return n / S['P'][p]
            return st[J_]['X'][k0]
        def wk(v, s=s, k0=k0):
            st = state_from_coords(S, coords, v)
            return st[(s[1], s[2])]['d'][k0] if s[0] == 'T' else -st[(s[1], s[2])]['b'][k0]
        dx, dw = lin_part(xk), lin_part(wk)
        si = r['si']
        _, lam, left, right, _ = eigdata(a, V, w)
        mt = Sys.meta[si]
        lr = [t / left[mt['rstar']][mt['jn'][mt['rstar']]] for t in left[mt['rstar']]]
        piv = dot(lr, xJ)
        for cj in range(len(coords)):
            pred = (dx[cj] + g * dw[cj]) / piv
            worst = max(worst, abs(J[ri, cj] - pred))
    check('L2 every trace row equals Lemma 4.3\'s prediction (dx_Jk - kappa dq)/pivot', worst < fd * scale,
          f'({mp.nstr(worst, 3)})')
    # L3 major block: dependence only on own column and ancestors; diagonal blocks nonsingular
    colsN = column_order(N)
    anc = lambda c, c2: c2[0] <= c[0] and c[1] <= c2[1]
    zt = mp.mpf(0)
    detmin = None
    for c in colsN:
        mr = [ri for ri, r in enumerate(rows) if r['class'] == 'maj' and r['col'] == c]
        mc = [cj for cj, cc in enumerate(coords) if cc[3] == 'maj' and cc[0] == c]
        for ri in mr:
            for cj, cc in enumerate(coords):
                if cc[3] == 'maj' and not anc(c, cc[0]):
                    zt = max(zt, abs(J[ri, cj]))
        assert len(mr) == len(mc), (c, len(mr), len(mc))
        B = mp.matrix(len(mr), len(mc))
        for x_, ri in enumerate(mr):
            for y_, cj in enumerate(mc):
                B[x_, y_] = J[ri, cj]
        sv = min(mp.svd_r(B, compute_uv=False))
        detmin = sv if detmin is None else min(detmin, sv)
    check('L3 major rows of c depend only on major coords of c and its ancestors; diagonal blocks nonsingular',
          zt < fd * scale and detmin > mp.mpf(10) ** -20, f'(off {mp.nstr(zt, 3)}, min sigma of blocks {mp.nstr(detmin, 4)})')
    # L4 trace blocks and A_k
    ok = True
    smin_T = None
    for k in range(1, N + 1):
        tr = [ri for ri, r in enumerate(rows) if r['class'] == k]
        tc = [cj for cj, cc in enumerate(coords) if cc[3] == k]
        if len(tr) != len(tc):
            ok = False
            continue
        if not tr:          # component k is major at every node (e.g. k = 2 for N = 3, pattern TB): empty block
            print(f'  (component {k}: no trace rows and no trace coordinates)')
            continue
        B = mp.matrix(len(tr), len(tc))
        for x_, ri in enumerate(tr):
            for y_, cj in enumerate(tc):
                B[x_, y_] = J[ri, cj]
        sv = min(mp.svd_r(B, compute_uv=False))
        smin_T = sv if smin_T is None else min(smin_T, sv)
        # A_k
        nodes = [J_ for J_ in [(i, m) for i in range(1, N + 1) for m in range(i, N + 1)] if node_comp(S, J_)[k - 1] == 0]
        lighter = sorted([J_ for J_ in nodes if J_[0] > k], key=h_)
        heavier = sorted([J_ for J_ in nodes if J_[1] < k], key=lambda J_: -h_(J_))
        ok &= len(lighter) + len(heavier) == len(nodes)
        order = lighter + heavier
        pos = {J_: n for n, J_ in enumerate(order)}
        A = mp.matrix(len(order), len(order))
        for J_ in order:
            if J_[0] == J_[1]:
                A[pos[J_], pos[J_]] += S['P'][J_[0]]
        for ri in tr:
            r = rows[ri]
            s = r['sec']
            V, L, w = section_flows(S, s)
            xJ = node_comp(S, r['closed'])
            kap = abs(Gx(a, xJ, a[k - 1], skip=k - 1)) / Psi(a, V, w, a[k - 1], skip=k - 1)
            ok &= kap > 0 and r['closed'] in pos
            A[pos[r['closed']], pos[r['closed']]] += 1 / kap
            if r['open'] in pos:
                A[pos[r['open']], pos[r['closed']]] -= 1 / kap
        for x_ in range(A.rows):
            ok &= A[x_, x_] > 0
            for y_ in range(x_):
                ok &= A[x_, y_] == 0
            for y_ in range(x_ + 1, A.cols):
                ok &= A[x_, y_] <= 0
    check('L4 trace blocks T_k square and nonsingular; A_k upper triangular, positive diagonal, Z-matrix', ok and smin_T > 0,
          f'(min sigma(T_k) {mp.nstr(smin_T, 4)})')
    sv = mp.svd_r(J, compute_uv=False)
    smin = min(sv)
    check('L5 DG_inf(v_inf) invertible', smin > mp.mpf(10) ** -25, f'(sigma_min {mp.nstr(smin, 4)}, sigma_max {mp.nstr(max(sv), 4)})')
    return J


def finite_n_one(S, Sys, Jinv, n, prev, ratio):
    N, a = Sys.N, Sys.a
    if True:
        nn = {s: n for s in Sys.secs}
        v = list(Sys.v0)
        it = 0
        res_prev = None
        nref = 0
        res = mp.inf
        try:
            for it in range(400):
                g = Sys.eval(v, nn)
                res = max(abs(t) for t in g)
                if res < mp.mpf(10) ** (-(mp.mp.dps - 10)):
                    break
                # chord step with DG_inf(v_inf)^{-1} (Corollary 4.7); if the contraction is slow (factor > 0.1),
                # refresh the matrix with the finite-difference Jacobian of G_n at the current iterate (Newton)
                if res_prev is not None and res > res_prev / 10 and nref < 12:
                    Jinv = Sys.jacobian(v, nn) ** -1
                    nref += 1
                res_prev = res
                dv = Jinv * mp.matrix(g)
                v = [v[j] - dv[j] for j in range(len(v))]
        except (ZeroDivisionError, ValueError, AssertionError) as exc:      # iteration left the domain of G_n
            print(f'  [INFO] F1 n = {n}: iteration left the domain ({type(exc).__name__}) after {it} it.: n below this '
                  f"case's n_0 (not a claim of the theorem)")
            return prev
        conv = res < mp.mpf(10) ** (-(mp.mp.dps - 12))
        dist = max(abs(v[j] - Sys.v0[j]) for j in range(len(v)))
        # reconstruct and check the physical equations
        _, st, X, x, cache = Sys.eval(v, nn, want=True)
        comp = lambda J_: x[J_[0]] if J_[0] == J_[1] else X[J_]
        y = lambda xx: [a[i] * xx[i] / mp.fsum(a[j] * xx[j] for j in range(N)) for i in range(N)]
        prof = {}
        minx, minS, topgap = mp.mpf(1), mp.mpf(10), mp.mpf(0)
        where = None
        for s in Sys.secs:
            M, lam, left, right, xb, xt, w, V, L = cache[s]
            xs = [xb]
            for j in range(n + 1):
                Sx = mp.fsum(a[i] * xs[-1][i] for i in range(N))
                minS = min(minS, Sx)
                nxt = M * mp.matrix(xs[-1])
                xs.append([nxt[i] / (L * Sx) for i in range(N)])
            topgap = max(topgap, max(abs(xs[-1][i] - xt[i]) for i in range(N)))
            xs[-1] = list(xt)                       # use the actual node composition at the top
            prof[s] = xs
            for j_, t_ in enumerate(xs):
                for k_ in range(N):
                    if t_[k_] < minx:
                        minx, where = t_[k_], (s, j_, k_ + 1)
        bal = mp.mpf(0)
        cols = S['cols']
        for s in Sys.secs:
            V, L, _ = section_flows(S, s)
            xs = prof[s]
            for j in range(1, n + 1):
                yin, yout = y(xs[j - 1]), y(xs[j])
                bal = max(bal, max(abs(V * yin[i] + L * xs[j + 1][i] - V * yout[i] - L * xs[j][i]) for i in range(N)))
        for (i, m), c in cols.items():
            Xc = X[(i, m)]
            ii = [c['LT'] * prof[('T', i, m)][1][k] + c['VB'] * y(prof[('B', i, m)][n])[k] for k in range(N)]
            if m < N:
                cp = cols[(i, m + 1)]
                ii = [ii[k] + cp['VT'] * y(prof[('T', i, m + 1)][n])[k] for k in range(N)]
            if i > 1:
                cp = cols[(i - 1, m)]
                ii = [ii[k] + cp['LB'] * prof[('B', i - 1, m)][1][k] for k in range(N)]
            if (i, m) == (1, N):
                ii = [ii[k] + S['D'].z[k] for k in range(N)]
            Lout = c['LB'] + (cols[(i, m + 1)]['LT'] if m < N else 0)
            Vout = c['VT'] + (cols[(i - 1, m)]['VB'] if i > 1 else 0)
            yX = y(Xc)
            bal = max(bal, max(abs(Lout * Xc[k] + Vout * yX[k] - ii[k]) for k in range(N)))
        for p in range(1, N + 1):
            xp = x[p]
            ii, Vout, Lout = vec(N), mp.mpf(0), S['P'][p]
            if p > 1:
                cB = cols[(p - 1, p)]
                ii = add(ii, scal(cB['LB'], prof[('B', p - 1, p)][1]))
                Vout += cB['VB']
            if p < N:
                cT = cols[(p, p + 1)]
                ii = add(ii, scal(cT['VT'], y(prof[('T', p, p + 1)][n])))
                Lout += cT['LT']
            yp = y(xp)
            bal = max(bal, max(abs(Lout * xp[k] + Vout * yp[k] - ii[k]) for k in range(N)))
            if p == 1:      # total condenser: no vapor outflow; composition rule x_1 = y(top stage of T[1,2])
                bal = max(bal, max(abs(xp[k] - y(prof[('T', 1, 2)][n])[k]) for k in range(N)))
        for p in range(1, N + 1):
            minx = min(minx, min(x[p]))
        for c in X:
            minx = min(minx, min(X[c]))
        imp = [1 - x[p][p - 1] for p in range(1, N + 1)]
        tol = mp.mpf(10) ** -40     # dps was raised by n*log10(max mode ratio) + 60 to absorb the forward amplification
        if not conv:
            # Corollary 4.7 claims a zero near v_inf only for n >= n_0 (not explicit): report, do not count as a FAIL
            print(f'  [INFO] F1 n = {n}: no convergence from v_inf ({it} it., residual {mp.nstr(res, 3)}): n below this '
                  f"case's n_0 (not a claim of the theorem)")
            return prev
        ok = bal < tol and minx > 0 and minS > 0 and topgap < tol
        if prev is not None:          # nlist is processed in decreasing n: |v_n - v_inf| must grow
            ok &= dist > prev
        if minx <= 0:
            s_, j_, k_ = where
            print(f'  [diag] smallest mole fraction {mp.nstr(minx, 5)} at section {s_}, stage {j_}, component {k_}; '
                  f'net flow of k there {mp.nstr(cache[s_][6][k_ - 1], 5)}; working precision 1e-{mp.mp.dps}')
        check(f'F1 n = {n}: chord/Newton ({it} it., {nref} Jacobian refreshes), unit balances, x_top reached, all x > 0', ok,
              f'(|v_n - v_inf| {mp.nstr(dist, 3)}, residual {mp.nstr(res, 3)}, max balance {mp.nstr(bal, 3)}, '
              f'min x {mp.nstr(minx, 3)}, min S {mp.nstr(minS, 3)}, top gap {mp.nstr(topgap, 3)}, dps {mp.mp.dps}, '
              f'impurities {[mp.nstr(t, 3) for t in imp]})')
        return dist


# ------------------------------------------------------------------------------------------------------------
def make_case(case, halvings=None):
    """returns (Data, delta, eta_nb, frac, pattern, eta_bin, halvings); all inputs are exact binary floats or
    integers, so the same case can be re-instantiated at a higher working precision."""
    parts = case.split(':')
    frac = mp.mpf('0.25')
    if parts[0] == 'n3':
        D = Data([4, 2, 1], [1, 1, 1], 1)
        pat, delta, eta_nb = parts[1], mp.mpf('0.05'), mp.mpf('0.025')
        eta_bin = {1: mp.mpf('0.003'), 2: mp.mpf('0.003')}
        return D, delta, eta_nb, frac, pat, eta_bin, 0
    if parts[0] == 'n4':
        D = Data([8, 4, 2, 1], [1, 1, 1, 1], 1)
        return (D, mp.mpf('0.05'), mp.mpf('0.005'), frac, parts[1],
                {1: mp.mpf('0.003'), 2: mp.mpf('0.003'), 3: mp.mpf('0.003')}, 0)
    if parts[0] in ('rand', 'near'):
        # 'near': as 'rand' but one consecutive volatility ratio is 1.05 (nearly equal volatilities)
        N, seed = int(parts[1]), int(parts[2])
        rnd = random.Random(seed)
        pat = parts[3] if len(parts) > 3 and parts[3] else ''.join(rnd.choice('BT') for _ in range(N - 1))
        if len(parts) > 4:
            frac = mp.mpf(parts[4])
        ratios = [rnd.uniform(1.3, 2.5) for _ in range(N - 1)]
        if parts[0] == 'near':
            ratios[rnd.randrange(N - 1)] = 1.05
        al = [mp.mpf(1)]
        for r_ in ratios:
            al.append(al[-1] * mp.mpf(r_))
        al = al[::-1]
        z = [mp.mpf(rnd.uniform(0.2, 1.0)) for _ in range(N)]
        q = mp.mpf(rnd.uniform(-0.5, 1.5))
        D = Data(al, z, q)
        delta = mp.mpf('0.05') * D.Hs
        eta = delta / 10
        if halvings is not None:
            return D, delta, eta / 2 ** halvings, frac, pat, None, halvings
        hv = 0
        while True:
            S = build(D, delta, eta, frac, pat)
            okv = all(all(c['d'][k] > 0 for k in range(c['i'] - 1, c['m'] - 1)) and
                      all(c['b'][k] > 0 for k in range(c['i'], c['m'])) for c in S['cols'].values() if c['kind'] == 'nb')
            okk = all(S['cols'][(p, p + 1)]['kappa'] >= D.e[p - 1] + delta / 2 for p in range(1, N))
            if okv and okk:
                break
            eta /= 2
            hv += 1
        return D, delta, eta, frac, pat, None, hv
    raise SystemExit('unknown case')


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else 'n3:BT'
    mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    nlist = sorted([int(t) for t in sys.argv[3].split(',')], reverse=True) if len(sys.argv) > 3 and sys.argv[3] else []
    extra = int(sys.argv[4]) if len(sys.argv) > 4 else 30      # extra digits for the finite-n reconstruction
    t0 = time.time()
    D, delta, eta_nb, frac, pat, eta_bin, hv = make_case(case)
    print(f'case {case}: N = {D.N}, alpha = {[mp.nstr(t, 5) for t in D.a]}, z = {[mp.nstr(t, 4) for t in D.z]}, '
          f'q = {mp.nstr(D.q, 4)}, pattern {pat}, delta = {mp.nstr(delta, 4)}, eta_nb = {mp.nstr(eta_nb, 4)}, dps {mp.mp.dps}')
    S = build(D, delta, eta_nb, frac, pat, eta_bin=eta_bin)
    limit_config(S)
    print(f'  H = {[mp.nstr(t, 6) for t in D.H]}, H* = {mp.nstr(D.Hs, 8)}, eta_bin = '
          f'{[mp.nstr(S["cols"][(p, p + 1)]["eta"], 4) for p in range(1, D.N)]}')
    run_design_checks(S, eta_nb)
    Sys = System(S)
    run_limit_checks(S, Sys)
    J = theorem_L_checks(S, Sys)
    if nlist:
        Jinv = J ** -1
        ratio = mp.mpf(1)
        for s in Sys.secs:
            V, L, w = section_flows(S, s)
            lam = eigdata(D.a, V, w)[1]
            ratio = max(ratio, lam[0] / lam[-1])
        prev = None
        base = mp.mp.dps
        for n in nlist:
            # the forward recursion amplifies rounding errors by at most ratio^(n+1): raise the precision accordingly
            dpsF = base + int((n + 1) * mp.log10(ratio)) + extra
            with mp.workdps(dpsF):
                D2, delta2, eta2, frac2, pat2, eb2, _ = make_case(case, halvings=hv)
                S2 = build(D2, delta2, eta2, frac2, pat2, eta_bin=eb2)
                limit_config(S2)
                Sys2 = System(S2)
                prev = finite_n_one(S2, Sys2, Jinv, n, prev, ratio)
    print(f'{"ALL PASS" if not FAIL else "FAIL: " + str(FAIL)}  ({time.time() - t0:.0f} s)')
    sys.exit(0 if not FAIL else 1)


if __name__ == '__main__':
    main()
