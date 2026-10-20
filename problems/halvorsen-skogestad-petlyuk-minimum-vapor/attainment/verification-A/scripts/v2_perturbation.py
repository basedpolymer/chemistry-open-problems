#!/usr/bin/env python3
"""V2 (independent verifier): end-to-end numerical check of the perturbation step and of the reconstruction
(attainment/theorem.md §§2-3, 6; attainment/analysis/analysis.md §§1-5).

NUMERICAL EVIDENCE ONLY.  The proofs are re-derived by hand in verification.md.  Written from scratch by V2: no code
of D, L, A or C is imported.  D's exported configurations (attainment/design/runs/config_*.json) are only READ, to
cross-check V2's own implementation of the design (theorem.md Definition 4.7 and §4.4).

What is checked, per case:
  (D)  V2's own implementation of F(delta, eta, pi) (Definition 4.7: preferred split relative to the actual feed,
       Astra's slack direction h, B/T pinning, vapors from the reboiler up) and of the limit configuration
       (Definition 4.12: interpolation vectors X_c, forced pinches) agrees with D's export (if one exists);
       the margins kappa_p, the validity (V), Lemma 4.9 (middle feed roots are global).
  (P)  the network map G_n of (3.3) built with EXPLICIT eigen-data: modes = roots of the polynomial P_w(t),
       Underwood left vectors l_phi = (alpha_i/(alpha_i - phi)) for major modes and the regularized trace vectors
       l~ of Lemma 3.4 / (3.1) for trace modes (not generic numerical eigenvectors); pinch = smallest major mode
       (top sections) or largest major mode (bottom sections), as in Table 1;
       G_inf(v_inf) = 0, pivots c_{r*} > 0 at both ends, half-profiles (H6) >= 0, sigma_min(DG_inf(v_inf)).
  (C1) the C^1 estimate of Lemma 5.2 / Lemma 3.5 (c): at random points v in a ball around v_inf,
       |G_n - G_inf| / rho_bar^(n+1) and ||D(G_n - G_inf)|| / ((n+2) rho_bar^(n+1)) stay bounded as n grows
       (n = 8 ... 128), rho_bar = max consecutive mode ratio over the sample.
  (N)  Newton zeros v_n of G_n for several n and unequal lengths; |v_n - v_inf| / rho_inf^(n+1) bounded;
       Lemma 6.1's bound |v_n - v_inf| <= 2 ||A^-1|| |G_n(v_inf)| (A = DG_inf(v_inf), sup norms);
       the CONTRACTION Phi(v) = v - A^-1 G_n(v) of Lemma 6.1, iterated from random starting points in a ball,
       converges to the same v_n (local uniqueness).
  (S)  physical reconstruction: every section is re-shot with the PHYSICAL cross-section balance
       x_{k+1} = (V y(x_k) - w)/L from its bottom-end composition (no eigen-data); the end x_{n+1} equals the top-end
       composition; every stage, junction, side-draw, condenser and reboiler component balance is recomputed from the
       actual liquid and vapour streams; S > 0 and every mole fraction > 0 at every unit; the purities; the duties.
Usage: python v2_perturbation.py [case ...]   cases: n3B n3BT n3q n3TT n4D n4Dlong n5D (default: n3B n3BT n3q n4D)
Environment: Python 3 + mpmath (WSL Fedora sci env), one thread.  Runtime: N = 3 about 2 min per case, N = 4 about
15 min, N = 5 about 1 h.  Prints PASS/FAIL per claim; exit code 1 on any FAIL.
"""
import sys
import os
import json
import time
import random

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
DRUNS = os.path.normpath(os.path.join(HERE, '..', '..', 'design', 'runs'))
FAILS = []


def check(name, ok, info=''):
    print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)
    if not ok:
        FAILS.append(name)


def nstr(x, k=4):
    return mp.nstr(x, k)


# ================================================================================================ design (V2's own)

def root_monotone(fun, a, b, increasing):
    """Root of a strictly monotone function on (a, b) which runs from -inf/+inf to +inf/-inf (bisection)."""
    lo, hi = mp.mpf(a), mp.mpf(b)
    for _ in range(int(mp.mp.prec * 1.05) + 20):
        mid = (lo + hi) / 2
        v = fun(mid)
        if (v < 0) == increasing:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def prod(it):
    r = mp.mpf(1)
    for x in it:
        r *= x
    return r


def Psi(V, w, al, t, skip=None):
    return V - mp.fsum(al[i] * w[i] / (al[i] - t) for i in range(len(al)) if i != skip and w[i] != 0)


def design(al, z, q, delta, eta_nb, eta_bin, pattern):
    """theorem.md Definition 4.7 and Definition 4.12, components 0-based, columns (i, m) 1-based."""
    N = len(al)
    A = lambda p: al[p - 1]                     # 1-based volatility
    omq = 1 - q
    theta = {}
    for j in range(1, N):
        fun = lambda t: omq - mp.fsum(al[p] * z[p] / (al[p] - t) for p in range(N))
        theta[j] = root_monotone(fun, A(j + 1), A(j), increasing=False)
    H = {j: mp.fsum(al[p - 1] * z[p - 1] / (theta[j] - al[p - 1]) for p in range(j + 1, N + 1)) for j in range(1, N)}
    Hs = max(H.values())
    R = Hs + delta
    col = {}

    def feed_of(c):
        i, m = c
        if c == (1, N):
            return [mp.mpf(x) for x in z], omq
        f = [mp.mpf(0)] * N
        U = mp.mpf(0)
        if m < N:
            f = [a + b for a, b in zip(f, col[(i, m + 1)]['d'])]
            U += col[(i, m + 1)]['VT']
        if i > 1:
            f = [a + b for a, b in zip(f, col[(i - 1, m)]['b'])]
            U -= col[(i - 1, m)]['VB']
        return f, U

    def column_roots(f, U, i, m):
        out = {}
        for p in range(i, m):
            fun = lambda t: U - mp.fsum(al[k - 1] * f[k - 1] / (al[k - 1] - t) for k in range(i, m + 1))
            out[p] = root_monotone(fun, A(p + 1), A(p), increasing=False)
        return out

    valid = True
    for Lc in range(N - 1, 1, -1):
        for i in range(1, N - Lc + 1):
            m = i + Lc
            c = (i, m)
            f, U = feed_of(c)
            th = column_roots(f, U, i, m)
            rho = {p: prod(A(p) - th[r] for r in range(i, m)) / prod(A(p) - A(pp) for pp in range(i, m) if pp != p)
                   for p in range(i, m)}
            VT0 = A(i) * f[i - 1] / rho[i]
            d0 = [mp.mpf(0)] * N
            for p in range(i, m):
                d0[p - 1] = VT0 * rho[p] / A(p)
            VB0 = VT0 - U
            h = [mp.mpf(0)] * N
            if m - i >= 3:
                beta = prod(A(p) for p in range(i + 1, m)) / prod(th[r] for r in range(i + 1, m - 1))
                for p in range(i + 1, m):
                    h[p - 1] = ((A(p) - beta) * prod(A(p) - th[r] for r in range(i + 1, m - 1)) /
                                (A(p) * prod(A(p) - A(pp) for pp in range(i + 1, m) if pp != p)))
            eta = eta_nb[c]
            d = [d0[k] + eta * h[k] for k in range(N)]
            b = [f[k] - d[k] for k in range(N)]
            VT, VB = VT0 + eta, VB0 + eta
            LT, LB = VT - mp.fsum(d), VB + mp.fsum(b)
            ok = all(d[p - 1] > 0 for p in range(i, m)) and all(b[p - 1] > 0 for p in range(i + 1, m + 1))
            valid &= ok and min(VT, VB, LT, LB) > 0
            col[c] = dict(kind='nb', f=f, U=U, th=th, d=d, b=b, VT=VT, VB=VB, LT=LT, LB=LB, eta=eta, h=h)
    # binary columns: feeds and U's
    for p in range(1, N):
        f, U = feed_of((p, p + 1))
        col[(p, p + 1)] = dict(kind='bin', f=f, U=U)
    # vapors from the reboiler up
    col[(N - 1, N)]['VB'] = R
    for p in range(N - 1, 0, -1):
        c = col[(p, p + 1)]
        c['VT'] = c['VB'] + c['U']
        if p >= 2:
            col[(p - 1, p)]['VB'] = c['VT']
    kappa = {}
    for p in range(1, N):
        c = col[(p, p + 1)]
        f, eta = c['f'], eta_bin[p]
        fp, fq = f[p - 1], f[p]
        assert all(f[k] == 0 for k in range(N) if k not in (p - 1, p)), 'binary feed not on {p, p+1}'
        thp = root_monotone(lambda t: c['U'] - al[p - 1] * fp / (al[p - 1] - t) - al[p] * fq / (al[p] - t),
                            A(p + 1), A(p), increasing=False)
        kappa[p] = c['VT'] - al[p - 1] * fp / (al[p - 1] - thp)
        d = [mp.mpf(0)] * N
        b = [mp.mpf(0)] * N
        if pattern[p - 1] == 'B':
            d[p - 1] = fp - eta
            b[p - 1], b[p] = eta, fq
        else:
            d[p - 1], d[p] = fp, eta
            b[p] = fq - eta
        c.update(d=d, b=b, LT=c['VT'] - mp.fsum(d), LB=c['VB'] + mp.fsum(b), eta=eta, pi=pattern[p - 1], thp=thp)
        valid &= min(c['VT'], c['VB'], c['LT'], c['LB']) > 0 and kappa[p] > 0 and 0 < eta < min(fp, fq)
    # products
    prods = {}
    for p in range(1, N + 1):
        n = [mp.mpf(0)] * N
        if p > 1:
            n = [a + b for a, b in zip(n, col[(p - 1, p)]['b'])]
        if p < N:
            n = [a + b for a, b in zip(n, col[(p, p + 1)]['d'])]
        P = mp.fsum(n)
        prods[p] = dict(n=n, P=P, x=[v / P for v in n])
    # limit junction compositions (Definition 4.12)
    X = {}
    for c, cd in col.items():
        i, m = c
        if cd['kind'] == 'nb':
            wT = cd['d']
            wB = [-x for x in cd['b']]
            rc = root_monotone(lambda t: Psi(cd['VT'], wT, al, t), A(i + 1), A(i), increasing=False)
            sc = root_monotone(lambda t: Psi(cd['VB'], wB, al, t), A(m), A(m - 1), increasing=True)
            chis = [rc] + [cd['th'][r] for r in range(i + 1, m - 1)] + [sc]
            x = [mp.mpf(0)] * N
            for p in range(i, m + 1):
                x[p - 1] = prod(A(p) - ch for ch in chis) / (A(p) * prod(A(p) - A(pp) for pp in range(i, m + 1) if pp != p))
            s = mp.fsum(x)
            X[c] = [v / s for v in x]
            cd.update(rc=rc, sc=sc)
        else:
            p = i
            if cd['pi'] == 'B':
                wB = [-x for x in cd['b']]
                hi = A(p) + 1
                while Psi(cd['VB'], wB, al, hi) < 0:
                    hi = A(p) + 2 * (hi - A(p))
                ph = root_monotone(lambda t: Psi(cd['VB'], wB, al, t), A(p), hi, increasing=True)
                x = [mp.mpf(0)] * N
                x[p - 1] = cd['b'][p - 1] / (ph - A(p))
                x[p] = cd['b'][p] / (ph - A(p + 1))
            else:
                wT = cd['d']
                ph = root_monotone(lambda t: Psi(cd['VT'], wT, al, t), mp.mpf(0), A(p + 1), increasing=False)
                x = [mp.mpf(0)] * N
                x[p - 1] = cd['d'][p - 1] / (A(p) - ph)
                x[p] = cd['d'][p] / (A(p + 1) - ph)
            s = mp.fsum(x)
            X[c] = [v / s for v in x]
            cd['forced_pinch'] = ph
    # sections
    secs = {}
    for c, cd in col.items():
        i, m = c
        secs[('T',) + c] = dict(V=cd['VT'], L=cd['LT'], w=list(cd['d']), bot=('J', i, m), top=('J', i, m - 1), col=c)
        secs[('B',) + c] = dict(V=cd['VB'], L=cd['LB'], w=[-x for x in cd['b']], bot=('J', i + 1, m), top=('J', i, m),
                                col=c)
    return dict(N=N, al=al, z=z, q=q, theta=theta, H=H, Hs=Hs, R=R, col=col, prods=prods, X=X, secs=secs,
                kappa=kappa, valid=valid, pattern=pattern)


def compare_with_D(Dz, fname):
    path = os.path.join(DRUNS, fname)
    if not os.path.exists(path):
        print('     (no export %s)' % fname)
        return
    with open(path) as fh:
        E = json.load(fh)
    worst = mp.mpf(0)
    for key, s in Dz['secs'].items():
        name = '%s[%d,%d]' % key
        es = E['sections'][name]
        worst = max(worst, abs(mp.mpf(es['V']) - s['V']), abs(mp.mpf(es['L']) - s['L']),
                    max(abs(mp.mpf(a) - b) for a, b in zip(es['w'], s['w'])))
    for c, x in Dz['X'].items():
        ex = E['nodes']['J[%d,%d]' % c]
        worst = max(worst, max(abs(mp.mpf(a) - b) for a, b in zip(ex, x)))
    for p, pr in Dz['prods'].items():
        worst = max(worst, abs(mp.mpf(E['products'][str(p)]['P']) - pr['P']),
                    max(abs(mp.mpf(a) - b) for a, b in zip(E['products'][str(p)]['x'], pr['x'])))
    check('(D) V2 design and limit configuration = D export %s' % fname, worst < mp.mpf('1e-45'),
          'max |difference| = %s' % nstr(worst, 2))


# ================================================================================================ network map

class Net:
    def __init__(self, Dz):
        self.D = Dz
        self.N = Dz['N']
        self.al = Dz['al']
        self.cols = sorted(Dz['col'].keys(), key=lambda c: (-(c[1] - c[0]), c[0]))
        self.dim = 2 * (self.N - 1) * len(self.cols)
        self.d0 = {c: list(Dz['col'][c]['d']) for c in self.cols}
        self.X0 = {c: list(Dz['X'][c]) for c in self.cols}
        self.DT = {c: Dz['col'][c]['VT'] - Dz['col'][c]['LT'] for c in self.cols}
        self.seckeys = sorted(Dz['secs'].keys(), key=lambda k: (-(k[2] - k[1]), k[1], k[0]))
        N = self.N
        self.w0, self.comp0 = self.derived(self.v0())
        # limit modes, trace labels, pinch index (Table 1 rule)
        self.lab, self.rstar = {}, {}
        for key in self.seckeys:
            s = Dz['secs'][key]
            t = self.modes(s['V'], self.w0[key])
            lab = []
            for r in range(N):
                k = [k for k in range(N) if abs(self.w0[key][k]) < mp.mpf(10) ** (-mp.mp.dps + 10) and abs(t[r] - self.al[k]) < mp.mpf(10) ** (-mp.mp.dps // 2)]
                lab.append(k[0] if k else None)
            self.lab[key] = lab
            major = [r for r in range(N) if lab[r] is None]
            self.rstar[key] = max(major) if key[0] == 'T' else min(major)   # top: smallest major; bottom: largest

    def v0(self):
        return [mp.mpf(0)] * self.dim

    def unpack(self, v):
        N = self.N
        d, X = {}, {}
        k = 0
        for c in self.cols:
            dd = v[k:k + N - 1]
            k += N - 1
            dx = v[k:k + N - 1]
            k += N - 1
            d[c] = [self.d0[c][j] + dd[j] for j in range(N - 1)]
            d[c].append(self.DT[c] - mp.fsum(d[c]))
            X[c] = [self.X0[c][j] + dx[j] for j in range(N - 1)]
            X[c].append(1 - mp.fsum(X[c]))
        return d, X

    def derived(self, v):
        N, Dz = self.N, self.D
        d, X = self.unpack(v)
        f, b = {}, {}
        for c in self.cols:
            i, m = c
            if c == (1, N):
                fc = [mp.mpf(x) for x in Dz['z']]
            else:
                fc = [mp.mpf(0)] * N
                if m < N:
                    fc = [a + x for a, x in zip(fc, d[(i, m + 1)])]
                if i > 1:
                    fc = [a + x for a, x in zip(fc, b[(i - 1, m)])]
            f[c] = fc
            b[c] = [fc[j] - d[c][j] for j in range(N)]
        comp = {('J',) + c: X[c] for c in self.cols}
        for p in range(1, N + 1):
            n = [mp.mpf(0)] * N
            if p > 1:
                n = [a + x for a, x in zip(n, b[(p - 1, p)])]
            if p < N:
                n = [a + x for a, x in zip(n, d[(p, p + 1)])]
            comp[('J', p, p)] = [x / Dz['prods'][p]['P'] for x in n]
        w = {}
        for key in self.seckeys:
            c = key[1:]
            w[key] = list(d[c]) if key[0] == 'T' else [-x for x in b[c]]
        return w, comp

    def modes(self, V, w):
        """roots of P_w(t) = V prod(al_i - t) - sum_i al_i w_i prod_{j != i}(al_j - t), decreasing."""
        N, al = self.N, self.al

        def pmul(p1, p2):
            out = [mp.mpf(0)] * (len(p1) + len(p2) - 1)
            for i, a in enumerate(p1):
                for j, b in enumerate(p2):
                    out[i + j] += a * b
            return out
        lin = [[al[i], mp.mpf(-1)] for i in range(N)]          # ascending coefficients of (al_i - t)
        P = [mp.mpf(V)]
        for i in range(N):
            P = pmul(P, lin[i])
        for i in range(N):
            Q = [-al[i] * w[i]]
            for j in range(N):
                if j != i:
                    Q = pmul(Q, lin[j])
            P = [a + b for a, b in zip(P, Q + [mp.mpf(0)] * (len(P) - len(Q)))]
        rts = mp.polyroots(list(reversed(P)), maxsteps=400, extraprec=mp.mp.prec)
        out = []
        for r in rts:
            if isinstance(r, mp.mpc):
                assert abs(r.imag) < mp.mpf(10) ** (-mp.mp.dps // 3), 'complex mode'
                r = r.real
            out.append(r)
        return sorted(out, reverse=True)

    def left(self, key, V, w, t, r):
        """explicit left eigenvector: Underwood for a major branch, (3.1) for a trace branch."""
        al, N = self.al, self.N
        k = self.lab[key][r]
        if k is None:
            return [al[i] / (al[i] - t[r]) for i in range(N)]
        Pk = Psi(V, w, al, t[r], skip=k)
        return [mp.mpf(1) if i == k else (w[k] / Pk) * al[i] / (al[i] - t[r]) for i in range(N)]

    def section_rows(self, key, w, xb, xt, n):
        """rows (3.2) of section key; n = None for E^(inf).  Returns (rows_n, rows_inf, modes)."""
        s = self.D['secs'][key]
        V = s['V']
        t = self.modes(V, w)
        rs = self.rstar[key]
        L = [self.left(key, V, w, t, r) for r in range(self.N)]
        dot = lambda a, b: mp.fsum(x * y for x, y in zip(a, b))
        pb, pt = dot(L[rs], xb), dot(L[rs], xt)
        rn, ri = [], []
        for r in range(self.N):
            if r == rs:
                continue
            if r < rs:
                mu = t[rs] / t[r]
                e_inf = dot(L[r], xb) * pt
                e_n = e_inf - (mu ** (n + 1) * dot(L[r], xt) * pb if n is not None else 0)
            else:
                mu = t[r] / t[rs]
                e_inf = dot(L[r], xt) * pb
                e_n = e_inf - (mu ** (n + 1) * dot(L[r], xb) * pt if n is not None else 0)
            rn.append(e_n)
            ri.append(e_inf)
        return rn, ri, t

    def G(self, v, ns):
        """ns: dict section -> n_sigma, or None for G_inf.  Returns (G_n, G_inf, max consecutive mode ratio)."""
        w, comp = self.derived(v)
        gn, gi = [], []
        rho = mp.mpf(0)
        for key in self.seckeys:
            s = self.D['secs'][key]
            rn, ri, t = self.section_rows(key, w[key], comp[s['bot']], comp[s['top']], None if ns is None else ns[key])
            gn += rn
            gi += ri
            rho = max(rho, max(t[r + 1] / t[r] for r in range(self.N - 1)))
        return gn, gi, rho

    def jac(self, v, ns, which):
        """central finite differences; which = 'n', 'inf' or 'diff' (G_n - G_inf)."""
        h = mp.mpf(10) ** (-(mp.mp.dps // 3))
        cols = []
        for j in range(self.dim):
            vp = list(v)
            vm = list(v)
            vp[j] += h
            vm[j] -= h
            gp_n, gp_i, _ = self.G(vp, ns)
            gm_n, gm_i, _ = self.G(vm, ns)
            if which == 'n':
                cols.append([(a - b) / (2 * h) for a, b in zip(gp_n, gm_n)])
            elif which == 'inf':
                cols.append([(a - b) / (2 * h) for a, b in zip(gp_i, gm_i)])
            else:
                cols.append([((a - c) - (b - e)) / (2 * h) for a, b, c, e in zip(gp_n, gm_n, gp_i, gm_i)])
        return mp.matrix([[cols[j][i] for j in range(self.dim)] for i in range(self.dim)])


def norm_inf(v):
    return max(abs(x) for x in v)


def matnorm_inf(M):
    return max(mp.fsum(abs(M[i, j]) for j in range(M.cols)) for i in range(M.rows))


def newton(net, ns, v, tol, maxit=60):
    """Damped Newton: the step is halved while the trial point leaves the domain of real simple modes or does not
    decrease the residual (needed only for small n, below the n_0 of Theorem 6.2)."""
    g, _, _ = net.G(v, ns)
    for it in range(maxit):
        res = norm_inf(g)
        if res < tol:
            return v, it, res
        J = net.jac(v, ns, 'n')
        dv = mp.lu_solve(J, mp.matrix(g))
        lam = mp.mpf(1)
        for _ in range(40):
            trial = [v[j] - lam * dv[j] for j in range(net.dim)]
            try:
                gt, _, _ = net.G(trial, ns)
                if norm_inf(gt) < res or lam < mp.mpf('1e-6'):
                    break
            except (AssertionError, ZeroDivisionError):
                pass
            lam /= 2
        v, g = trial, gt
    return v, maxit, norm_inf(g)


# ================================================================================================ physics

def y_of(al, x):
    S = mp.fsum(a * b for a, b in zip(al, x))
    return [a * b / S for a, b in zip(al, x)], S


def reconstruct_and_check(net, v, ns, tag):
    """Shoot every section with the physical cross-section balance and check every unit balance."""
    Dz, al, N = net.D, net.al, net.N
    w, comp = net.derived(v)
    prof = {}
    worst_end, minS, minx = mp.mpf(0), mp.mpf(10), mp.mpf(10)
    for key in net.seckeys:
        s = Dz['secs'][key]
        V, L = s['V'], s['L']
        x = list(comp[s['bot']])
        xs = [x]
        for k in range(ns[key] + 1):
            y, S = y_of(al, x)
            minS = min(minS, S)
            x = [(V * a - b) / L for a, b in zip(y, w[key])]
            xs.append(x)
        worst_end = max(worst_end, norm_inf([a - b for a, b in zip(xs[-1], comp[s['top']])]))
        prof[key] = xs
        for xx in xs[1:-1]:
            minx = min(minx, min(xx))
    # unit balances with actual streams: stages
    worst_stage = mp.mpf(0)
    for key, xs in prof.items():
        s = Dz['secs'][key]
        V, L = s['V'], s['L']
        n = ns[key]
        for k in range(1, n + 1):
            yb, _ = y_of(al, xs[k - 1])
            yk, _ = y_of(al, xs[k])
            above = xs[k + 1] if k < n else comp[s['top']]      # stage n receives the top node's liquid
            r = [L * above[i] + V * yb[i] - L * xs[k][i] - V * yk[i] for i in range(N)]
            worst_stage = max(worst_stage, norm_inf(r))
    # nodes
    worst_node = mp.mpf(0)
    nodes = set([s['bot'] for s in Dz['secs'].values()] + [s['top'] for s in Dz['secs'].values()])
    for J in nodes:
        xJ = comp[J]
        bal = [mp.mpf(0)] * N
        cond = (J == ('J', 1, 1))
        if not cond:
            yJ, SJ = y_of(al, xJ)
            minS = min(minS, SJ)
        minx = min(minx, min(xJ))
        for key, s in Dz['secs'].items():
            xs = prof[key]
            if s['bot'] == J:     # J sends vapour V y(x_J) up, receives liquid L x_1
                bal = [bb + s['L'] * xs[1][i] - s['V'] * yJ[i] for i, bb in enumerate(bal)]
            if s['top'] == J:     # J receives vapour V y(x_n), sends liquid L x_J down
                yn, _ = y_of(al, xs[-2])
                bal = [bb + s['V'] * yn[i] - s['L'] * xJ[i] for i, bb in enumerate(bal)]
        if J == ('J', 1, N):
            bal = [bb + mp.mpf(Dz['z'][i]) for i, bb in enumerate(bal)]
        if J[1] == J[2]:
            P = Dz['prods'][J[1]]['P']
            bal = [bb - P * xJ[i] for i, bb in enumerate(bal)]
        worst_node = max(worst_node, norm_inf(bal))
    # condenser composition rule: condensed vapour of the top stage of T[1,2] = x_1
    yn, _ = y_of(al, prof[('T', 1, 2)][-2])
    cond_rule = norm_inf([a - b for a, b in zip(yn, comp[('J', 1, 1)])])
    # vapour compositions positive (y = K x with S > 0) follows from x > 0; purities
    imp = [1 - comp[('J', p, p)][p - 1] for p in range(1, N + 1)]
    return dict(end=worst_end, stage=worst_stage, node=worst_node, cond=cond_rule, minS=minS, minx=minx, imp=imp)


def flows_check(Dz):
    """Total and CMO balances and duties of the specification at every node (theorem.md Lemma 4.8)."""
    N = Dz['N']
    worst_tot, worst_duty = mp.mpf(0), mp.mpf(0)
    duties = {}
    nodes = set([s['bot'] for s in Dz['secs'].values()] + [s['top'] for s in Dz['secs'].values()])
    for J in nodes:
        Vout = Vin = Lout = Lin = mp.mpf(0)
        for s in Dz['secs'].values():
            if s['bot'] == J:
                Vout += s['V']
                Lin += s['L']
            if s['top'] == J:
                Vin += s['V']
                Lout += s['L']
        feedV = feedL = mp.mpf(0)
        if J == ('J', 1, N):
            feedV, feedL = 1 - Dz['q'], Dz['q']
        P = Dz['prods'][J[1]]['P'] if J[1] == J[2] else mp.mpf(0)
        if J == ('J', 1, 1):        # total condenser: vapour in -> liquid out
            worst_tot = max(worst_tot, abs(Vin - (Lout + P)))
            duties[J] = -Vin
            continue
        tot = (Vout + Lout + P) - (Vin + Lin + feedV + feedL)
        worst_tot = max(worst_tot, abs(tot))
        duties[J] = Vout - Vin - feedV
    for J, Q in duties.items():
        if J not in (('J', 1, 1), ('J', N, N)):
            worst_duty = max(worst_duty, abs(Q))
    return worst_tot, worst_duty, duties


# ================================================================================================ the case driver

def half_profiles_ok(net, K=300):
    """(H6)/(A5) with explicit right eigenvectors, 1^T e = 1, l.e = 1."""
    Dz, al, N = net.D, net.al, net.N
    worst_neg, c0, minpiv = mp.mpf(0), mp.mpf(10), mp.mpf(10)
    for key in net.seckeys:
        s = Dz['secs'][key]
        V, L, w = s['V'], s['L'], net.w0[key]
        t = net.modes(V, w)
        E, Lv = [], []
        for r in range(N):
            k = net.lab[key][r]
            if k is None:
                e = [w[i] / (al[i] - t[r]) for i in range(N)]
                l = [al[i] / (al[i] - t[r]) for i in range(N)]
            else:
                e = [w[i] / (al[i] - al[k]) if i != k else Psi(V, w, al, al[k], skip=k) / al[k] for i in range(N)]
                l = [mp.mpf(1) if i == k else mp.mpf(0) for i in range(N)]
            se = mp.fsum(e)
            e = [x / se for x in e]
            le = mp.fsum(a * b for a, b in zip(l, e))
            l = [x / le for x in l]
            E.append(e)
            Lv.append(l)
        rs = net.rstar[key]
        xb, xt = net.comp0[s['bot']], net.comp0[s['top']]
        cb = [mp.fsum(a * b for a, b in zip(Lv[r], xb)) for r in range(N)]
        ct = [mp.fsum(a * b for a, b in zip(Lv[r], xt)) for r in range(N)]
        minpiv = min(minpiv, cb[rs], ct[rs])
        for k in range(K):
            vb = [E[rs][i] + mp.fsum((cb[r] / cb[rs]) * (t[r] / t[rs]) ** k * E[r][i] for r in range(rs + 1, N))
                  for i in range(N)]
            vt = [E[rs][i] + mp.fsum((ct[r] / ct[rs]) * (t[rs] / t[r]) ** k * E[r][i] for r in range(rs))
                  for i in range(N)]
            worst_neg = min(worst_neg, min(vb), min(vt))
            c0 = min(c0, mp.fsum(vb), mp.fsum(vt))
    return worst_neg, c0, minpiv


CASES = {
    'n3B': dict(al=['4', '2', '1'], z=['1/3', '1/3', '1/3'], q='1', delta='0.05', eta_nb='0.025', eta_bin='0.003',
                pattern='BB', export='config_n3_B.json', ns=[64, 48, 32, 16, 8], nmax=128),
    'n3BT': dict(al=['4', '2', '1'], z=['1/3', '1/3', '1/3'], q='1', delta='0.05', eta_nb='0.025', eta_bin='0.003',
                 pattern='BT', export='config_n3_N.json', ns=[64, 48, 32, 16, 8], nmax=128),
    'n3q': dict(al=['5', '2.2', '1'], z=['0.2', '0.5', '0.3'], q='1.4', delta='0.02', eta_nb='0.01', eta_bin='0.002',
                pattern='TB', export=None, ns=[64, 48, 32, 16], nmax=128),
    'n3TT': dict(al=['3', '1.8', '1'], z=['0.3', '0.3', '0.4'], q='-0.5', delta='0.03', eta_nb='0.01', eta_bin='0.002',
                 pattern='TT', export=None, ns=[64, 48, 32, 16], nmax=128),
    'n4D': dict(al=['8', '4', '2', '1'], z=['0.25'] * 4, q='1', delta='0.05', eta_nb='0.005', eta_bin='0.003',
                pattern='BBB', export='config_n4_1.json', ns=[48, 32, 16], nmax=64),
    'n4Dlong': dict(al=['8', '4', '2', '1'], z=['0.25'] * 4, q='1', delta='0.05', eta_nb='0.005', eta_bin='0.003',
                    pattern='BBB', export='config_n4_1.json', ns=[96, 80, 64], nmax=96),
    'n5D': dict(al=['16', '8', '4', '2', '1'], z=['0.2'] * 5, q='0.5', delta='0.05', eta_nb='0.002', eta_bin='0.002',
                pattern='BBBB', export='config_n5_1.json', ns=[36, 24], nmax=36),
}


def frac(s):
    if '/' in s:
        a, b = s.split('/')
        return mp.mpf(a) / mp.mpf(b)
    return mp.mpf(s)


def run_case(name):
    C = CASES[name]
    t0 = time.time()
    N = len(C['al'])
    # precision: forward reconstruction amplifies by (lambda_max/lambda_min)^(n+1)
    al_ratio = float(frac(C['al'][0]) / frac(C['al'][-1]))
    mp.mp.dps = int(40 + (max(C['ns']) + 2) * mp.log10(3 * al_ratio))
    al = [frac(a) for a in C['al']]
    z = [frac(a) for a in C['z']]
    q, delta = frac(C['q']), frac(C['delta'])
    eta_nb = {}
    for Lc in range(N - 1, 1, -1):
        for i in range(1, N - Lc + 1):
            eta_nb[(i, i + Lc)] = frac(C['eta_nb'])
    eta_bin = {p: frac(C['eta_bin']) for p in range(1, N)}
    print('=' * 110)
    print('case %s: alpha=%s z=%s q=%s delta=%s eta_nb=%s eta_bin=%s pattern=%s | %d digits' %
          (name, C['al'], C['z'], C['q'], C['delta'], C['eta_nb'], C['eta_bin'], C['pattern'], mp.mp.dps), flush=True)
    Dz = design(al, z, q, delta, eta_nb, eta_bin, C['pattern'])
    print('     H =', [nstr(Dz['H'][j], 8) for j in sorted(Dz['H'])], ' H* + delta = R =', nstr(Dz['R'], 10))
    # Lemma 4.9: middle feed roots global
    worst49 = mp.mpf(0)
    for c, cd in Dz['col'].items():
        if cd['kind'] == 'nb':
            for r in range(c[0] + 1, c[1] - 1):
                worst49 = max(worst49, abs(cd['th'][r] - Dz['theta'][r]))
    e = {p: Dz['Hs'] - Dz['H'][p] for p in Dz['H']}
    check('(D) admissible: valid non-binary columns, kappa_p > 0, forced flows in range', Dz['valid'],
          'kappa - e_p - delta = %s' % [nstr(Dz['kappa'][p] - e[p] - delta, 3) for p in sorted(Dz['kappa'])])
    check('(D) Lemma 4.9: middle feed roots of non-binary columns are global', worst49 < mp.mpf(10) ** (-mp.mp.dps + 10),
          'max dev %s' % nstr(worst49, 2))
    if C['export']:
        compare_with_D(Dz, C['export'])
    wt, wd, duties = flows_check(Dz)
    check('(A1) total balances; zero duty at junctions and side draws; reboiler duty = R; condenser removes V_T[1,2]',
          wt < mp.mpf(10) ** (-mp.mp.dps + 10) and wd < mp.mpf(10) ** (-mp.mp.dps + 10)
          and abs(duties[('J', N, N)] - Dz['R']) < mp.mpf(10) ** (-mp.mp.dps + 10) and duties[('J', 1, 1)] < 0,
          'reboiler %s, condenser %s' % (nstr(duties[('J', N, N)], 10), nstr(duties[('J', 1, 1)], 10)))
    net = Net(Dz)
    v0 = net.v0()
    g_n, g_i, rho0 = net.G(v0, None)
    check('(P) G_inf(v_inf) = 0 with explicit Underwood / trace eigen-data', norm_inf(g_i) < mp.mpf(10) ** (-mp.mp.dps + 15),
          'max |row| = %s' % nstr(norm_inf(g_i), 2))
    neg, c0, minpiv = half_profiles_ok(net)
    check('(P) pivots c_r*(x_bot), c_r*(x_top) > 0; half-profiles (H6) >= 0 (k < 300), inf 1^T v_k = c_0 > 0',
          minpiv > 0 and neg > -mp.mpf(10) ** (-mp.mp.dps + 15) and c0 > 0,
          'min pivot %s, most negative entry %s, c_0 %s' % (nstr(minpiv, 4), nstr(neg, 2), nstr(c0, 4)))
    # rho_inf (adjacent to the pinch) and rho_bar (all consecutive) at v_inf
    rho_inf = mp.mpf(0)
    for key in net.seckeys:
        t = net.modes(Dz['secs'][key]['V'], net.w0[key])
        rs = net.rstar[key]
        if rs + 1 < N:
            rho_inf = max(rho_inf, t[rs + 1] / t[rs])
        if rs - 1 >= 0:
            rho_inf = max(rho_inf, t[rs] / t[rs - 1])
    A = net.jac(v0, None, 'inf')
    sv = mp.svd_r(A, compute_uv=False)
    smin = min(sv)
    Ainv = mp.inverse(A)
    nAinv = matnorm_inf(Ainv)
    check('(P) DG_inf(v_inf) invertible', smin > mp.mpf('1e-12'),
          'dim %d, sigma_min %s, sigma_max %s, ||A^-1||_inf %s, rho_inf %s, rho_bar(v_inf) %s' %
          (net.dim, nstr(smin, 4), nstr(max(sv), 4), nstr(nAinv, 4), nstr(rho_inf, 6), nstr(rho0, 6)))
    # ---------------------------------------------------------------- (C1) the C^1 estimate at random points
    rnd = random.Random(12345)
    r_dom, in_domain = mp.mpf(0), True     # largest tested radius whose samples all have real simple modes
    for rad in ['3e-4', '1e-3', '1e-2']:
        rad = mp.mpf(rad)
        pts = [[rad * mp.mpf(rnd.uniform(-1, 1)) for _ in range(net.dim)] for _ in range(4)]
        try:
            rho_bar = max(net.G(p, None)[2] for p in pts)
            r0s, r1s = [], []
            for n in [8, 16, 32, 64, C['nmax']]:
                ns = {key: n for key in net.seckeys}
                m0 = m1 = mp.mpf(0)
                for p in pts:
                    gn, gi, _ = net.G(p, ns)
                    m0 = max(m0, norm_inf([a - b for a, b in zip(gn, gi)]) / rho_bar ** (n + 1))
                    Dd = net.jac(p, ns, 'diff')
                    m1 = max(m1, matnorm_inf(Dd) / ((n + 2) * rho_bar ** (n + 1)))
                r0s.append(m0)
                r1s.append(m1)
        except (AssertionError, ZeroDivisionError):
            print('INFO (C1) radius %s: a sample point lies outside the domain of real simple modes (the neighborhood '
                  'of Lemma 5.2 is smaller than this radius)' % nstr(rad, 2), flush=True)
            in_domain = False
            continue
        if in_domain:
            r_dom = rad
        ok = max(r0s[1:]) <= 2 * r0s[0] + mp.mpf('1e-30') and max(r1s[1:]) <= 2 * r1s[0] + mp.mpf('1e-30')
        check('(C1) radius %s: |G_n-G_inf|/rho^(n+1) and ||D(G_n-G_inf)||/((n+2)rho^(n+1)) bounded in n' % nstr(rad, 2),
              ok, 'rho_bar %s | n = 8,16,32,64,%d: %s | %s' % (nstr(rho_bar, 5), C['nmax'],
                                                            [nstr(x, 3) for x in r0s], [nstr(x, 3) for x in r1s]))
    # ---------------------------------------------------------------- (N) Newton zeros, Lemma 6.1, contraction
    v = list(v0)
    tol = mp.mpf(10) ** (-mp.mp.dps + 20)
    configs = [('n=%d' % n, {key: n for key in net.seckeys}) for n in C['ns']]
    rnd2 = random.Random(777)
    mixed = {key: rnd2.randint(min(C['ns'][:2]) // 2, max(C['ns'])) for key in net.seckeys}
    configs.insert(1, ('mixed(min %d)' % min(mixed.values()), mixed))
    for label, ns in configs:
        nmin = min(ns.values())
        v, it, res = newton(net, ns, v, tol)
        dist = norm_inf(v)
        g0, _, _ = net.G(v0, ns)
        bound = 2 * nAinv * norm_inf(g0)
        rec = reconstruct_and_check(net, v, ns, label)
        okN = res < tol * 10
        # perturbative regime: Lemma 6.1's hypothesis ||A^-1|| |G_n(v_inf)| <= r_1/2 is plausible with r_1 = r_check :=
        # min(1e-3, r_dom), r_dom the largest tested radius whose random samples all have real simple modes (the ball of
        # Lemma 6.1 must lie in that domain); outside the regime the theorem says nothing and the results are information
        # only (e.g. n = 16 in case n3TT: a NON-physical zero of G_n, cf. theorem.md §8.5).
        r_check = min(mp.mpf('1e-3'), r_dom)
        regime = bound <= r_check
        tagN = '(N) %s: zero v_n of G_n (Newton it %d, |G_n| %s)' % (label, it, nstr(res, 2))
        infoN = '|v_n - v_inf| = %s, /rho_inf^(n+1) = %s, Lemma 6.1 bound 2||A^-1|| |G_n(v_inf)| = %s' % (
            nstr(dist, 4), nstr(dist / rho_inf ** (nmin + 1), 4), nstr(bound, 4))
        okS = (rec['end'] < mp.mpf(10) ** (-30) and rec['stage'] < mp.mpf(10) ** (-30) and rec['node'] < mp.mpf(10) ** (-30)
               and rec['cond'] < mp.mpf(10) ** (-30) and rec['minS'] > 0 and rec['minx'] > 0)
        tagS = '(S) %s: physical re-shooting reaches x_top; all stage/node/condenser balances; S > 0, all x > 0' % label
        infoS = 'end %s, stage %s, node %s, cond %s | min S %s, min x %s | impurities %s' % (
            nstr(rec['end'], 2), nstr(rec['stage'], 2), nstr(rec['node'], 2), nstr(rec['cond'], 2),
            nstr(rec['minS'], 4), nstr(rec['minx'], 3), [nstr(x, 3) for x in rec['imp']])
        if regime:
            check(tagN, okN, infoN)
            check('(N) %s: Lemma 6.1 displacement bound |v_n - v_inf| <= 2||A^-1|| |G_n(v_inf)|' % label, dist <= bound, '')
            check(tagS, okS, infoS)
        else:
            print('INFO (outside the perturbative regime, 2||A^-1|| |G_n(v_inf)| > r_check = %s) %s [%s] %s' %
                  (nstr(r_check, 2), tagN, 'ok' if okN else 'not converged', infoN))
            print('INFO (outside the perturbative regime) %s [%s] %s' % (tagS, 'physical' if okS else 'NOT physical', infoS),
                  flush=True)
        if label in ('n=%d' % C['ns'][0], 'n=%d' % C['ns'][1]):
            # contraction Phi from random starts at radius r around v_inf
            for ir, rr in enumerate([r_check, 3 * r_check]):
                rr = mp.mpf(rr)
                conv = []
                for trial in range(3):
                    u = [v0[j] + rr * mp.mpf(rnd2.uniform(-1, 1)) for j in range(net.dim)]
                    ok, factor, prev = False, mp.mpf(0), None
                    try:
                        for k in range(400):
                            g, _, _ = net.G(u, ns)
                            du = Ainv * mp.matrix(g)
                            u = [u[j] - du[j] for j in range(net.dim)]
                            dcur = norm_inf([a - b for a, b in zip(u, v)])
                            if prev is not None and prev > mp.mpf('1e-35'):
                                factor = max(factor, dcur / prev)
                            prev = dcur
                            if dcur < mp.mpf('1e-40'):
                                ok = True
                                break
                    except (AssertionError, ZeroDivisionError):
                        ok = None            # an iterate left the domain of the eigen-data (outside the ball B(v_inf, r_1))
                    conv.append((ok, k, factor))
                info = 'iterations %s, max observed step ratio %s' % ([c[1] for c in conv], [nstr(c[2], 3) for c in conv])
                if ir == 0 and 2 * bound <= rr:     # Lemma 6.1 hypotheses plausible with r_1 = r_check
                    # a start whose orbit leaves the domain of real simple modes (None) was not in the ball of Lemma 6.1
                    # (that ball lies in the domain); a start that stays in the domain and does not converge (False), or a
                    # step ratio > 1/2, would contradict the contraction
                    okc = (all(c[0] is not False for c in conv) and any(c[0] is True for c in conv)
                           and all(c[2] <= mp.mpf('0.5') for c in conv if c[0] is True))
                    check('(N) %s: Phi = v - A^-1 G_n contracts to v_n (step ratio <= 1/2) from random starts at radius %s'
                          % (label, nstr(rr, 2)), okc, info + ' | starts that left the domain: %d of 3'
                          % sum(1 for c in conv if c[0] is None))
                else:
                    print('INFO (N) %s: Phi from 3 random starts at radius %s: converged to v_n %s; %s' %
                          (label, nstr(rr, 2), [c[0] for c in conv], info), flush=True)
    print('case %s done in %.1f s' % (name, time.time() - t0), flush=True)


if __name__ == '__main__':
    names = sys.argv[1:] or ['n3B', 'n3BT', 'n3q', 'n4D']
    for nm in names:
        run_case(nm)
    print('=' * 110)
    print('ALL PASS' if not FAILS else 'FAILURES: %s' % FAILS)
    sys.exit(1 if FAILS else 0)
