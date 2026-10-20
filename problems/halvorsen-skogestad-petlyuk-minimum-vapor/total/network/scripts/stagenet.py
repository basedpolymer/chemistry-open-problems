#!/usr/bin/env python3
"""Small steady-state solver for networks of equilibrium stages (Agent N, used by the check scripts).

Model: constant relative volatilities alpha, constant molar overflow (CMO), equilibrium stages
y = alpha x / (alpha . x).  A network is a list of stages; every stage k sends its liquid and its vapor to
destinations (other stages or product sinks) with fixed flows (this covers splitters: one outlet divided among
several destinations with the same composition), receives any number of streams (mixing), may receive feed portions,
and has a net vaporization duty Delta_k = V_out - V_in - Fv (reboiler/heater if > 0, condenser/cooler if < 0).
Total condensers / total vaporizers are stages with zero vapor / zero liquid outflow.  Flows are given;
compositions are solved.

Composition solve: successive substitution (theta method, sparse linear solves per component) followed by Newton's
method with a sparse analytic Jacobian.  Floating point; used for checks only.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class Net:
    def __init__(self, alpha, z, q):
        self.alpha = np.asarray(alpha, float)
        self.N = len(alpha)
        self.z = np.asarray(z, float)
        self.q = float(q)
        self.nst = 0
        self.liq = []   # liq[k] = list of (dest, flow); dest = stage index or ('P', name)
        self.vap = []
        self.feed = {}  # stage -> (F_liq, F_vap, F_total) with composition z
        self.names = []
        self._mats = None

    def add_stage(self, name=None):
        self.liq.append([])
        self.vap.append([])
        self.names.append(name if name is not None else f"s{self.nst}")
        self.nst += 1
        self._mats = None
        return self.nst - 1

    def add_liq(self, k, dest, flow):
        self.liq[k].append((dest, float(flow)))
        self._mats = None

    def add_vap(self, k, dest, flow):
        self.vap[k].append((dest, float(flow)))
        self._mats = None

    def add_feed(self, k, F=1.0):
        fl, fv, ft = self.feed.get(k, (0.0, 0.0, 0.0))
        self.feed[k] = (fl + F * self.q, fv + F * (1 - self.q), ft + F)
        self._mats = None

    # ------------------------------------------------------------------ flows
    def flows(self):
        n = self.nst
        Lout = np.array([sum(f for _, f in self.liq[k]) for k in range(n)])
        Vout = np.array([sum(f for _, f in self.vap[k]) for k in range(n)])
        Lin = np.zeros(n)
        Vin = np.zeros(n)
        for k in range(n):
            for d, f in self.liq[k]:
                if not isinstance(d, tuple):
                    Lin[d] += f
            for d, f in self.vap[k]:
                if not isinstance(d, tuple):
                    Vin[d] += f
        Fl = np.zeros(n)
        Fv = np.zeros(n)
        for k, (fl, fv, _) in self.feed.items():
            Fl[k] += fl
            Fv[k] += fv
        delta = Vout - Vin - Fv            # net vaporization
        massres = (Lin + Vin + Fl + Fv) - (Lout + Vout)
        return Lout, Vout, Lin, Vin, Fl, Fv, delta, massres

    def mats(self):
        """Sparse matrices: Mliq[d,k] = liquid flow k->d, Mvap[d,k] = vapor flow k->d (internal streams);
        Lout, Vout totals (including products); feed vector b[k] = total feed at k."""
        if self._mats is not None:
            return self._mats
        n = self.nst
        rl, cl, vl, rv, cv, vv = [], [], [], [], [], []
        for k in range(n):
            for d, f in self.liq[k]:
                if not isinstance(d, tuple):
                    rl.append(d); cl.append(k); vl.append(f)
            for d, f in self.vap[k]:
                if not isinstance(d, tuple):
                    rv.append(d); cv.append(k); vv.append(f)
        Ml = sp.csr_matrix((vl, (rl, cl)), shape=(n, n))
        Mv = sp.csr_matrix((vv, (rv, cv)), shape=(n, n))
        Lout = np.array([sum(f for _, f in self.liq[k]) for k in range(n)])
        Vout = np.array([sum(f for _, f in self.vap[k]) for k in range(n)])
        bF = np.zeros(n)
        for k, (_, _, F) in self.feed.items():
            bF[k] = F
        self._mats = (Ml, Mv, Lout, Vout, bF)
        return self._mats

    # ------------------------------------------------------------ compositions
    def residual(self, X):
        a = self.alpha
        n = self.nst
        Ml, Mv, Lout, Vout, bF = self.mats()
        X = X.reshape(n, self.N)
        S = X @ a
        Y = X * a / S[:, None]
        res = Ml @ X + Mv @ Y - Lout[:, None] * X - Vout[:, None] * Y + bF[:, None] * self.z[None, :]
        return res.ravel()

    def jacobian(self, X):
        a = self.alpha
        n, N = self.nst, self.N
        Ml, Mv, Lout, Vout, bF = self.mats()
        X = X.reshape(n, N)
        S = X @ a
        # dY_k/dX_k = diag(a)/S - (a x)(a)^T / S^2   (N x N block per stage)
        blocks = []
        for k in range(n):
            blocks.append(np.diag(a) / S[k] - np.outer(a * X[k], a) / S[k] ** 2)
        DY = sp.block_diag(blocks, format='csr')
        I_N = sp.identity(N, format='csr')
        J = sp.kron(Ml, I_N) + sp.kron(Mv, I_N) @ DY - sp.kron(sp.diags(Lout), I_N) \
            - sp.kron(sp.diags(Vout), I_N) @ DY
        return J.tocsc()

    def solve(self, tol=1e-13, maxit=100, X0=None):
        n, N, a = self.nst, self.N, self.alpha
        Ml, Mv, Lout, Vout, bF = self.mats()
        X = np.tile(self.z, (n, 1)) if X0 is None else X0.copy()
        # successive substitution (a few sweeps only when a good initial guess is supplied)
        for it in range(200 if X0 is None else 3):
            S = X @ a
            Xn = np.zeros_like(X)
            ok = True
            for p in range(N):
                K = a[p] / S
                A = Ml + Mv @ sp.diags(K) - sp.diags(Lout + Vout * K)
                try:
                    Xn[:, p] = spla.spsolve(A.tocsc(), -bF * self.z[p])
                except Exception:
                    ok = False
                    break
            if not ok or not np.all(np.isfinite(Xn)):
                break
            Xn = np.maximum(Xn, 1e-300)
            Xn = Xn / Xn.sum(axis=1, keepdims=True)
            dif = np.max(np.abs(Xn - X))
            X = Xn
            if dif < 1e-12:
                break
        x = X.ravel().copy()
        r = self.residual(x)
        for it in range(maxit):
            nr = np.max(np.abs(r))
            if nr < tol:
                break
            J = self.jacobian(x)
            try:
                dx = spla.spsolve(J, -r)
            except Exception:
                break
            if not np.all(np.isfinite(dx)):
                break
            # fraction-to-boundary rule: keep every mole fraction >= 1% of its current value
            neg = dx < 0
            lam = 1.0
            if np.any(neg):
                lam = min(1.0, 0.99 * np.min(-x[neg] / dx[neg]))
            accepted = False
            while lam > 1e-10:
                xn = x + lam * dx
                rn = self.residual(xn)
                if np.max(np.abs(rn)) < nr * (1 - 1e-4 * lam):
                    accepted = True
                    break
                lam /= 2
            if not accepted:
                break
            x, r = xn, rn
        X = x.reshape(n, N)
        return X, np.max(np.abs(self.residual(x)))


def solve_log(net, X0, tol=1e-12, maxit=200):
    """Newton's method in log mole fractions on the relative component balances
    g_{k,i} = log(inflow_{k,i}) - log(outflow_{k,i}) (robust for tiny compositions near pinches)."""
    a = net.alpha
    n, N = net.nst, net.N
    Ml, Mv, Lout, Vout, bF = net.mats()
    Mlc, Mvc = Ml.tocoo(), Mv.tocoo()
    xi = np.log(np.maximum(X0, 1e-300)).ravel()

    def parts(xi):
        X = np.exp(xi.reshape(n, N))
        S = X @ a
        K = a[None, :] / S[:, None]
        I = Ml @ X + Mv @ (K * X) + bF[:, None] * net.z[None, :]
        O = (Lout[:, None] + Vout[:, None] * K) * X
        return X, S, K, I, O

    def g_of(xi):
        X, S, K, I, O = parts(xi)
        with np.errstate(divide='ignore', invalid='ignore'):
            g = np.log(I) - np.log(O)
        return g.ravel(), (X, S, K, I, O)

    g, P = g_of(xi)
    for it in range(maxit):
        ng = np.max(np.abs(g)); m2 = float(g @ g)
        if not np.isfinite(ng) or ng < tol:
            break
        X, S, K, I, O = P
        rows, cols, vals = [], [], []
        # -d log O_{k,i} / d xi_{k,j}
        for k in range(n):
            base = k * N
            blk = -(np.eye(N) - Vout[k] * np.outer(a * X[k], a * X[k]) / (S[k] ** 2 * O[k][:, None]))
            r_, c_ = np.meshgrid(range(N), range(N), indexing='ij')
            rows.append((base + r_).ravel()); cols.append((base + c_).ravel()); vals.append(blk.ravel())
        # d log I_{k,i} / d xi_{k',j}
        for (kk, kp, fl, isv) in [(r, c, v, False) for r, c, v in zip(Mlc.row, Mlc.col, Mlc.data)] + \
                                 [(r, c, v, True) for r, c, v in zip(Mvc.row, Mvc.col, Mvc.data)]:
            if isv:
                blk = fl * (np.diag(K[kp] * X[kp]) - np.outer(a * X[kp], a * X[kp]) / S[kp] ** 2) / I[kk][:, None]
            else:
                blk = fl * np.diag(X[kp]) / I[kk][:, None]
            r_, c_ = np.meshgrid(range(N), range(N), indexing='ij')
            rows.append((kk * N + r_).ravel()); cols.append((kp * N + c_).ravel()); vals.append(blk.ravel())
        Jm = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n * N, n * N))
        try:
            dx = spla.spsolve(Jm, -g)
        except Exception:
            break
        if not np.all(np.isfinite(dx)):
            break
        lam = min(1.0, 2.0 / max(1e-300, np.max(np.abs(dx))))  # at most a factor e^2 per step
        accepted = False
        while lam > 1e-10:
            xin = xi + lam * dx
            gn, Pn = g_of(xin)
            if np.all(np.isfinite(gn)) and float(gn @ gn) < m2 * (1 - 1e-4 * lam):
                accepted = True
                break
            lam /= 2
        if not accepted:
            break
        xi, g, P = xin, gn, Pn
    X = np.exp(xi.reshape(n, N))
    return X, np.max(np.abs(net.residual(X.ravel()))), np.max(np.abs(g))


def feed_roots(alpha, z, q):
    """Roots theta_j in (alpha_{j+1}, alpha_j) of sum alpha z/(alpha - theta) = 1 - q (bisection)."""
    alpha = np.asarray(alpha, float)
    z = np.asarray(z, float)
    f = lambda t: np.sum(alpha * z / (alpha - t)) - (1 - q)
    roots = []
    for j in range(len(alpha) - 1):
        lo, hi = alpha[j + 1], alpha[j]
        a_, b_ = lo, hi
        for _ in range(200):
            m = 0.5 * (a_ + b_)
            if m <= lo or m >= hi:
                break
            if f(m) > 0:
                b_ = m
            else:
                a_ = m
        roots.append(0.5 * (a_ + b_))
    return np.array(roots)


def cut_check(net, X, theta):
    """Terms of the general cut identity at theta with U = {stages with u(x_k; theta) >= 0}."""
    a = net.alpha
    n = net.nst
    S = X @ a
    Y = X * a / S[:, None]
    u = (X * a / (a - theta)).sum(axis=1)
    inU = u >= 0
    Lout, Vout, Lin, Vin, Fl, Fv, delta, massres = net.flows()

    def phi(comp_flows, vapor):
        tot = comp_flows.sum()
        return (tot if vapor else 0.0) - np.sum(a * comp_flows / (a - theta))

    cross = 0.0          # Phi(U^c -> U) - Phi(U -> U^c), streams between stages
    min_in = np.inf      # min flux of streams U^c -> U (should be >= 0)
    max_out = -np.inf    # max flux of streams U -> U^c (should be <= 0)
    Hprime = 0.0         # flux of product streams leaving U^c units
    prods = []
    for k in range(n):
        for lst, comp, vapor in ((net.liq[k], X[k], False), (net.vap[k], Y[k], True)):
            for d, f in lst:
                fl = phi(f * comp, vapor)
                if isinstance(d, tuple):
                    prods.append((d[1], k, f, vapor, fl, inU[k]))
                    if not inU[k]:
                        Hprime += fl
                    continue
                if inU[k] and not inU[d]:
                    cross -= fl
                    max_out = max(max_out, fl)
                elif (not inU[k]) and inU[d]:
                    cross += fl
                    min_in = min(min_in, fl)
    feedflux = 0.0
    for k, (fl_, fv_, F) in net.feed.items():
        ff = fv_ - np.sum(a * F * net.z / (a - theta))
        if not inU[k]:
            feedflux += ff
    R = np.sum(np.maximum(delta, 0))
    RUc = np.sum(np.maximum(delta, 0)[~inU])
    QUc = np.sum(np.maximum(-delta, 0)[~inU])
    RU = np.sum(np.maximum(delta, 0)[inU])
    # sum of the balances of the U^c units: R(U^c) - Q(U^c) + feedflux(U^c) = H' + cross
    lhs = RUc - QUc + feedflux
    rhs = Hprime + cross
    return dict(R=R, RU=RU, RUc=RUc, QUc=QUc, Hprime=Hprime, cross=cross, min_in=min_in, max_out=max_out,
                identity_residual=lhs - rhs, feedflux=feedflux, u=u, inU=inU, prods=prods,
                massres=np.max(np.abs(massres)))


def solve_lse(net, X0, tol=1e-11, maxit=300, verbose=False):
    """Newton's method in log mole fractions with log-sum-exp inflows (no underflow for compositions far below
    1e-300): g_{k,i} = log(sum of inflows of i into k) - log((L_k + V_k K_{k,i}) x_{k,i}).  Returns
    (X as log mole fractions xi, max relative residual)."""
    a = net.alpha
    la = np.log(a)
    n, N = net.nst, net.N
    Ml, Mv, Lout, Vout, bF = net.mats()
    src, dst, lfl, isv = [], [], [], []
    for M_, v_ in ((Ml.tocoo(), False), (Mv.tocoo(), True)):
        for r, c, f in zip(M_.row, M_.col, M_.data):
            if f > 0:
                src.append(c); dst.append(r); lfl.append(np.log(f)); isv.append(v_)
    feeds = [(k, F) for k, F in enumerate(bF) if F > 0]
    for k, F in feeds:   # feed as a pseudo-stream with constant log term
        src.append(-1); dst.append(k); lfl.append(np.log(F)); isv.append(False)
    src = np.array(src); dst = np.array(dst); lfl = np.array(lfl); isv = np.array(isv)
    order = np.argsort(dst, kind='stable')
    src, dst, lfl, isv = src[order], dst[order], lfl[order], isv[order]
    starts = np.r_[0, np.flatnonzero(np.diff(dst)) + 1]
    dsts = dst[starts]
    assert len(dsts) == n, "every stage must receive a stream"
    lz = np.log(net.z)
    xi = np.asarray(X0, float)   # X0 = initial LOG mole fractions
    xi = xi.reshape(n, N).copy()

    def evalg(xi):
        m = xi.max(axis=1, keepdims=True)
        lS = (m[:, 0] + np.log(np.exp(xi - m) @ a))             # log S_k
        lK = la[None, :] - lS[:, None]
        T = np.empty((len(src), N))
        real = src >= 0
        T[real] = lfl[real, None] + xi[src[real]] + np.where(isv[real, None], lK[src[real]], 0.0)
        T[~real] = lfl[~real, None] + lz[None, :]
        Tm = np.maximum.reduceat(T, starts, axis=0)
        W = np.exp(T - Tm[np.searchsorted(starts, np.arange(len(src)), side='right') - 1])
        sW = np.add.reduceat(W, starts, axis=0)
        lI = Tm + np.log(sW)
        lO = np.log(Lout[:, None] + Vout[:, None] * np.exp(lK)) + xi
        g = lI - lO
        return g, (lS, lK, T, lI)

    g, P = evalg(xi)
    seg = np.searchsorted(starts, np.arange(len(src)), side='right') - 1
    for it in range(maxit):
        ng = np.max(np.abs(g))
        m2 = float(np.sum(g * g))
        if verbose:
            print("  it", it, ng)
        if not np.isfinite(ng) or ng < tol:
            break
        lS, lK, T, lI = P
        X = np.exp(xi)
        rows, cols, vals = [], [], []
        RR, CC = np.meshgrid(range(N), range(N), indexing='ij')
        # outflow part
        K = np.exp(lK)
        S = np.exp(lS)
        for k in range(n):
            coef = Vout[k] * K[k] / (Lout[k] + Vout[k] * K[k])          # per i
            blk = -(np.eye(N) - np.outer(coef, a * X[k] / S[k]))
            rows.append((k * N + RR).ravel()); cols.append((k * N + CC).ravel()); vals.append(blk.ravel())
        # inflow part
        Wn = np.exp(T - lI[dst])                                          # softmax weights per stream, comp
        for s_ in range(len(src)):
            kp = src[s_]
            if kp < 0:
                continue
            k = dst[s_]
            w = Wn[s_]
            if isv[s_]:
                blk = w[:, None] * (np.eye(N) - (a * X[kp] / S[kp])[None, :])
            else:
                blk = np.diag(w)
            rows.append((k * N + RR).ravel()); cols.append((kp * N + CC).ravel()); vals.append(blk.ravel())
        Jm = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n * N, n * N))
        try:
            dx = spla.spsolve(Jm, -g.ravel()).reshape(n, N)
        except Exception:
            break
        if not np.all(np.isfinite(dx)):
            break
        lam = min(1.0, 5.0 / max(1e-300, np.max(np.abs(dx))))
        accepted = False
        while lam > 1e-12:
            xin = xi + lam * dx
            gn, Pn = evalg(xin)
            if np.all(np.isfinite(gn)) and float(np.sum(gn * gn)) < m2 * (1 - 1e-4 * lam):
                accepted = True
                break
            lam /= 2
        if not accepted:
            break
        xi, g, P = xin, gn, Pn
    return xi, float(np.max(np.abs(g)))


def solve_ptc(net, X0, tol=1e-12, dt0=0.1, maxit=2000, verbose=False):
    """Pseudo-transient continuation: implicit Euler on h dx_k/dt = (component balances of stage k), unit liquid
    holdup h = 1 per stage, time step increased while the Newton corrector converges (SER rule); ends with plain
    Newton steps.  Robust for the ill-conditioned (near-total-reflux, pinched) steady states."""
    n, N = net.nst, net.N
    x = np.asarray(X0, float).ravel().copy()
    r = net.residual(x)
    nr0 = np.max(np.abs(r))
    dt = dt0
    I = sp.identity(n * N, format='csc')
    for it in range(maxit):
        nr = np.max(np.abs(r))
        if verbose and it % 20 == 0:
            print("  ptc", it, nr, dt)
        if nr < tol:
            break
        J = net.jacobian(x)
        A = (I / dt - J).tocsc()
        try:
            dx = spla.spsolve(A, r)
        except Exception:
            dt /= 4
            continue
        neg = dx < 0
        lam = 1.0
        if np.any(neg):
            lam = min(1.0, 0.95 * np.min(-x[neg] / dx[neg]))
        xn = x + lam * dx
        rn = net.residual(xn)
        nrn = np.max(np.abs(rn))
        if np.isfinite(nrn) and nrn < 2 * nr:
            x, r = xn, rn
            dt = min(dt * max(1.0, min(4.0, nr / max(nrn, 1e-300))), 1e12)
        else:
            dt /= 4
        if dt < 1e-8:
            break
    return x.reshape(n, N), float(np.max(np.abs(net.residual(x))))
