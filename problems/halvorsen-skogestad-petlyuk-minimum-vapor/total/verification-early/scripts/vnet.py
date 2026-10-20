#!/usr/bin/env python3
"""Independent CMO / constant-alpha steady-state network solver (early verifier, wave 2).

Written from scratch by the verifier; it does not import or copy any code of the solving agents.

Model.  Units:
  'E'  equilibrium stage: any inflows, duty Q (external), outflows = saturated liquid (flow Lo, composition x)
       and saturated vapor (flow Vo, composition y_i = alpha_i x_i / S(x)), each phase split by fixed fractions
       among destination units / product outlets.
  'ML' / 'MV'  single-phase mixer-splitter (liquid / vapor), Q = 0.
  'C'  total condenser (vapor in, liquid out, Q = -vapor in).
  'R'  total vaporizer (liquid in, vapor out, Q = +liquid in).
Feeds: list of (unit, amount, z, q) (several feeds allowed; the theory needs one composition/thermal state).
Energy (CMO, heat of vaporization 1): sum_out v - sum_in v = Q at E units; a feed of amount f and liquid fraction q
brings vapor-equivalent (1 - q) f.
Given split fractions and E duties, total flows solve a linear system; given S_u (one scalar per two-phase E unit),
each component's flows solve a linear system; the S_u are found by Newton's method on sum_i y_{u,i} = 1 with an
analytic Jacobian.

Environment: Python 3, numpy, scipy (WSL Fedora sci env).  Library module, no output of its own.
"""
import numpy as np

ACCEPTS = {'E': ('L', 'V'), 'ML': ('L',), 'MV': ('V',), 'C': ('V',), 'R': ('L',)}
OUTPH = {'E': ('L', 'V'), 'ML': ('L',), 'MV': ('V',), 'C': ('L',), 'R': ('V',)}


class Net:
    def __init__(self, alpha):
        self.a = np.asarray(alpha, float)
        self.N = len(self.a)
        self.types = []
        self.Q = []
        self.streams = []      # (src, phase, dst, frac); dst = int unit or ('P', label)
        self.feeds = []        # (unit, f, z(np), q)

    def add(self, typ, Q=0.0):
        self.types.append(typ)
        self.Q.append(float(Q))
        return len(self.types) - 1

    def link(self, src, phase, dst, frac):
        assert phase in OUTPH[self.types[src]], (src, phase)
        if not isinstance(dst, tuple):
            assert phase in ACCEPTS[self.types[dst]], (src, phase, dst)
        self.streams.append((src, phase, dst, float(frac)))

    def feed(self, unit, f, z, q):
        assert self.types[unit] == 'E'
        self.feeds.append((unit, float(f), np.asarray(z, float), float(q)))

    # ------------------------------------------------------------------ total flows
    def total_flows(self):
        K = len(self.types)
        # unknown vector: [Lo_0..Lo_{K-1}, Vo_0..Vo_{K-1}]
        A = np.zeros((2 * K, 2 * K))
        b = np.zeros(2 * K)
        for u in range(K):
            A[u, u] = 1.0
            A[K + u, K + u] = 1.0
        fq = np.zeros(K)
        fv = np.zeros(K)
        for (u, f, z, q) in self.feeds:
            fq[u] += q * f
            fv[u] += (1 - q) * f
        for (s, ph, d, fr) in self.streams:
            if isinstance(d, tuple):
                continue
            col = s if ph == 'L' else K + s
            t = self.types[d]
            if t == 'E':
                row = d if ph == 'L' else K + d
            elif t == 'ML':
                row = d
            elif t == 'MV':
                row = K + d
            elif t == 'C':
                row = d          # vapor in -> liquid out
            elif t == 'R':
                row = K + d      # liquid in -> vapor out
            A[row, col] -= fr
        for u, t in enumerate(self.types):
            if t == 'E':
                b[u] = fq[u] - self.Q[u]
                b[K + u] = fv[u] + self.Q[u]
            else:
                # zero-flow phases: force Lo or Vo = 0
                if 'L' not in OUTPH[t]:
                    A[u, :] = 0.0
                    A[u, u] = 1.0
                    b[u] = 0.0
                if 'V' not in OUTPH[t]:
                    A[K + u, :] = 0.0
                    A[K + u, K + u] = 1.0
                    b[K + u] = 0.0
        sol = np.linalg.solve(A, b)
        self.Lo = sol[:K].copy()
        self.Vo = sol[K:].copy()
        # duties of C and R units
        vin = np.zeros(K)
        lin = np.zeros(K)
        for (s, ph, d, fr) in self.streams:
            if isinstance(d, tuple):
                continue
            fl = fr * (self.Lo[s] if ph == 'L' else self.Vo[s])
            if ph == 'V':
                vin[d] += fl
            else:
                lin[d] += fl
        self.Qeff = np.array(self.Q, float)
        for u, t in enumerate(self.types):
            if t == 'C':
                self.Qeff[u] = -vin[u]
            elif t == 'R':
                self.Qeff[u] = lin[u]
            elif t in ('ML', 'MV'):
                self.Qeff[u] = 0.0
        return self.Lo, self.Vo

    # ------------------------------------------------------------------ compositions
    def _setup(self):
        K = len(self.types)
        self.two = [u for u, t in enumerate(self.types) if t == 'E' and self.Lo[u] > 1e-14 and self.Vo[u] > 1e-14]
        self.idx = {u: k for k, u in enumerate(self.two)}
        cV = np.zeros((K, K))   # cV[d, s] = sum of fractions of vapor of s going to unit d
        cL = np.zeros((K, K))
        for (s, ph, d, fr) in self.streams:
            if isinstance(d, tuple):
                continue
            if ph == 'V':
                cV[d, s] += fr
            else:
                cL[d, s] += fr
        self.cV, self.cL = cV, cL
        self.bfeed = np.zeros((K, self.N))
        for (u, f, z, q) in self.feeds:
            self.bfeed[u] += f * z

    def shares(self, S):
        """vapor share beta[u, i] of component i leaving unit u."""
        K = len(self.types)
        beta = np.zeros((K, self.N))
        dbeta = np.zeros((K, self.N))
        for u, t in enumerate(self.types):
            if t in ('MV', 'R'):
                beta[u] = 1.0
            elif t == 'E':
                if u in self.idx:
                    Su = S[self.idx[u]]
                    den = self.Vo[u] * self.a + self.Lo[u] * Su
                    beta[u] = self.Vo[u] * self.a / den
                    dbeta[u] = -self.Vo[u] * self.a * self.Lo[u] / den ** 2
                else:
                    beta[u] = 1.0 if self.Vo[u] > 1e-14 else 0.0
        return beta, dbeta

    def comp_flows(self, S, want_jac=False):
        K = len(self.types)
        beta, dbeta = self.shares(S)
        n = np.zeros((K, self.N))
        Ms = []
        for i in range(self.N):
            Ai = self.cV * beta[:, i][None, :] + self.cL * (1 - beta[:, i])[None, :]
            Mi = np.eye(K) - Ai
            n[:, i] = np.linalg.solve(Mi, self.bfeed[:, i])
            Ms.append(Mi)
        r = np.array([np.dot(beta[u], n[u]) / self.Vo[u] - 1.0 for u in self.two])
        if not want_jac:
            return n, r, beta
        m = len(self.two)
        J = np.zeros((m, m))
        Dm = (self.cV - self.cL)[:, self.two]            # K x m
        for i in range(self.N):
            Y = np.linalg.solve(Ms[i], Dm)                # dn_i/dS_u = Y[:, u] * dbeta[u,i] * n[u,i]
            coef = np.array([dbeta[u, i] * n[u, i] for u in self.two])
            for w_k, w in enumerate(self.two):
                J[w_k, :] += beta[w, i] * Y[w, :] * coef / self.Vo[w]
                J[w_k, w_k] += dbeta[w, i] * n[w, i] / self.Vo[w]
        J = J * S[None, :]                                # derivative w.r.t. log S
        return n, r, beta, J

    def solve(self, tol=1e-12, maxit=200, S0=None, verbose=False, bp_iter=0):
        self.resid = np.inf
        self.total_flows()
        K = len(self.types)
        if np.any(self.Lo < -1e-12) or np.any(self.Vo < -1e-12):
            return False
        self._setup()
        m = len(self.two)
        if S0 is None:
            zbar = sum(f * z for (_, f, z, _) in self.feeds) / sum(f for (_, f, _, _) in self.feeds)
            S0 = float(np.dot(self.a, zbar))
        s = np.full(m, np.log(S0)) if np.isscalar(S0) else np.log(np.asarray(S0, float))
        ok = False
        # phase 1: bubble-point type fixed-point iteration S_u <- sum_i alpha_i x_i / sum_i x_i (robust, slow)
        if m:
            two = np.array(self.two)
            for it in range(bp_iter):
                S = np.exp(s)
                n, r, beta = self.comp_flows(S)
                if np.max(np.abs(r)) < 1e-4:
                    break
                x = (1 - beta[two]) * n[two] / self.Lo[two][:, None]
                Snew = (x @ self.a) / np.maximum(x.sum(axis=1), 1e-300)
                s = np.log(np.maximum(Snew, 1e-300))
        for it in range(maxit):
            S = np.exp(s)
            n, r, beta, J = self.comp_flows(S, want_jac=True)
            nr = np.max(np.abs(r)) if m else 0.0
            if verbose:
                print(it, nr)
            if nr < tol:
                ok = True
                break
            try:
                ds = np.linalg.solve(J, -r)
            except np.linalg.LinAlgError:
                ds = np.linalg.lstsq(J, -r, rcond=None)[0]
            mx = np.max(np.abs(ds))
            if mx > 0.5:
                ds *= 0.5 / mx
            lam = 1.0
            for _ in range(30):
                st = s + lam * ds
                _, r2, _ = self.comp_flows(np.exp(st))
                if np.max(np.abs(r2)) < (1 - 1e-4 * lam) * nr or lam < 1e-6:
                    break
                lam *= 0.5
            s = st
        self.S = np.exp(s)
        self.n, r, self.beta = self.comp_flows(self.S)
        self.resid = np.max(np.abs(r)) if m else 0.0
        if np.any(self.n < -1e-9):
            return False
        return ok or self.resid < 1e-9

    # ------------------------------------------------------------------ stream data
    def stream_data(self):
        """list of dicts: src, dst, phase, comp flows nvec, vapor v."""
        out = []
        for (s, ph, d, fr) in self.streams:
            sh = self.beta[s] if ph == 'V' else (1 - self.beta[s])
            nv = fr * sh * self.n[s]
            nv = np.maximum(nv, 0.0)
            out.append(dict(src=s, dst=d, phase=ph, n=nv, v=(nv.sum() if ph == 'V' else 0.0)))
        return out

    def R(self):
        return float(sum(max(q, 0.0) for q in self.Qeff))


def phi_stream(a, nv, v, th):
    return v - float(np.sum(a * nv / (a - th)))


def phi_feed(a, f, z, q, th):
    return f * ((1 - q) - float(np.sum(a * z / (a - th))))


def feed_roots(a, z, q):
    from scipy.optimize import brentq
    g = lambda t: float(np.sum(a * z / (a - t))) - (1 - q)
    roots = []
    for j in range(len(a) - 1):
        lo, hi = float(a[j + 1]), float(a[j])   # g increases from -inf to +inf on (lo, hi): bisection
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if mid <= lo or mid >= hi:
                break
            if g(mid) > 0:
                hi = mid
            else:
                lo = mid
        roots.append(0.5 * (lo + hi))
    return np.array(roots)


def peaks(a, z, th):
    N = len(a)
    return np.array([sum(a[p] * z[p] / (th[j] - a[p]) for p in range(j + 1, N)) for j in range(N - 1)])


def check_bound(net, thetas, tol=1e-8):
    """Return dict with node-balance residual, sign uniformity violations, bound slack (R - B_j), cut identity
    residual, for each theta in thetas.  Works with several feeds (feed term kept)."""
    a = net.a
    sd = net.stream_data()
    K = len(net.types)
    R = net.R()
    res = dict(R=R, node=0.0, signbad=0, slack=[], B=[], cutres=0.0, ffeed=[])
    flows = [s['n'].sum() for s in sd]
    scale = max(1.0, max(flows) if flows else 1.0)
    for th in thetas:
        ph = [phi_stream(a, s['n'], s['v'], th) for s in sd]
        pf = np.zeros(K)
        for (u, f, z, q) in net.feeds:
            pf[u] += phi_feed(a, f, z, q, th)
        bal = -np.array(net.Qeff) - pf
        for k, s in enumerate(sd):
            bal[s['src']] += ph[k]
            if not isinstance(s['dst'], tuple):
                bal[s['dst']] -= ph[k]
        res['node'] = max(res['node'], float(np.max(np.abs(bal))) / scale)
        # sign uniformity: classify units
        sgn = np.zeros(K)    # +1 heavy-ish (phi >= 0 for all outflows), -1 light-ish
        for u in range(K):
            vals = [ph[k] for k, s in enumerate(sd) if s['src'] == u and flows[k] > 1e-12 * scale]
            pos = any(v > tol * scale for v in vals)
            neg = any(v < -tol * scale for v in vals)
            if pos and neg:
                res['signbad'] += 1
            sgn[u] = -1 if (neg and not pos) else +1
        Up = sgn > 0
        B = sum(max(0.0, ph[k]) for k, s in enumerate(sd) if isinstance(s['dst'], tuple))
        Fterm = sum(max(0.0, phi_feed(a, f, z, q, th)) for (u, f, z, q) in net.feeds)
        # exact defect identity: R - B + (feed part entering U+) = R(U-) + Q(U+) + sum crossing |phi|
        RUm = sum(max(net.Qeff[u], 0.0) for u in range(K) if not Up[u])
        QUp = sum(max(-net.Qeff[u], 0.0) for u in range(K) if Up[u])
        cross = 0.0
        for k, s in enumerate(sd):
            d = s['dst']
            if isinstance(d, tuple):
                continue
            if Up[s['src']] != Up[d]:
                cross += abs(ph[k])
        fUp = sum(pf[u] for u in range(K) if Up[u])
        lhs = R - B + fUp
        rhs = RUm + QUp + cross
        res['cutres'] = max(res['cutres'], abs(lhs - rhs) / scale)
        res['slack'].append(R - B + Fterm)
        res['B'].append(B)
        res['ffeed'].append(Fterm)
    return res
