#!/usr/bin/env python3
"""
stagenet.py (Agent S) -- minimal steady-state solver for networks of equilibrium stages with constant
molar overflow (CMO) and constant relative volatilities, used to test the section/junction lemmas of
sections.md on genuine finite-stage physics (no Underwood formulation inside).

Network = nodes + streams.
  node types: 'stage' -- equilibrium stage: all outgoing liquid has composition x, all outgoing vapor has
                         y_i = alpha_i x_i / sum_k alpha_k x_k;
              'mix'   -- perfect mixer/splitter/total condenser/total reboiler: every outgoing stream (any
                         phase) has the composition of the mixed inflow.
  streams: (src, dst, phase, flow); src = node index or 'F' (a piece of the process feed, phase 'F',
           flow = fraction of the feed, same composition z and quality q); dst = node index or ('P', name).
The flows are inputs (CMO: they do not depend on compositions); the heat duty of a node is whatever closes its
vapor balance: netvap(node) = V_out - V_in - (1-q)*(feed pieces in).  R = sum of positive netvap (vapor
generated), Q = sum of negative netvap (vapor condensed).
Compositions: Gauss-Seidel sweeps over the nodes (exact flash at each stage), then damped Newton with the
analytic Jacobian.  Environment: Python 3, numpy, scipy.
"""
import numpy as np
from scipy.optimize import brentq


class Net:
    def __init__(self, alpha, z, q):
        self.alpha = np.asarray(alpha, float)
        self.z = np.asarray(z, float)
        self.q = float(q)
        self.N = len(self.alpha)
        self.nodes = []
        self.streams = []

    def add_node(self, typ, name=None):
        self.nodes.append({'type': typ, 'name': name or f'n{len(self.nodes)}'})
        return len(self.nodes) - 1

    def add_stream(self, src, dst, phase, flow):
        assert flow >= 0
        self.streams.append((src, dst, phase, float(flow)))

    # ---------------------------------------------------------------- flows and duties
    def balances(self):
        n = len(self.nodes)
        Lin = np.zeros(n); Vin = np.zeros(n); Lout = np.zeros(n); Vout = np.zeros(n); Fin = np.zeros(n)
        for (s, d, ph, f) in self.streams:
            if s == 'F':
                Fin[d] += f
                continue
            if ph == 'L':
                Lout[s] += f
            else:
                Vout[s] += f
            if isinstance(d, int):
                if ph == 'L':
                    Lin[d] += f
                else:
                    Vin[d] += f
        tot = Lin + Vin + Fin - Lout - Vout
        netvap = Vout - Vin - (1 - self.q) * Fin
        return dict(Lin=Lin, Vin=Vin, Lout=Lout, Vout=Vout, Fin=Fin, tot=tot, netvap=netvap)

    def duties(self):
        nv = self.balances()['netvap']
        return float(np.sum(np.maximum(nv, 0))), float(np.sum(np.maximum(-nv, 0))), nv

    # ---------------------------------------------------------------- composition model
    def outcomp(self, k, c, ph):
        """composition of a stream of phase ph leaving node k whose state is c"""
        if self.nodes[k]['type'] == 'stage' and ph == 'V':
            y = self.alpha * c
            return y / y.sum()
        return c

    def _prep(self):
        n = len(self.nodes)
        self._in = [[] for _ in range(n)]
        self._out = [[] for _ in range(n)]
        for (s, d, ph, f) in self.streams:
            if isinstance(d, int):
                self._in[d].append((s, ph, f))
            if s != 'F':
                self._out[s].append((d, ph, f))

    def inflow(self, k, C):
        m = np.zeros(self.N)
        for (s, ph, f) in self._in[k]:
            if s == 'F':
                m += f * self.z
            else:
                m += f * self.outcomp(s, C[s], ph)
        return m

    def node_update(self, k, C, b):
        m = self.inflow(k, C)
        M = m.sum()
        if M <= 0:
            return C[k]
        if self.nodes[k]['type'] == 'mix' or b['Vout'][k] == 0 or b['Lout'][k] == 0:
            if self.nodes[k]['type'] == 'stage' and b['Lout'][k] == 0 and b['Vout'][k] > 0:
                y = m / M
                x = y / self.alpha
                return x / x.sum()
            return m / M
        L, V = b['Lout'][k], b['Vout'][k]
        a = self.alpha
        g = lambda s: np.sum(a * m / (L * s + V * a)) - 1.0
        s = brentq(g, 1e-14, a.max() * 10, xtol=1e-15, rtol=1e-15)
        x = m / (L + V * a / s)
        return x / x.sum()

    def res_jac(self, C):
        n, N, a = len(self.nodes), self.N, self.alpha
        res = np.zeros((n, N))
        J = np.zeros((n * N, n * N))
        I = np.eye(N)
        dy = {}
        for k in range(n):
            if self.nodes[k]['type'] == 'stage':
                x = C[k]; s = a @ x
                dy[k] = (np.diag(a) * s - np.outer(a * x, a)) / s ** 2

        def cj(src, ph):
            if self.nodes[src]['type'] == 'stage' and ph == 'V':
                y = a * C[src]
                return y / y.sum(), dy[src]
            return C[src], I
        for k in range(n):
            for (s, ph, f) in self._in[k]:
                if s == 'F':
                    res[k] += f * self.z
                else:
                    c, dc = cj(s, ph)
                    res[k] += f * c
                    J[k * N:(k + 1) * N, s * N:(s + 1) * N] += f * dc
            for (d, ph, f) in self._out[k]:
                c, dc = cj(k, ph)
                res[k] -= f * c
                J[k * N:(k + 1) * N, k * N:(k + 1) * N] -= f * dc
        # the component balances of a node sum to its (already closed) molar balance: replace the last one
        # by the normalisation sum c = 1
        for k in range(n):
            res[k, -1] = C[k].sum() - 1.0
            J[k * N + N - 1, :] = 0.0
            J[k * N + N - 1, k * N:(k + 1) * N] = 1.0
        return res.ravel(), J

    def residual(self, Cflat):
        C = Cflat.reshape(len(self.nodes), self.N)
        res = np.zeros_like(C)
        for k in range(len(self.nodes)):
            res[k] = self.inflow(k, C)
            for (d, ph, f) in self._out[k]:
                res[k] -= f * self.outcomp(k, C[k], ph)
        return res.ravel()

    def solve(self, maxsweep=3000, tol=1e-14, C0=None):
        b = self.balances()
        self._b = b
        self._prep()
        assert np.max(np.abs(b['tot'])) < 1e-9, 'total molar balance violated'
        n = len(self.nodes)
        C = np.tile(self.z, (n, 1)) if C0 is None else np.array(C0, float)
        for it in range(maxsweep):
            Cold = C.copy()
            for k in range(n):
                C[k] = self.node_update(k, C, b)
            if np.max(np.abs(C - Cold)) < 1e-10:
                break
        # damped Newton in logarithmic variables u = log c (keeps every composition positive, handles traces)
        C = np.maximum(C, 1e-300)
        U = np.log(C)
        for it in range(400):
            r, J = self.res_jac(np.exp(U))
            nr = np.max(np.abs(r))
            if nr < tol:
                break
            Ju = J * np.exp(U).ravel()[None, :]
            try:
                du = np.linalg.solve(Ju, -r).reshape(n, self.N)
            except np.linalg.LinAlgError:
                break
            du = np.clip(du, -30.0, 30.0)
            t = 1.0
            ok = False
            while t > 1e-10:
                Un = U + t * du
                rn, _ = self.res_jac(np.exp(Un))
                if np.max(np.abs(rn)) < nr:
                    ok = True
                    break
                t /= 2
            if not ok:
                break
            U = Un
        C = np.exp(U)
        self.C = C
        self.res = np.max(np.abs(self.residual(C.ravel())))
        return C

    # ---------------------------------------------------------------- Underwood quantities
    def feed_roots(self):
        a, z, q = self.alpha, self.z, self.q
        f = lambda t: np.sum(a * z / (a - t)) - (1 - q)
        roots = []
        for j in range(self.N - 1):
            lo, hi = a[j + 1], a[j]
            eps = 1e-13 * (hi - lo)
            roots.append(brentq(f, lo + eps, hi - eps, xtol=1e-15, rtol=1e-15, maxiter=500))
        return roots

    def T(self, comp_flows, ph, th):
        a = self.alpha
        if ph == 'L':
            return float(np.sum(a * comp_flows / (a - th)))
        return float(np.sum(th * comp_flows / (a - th)))

    def stream_T(self, st, th):
        s, d, ph, f = st
        if s == 'F':
            return f * (float(np.sum(self.alpha * self.z / (self.alpha - th))) - (1 - self.q))
        return self.T(f * self.outcomp(s, self.C[s], ph), ph, th)

    def products(self):
        return [st for st in self.streams if not isinstance(st[1], int)]
