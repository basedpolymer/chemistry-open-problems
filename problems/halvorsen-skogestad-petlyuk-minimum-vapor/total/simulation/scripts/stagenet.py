"""
stagenet.py -- independent stage-by-stage (MESH) model of arbitrary adiabatic distillation networks.

Agent P (Claude Opus 5.5, effort xhigh).  No Underwood quantity is used anywhere in this module.

Physical model (ideal mixture, constant relative volatilities, constant molar overflow, constant pressure):
  * every stage is an equilibrium stage: y_p = alpha_p x_p / sum_k alpha_k x_k;
  * a NODE is one equilibrium stage that may receive any set of liquid / vapor streams, the process feed,
    and an external heat duty Q (Q > 0 vapor generated = reboiler/heater, Q < 0 vapor condensed =
    condenser/cooler); its liquid (composition x) and vapor (composition y) outflows may be divided among
    any number of destinations (composition-preserving splitters are therefore automatic);
  * a SECTION is a chain of n adiabatic equilibrium stages with constant vapor V and liquid L (CMO); its top
    stage receives liquid from its top node and sends vapor to it, its bottom stage receives vapor from its
    bottom node and sends liquid to it.  n = 0 is allowed (direct two-way link between the two nodes);
  * a STREAM carries liquid or vapor of a node to another node (bypasses, transfer streams, recycles);
  * a PRODUCT withdraws liquid or vapor from a node.
Flows (section V and L, stream flows, product flows, node duties Q) must satisfy at every node
  V_out - V_in = Q,   L_out - L_in = -Q   (total moles conserved; Q converts liquid into vapor).
The objective of the conjecture is R = sum of positive duties (total vapor generated).

Given the flows, the steady-state compositions solve the component balances at every stage
  sum_in f * c_src + F z [node with feed] - L_out x - V_out y(x) = 0,
solved here by successive substitution on the K-values (bubble-point method; linear systems per component)
followed by Newton's method on the full sparse system.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class Network:
    def __init__(self, alpha):
        self.alpha = np.asarray(alpha, dtype=float)
        self.N = len(self.alpha)
        self.nodes = []          # node names
        self.node_kind = []      # 'adiabatic' | 'heater' (Q>=0) | 'cooler' (Q<=0) | 'free'
        self.sections = []       # dicts: name, n, top, bot
        self.streams = []        # dicts: name, src, phase, dst
        self.feeds = []          # dicts: node, F, q, z
        self.products = []       # dicts: name, node, phase, comp
        self._built = False

    # ---------------------------------------------------------------- construction
    def node(self, name, kind='adiabatic'):
        assert name not in self.nodes, name
        self.nodes.append(name)
        self.node_kind.append(kind)
        return len(self.nodes) - 1

    def nid(self, name):
        return self.nodes.index(name) if isinstance(name, str) else name

    def section(self, name, n, top, bot):
        self.sections.append(dict(name=name, n=int(n), top=self.nid(top), bot=self.nid(bot)))
        self._built = False

    def stream(self, name, src, phase, dst):
        assert phase in 'LV'
        self.streams.append(dict(name=name, src=self.nid(src), phase=phase, dst=self.nid(dst)))
        self._built = False

    def feed(self, node, F, q, z):
        self.feeds.append(dict(node=self.nid(node), F=float(F), q=float(q), z=np.asarray(z, float)))
        self._built = False

    def product(self, name, node, phase, comp):
        self.products.append(dict(name=name, node=self.nid(node), phase=phase, comp=int(comp)))
        self._built = False

    # ---------------------------------------------------------------- flow vector layout
    def build(self):
        """Index the flow vector and the stages; build the static connection lists."""
        names = []
        for s in self.sections:
            names += ['V:' + s['name'], 'L:' + s['name']]
        for e in self.streams:
            names.append('S:' + e['name'])
        for p in self.products:
            names.append('P:' + p['name'])
        for k, nm in enumerate(self.nodes):
            names.append('Q:' + nm)
        self.fnames = names
        self.fidx = {nm: i for i, nm in enumerate(names)}
        self.nf = len(names)
        # stages: nodes first, then section stages
        K = len(self.nodes)
        self.nst = K
        for s in self.sections:
            s['first'] = self.nst
            self.nst += s['n']
        # connections: (dst stage, src stage, phase 0=L 1=V, flow index)
        dst, src, ph, fi = [], [], [], []
        for s in self.sections:
            iV, iL = self.fidx['V:' + s['name']], self.fidx['L:' + s['name']]
            n, a = s['n'], s['first']
            top, bot = s['top'], s['bot']
            if n == 0:
                dst += [top, bot]; src += [bot, top]; ph += [1, 0]; fi += [iV, iL]
                continue
            for i in range(n):
                st = a + i
                up = top if i == 0 else st - 1          # liquid comes from above
                dn = bot if i == n - 1 else st + 1      # vapor comes from below
                dst += [st, st]; src += [up, dn]; ph += [0, 1]; fi += [iL, iV]
            dst += [top, bot]; src += [a, a + n - 1]; ph += [1, 0]; fi += [iV, iL]
        for e in self.streams:
            dst.append(e['dst']); src.append(e['src']); ph.append(0 if e['phase'] == 'L' else 1)
            fi.append(self.fidx['S:' + e['name']])
        self.c_dst = np.array(dst, int); self.c_src = np.array(src, int)
        self.c_ph = np.array(ph, int); self.c_fi = np.array(fi, int)
        # outflow coefficient lists: Lout[stage] = sum f[idx]; Vout[stage] = sum f[idx]
        Lo, Vo = [], []    # (stage, flow index)
        for s in self.sections:
            iV, iL = self.fidx['V:' + s['name']], self.fidx['L:' + s['name']]
            for i in range(s['n']):
                Lo.append((s['first'] + i, iL)); Vo.append((s['first'] + i, iV))
            Lo.append((s['top'], iL)); Vo.append((s['bot'], iV))
        for e in self.streams:
            (Lo if e['phase'] == 'L' else Vo).append((e['src'], self.fidx['S:' + e['name']]))
        for p in self.products:
            (Lo if p['phase'] == 'L' else Vo).append((p['node'], self.fidx['P:' + p['name']]))
        self.Lo = np.array(Lo, int).reshape(-1, 2); self.Vo = np.array(Vo, int).reshape(-1, 2)
        # feed vector per stage
        self.feedvec = np.zeros((self.nst, self.N))
        for fd in self.feeds:
            self.feedvec[fd['node']] += fd['F'] * fd['z']
        self._built = True

    # ---------------------------------------------------------------- linear flow balances
    def balance_matrix(self):
        """Rows: for each node, vapor balance V_out - V_in - Q = Fvap_in, liquid balance L_out - L_in + Q = Fliq_in.
        Returns (A, b) with A f = b."""
        if not self._built:
            self.build()
        K = len(self.nodes)
        A = np.zeros((2 * K, self.nf)); b = np.zeros(2 * K)
        for s in self.sections:
            iV, iL = self.fidx['V:' + s['name']], self.fidx['L:' + s['name']]
            t, bo = s['top'], s['bot']
            A[2 * t, iV] -= 1      # vapor enters top node
            A[2 * t + 1, iL] += 1  # liquid leaves top node
            A[2 * bo, iV] += 1     # vapor leaves bottom node
            A[2 * bo + 1, iL] -= 1  # liquid enters bottom node
        for e in self.streams:
            i = self.fidx['S:' + e['name']]
            r = 0 if e['phase'] == 'V' else 1
            A[2 * e['src'] + r, i] += 1
            A[2 * e['dst'] + r, i] -= 1
        for p in self.products:
            i = self.fidx['P:' + p['name']]
            r = 0 if p['phase'] == 'V' else 1
            A[2 * p['node'] + r, i] += 1
        for k in range(K):
            i = self.fidx['Q:' + self.nodes[k]]
            A[2 * k, i] -= 1
            A[2 * k + 1, i] += 1
        for fd in self.feeds:
            b[2 * fd['node']] += fd['F'] * (1 - fd['q'])
            b[2 * fd['node'] + 1] += fd['F'] * fd['q']
        return A, b

    def flows_from(self, d):
        f = np.zeros(self.nf)
        for k, v in d.items():
            f[self.fidx[k]] = v
        return f

    def R(self, f):
        return sum(max(f[self.fidx['Q:' + nm]], 0.0) for nm in self.nodes)

    def check_flows(self, f, tol=1e-9):
        A, b = self.balance_matrix()
        res = A @ f - b
        return np.max(np.abs(res)) if len(res) else 0.0

    # ---------------------------------------------------------------- composition solve
    def outflows(self, f):
        Lout = np.zeros(self.nst); Vout = np.zeros(self.nst)
        np.add.at(Lout, self.Lo[:, 0], f[self.Lo[:, 1]])
        np.add.at(Vout, self.Vo[:, 0], f[self.Vo[:, 1]])
        return Lout, Vout

    def yof(self, x):
        ax = x * self.alpha
        return ax / ax.sum(axis=1, keepdims=True)

    def residual(self, x, f, Lout=None, Vout=None):
        if Lout is None:
            Lout, Vout = self.outflows(f)
        y = self.yof(x)
        c = np.where(self.c_ph[:, None] == 0, x[self.c_src], y[self.c_src])
        r = self.feedvec - Lout[:, None] * x - Vout[:, None] * y
        np.add.at(r, self.c_dst, f[self.c_fi][:, None] * c)
        dead = (Lout + Vout) <= 1e-14
        if np.any(dead):     # stages without any flow: composition is immaterial, pin it
            r[dead] = 1.0 / self.N - x[dead]
        return r

    def jacobian(self, x, f, Lout=None, Vout=None):
        if Lout is None:
            Lout, Vout = self.outflows(f)
        N = self.N; a = self.alpha
        S = (x * a).sum(axis=1)
        y = (x * a) / S[:, None]
        # Jy[s] = (diag(a) - y a^T)/S
        Jy = (np.diag(a)[None, :, :] - y[:, :, None] * a[None, None, :]) / S[:, None, None]
        rows, cols, vals = [], [], []
        I = np.eye(N)
        # diagonal blocks
        D = -Lout[:, None, None] * I[None] - Vout[:, None, None] * Jy
        dead = (Lout + Vout) <= 1e-14
        if np.any(dead):
            D[dead] = -I[None]
        st = np.arange(self.nst)
        ii, jj = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        rows.append((st[:, None, None] * N + ii[None]).ravel()); cols.append((st[:, None, None] * N + jj[None]).ravel())
        vals.append(D.ravel())
        # connections
        fc = f[self.c_fi]
        B = np.where(self.c_ph[:, None, None] == 0, I[None], Jy[self.c_src]) * fc[:, None, None]
        rows.append((self.c_dst[:, None, None] * N + ii[None]).ravel())
        cols.append((self.c_src[:, None, None] * N + jj[None]).ravel())
        vals.append(B.ravel())
        M = self.nst * N
        J = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(M, M))
        return J

    def bp_sweep(self, x, f, Lout, Vout, iters=50, tol=1e-10):
        """Successive substitution on K_s = alpha/S_s: for fixed S the balances are linear per component."""
        N = self.N; a = self.alpha
        nst = self.nst
        for it in range(iters):
            S = (x * a).sum(axis=1)
            xn = np.zeros_like(x)
            fc = f[self.c_fi]
            for p in range(N):
                Kp = a[p] / S
                diag = -(Lout + Vout * Kp)
                diag = np.where((Lout + Vout) <= 1e-14, -1.0, diag)
                off = fc * np.where(self.c_ph == 0, 1.0, Kp[self.c_src])
                M = sp.csc_matrix((np.concatenate([diag, off]),
                                   (np.concatenate([np.arange(nst), self.c_dst]),
                                    np.concatenate([np.arange(nst), self.c_src]))), shape=(nst, nst))
                rhs = -self.feedvec[:, p].copy()
                rhs[(Lout + Vout) <= 1e-14] = -1.0 / N
                xn[:, p] = spla.spsolve(M, rhs)
            xn = np.maximum(xn, 1e-300)
            xn /= xn.sum(axis=1, keepdims=True)
            ch = np.max(np.abs(xn - x))
            x = xn
            if ch < tol:
                break
        return x

    def weights(self, x, f, Lout, Vout):
        """Gross throughput of each component at each stage (row scaling of the balances)."""
        y = self.yof(x)
        w = Lout[:, None] * x + Vout[:, None] * y + self.feedvec
        c = np.where(self.c_ph[:, None] == 0, x[self.c_src], y[self.c_src])
        np.add.at(w, self.c_dst, f[self.c_fi][:, None] * c)
        dead = (Lout + Vout) <= 1e-14
        if np.any(dead):
            w[dead] = 1.0
        # components below TAU (mole fraction) are resolved in absolute terms only
        return np.maximum(w, self.TAU * np.maximum(Lout + Vout, 1e-300)[:, None])

    def scaled_jacobian(self, x, f, Lout, Vout, w):
        """Js = diag(1/w) J diag(x): Newton in u = log x with relative residuals (high-purity safe)."""
        J = self.jacobian(x, f, Lout, Vout)
        return sp.csc_matrix(sp.diags(1.0 / w.ravel()) @ J @ sp.diags(np.maximum(x, self.TAU).ravel()))

    XMIN = 1e-250
    TAU = 1e-30

    def _update(self, x, du, t=1.0):
        """x_new from a step du in the scaled variables (dx = max(x, TAU) du), positivity-preserving."""
        rel = t * du * np.maximum(x, self.TAU) / x
        return np.clip(x * np.exp(np.clip(rel, -30, 3)), self.XMIN, 1.5)

    def ptc(self, x, f, Lout, Vout, tol=1e-12, maxit=400, dt0=0.5, verbose=False):
        """Pseudo-transient continuation: implicit Euler on the composition dynamics  dx/dt = r(x)
        (unit holdup per stage, constant flows), in log variables, time step grown by switched evolution
        relaxation (-> Newton)."""
        w = self.weights(x, f, Lout, Vout)
        r = self.residual(x, f, Lout, Vout) / w
        nr = np.linalg.norm(r)
        dt = dt0
        for it in range(maxit):
            if np.max(np.abs(r)) < tol:
                break
            Js = self.scaled_jacobian(x, f, Lout, Vout, w)
            A = sp.csc_matrix(sp.diags((np.maximum(x, self.TAU) / (w * dt)).ravel()) - Js)
            try:
                du = spla.spsolve(A, r.ravel()).reshape(x.shape)
            except Exception:
                dt /= 4
                continue
            if not np.all(np.isfinite(du)):
                dt /= 4
                continue
            xn = self._update(x, du)
            wn = self.weights(xn, f, Lout, Vout)
            rn = self.residual(xn, f, Lout, Vout) / wn
            nn = np.linalg.norm(rn)
            if nn < 5 * nr or np.max(np.abs(rn)) < tol:
                dt = min(1e15, dt * max(1.5 if nn < nr else 0.7, min(10.0, nr / max(nn, 1e-300))))
                x, r, nr, w = xn, rn, nn, wn
            else:
                dt /= 4
                if dt < 1e-12:
                    break
            if verbose:
                print('ptc', it, nr, dt)
        return x

    def solve(self, f, x0=None, tol=1e-11, maxit=60, bp_iters=30, verbose=False, ptc_fallback=2000):
        """Steady-state compositions for the given flows (relative residual tol on every component of every
        stage).  Returns (x, info)."""
        if not self._built:
            self.build()
        Lout, Vout = self.outflows(f)
        if x0 is None:
            x = np.full((self.nst, self.N), 1.0 / self.N)
            try:
                xb = self.bp_sweep(x, f, Lout, Vout, iters=bp_iters)
                xb = np.maximum(xb, 1e-30); xb /= xb.sum(axis=1, keepdims=True)
                if not np.all(np.isfinite(xb)):
                    raise ValueError
                x = xb
            except Exception:
                pass
            x = self.ptc(x, f, Lout, Vout, tol=tol, maxit=400, verbose=verbose)
            w = self.weights(x, f, Lout, Vout)
            if not np.max(np.abs(self.residual(x, f, Lout, Vout) / w)) < 1e-6:
                x = self.ptc(np.full((self.nst, self.N), 1.0 / self.N), f, Lout, Vout, tol=tol, maxit=3000)
        else:
            x = np.maximum(np.array(x0, float), self.XMIN)
        bad = ~np.all(np.isfinite(x), axis=1)
        if np.any(bad):
            x[bad] = 1.0 / self.N
        w = self.weights(x, f, Lout, Vout)
        r = self.residual(x, f, Lout, Vout) / w
        nr = np.max(np.abs(r))
        atol = 1e-13 * max(1.0, np.max(np.abs(self.feedvec)))
        na = np.max(np.abs(r * w))
        it = 0
        t = 1.0
        nfail = 0
        slow = 0
        for it in range(maxit):
            if nr < tol and na < atol:
                break
            Js = self.scaled_jacobian(x, f, Lout, Vout, w)
            try:
                lu = spla.splu(Js)
                du = lu.solve(-r.ravel()).reshape(x.shape)
                okstep = np.all(np.isfinite(du))
            except Exception:
                okstep = False
            if not okstep:
                x = self.ptc(x, f, Lout, Vout, tol=tol, maxit=ptc_fallback)
                w = self.weights(x, f, Lout, Vout)
                r = self.residual(x, f, Lout, Vout) / w; nr = np.max(np.abs(r)); na = np.max(np.abs(r * w))
                nfail += 1
                if nfail > 3:
                    break
                continue
            ndu = np.max(np.abs(du))
            if ndu < 1e-15:
                break
            # affine-invariant damping (Deuflhard): accept t if the simplified Newton correction shrinks
            t = min(1.0, 4 * t)
            ok = False
            for _ in range(40):
                xn = self._update(x, du, t)
                rn = self.residual(xn, f, Lout, Vout) / w
                dub = lu.solve(-rn.ravel())
                if np.all(np.isfinite(dub)) and np.max(np.abs(dub)) <= (1 - t / 4) * ndu:
                    ok = True
                    break
                t *= 0.5
            if not ok:
                x = self.ptc(x, f, Lout, Vout, tol=tol, maxit=ptc_fallback)
                w = self.weights(x, f, Lout, Vout)
                r = self.residual(x, f, Lout, Vout) / w; nr = np.max(np.abs(r)); na = np.max(np.abs(r * w))
                t = 1.0
                nfail += 1
                if nfail > 3:
                    break
                continue
            x = xn
            w = self.weights(x, f, Lout, Vout)
            r = self.residual(x, f, Lout, Vout) / w
            nr = np.max(np.abs(r)); na = np.max(np.abs(r * w))
            slow = slow + 1 if t < 1e-2 else 0
            if slow >= 5:
                x = self.ptc(x, f, Lout, Vout, tol=tol, maxit=ptc_fallback)
                w = self.weights(x, f, Lout, Vout)
                r = self.residual(x, f, Lout, Vout) / w; nr = np.max(np.abs(r)); na = np.max(np.abs(r * w))
                slow = 0; t = 1.0
                nfail += 1
                if nfail > 3:
                    break
            if verbose:
                print(it, nr, t, ndu)
        return x, dict(res=nr, it=it, ok=bool(nr < max(tol * 100, 1e-9)))

    def interpolate_from(self, other, x_other):
        """Initial guess for this network from a solution of a network with the same nodes and sections but other
        stage counts (profiles stretched by relative position)."""
        x = np.full((self.nst, self.N), 1.0 / self.N)
        K = len(self.nodes)
        x[:K] = x_other[:K]
        osec = {s['name']: s for s in other.sections}
        for s in self.sections:
            o = osec[s['name']]
            if s['n'] == 0:
                continue
            if o['n'] == 0:
                for i in range(s['n']):
                    t = (i + 1) / (s['n'] + 1)
                    x[s['first'] + i] = (1 - t) * x_other[o['top']] + t * x_other[o['bot']]
                continue
            for i in range(s['n']):
                j = min(o['n'] - 1, int((i + 0.5) * o['n'] / s['n']))
                x[s['first'] + i] = x_other[o['first'] + j]
        return x

    # ---------------------------------------------------------------- sensitivities
    def dres_dflow(self, x, f):
        """Sparse matrix G = d r / d f  (rows: stage*N + comp, cols: flow index)."""
        N = self.N
        y = self.yof(x)
        rows, cols, vals = [], [], []
        comp = np.arange(N)
        c = np.where(self.c_ph[:, None] == 0, x[self.c_src], y[self.c_src])
        rows.append((self.c_dst[:, None] * N + comp[None]).ravel())
        cols.append(np.repeat(self.c_fi, N)); vals.append(c.ravel())
        rows.append((self.Lo[:, 0][:, None] * N + comp[None]).ravel())
        cols.append(np.repeat(self.Lo[:, 1], N)); vals.append(-x[self.Lo[:, 0]].ravel())
        rows.append((self.Vo[:, 0][:, None] * N + comp[None]).ravel())
        cols.append(np.repeat(self.Vo[:, 1], N)); vals.append(-y[self.Vo[:, 0]].ravel())
        return sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.nst * N, self.nf))

    def log_impurity_and_grad_x(self, x, pr):
        """log(impurity) of a product and its gradient w.r.t. the node composition (length N)."""
        k, p = pr['node'], pr['comp']
        a = self.alpha
        if pr['phase'] == 'L':
            c = x[k]
            imp = np.sum(np.delete(c, p))
            g = np.ones(self.N); g[p] = 0.0
            # make it scale-invariant in x (sum x = 1 at the solution): imp/sum(x)
            s = np.sum(c)
            gi = g / s - imp / s ** 2
            return np.log(imp / s), gi / (imp / s)
        else:
            S = np.dot(a, x[k]); y = a * x[k] / S
            imp = 1.0 - y[p]
            imp = np.sum(np.delete(y, p))
            Jy = (np.diag(a) - np.outer(y, a)) / S
            gi = -Jy[p]
            return np.log(imp), gi / imp

    # ---------------------------------------------------------------- helpers
    def node_x(self, x, name):
        return x[self.nid(name)]

    def node_y(self, x, name):
        return self.yof(x[[self.nid(name)]])[0]

    def product_comp(self, x, pr):
        k = pr['node']
        return x[k] if pr['phase'] == 'L' else self.yof(x[[k]])[0]

    def impurities(self, x):
        out = {}
        for pr in self.products:
            c = self.product_comp(x, pr)
            out[pr['name']] = float(np.sum(np.delete(c, pr['comp'])))
        return out

    def section_netflows(self, x, f):
        """Upward net component flows w = V y(stage below top end) - L x(stage) ... evaluated at the top end:
        w = V y_1 - L x_top_node (vapor leaving the top stage minus liquid entering it)."""
        y = self.yof(x)
        out = {}
        for s in self.sections:
            V = f[self.fidx['V:' + s['name']]]; L = f[self.fidx['L:' + s['name']]]
            ytop = y[s['first']] if s['n'] > 0 else y[s['bot']]
            out[s['name']] = V * ytop - L * x[s['top']]
        return out
