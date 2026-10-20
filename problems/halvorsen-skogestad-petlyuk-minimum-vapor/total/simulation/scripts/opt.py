"""
opt.py -- minimum total vaporization of a stage-by-stage network at given product purities.

Agent P (Claude Opus 5.5, effort xhigh).  Reduced-space SQP: the flows are parametrized on the affine space of
the node balances, f = f0 + Z u (exact balances at every iterate), every evaluation is a converged steady state
of the MESH model of stagenet.py, and the gradients of the purity constraints come from the adjoint of the
composition equations (implicit function theorem).  No Underwood quantity is used.

  minimize   R(f) = sum of heater/reboiler duties
  subject to node balances (exact), bounds on flows, impurity_k(x(f)) <= eps_k for every product k.
"""
import numpy as np
import scipy.optimize as so
import scipy.linalg as sl
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def flow_bounds(net, vmin=1e-6, pmin=0.0):
    nf = net.nf
    lb = np.zeros(nf); ub = np.full(nf, np.inf)
    for i, nm in enumerate(net.fnames):
        if nm.startswith('V:') or nm.startswith('L:'):
            lb[i] = vmin
        elif nm.startswith('S:'):
            lb[i] = 0.0
        elif nm.startswith('P:'):
            lb[i] = pmin
    fixed = {}
    for k, nm in enumerate(net.nodes):
        i = net.fidx['Q:' + nm]
        kind = net.node_kind[k]
        if kind == 'adiabatic':
            fixed[i] = 0.0
            lb[i] = -np.inf
        elif kind == 'heater':
            lb[i] = 0.0
        elif kind == 'cooler':
            lb[i] = -np.inf; ub[i] = 0.0
        else:
            lb[i] = -np.inf
    return lb, ub, fixed


def equality_system(net, fixed_idx, extra=()):
    A, b = net.balance_matrix()
    rows = [A]; rhs = [b]
    for i, v in fixed_idx.items():
        r = np.zeros(net.nf); r[i] = 1.0
        rows.append(r[None]); rhs.append(np.array([v]))
    for coef, v in extra:
        r = np.zeros(net.nf)
        for nm, c in coef.items():
            r[net.fidx[nm]] += c
        rows.append(r[None]); rhs.append(np.array([v]))
    return np.vstack(rows), np.concatenate(rhs)


def heater_vector(net):
    c = np.zeros(net.nf)
    for k, nm in enumerate(net.nodes):
        if net.node_kind[k] == 'heater':
            c[net.fidx['Q:' + nm]] = 1.0
    return c


def lp_init(net, extra=(), vmin=0.3, lmin=0.3, objective_weightV=1e-3, fixed_names=None):
    """A flow vector satisfying all balances, bounds V,L >= vmin/lmin and the extra equalities
    (list of ({flowname: coef}, value)); minimizes R + w * sum V."""
    lb, ub, fixed = flow_bounds(net)
    for i, nm in enumerate(net.fnames):
        if nm.startswith('V:'):
            lb[i] = max(lb[i], vmin)
        if nm.startswith('L:'):
            lb[i] = max(lb[i], lmin)
    if fixed_names:
        for nm, v in fixed_names.items():
            fixed[net.fidx[nm]] = v
    Aeq, beq = equality_system(net, fixed, extra)
    c = heater_vector(net).copy()
    for i, nm in enumerate(net.fnames):
        if nm.startswith('V:'):
            c[i] += objective_weightV
    bounds = [(None if not np.isfinite(l) else l, None if not np.isfinite(u) else u) for l, u in zip(lb, ub)]
    res = so.linprog(c, A_eq=Aeq, b_eq=beq, bounds=bounds, method='highs')
    if res.status != 0:
        raise RuntimeError('LP infeasible: ' + res.message)
    return res.x


class MinR:
    def __init__(self, net, eps, vmin=1e-5, extra=(), fixed_names=None, products=None, verbose=False):
        """eps: float or dict product name -> max impurity."""
        self.net = net
        net.build()
        lb, ub, fixed = flow_bounds(net, vmin=vmin)
        if fixed_names:
            for nm, v in fixed_names.items():
                fixed[net.fidx[nm]] = v
        self.lb, self.ub = lb, ub
        Aeq, beq = equality_system(net, fixed, extra)
        self.f0 = np.linalg.lstsq(Aeq, beq, rcond=None)[0]
        self.Z = sl.null_space(Aeq)
        self.nu = self.Z.shape[1]
        self.c = heater_vector(net)
        self.prods = [p for p in net.products if (products is None or p['name'] in products)]
        if isinstance(eps, dict):
            self.leps = np.array([np.log(eps[p['name']]) for p in self.prods])
        else:
            self.leps = np.full(len(self.prods), np.log(eps))
        self.x = None
        self.last = None
        self.best = None
        self.newton_maxit = 25
        self.ptc_fallback = 100
        self._cache = {}
        self.verbose = verbose
        self.nfev = 0
        self.ibl = np.where(np.isfinite(lb))[0]
        self.ibu = np.where(np.isfinite(ub))[0]

    def f_of(self, u):
        return self.f0 + self.Z @ u

    def u_of(self, f):
        return self.Z.T @ (f - self.f0)

    def evaluate(self, u):
        key = u.tobytes()
        if key in self._cache:
            return self._cache[key]
        net = self.net
        f = self.f_of(u)
        out = None
        viol = np.max(self.lb[self.ibl] - f[self.ibl]) if len(self.ibl) else -1.0
        if len(self.ibu):
            viol = max(viol, np.max(f[self.ibu] - self.ub[self.ibu]))
        if viol <= 1e-10:
            fc = np.maximum(f, np.where(np.isfinite(self.lb), self.lb, -np.inf))
            try:
                x, info = net.solve(fc, x0=self.x, maxit=self.newton_maxit, ptc_fallback=self.ptc_fallback)
                if not info['ok'] and self.x is None:
                    x, info = net.solve(fc, x0=None)
                self.nfev += 1
                if info['ok']:
                    g = np.zeros(len(self.prods)); dg = np.zeros((len(self.prods), self.nu))
                    Lout, Vout = net.outflows(fc)
                    w = net.weights(x, fc, Lout, Vout)
                    Js = net.scaled_jacobian(x, fc, Lout, Vout, w)     # diag(1/w) J diag(x)
                    lu = spla.splu(Js)
                    Gs = sp.diags(1.0 / w.ravel()) @ net.dres_dflow(x, fc)
                    xr = x.ravel()
                    for i, pr in enumerate(self.prods):
                        li, gx = net.log_impurity_and_grad_x(x, pr)
                        v = np.zeros(Js.shape[0]); sl_ = slice(pr['node'] * net.N, (pr['node'] + 1) * net.N)
                        v[sl_] = gx * xr[sl_]
                        lam = lu.solve(v, trans='T')
                        dgdf = -(Gs.T @ lam)
                        g[i] = self.leps[i] - li
                        dg[i] = -(dgdf @ self.Z)
                    if np.all(np.isfinite(g)) and np.all(np.isfinite(dg)):
                        self.x = x
                        out = dict(f=f, x=x, ok=True, g=g, dg=dg, R=float(self.c @ f))
                        self.last = out
                        if np.min(g) > -1e-7 and (self.best is None or out['R'] < self.best['R']):
                            self.best = out
            except Exception:
                out = None
        if out is None:
            # failure (bound violation or no converged steady state): report a strongly infeasible point
            lg = self.last
            g = (lg['g'] if lg is not None else np.zeros(len(self.prods))) - 20.0
            dg = lg['dg'] if lg is not None else np.zeros((len(self.prods), self.nu))
            out = dict(f=f, x=None, ok=False, g=g, dg=dg, R=float(self.c @ f))
        self._cache = {key: out}
        return out

    def run(self, f_start, maxiter=300, ftol=1e-10):
        self.last = None; self.best = None
        u0 = self.u_of(f_start)
        c = self.c @ self.Z

        def fun(u):
            return float(self.c @ self.f_of(u)), c

        cons = [
            dict(type='ineq', fun=lambda u: self.evaluate(u)['g'], jac=lambda u: self.evaluate(u)['dg']),
            dict(type='ineq', fun=lambda u: self.f_of(u)[self.ibl] - self.lb[self.ibl],
                 jac=lambda u: self.Z[self.ibl]),
        ]
        if len(self.ibu):
            cons.append(dict(type='ineq', fun=lambda u: self.ub[self.ibu] - self.f_of(u)[self.ibu],
                             jac=lambda u: -self.Z[self.ibu]))
        res = so.minimize(fun, u0, jac=True, constraints=cons, method='SLSQP',
                          options=dict(maxiter=maxiter, ftol=ftol, disp=self.verbose))
        for _ in range(4):        # restarts after a line-search breakdown
            if res.status == 0:
                break
            ub_ = self.best['f'] if self.best is not None else self.f_of(res.x)
            u1 = self.u_of(ub_)
            self._cache = {}
            res2 = so.minimize(fun, u1, jac=True, constraints=cons, method='SLSQP',
                               options=dict(maxiter=maxiter, ftol=ftol, disp=self.verbose))
            improved = self.best is not None and self.best['R'] < float(self.c @ ub_) - 1e-12
            res = res2
            if not improved and res2.status != 0:
                break
        ev = self.evaluate(res.x)
        if not ev['ok'] or np.min(ev['g']) < -1e-7:
            if self.best is not None:
                ev = self.best
        if ev['x'] is None:
            return dict(R=np.inf, f=ev['f'], x=None, ok=False, g=ev['g'], min_g=-np.inf, viol_bounds=np.inf,
                        status=res.status, message=res.message, nit=res.nit, nfev=self.nfev)
        f = ev['f']
        viol_bounds = max(0.0, np.max(self.lb[self.ibl] - f[self.ibl])) if len(self.ibl) else 0.0
        return dict(R=ev['R'], f=f, x=ev['x'], ok=ev['ok'], g=ev['g'], min_g=float(np.min(ev['g'])),
                    viol_bounds=viol_bounds, status=res.status, message=res.message, nit=res.nit,
                    nfev=self.nfev)


def report(net, f, x):
    imp = net.impurities(x)
    flows = {nm: float(f[i]) for i, nm in enumerate(net.fnames) if abs(f[i]) > 0}
    return dict(R=float(net.R(f)), impurities=imp, flows=flows)
