"""
ipopt_net.py -- full-space minimum-vaporization NLP (pyomo + ipopt) for a stagenet.Network.

Agent P (Claude Opus 5.5, effort xhigh).  Same MESH/CMO model as stagenet.py (no Underwood quantity):
variables x, y (every stage), flows f (sections, streams, products, duties); constraints: component balances,
equilibrium y_p * sum_k alpha_k x_k = alpha_p x_p, node balances (linear), purity of every product.
Every optimum is afterwards re-solved by stagenet.Network.solve at the optimal flows (independent Newton solver,
relative residual 1e-11 per component and stage) and the impurities reported are those of that re-solve.
"""
import numpy as np
import pyomo.environ as pe
from opt import flow_bounds


def build_model(net, eps, x0, f0, extra=(), fixed_names=None, purity='mole'):
    net.build()
    N, a = net.N, net.alpha
    m = pe.ConcreteModel()
    S = range(net.nst); P = range(N); F = range(net.nf)
    lb, ub, fixed = flow_bounds(net, vmin=1e-7)
    if fixed_names:
        for nm, v in fixed_names.items():
            fixed[net.fidx[nm]] = v
    y0 = net.yof(x0)
    m.x = pe.Var(S, P, bounds=(0.0, 1.0), initialize=lambda m, s, p: float(min(1.0, max(x0[s, p], 1e-20))))
    m.y = pe.Var(S, P, bounds=(0.0, 1.0), initialize=lambda m, s, p: float(min(1.0, max(y0[s, p], 1e-20))))

    def fb(m, i):
        return (None if not np.isfinite(lb[i]) else float(lb[i]), None if not np.isfinite(ub[i]) else float(ub[i]))
    m.f = pe.Var(F, bounds=fb, initialize=lambda m, i: float(f0[i]))
    for i, v in fixed.items():
        m.f[i].fix(float(v))
    # outflow expressions
    Lo = {s: [] for s in S}; Vo = {s: [] for s in S}
    for s, i in net.Lo:
        Lo[int(s)].append(int(i))
    for s, i in net.Vo:
        Vo[int(s)].append(int(i))
    conns = {s: [] for s in S}
    for d, sr, ph, fi in zip(net.c_dst, net.c_src, net.c_ph, net.c_fi):
        conns[int(d)].append((int(sr), int(ph), int(fi)))

    def bal(m, s, p):
        e = float(net.feedvec[s, p])
        for sr, ph, fi in conns[s]:
            e += m.f[fi] * (m.x[sr, p] if ph == 0 else m.y[sr, p])
        e -= sum(m.f[i] for i in Lo[s]) * m.x[s, p]
        e -= sum(m.f[i] for i in Vo[s]) * m.y[s, p]
        return e == 0
    m.bal = pe.Constraint(S, P, rule=bal)

    def eq(m, s, p):
        return m.y[s, p] * sum(float(a[k]) * m.x[s, k] for k in P) == float(a[p]) * m.x[s, p]
    m.eq = pe.Constraint(S, P, rule=eq)
    A, b = net.balance_matrix()
    rows = list(range(A.shape[0]))

    def nb(m, r):
        nz = np.nonzero(A[r])[0]
        if len(nz) == 0:
            return pe.Constraint.Skip
        return sum(float(A[r, i]) * m.f[int(i)] for i in nz) == float(b[r])
    m.nb = pe.Constraint(rows, rule=nb)
    m.extra = pe.ConstraintList()
    for coef, v in extra:
        m.extra.add(sum(c * m.f[net.fidx[nm]] for nm, c in coef.items()) == v)
    m.pur = pe.ConstraintList()
    for pr in net.products:
        k, p = pr['node'], pr['comp']
        c = m.x if pr['phase'] == 'L' else m.y
        e = eps[pr['name']] if isinstance(eps, dict) else eps
        m.pur.add(sum(c[k, q] for q in P if q != p) <= e * sum(c[k, q] for q in P))
    heat = [net.fidx['Q:' + nm] for k, nm in enumerate(net.nodes) if net.node_kind[k] == 'heater']
    m.obj = pe.Objective(expr=sum(m.f[i] for i in heat))
    return m


def solve_min_R(net, eps, x0, f0, extra=(), fixed_names=None, tol=1e-10, max_iter=3000, verbose=False,
                options=None):
    m = build_model(net, eps, x0, f0, extra, fixed_names)
    opt = pe.SolverFactory('ipopt')
    opt.options['tol'] = tol
    opt.options['constr_viol_tol'] = 1e-12
    opt.options['max_iter'] = max_iter
    opt.options['mu_strategy'] = 'adaptive'
    opt.options['bound_push'] = 1e-12
    opt.options['bound_frac'] = 1e-12
    opt.options['mu_init'] = 1e-4
    opt.options['bound_relax_factor'] = 0.0
    opt.options['linear_solver'] = 'mumps'
    if options:
        for k, v in options.items():
            opt.options[k] = v
    res = opt.solve(m, tee=verbose)
    f = np.array([pe.value(m.f[i]) for i in range(net.nf)])
    x = np.array([[pe.value(m.x[s, p]) for p in range(net.N)] for s in range(net.nst)])
    status = str(res.solver.termination_condition)
    # independent re-solve of the steady state at the optimal flows
    xs, info = net.solve(f, x0=np.maximum(x, 1e-300))
    if not info['ok']:
        xs, info = net.solve(f)
    return dict(status=status, R=float(net.R(f)), f=f, x_ipopt=x, x=xs, resolve_ok=info['ok'],
                impurities=net.impurities(xs), max_imp=max(net.impurities(xs).values()))
