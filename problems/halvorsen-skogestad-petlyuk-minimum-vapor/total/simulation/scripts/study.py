"""
study.py -- helpers: minimum R of an arrangement at (n stages, purity 1-eps), stage-level bound at the optimum.
Agent P (Claude Opus 5.5, effort xhigh).
"""
import numpy as np, time
from arrangements import regular_network, init_extra
from opt import lp_init
from ipopt_net import solve_min_R
from refvals import feed_roots


def stage_bound(net, f, x, alpha, z, q):
    """max_j sum_P max(0, phi_P(theta_j)) at the steady state (S/F stage-level lower bound on R)."""
    a = np.asarray(alpha, float)
    th = [float(t) for t in feed_roots(alpha, z, q)]
    y = net.yof(x)
    best = -np.inf; per = []
    for t in th:
        b = 0.0
        for pr in net.products:
            fl = f[net.fidx['P:' + pr['name']]]
            c = x[pr['node']] if pr['phase'] == 'L' else y[pr['node']]
            ph = -np.sum(a * c / (a - t)) if pr['phase'] == 'L' else 1.0 - np.sum(a * c / (a - t))
            b += max(0.0, fl * ph)
        per.append(b); best = max(best, b)
    return best, per


def run_config(alpha, z, q, cols, transfer, n, eps, nT=None, vmin=None, lmin=None, init=None, verbose=False,
               ipopt_options=None):
    N = len(alpha)
    net = regular_network(alpha, z, q, cols, transfer, n=n, nT=nT)
    t = time.time()
    if init is None:
        vm = vmin if vmin is not None else 1.5 + (1 - q if q < 1 else 0)
        f0 = lp_init(net, extra=init_extra(net, N, z, cols), vmin=vm, lmin=lmin if lmin is not None else 1.0)
        x0, info = net.solve(f0)
    else:
        f0, x0 = init
    r = solve_min_R(net, eps, x0, f0, verbose=verbose, options=ipopt_options)
    b, per = stage_bound(net, r['f'], r['x'], alpha, z, q)
    r.update(bound=b, bound_per_j=per, sec=time.time() - t, net=net)
    return r
