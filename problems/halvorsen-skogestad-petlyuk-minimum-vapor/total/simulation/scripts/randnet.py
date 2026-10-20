"""
randnet.py -- random networks of equilibrium stages and random consistent flows (for lemma tests).
Agent P (Claude Opus 5.5, effort xhigh).

random_network(rng, N): 3-7 equilibrium-stage nodes (adiabatic or with a free-sign duty), 2-6 CMO sections of
0-11 stages between random node pairs, 0-4 liquid/vapor streams between random nodes (bypasses, recycles), the
process feed (possibly split in two portions), 2-4 products of either phase at random nodes (+ one product at every
node not touched by a section).  random_flows(rng, net): a random point of the flow polytope (convex combination
of 6 random LP vertices), all section flows >= a random positive lower bound.
"""
import numpy as np
import scipy.optimize as so
from stagenet import Network


def phi_liq(x, th, a):
    """stream flux per unit flow of a liquid of composition x: -sum alpha x/(alpha - theta)."""
    return -np.sum(a * x / (a - th), axis=-1)


def phi_vap(y, th, a):
    """stream flux per unit flow of a vapor of composition y: 1 - sum alpha y/(alpha - theta)."""
    return 1.0 - np.sum(a * y / (a - th), axis=-1)


def random_network(rng, N):
    a = np.sort(np.exp(rng.uniform(0, np.log(12), N)))[::-1]
    a = a / a[-1]
    z = rng.dirichlet(np.ones(N)); z = np.maximum(z, 0.02); z /= z.sum()
    q = rng.choice([rng.uniform(-0.3, 1.3), 0.0, 1.0, 0.5])
    net = Network(a)
    K = int(rng.integers(3, 8))
    for k in range(K):
        net.node('n%d' % k, rng.choice(['adiabatic', 'adiabatic', 'free', 'free', 'free']))
    used = set()
    for s in range(int(rng.integers(2, 7))):
        t, b = rng.choice(K, 2, replace=False)
        used |= {int(t), int(b)}
        net.section('s%d' % s, int(rng.integers(0, 12)), int(t), int(b))
    for e in range(int(rng.integers(0, 5))):
        u, v = rng.choice(K, 2, replace=False)
        net.stream('e%d' % e, int(u), rng.choice(['L', 'V']), int(v))
    if rng.random() < 0.3:
        k1, k2 = rng.choice(K, 2, replace=False); w = rng.uniform(0.2, 0.8)
        net.feed(int(k1), w, q, z); net.feed(int(k2), 1 - w, q, z)
    else:
        net.feed(int(rng.integers(K)), 1.0, q, z)
    for p in range(int(rng.integers(2, 5))):
        net.product('p%d' % p, int(rng.integers(K)), rng.choice(['L', 'V']), int(rng.integers(N)))
    for k in range(K):
        if k not in used:
            net.product('px%d' % k, k, rng.choice(['L', 'V']), int(rng.integers(N)))
    net.build()
    return net, a, z, q


def random_flows(rng, net, nvert=6):
    A, b = net.balance_matrix()
    nf = net.nf
    lb = np.zeros(nf); ub = np.full(nf, 5.0)
    for k, nm in enumerate(net.nodes):
        i = net.fidx['Q:' + nm]
        kind = net.node_kind[k]
        if kind == 'adiabatic':
            lb[i] = ub[i] = 0.0
        elif kind == 'heater':
            lb[i], ub[i] = 0.0, 5.0
        elif kind == 'cooler':
            lb[i], ub[i] = -5.0, 0.0
        else:
            lb[i], ub[i] = -5.0, 5.0
    for i, nm in enumerate(net.fnames):
        if nm[0] in 'VL':
            lb[i] = rng.uniform(0.005, 0.1)
    pts = []
    for _ in range(nvert):
        res = so.linprog(rng.normal(size=nf), A_eq=A, b_eq=b, bounds=list(zip(lb, ub)), method='highs')
        if res.status != 0:
            return None
        pts.append(res.x)
    w = rng.dirichlet(np.ones(len(pts)))
    return np.sum([wi * p for wi, p in zip(w, pts)], axis=0)
