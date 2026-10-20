"""
test_stage_bound.py -- numerical test of the stage-level bound of Agents S and F on MESH steady states.

Agent P (Claude Opus 5.5, effort xhigh).  For random networks (random equilibrium-stage nodes, sections between
random node pairs, liquid/vapor bypass and recycle streams, split feed, heaters/coolers, products of either phase)
and random flows satisfying the CMO balances, the steady state is computed with stagenet.py; then, at every feed
root theta_j (and at random theta > 0), we check
  (B) unit balance  sum_out phi - sum_in phi = Q_u   at every stage (nodes and section stages),
  (S) sign uniformity: the liquid and vapor outflows of every stage have phi of the same sign,
  (C) the bound  R >= sum_{product streams} max(0, phi_P(theta_j)),
with phi_e(theta) = v_e - sum_i alpha_i n_{e,i}/(alpha_i - theta) (stream flux).
Also records the reverse-flow content of the sections (components with net flow against the majority direction)
and the minimum slack R - bound.

Run: python test_stage_bound.py [ntrials] [seed]   (default 300 trials, ~2-5 min, 1 thread).
Prints PASS/FAIL per claim; exit code 1 on failure.
"""
import sys, json, time
import numpy as np
import scipy.optimize as so
from stagenet import Network
from refvals import feed_roots

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 12345)
ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 300


def phi_liq(x, th, a):
    return -np.sum(a * x / (a - th), axis=-1)


def phi_vap(y, th, a):
    return 1.0 - np.sum(a * y / (a - th), axis=-1)


def random_network(N):
    a = np.sort(np.exp(rng.uniform(0, np.log(12), N)))[::-1]
    a = a / a[-1]
    z = rng.dirichlet(np.ones(N)); z = np.maximum(z, 0.02); z /= z.sum()
    q = rng.choice([rng.uniform(-0.3, 1.3), 0.0, 1.0, 0.5])
    net = Network(a)
    K = rng.integers(3, 8)
    kinds = []
    for k in range(K):
        kind = rng.choice(['adiabatic', 'adiabatic', 'free', 'free', 'free'])
        net.node('n%d' % k, kind); kinds.append(kind)
    ns = rng.integers(2, 7)
    used = set()
    for s in range(ns):
        t, b = rng.choice(K, 2, replace=False)
        used |= {int(t), int(b)}
        net.section('s%d' % s, int(rng.integers(0, 12)), int(t), int(b))
    nstr = rng.integers(0, 5)
    for e in range(nstr):
        u, v = rng.choice(K, 2, replace=False)
        net.stream('e%d' % e, int(u), rng.choice(['L', 'V']), int(v))
    # feed, possibly split in two portions
    if rng.random() < 0.3:
        k1, k2 = rng.choice(K, 2, replace=False); w = rng.uniform(0.2, 0.8)
        net.feed(int(k1), w, q, z); net.feed(int(k2), 1 - w, q, z)
    else:
        net.feed(int(rng.integers(K)), 1.0, q, z)
    npr = rng.integers(2, 5)
    for p in range(npr):
        net.product('p%d' % p, int(rng.integers(K)), rng.choice(['L', 'V']), int(rng.integers(N)))
    for k in range(K):
        if k not in used:
            net.product('px%d' % k, k, rng.choice(['L', 'V']), int(rng.integers(N)))
    net.build()
    return net, a, z, q


def random_flows(net):
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
        if nm[0] == 'P':
            lb[i] = 0.0
    pts = []
    for _ in range(6):
        c = rng.normal(size=nf)
        res = so.linprog(c, A_eq=A, b_eq=b, bounds=list(zip(lb, ub)), method='highs')
        if res.status != 0:
            return None
        pts.append(res.x)
    w = rng.dirichlet(np.ones(len(pts)))
    return np.sum([wi * p for wi, p in zip(w, pts)], axis=0)


def check(net, a, z, q, f, x):
    y = net.yof(x)
    Lout, Vout = net.outflows(f)
    thetas = list(map(float, feed_roots(a, z, q)))
    worst_B = 0.0; sign_bad = 0; slack = []
    for th in thetas + list(rng.uniform(0.05, a[0] * 1.5, 2)):
        if np.min(np.abs(a - th)) < 1e-6:
            continue
        pl = phi_liq(x, th, a); pv = phi_vap(y, th, a)
        phiF = (1 - q) - np.sum(a * z / (a - th))
        # unit balances
        out = Lout * pl + Vout * pv
        inn = np.zeros(net.nst)
        cphi = np.where(net.c_ph == 0, pl[net.c_src], pv[net.c_src]) * f[net.c_fi]
        np.add.at(inn, net.c_dst, cphi)
        for fd in net.feeds:
            inn[fd['node']] += fd['F'] * phiF
        Q = np.zeros(net.nst)
        for k, nm in enumerate(net.nodes):
            Q[k] = f[net.fidx['Q:' + nm]]
        # scale: sum of absolute values of every term (phi itself is a sum that can cancel near poles)
        al = np.sum(np.abs(a * x / (a - th)), axis=-1); av = 1 + np.sum(np.abs(a * y / (a - th)), axis=-1)
        scale = 1 + Lout * al + Vout * av + np.abs(Q)
        np.add.at(scale, net.c_dst, np.where(net.c_ph == 0, al[net.c_src], av[net.c_src]) * f[net.c_fi])
        for fd in net.feeds:
            scale[fd['node']] += fd['F'] * (abs(1 - q) + np.sum(np.abs(a * z / (a - th))))
        worst_B = max(worst_B, np.max(np.abs(out - inn - Q) / scale))
        # sign uniformity (only meaningful where both phases leave)
        both = (Lout > 1e-12) & (Vout > 1e-12)
        sign_bad += int(np.sum(both & (pl * pv < 0) & (np.abs(pl) > 1e-12) & (np.abs(pv) > 1e-12)))
        if th in thetas:
            bound = 0.0
            for pr in net.products:
                fl = f[net.fidx['P:' + pr['name']]]
                ph = pl[pr['node']] if pr['phase'] == 'L' else pv[pr['node']]
                bound += max(0.0, fl * ph)
            slack.append(net.R(f) - bound)
    return worst_B, sign_bad, min(slack)


t0 = time.time()
stats = dict(trials=0, solved=0, worst_B=0.0, sign_bad=0, min_slack=np.inf, reverse_sections=0, recycles=0)
fails = []
for tr in range(ntr):
    N = int(rng.integers(2, 6))
    net, a, z, q = random_network(N)
    f = random_flows(net)
    if f is None:
        continue
    stats['trials'] += 1
    x, info = net.solve(f)
    if not info['ok'] or info['res'] > 1e-11:
        continue
    stats['solved'] += 1
    wB, sb, sl = check(net, a, z, q, f, x)
    stats['worst_B'] = max(stats['worst_B'], wB); stats['sign_bad'] += sb
    stats['min_slack'] = min(stats['min_slack'], sl)
    # reverse flows: a section whose net flows are not "lights up, heavies down" (sign pattern not monotone)
    for nm, w in net.section_netflows(x, f).items():
        sg = np.sign(np.where(np.abs(w) > 1e-9, w, 0))
        nz = sg[sg != 0]
        if len(nz) > 1 and np.any(np.diff(nz) > 0):
            stats['reverse_sections'] += 1
    if sl < -1e-9:
        fails.append(dict(trial=tr, slack=sl))
stats['seconds'] = time.time() - t0
print(json.dumps(stats))
ok = stats['worst_B'] < 1e-9 and stats['sign_bad'] == 0 and stats['min_slack'] > -1e-9 and stats['solved'] > 0
print('(B) unit balances      :', 'PASS' if stats['worst_B'] < 1e-9 else 'FAIL', stats['worst_B'])
print('(S) sign uniformity    :', 'PASS' if stats['sign_bad'] == 0 else 'FAIL', stats['sign_bad'])
print('(C) R >= sum max(0,phi):', 'PASS' if stats['min_slack'] > -1e-9 else 'FAIL', stats['min_slack'])
json.dump(dict(stats=stats, fails=fails), open('../runs/test_stage_bound_%d.json' % (int(sys.argv[2]) if len(sys.argv) > 2 else 12345), 'w'), indent=1)
sys.exit(0 if ok else 1)
