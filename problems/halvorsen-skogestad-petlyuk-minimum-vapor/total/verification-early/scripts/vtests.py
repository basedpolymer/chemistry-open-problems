#!/usr/bin/env python3
"""Early verifier's numerical tests of the stage-level bound R >= sum_P max(0, phi_P(theta_j)) (claims 1-2).

Uses only the verifier's own solver vnet.py.  Tests:
  T1  binary columns near minimum reflux (several q, including subcooled and superheated feeds), many stages;
  T2  ternary Petlyuk (fully thermally coupled) columns with HS's preferred-split flows, increasing stage numbers;
  T3  random networks: E stages, liquid/vapor mixer-splitters, total condensers and vaporizers, recycles, bypasses,
      split feeds, vapor and liquid products, heaters/coolers anywhere, N = 2..5, q in [-0.6, 1.6];
  T4  several feeds of different compositions / thermal states, random theta (bound with the feed term kept);
  T5  adversarial local minimisation of the slack R - B_j over the parameters of small random networks;
  T6  binary column with Murphree trays (E unit + vapor bypass + vapor mixer), near minimum reflux.
At every converged steady state and every theta: node balance residual, sign uniformity of every unit, exact
defect identity, bound.  Prints PASS/FAIL per test, exits 1 on any failure.

Environment: Python 3, numpy, scipy (WSL Fedora sci env), 1 thread.  Runtime: about 5-15 minutes.
Usage: python vtests.py [seed] [scale]   (scale multiplies the number of random trials, default 1)
"""
import sys
import time
import numpy as np
from scipy.optimize import minimize
from vnet import Net, feed_roots, peaks, check_bound, phi_stream

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1
SCALE = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
rng = np.random.default_rng(SEED)
FAIL = []
TOL = 1e-8


def chain(net, n, sec='s'):
    st = [net.add('E') for _ in range(n)]
    if not hasattr(net, 'coord'):
        net.coord = {}
    for k, u in enumerate(st):
        net.coord[u] = (sec, (k + 0.5) / n)
    for k in range(n - 1):
        net.link(st[k], 'V', st[k + 1], 1.0)
        net.link(st[k + 1], 'L', st[k], 1.0)
    return st


def single(net, sec, Q=0.0):
    u = net.add('E', Q=Q)
    if not hasattr(net, 'coord'):
        net.coord = {}
    net.coord[u] = (sec, 0.5)
    return u


def interp_S(prev, net):
    """initial S for net from the solution of prev (same builder, other stage numbers), by section position."""
    table = {}
    for u in prev.two:
        sec, t = prev.coord[u]
        table.setdefault(sec, []).append((t, prev.S[prev.idx[u]]))
    S0 = []
    net.total_flows(); net._setup()
    for u in net.two:
        sec, t = net.coord[u]
        pts = sorted(table.get(sec, []))
        if not pts:
            S0.append(float(np.mean(prev.S)))
        else:
            S0.append(float(np.interp(t, [p[0] for p in pts], [p[1] for p in pts])))
    return np.array(S0)


def solve_cont(builder, n_list):
    """Solve builder(n) for increasing n, warm-starting each from the previous solution."""
    prev = None
    for n in n_list:
        net = builder(n)
        S0 = interp_S(prev, net) if prev is not None else None
        ok = net.solve(S0=S0, bp_iter=0)
        yield n, net, ok
        if ok:
            prev = net


def purity(net, labels):
    """labels: dict product label -> component index. Returns (min purity, product flows dict)."""
    sd = net.stream_data()
    tot = {}
    worst = 1.0
    for s in sd:
        if isinstance(s['dst'], tuple):
            lab = s['dst'][1]
            f = s['n'].sum()
            if f > 1e-14:
                worst = min(worst, s['n'][labels[lab]] / f)
            tot[lab] = tot.get(lab, 0) + f
    return worst, tot


# ---------------------------------------------------------------------------------------------- T1
def binary_column(alpha, z, q, VB, nb, nt, murphree=None):
    a = np.array([alpha, 1.0])
    net = Net(a)
    D = z
    VT = VB + (1 - q)
    LT = VT - D
    reb = single(net, 'reb', Q=VB)
    if murphree is None:
        bot = chain(net, nb, 'b'); fs = single(net, 'f'); top = chain(net, nt, 't')
        stages = bot + [fs] + top
        if nb:
            net.link(bot[-1], 'V', fs, 1.0); net.link(fs, 'L', bot[-1], 1.0)
        if nt:
            net.link(fs, 'V', top[0], 1.0); net.link(top[0], 'L', fs, 1.0)
    else:
        # Murphree trays: tray k = E unit receiving fraction E of vapor from below + MV mixer of its vapor with bypass
        E = murphree
        stages = []
        mixers = []
        nst = nb + 1 + nt
        for k in range(nst):
            u = net.add('E')
            net.coord[u] = ('m', (k + 0.5) / nst)
            stages.append(u)
        # vapor from below of tray k comes from mixer of tray k-1 (or reboiler); split by an MV splitter
        prev_v = reb
        for k in range(nst):
            sp = net.add('MV')      # splitter of the incoming vapor
            mx = net.add('MV')      # mixer of equilibrium vapor and bypass
            if prev_v == reb:
                net.link(reb, 'V', sp, 1.0)
            else:
                net.link(prev_v, 'V', sp, 1.0)
            net.link(sp, 'V', stages[k], E)
            net.link(sp, 'V', mx, 1 - E)
            net.link(stages[k], 'V', mx, 1.0)
            prev_v = mx
            mixers.append(mx)
        for k in range(1, nst):
            net.link(stages[k], 'L', stages[k - 1], 1.0)
        fs = stages[nb]
    net.link(reb, 'L', ('P', 'B'), 1.0)
    if murphree is None:
        net.link(reb, 'V', stages[0], 1.0)
        top_v = stages[-1]
    else:
        top_v = mixers[-1]
    net.link(stages[0], 'L', reb, 1.0)
    cd = net.add('C')
    net.link(top_v, 'V', cd, 1.0)
    net.link(cd, 'L', stages[-1], LT / VT)
    net.link(cd, 'L', ('P', 'D'), D / VT)
    net.feed(fs, 1.0, [z, 1 - z], q)
    return net


def T1():
    worst_slack, worst_res, nconv = np.inf, 0.0, 0
    rows = []
    for (alpha, z, q) in [(2.0, 0.5, 1.0), (1.5, 0.3, 0.0), (3.0, 0.7, 0.5), (1.4, 0.5, 1.4), (2.5, 0.4, -0.3)]:
        a = np.array([alpha, 1.0]); zz = np.array([z, 1 - z])
        th = feed_roots(a, zz, q)
        H = peaks(a, zz, th)[0]
        for delta in (0.05, 0.01):
            for n, net, ok in solve_cont(lambda n: binary_column(alpha, z, q, H * (1 + delta), n, n),
                                         (10, 20, 30, 45, 60, 90, 120)):
                if not ok:
                    print('  T1 nonconv', alpha, z, q, delta, n, net.resid)
                    continue
                nconv += 1
                r = check_bound(net, th)
                pur, _ = purity(net, {'B': 1, 'D': 0})
                eps = 1 - pur
                c = max(a[0] / (a[0] - th[0]), a[1] / (th[0] - a[1]))
                sl = r['slack'][0]
                worst_slack = min(worst_slack, sl)
                worst_res = max(worst_res, r['node'], r['cutres'])
                if r['signbad'] or sl < -TOL or r['R'] < H - c * eps - TOL:
                    FAIL.append(('T1', alpha, z, q, delta, n, r))
                rows.append((alpha, z, q, delta, n, r['R'] / H, eps, sl, r['R'] - (H - c * eps)))
    for row in rows:
        if row[4] in (10, 120):
            print('  T1 alpha=%.1f z=%.1f q=%+.1f delta=%.2f n=%3d  R/H1=%.4f  impurity=%.2e  R-B1=%.3e  R-(H1-c eps)=%.3e' % row)
    print('T1 binary columns: converged %d, min slack %.3e, max residual %.1e -> %s'
          % (nconv, worst_slack, worst_res, 'PASS' if not [f for f in FAIL if f[0] == 'T1'] else 'FAIL'))


# ---------------------------------------------------------------------------------------------- T2
def petlyuk3(a, z, q, VB, n, npre=None):
    """Ternary Petlyuk: prefractionator at HS preferred split, main column with liquid side draw."""
    th = feed_roots(a, z, q)
    tA, tB = th
    # preferred split of the prefractionator: top net flows (zA, beta zB, 0)
    num = a[0] * z[0] / (a[0] - tB) - a[0] * z[0] / (a[0] - tA)
    den = a[1] * z[1] / (a[1] - tA) - a[1] * z[1] / (a[1] - tB)
    beta = num / den
    VT1 = a[0] * z[0] / (a[0] - tA) + a[1] * beta * z[1] / (a[1] - tA)
    D1 = z[0] + beta * z[1]
    LT1 = VT1 - D1
    VB1 = VT1 - (1 - q)
    LB1 = LT1 + q
    Lbot = VB + z[2]
    L2 = Lbot - LB1
    V2 = VB - VB1
    L3 = L2 + z[1]
    L4 = L3 + LT1
    V4 = V2 + VT1
    D = z[0]
    assert min(L2, V2, LT1, VB1, L4 - D) > 0, (L2, V2, LT1, VB1)
    npre = n if npre is None else npre
    net = Net(a)
    reb = single(net, 'reb', Q=VB)
    MB = chain(net, n, 'MB'); J2 = single(net, 'J2'); M2 = chain(net, n, 'M2'); S = single(net, 'S')
    M3 = chain(net, n, 'M3'); J1 = single(net, 'J1'); MT = chain(net, n, 'MT'); cd = net.add('C')
    PB = chain(net, npre, 'PB'); PF = single(net, 'PF'); PT = chain(net, npre, 'PT')
    net.link(reb, 'L', ('P', 'C'), 1.0); net.link(reb, 'V', MB[0], 1.0); net.link(MB[0], 'L', reb, 1.0)
    net.link(MB[-1], 'V', J2, 1.0); net.link(J2, 'L', MB[-1], 1.0)
    net.link(J2, 'V', PB[0], VB1 / VB); net.link(J2, 'V', M2[0], V2 / VB)
    net.link(M2[0], 'L', J2, 1.0); net.link(M2[-1], 'V', S, 1.0)
    net.link(S, 'L', M2[-1], L2 / L3); net.link(S, 'L', ('P', 'B'), z[1] / L3)
    net.link(S, 'V', M3[0], 1.0); net.link(M3[0], 'L', S, 1.0); net.link(M3[-1], 'V', J1, 1.0)
    net.link(J1, 'L', M3[-1], L3 / L4); net.link(J1, 'L', PT[-1], LT1 / L4)
    net.link(J1, 'V', MT[0], 1.0); net.link(MT[0], 'L', J1, 1.0); net.link(MT[-1], 'V', cd, 1.0)
    net.link(cd, 'L', MT[-1], (V4 - D) / V4); net.link(cd, 'L', ('P', 'A'), D / V4)   # reflux L4 = V4 - D
    net.link(PB[0], 'L', J2, 1.0); net.link(PB[-1], 'V', PF, 1.0); net.link(PF, 'L', PB[-1], 1.0)
    net.link(PF, 'V', PT[0], 1.0); net.link(PT[0], 'L', PF, 1.0); net.link(PT[-1], 'V', J1, 1.0)
    net.feed(PF, 1.0, z, q)
    return net, th


def T2():
    bad = 0
    for (a, z, q) in [((4.0, 2.0, 1.0), (1 / 3, 1 / 3, 1 / 3), 1.0), ((6.0, 2.5, 1.0), (0.2, 0.5, 0.3), 1.0),
                      ((4.0, 2.0, 1.0), (0.3, 0.3, 0.4), 0.7)]:
        a = np.array(a); z = np.array(z)
        th = feed_roots(a, z, q)
        H = peaks(a, z, th)
        Vp = H.max()
        for delta in (0.05, 0.01):
            for n, net, ok in solve_cont(lambda n: petlyuk3(a, z, q, Vp * (1 + delta), n)[0],
                                         (8, 12, 16, 24, 32, 40, 48)):
                if not ok:
                    print('  T2 nonconv', a, z, q, delta, n, net.resid)
                    continue
                r = check_bound(net, th)
                pur, _ = purity(net, {'A': 0, 'B': 1, 'C': 2})
                eps = 1 - pur
                c = [max(a[j] / (a[j] - th[j]), a[j + 1] / (th[j] - a[j + 1])) for j in range(2)]
                lb = max(H[j] - c[j] * eps for j in range(2))
                okk = (r['signbad'] == 0 and min(r['slack']) > -TOL and r['R'] >= lb - TOL
                       and r['node'] < 1e-9 and r['cutres'] < 1e-9)
                bad += (not okk)
                if not okk:
                    FAIL.append(('T2', a, z, q, delta, n, r))
                if n in (8, 48):
                    print('  T2 a=%s z=%s q=%.1f delta=%.2f n=%2d R/Vpet=%.4f impurity=%.2e  R-B_j=%s  maxres=%.1e'
                          % (a, np.round(z, 2), q, delta, n, r['R'] / Vp, eps,
                             np.array2string(np.array(r['slack']), precision=3), max(r['node'], r['cutres'])))
    print('T2 ternary Petlyuk: %s' % ('PASS' if bad == 0 else 'FAIL'))


# ---------------------------------------------------------------------------------------------- T3/T4
def random_net(N, K, multi_feed=False):
    a = np.sort(rng.uniform(1.0, 6.0, N))[::-1]
    a[-1] = 1.0
    if len(set(np.round(a, 6))) < N:
        a = np.linspace(4, 1, N)
    z = rng.dirichlet(np.ones(N))
    q = rng.uniform(-0.6, 1.6)
    types = list(rng.choice(['E', 'E', 'E', 'E', 'ML', 'MV', 'C', 'R'], size=K))
    types[0] = 'E'
    net = Net(a)
    for t in types:
        net.add(t)
    accepts = {'E': 'LV', 'ML': 'L', 'MV': 'V', 'C': 'V', 'R': 'L'}
    acc = {ph: [u for u, t in enumerate(types) if ph in accepts[t]] for ph in 'LV'}
    plab = 0
    for u, t in enumerate(types):
        for ph in {'E': 'LV', 'ML': 'L', 'MV': 'V', 'C': 'L', 'R': 'V'}[t]:
            ndst = rng.integers(1, 4)
            dsts = []
            for _ in range(ndst):
                if rng.random() < 0.3 or not acc[ph]:
                    dsts.append(('P', plab)); plab += 1
                else:
                    dsts.append(int(rng.choice(acc[ph])))
            fr = rng.dirichlet(np.ones(len(dsts)))
            for d, f in zip(dsts, fr):
                if d == u:            # self loop: allowed for E (recycle of own outflow)? keep it simple: product
                    d = ('P', plab); plab += 1
                net.link(u, ph, d, f)
    Eunits = [u for u, t in enumerate(types) if t == 'E']
    nf = rng.integers(1, min(3, len(Eunits)) + 1)
    fu = rng.choice(Eunits, size=nf, replace=False)
    ff = rng.dirichlet(np.ones(nf))
    feeds = []
    for k, (u, f) in enumerate(zip(fu, ff)):
        if multi_feed and k > 0:
            zz = rng.dirichlet(np.ones(N)); qq = rng.uniform(-0.6, 1.6)
        else:
            zz, qq = z, q
        net.feed(int(u), f, zz, qq)
        feeds.append((zz, qq))
    for u in Eunits:
        net.Q[u] = rng.normal(0, 0.5) if rng.random() < 0.6 else 0.0
    return net, a, z, q


def connected_ok(net):
    K = len(net.types)
    succ = {u: set() for u in range(K)}
    outp = set()
    for (s, ph, d, fr) in net.streams:
        if isinstance(d, tuple):
            outp.add(s)
        else:
            succ[s].add(d)
    # every unit reachable from a feed and reaching a product
    seen = set(u for (u, _, _, _) in net.feeds)
    stack = list(seen)
    while stack:
        u = stack.pop()
        for v in succ[u]:
            if v not in seen:
                seen.add(v); stack.append(v)
    if len(seen) < K:
        return False
    pred = {u: set() for u in range(K)}
    for u in range(K):
        for v in succ[u]:
            pred[v].add(u)
    good = set(outp); stack = list(good)
    while stack:
        u = stack.pop()
        for v in pred[u]:
            if v not in good:
                good.add(v); stack.append(v)
    return len(good) == K


def T3(ntr, multi=False):
    tag = 'T4' if multi else 'T3'
    stats = dict(tried=0, solved=0, signbad=0, minslack=np.inf, maxres=0.0, recycles=0, reverse=0)
    t0 = time.time()
    while stats['tried'] < ntr:
        N = int(rng.integers(2, 6)); K = int(rng.integers(2, 13))
        net, a, z, q = random_net(N, K, multi_feed=multi)
        if not connected_ok(net):
            continue
        # scale duties down until total flows are nonnegative
        okflows = False
        for _ in range(12):
            try:
                Lo, Vo = net.total_flows()
            except np.linalg.LinAlgError:
                break
            if np.all(Lo > 1e-6 * (Lo.max() + 1) - 1e-12) and np.all(Vo >= -1e-12) and max(Lo.max(), Vo.max()) < 1e4:
                # E units need positive flow of each phase that has streams (all E have both)
                Eok = all(Lo[u] > 1e-7 and Vo[u] > 1e-7 for u, t in enumerate(net.types) if t == 'E')
                if Eok:
                    okflows = True
                    break
            net.Q = [0.6 * x for x in net.Q]
            # shift: give E units a positive duty to create vapor if vapor flows are negative
            for u, t in enumerate(net.types):
                if t == 'E' and Vo[u] <= 1e-7:
                    net.Q[u] += abs(Vo[u]) + rng.uniform(0.05, 0.5)
                if t == 'E' and Lo[u] <= 1e-7:
                    net.Q[u] -= abs(Lo[u]) + rng.uniform(0.01, 0.1)
        if not okflows:
            continue
        stats['tried'] += 1
        if not net.solve():
            continue
        stats['solved'] += 1
        # recycle detection
        K = len(net.types)
        adj = np.zeros((K, K), bool)
        for (s, ph, d, fr) in net.streams:
            if not isinstance(d, tuple):
                adj[s, d] = True
        reach = adj.copy()
        for k in range(K):
            reach |= reach[:, [k]] & reach[[k], :]
        if np.any(np.diag(reach)):
            stats['recycles'] += 1
        if multi:
            thetas = list(rng.uniform(0.2, a[0] * 1.5, 6))
            thetas = [t for t in thetas if np.min(np.abs(a - t)) > 1e-3]
        else:
            thetas = list(feed_roots(a, z, q))
        r = check_bound(net, thetas)
        stats['signbad'] += r['signbad']
        sl = min(r['slack'])
        stats['minslack'] = min(stats['minslack'], sl / max(1.0, r['R']))
        stats['maxres'] = max(stats['maxres'], r['node'], r['cutres'])
        if r['signbad'] or sl < -TOL * max(1.0, r['R']) or r['node'] > 1e-8 or r['cutres'] > 1e-8:
            FAIL.append((tag, r))
    ok = not [f for f in FAIL if f[0] == tag]
    print('%s %s: %s  (%.0f s)' % (tag, 'several feeds, random theta' if multi else 'random networks, one feed, theta_j',
                                   stats, time.time() - t0) + ' -> ' + ('PASS' if ok else 'FAIL'))


# ---------------------------------------------------------------------------------------------- T5
def T5(nstarts):
    """Adversarial: fix a random topology, minimise min_j (R - B_j)/(1+R) over fractions and duties."""
    best = np.inf
    nopt = 0
    t0 = time.time()
    for st in range(nstarts):
        N = int(rng.integers(2, 5)); K = int(rng.integers(3, 8))
        net0, a, z, q = random_net(N, K)
        if not connected_ok(net0):
            continue
        th = feed_roots(a, z, q)
        # parameters: logits of fractions per (unit, phase) group, duties of E units
        groups = {}
        for k, (s, ph, d, fr) in enumerate(net0.streams):
            groups.setdefault((s, ph), []).append(k)
        Eu = [u for u, t in enumerate(net0.types) if t == 'E']
        nfr = len(net0.streams)

        def build(p):
            net = Net(a)
            net.types = list(net0.types); net.Q = [0.0] * len(net0.types)
            for u_k, u in enumerate(Eu):
                net.Q[u] = p[nfr + u_k]
            fr = np.zeros(nfr)
            for g, ks in groups.items():
                w = np.exp(np.clip(p[ks], -20, 20)); fr[ks] = w / w.sum()
            net.streams = [(s, ph, d, fr[k]) for k, (s, ph, d, _) in enumerate(net0.streams)]
            net.feeds = list(net0.feeds)
            return net

        def obj(p):
            net = build(p)
            try:
                Lo, Vo = net.total_flows()
            except np.linalg.LinAlgError:
                return 10.0
            if any((Lo[u] <= 1e-7 or Vo[u] <= 1e-7) for u, t in enumerate(net.types) if t == 'E') or \
               np.any(Lo < -1e-12) or np.any(Vo < -1e-12) or max(Lo.max(), Vo.max()) > 1e3:
                return 10.0
            if not net.solve(tol=1e-11, maxit=60):
                return 10.0
            r = check_bound(net, th)
            if r['signbad']:
                FAIL.append(('T5-sign', r))
            return min(r['slack']) / (1.0 + r['R'])

        p0 = np.concatenate([rng.normal(0, 1, nfr), [net0.Q[u] for u in Eu]])
        if obj(p0) >= 10:
            continue
        nopt += 1
        res = minimize(obj, p0, method='Nelder-Mead', options=dict(maxfev=600, xatol=1e-6, fatol=1e-10))
        best = min(best, res.fun)
        if res.fun < -1e-8:
            FAIL.append(('T5', res.fun, res.x))
    ok = not [f for f in FAIL if f[0].startswith('T5')]
    print('T5 adversarial minimisation: %d optimised topologies, smallest relative slack found %.3e (%.0f s) -> %s'
          % (nopt, best, time.time() - t0, 'PASS' if ok else 'FAIL'))


def T6():
    bad = 0
    alpha, z, q = 2.0, 0.5, 1.0
    a = np.array([alpha, 1.0]); zz = np.array([z, 1 - z])
    th = feed_roots(a, zz, q); H = peaks(a, zz, th)[0]
    for E in (1.0, 0.7, 0.4):
        for n, net, ok in solve_cont(lambda n: binary_column(alpha, z, q, H * 1.02, n, n, murphree=E),
                                     (10, 20, 40, 60, 90, 120)):
            if not ok:
                print('  T6 nonconv', E, n); continue
            if n not in (20, 60, 120):
                continue
            r = check_bound(net, th)
            pur, _ = purity(net, {'B': 1, 'D': 0})
            okk = r['signbad'] == 0 and r['slack'][0] > -TOL
            bad += not okk
            if not okk:
                FAIL.append(('T6', E, n, r))
            print('  T6 Murphree E=%.1f n=%3d  R/H1=%.4f impurity=%.2e  R-B1=%.3e' % (E, n, r['R'] / H, 1 - pur, r['slack'][0]))
    print('T6 Murphree trays: %s' % ('PASS' if bad == 0 else 'FAIL'))


if __name__ == '__main__':
    t0 = time.time()
    T1()
    T2()
    T6()
    T3(int(400 * SCALE))
    T3(int(150 * SCALE), multi=True)
    T5(int(25 * SCALE))
    print('total time %.0f s' % (time.time() - t0))
    if FAIL:
        print('FAILURES:', len(FAIL))
        for f in FAIL[:5]:
            print(f)
        print('FAIL')
        sys.exit(1)
    print('ALL PASS')
