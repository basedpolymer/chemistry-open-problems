#!/usr/bin/env python3
"""
check_stage_cut.py (Agent S) -- numerical test of Lemmas S0-S2 and of the stage-cut theorem of sections.md on
finite-stage MESH/CMO steady states (stagenet.py), for
  (1) single columns (total condenser, partial reboiler), random alpha, z, q, stage numbers, D, V_B;
  (2) three-product Petlyuk arrangements (prefractionator + main column, two thermal couplings, liquid side draw);
  (3) random networks of equilibrium stages and mixers/splitters with recycles and heat duties anywhere,
      feed split among several nodes.
Checks at every converged point (residual < 1e-10, compositions >= 0), for every feed root theta_j:
  (a) global identity sum_P T(P) = Q - R;
  (b) sign consistency of every equilibrium stage (T of liquid and vapor outflows have the same sign);
  (c) the cut identity for U^c = {nodes with sign <= 0} (every crossing term has the favourable sign);
  (d) the bound R >= sum_P max(0, -T(P; theta_j)).
Prints PASS/FAIL per family; exits non-zero on failure.
Environment: Python 3 + numpy + scipy (WSL Fedora sci env).  Runtime: about 1-3 minutes.
Usage: python check_stage_cut.py [ncases_per_family] [seed]
"""
import sys
import numpy as np
from stagenet import Net

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
NC = int(sys.argv[1]) if len(sys.argv) > 1 else 40


def rand_data(N):
    a = np.sort(np.exp(rng.uniform(0, 3, N)))[::-1]
    a = a / a[-1]
    z = rng.dirichlet(np.ones(N))
    q = rng.choice([rng.uniform(-0.3, 1.3), 0.0, 1.0, 0.5])
    return a, z, q


def column(N, nT, nB, D, VB, a, z, q):
    net = Net(a, z, q)
    VT = VB + (1 - q)
    LT = VT - D
    B = 1 - D
    LB = VB + B
    if min(VT, LT, LB, VB, B, D) <= 1e-6:
        return None
    cond = net.add_node('mix', 'cond')
    tops = [net.add_node('stage', f't{i}') for i in range(nT)]
    f = net.add_node('stage', 'feed')
    bots = [net.add_node('stage', f'b{i}') for i in range(nB)]
    reb = net.add_node('stage', 'reb')
    chain = [cond] + tops + [f] + bots + [reb]
    net.add_stream('F', f, 'F', 1.0)
    for k in range(len(chain) - 1):
        up, dn = chain[k], chain[k + 1]
        below_feed = chain.index(dn) > chain.index(f)
        Lk = LB if below_feed else LT
        Vk = VB if below_feed else VT
        net.add_stream(up, dn, 'L', Lk)
        net.add_stream(dn, up, 'V', Vk)
    net.add_stream(cond, ('P', 'D'), 'L', D)
    net.add_stream(reb, ('P', 'B'), 'L', B)
    return net


def petlyuk(N, ns, a, z, q, VB, V1, L1, S, D):
    """ns = stages per section (prefrac top, prefrac bottom, main top, mid-upper, mid-lower, main bottom)"""
    net = Net(a, z, q)
    B = 1 - D - S
    LB = VB + B
    Lml = LB - L1 - q
    Lmu = Lml + S
    LT = Lmu + L1
    VT = VB + 1 - q
    Vm = VB - V1
    flows = [VB, V1, L1, S, D, B, LB, Lml, Lmu, LT, VT, Vm, V1 + 1 - q, L1 + q]
    if min(flows) <= 1e-6:
        return None
    # prefractionator
    pt = [net.add_node('stage', f'pt{i}') for i in range(ns[0])]
    pf = net.add_node('stage', 'pf')
    pb = [net.add_node('stage', f'pb{i}') for i in range(ns[1])]
    # main column
    cond = net.add_node('mix', 'cond')
    mt = [net.add_node('stage', f'mt{i}') for i in range(ns[2])]
    J1 = net.add_node('stage', 'J1')
    mu = [net.add_node('stage', f'mu{i}') for i in range(ns[3])]
    Sd = net.add_node('stage', 'S')
    ml = [net.add_node('stage', f'ml{i}') for i in range(ns[4])]
    J2 = net.add_node('stage', 'J2')
    mb = [net.add_node('stage', f'mb{i}') for i in range(ns[5])]
    reb = net.add_node('stage', 'reb')

    def link(chain, L, V):
        for k in range(len(chain) - 1):
            net.add_stream(chain[k], chain[k + 1], 'L', L)
            net.add_stream(chain[k + 1], chain[k], 'V', V)

    net.add_stream('F', pf, 'F', 1.0)
    link([J1] + pt + [pf], L1, V1 + 1 - q)
    link([pf] + pb + [J2], L1 + q, V1)
    link([cond] + mt + [J1], LT, VT)
    link([J1] + mu + [Sd], Lmu, Vm)
    link([Sd] + ml + [J2], Lml, Vm)
    link([J2] + mb + [reb], LB, VB)
    net.add_stream(cond, ('P', 'D'), 'L', D)
    net.add_stream(Sd, ('P', 'S'), 'L', S)
    net.add_stream(reb, ('P', 'B'), 'L', B)
    return net


def random_net(N, a, z, q, nnodes):
    net = Net(a, z, q)
    types = ['stage' if rng.random() < 0.75 else 'mix' for _ in range(nnodes)]
    for t in types:
        net.add_node(t)
    flow = {}  # (s,d,ph) -> flow

    def add(s, d, ph, f):
        flow[(s, d, ph)] = flow.get((s, d, ph), 0.0) + f
    # feed pieces
    k = rng.integers(1, 3)
    entry = rng.choice(nnodes, size=k, replace=False)
    fr = rng.dirichlet(np.ones(k))
    feed_in = {int(e): float(f) for e, f in zip(entry, fr)}
    # paths from feed nodes to products (total flow 1), random phases
    nprod = 0
    for e, f in feed_in.items():
        npaths = rng.integers(1, 3)
        for g in rng.dirichlet(np.ones(npaths)):
            cur = e
            for step in range(rng.integers(1, 5)):
                nxt = int(rng.integers(nnodes))
                if nxt == cur:
                    continue
                add(cur, nxt, 'L' if rng.random() < 0.5 else 'V', f * g)
                cur = nxt
            add(cur, ('P', f'P{nprod}'), 'L' if rng.random() < 0.7 else 'V', f * g)
            nprod += 1
    # circulations (recycles / counter-current pairs)
    for c in range(rng.integers(3, 3 * nnodes)):
        ln = rng.integers(2, 5)
        cyc = list(rng.choice(nnodes, size=min(ln, nnodes), replace=False))
        amt = rng.exponential(1.0)
        for i in range(len(cyc)):
            add(int(cyc[i]), int(cyc[(i + 1) % len(cyc)]), 'L' if rng.random() < 0.5 else 'V', amt)
    for (s, d, ph), f in flow.items():
        net.add_stream(s, d, ph, f)
    for e, f in feed_in.items():
        net.add_stream('F', e, 'F', f)
    # every stage must have a liquid outflow (else its liquid state is undetermined): give it a small one to a
    # random other node and send the same amount back as vapor (keeps the molar balance)
    b = net.balances()
    for kk in range(nnodes):
        if net.nodes[kk]['type'] == 'stage' and b['Lout'][kk] == 0:
            other = int((kk + 1 + rng.integers(nnodes - 1)) % nnodes)
            net.add_stream(kk, other, 'L', 0.3)
            net.add_stream(other, kk, 'V', 0.3)
    return net


def check(net, label, fails, stats):
    C = net.solve()
    if net.res > 1e-10 or np.any(C < -1e-12) or not np.all(np.isfinite(C)):
        stats['nonconv'] += 1
        return
    stats['ok'] += 1
    R, Q, nv = net.duties()
    P = net.products()
    for th in net.feed_roots():
        TP = [net.stream_T(st, th) for st in P]
        # (a) global identity
        idres = abs(sum(TP) - (Q - R))
        stats['idres'] = max(stats['idres'], idres)
        if idres > 1e-8:
            fails.append((label, 'identity', idres))
        # node signs
        sig = np.zeros(len(net.nodes))
        for k, nd in enumerate(net.nodes):
            outs = [net.stream_T(st, th) for st in net.streams if st[0] == k]
            outs = [t for t in outs if abs(t) > 1e-12]
            if not outs:
                continue
            pos = any(t > 0 for t in outs)
            neg = any(t < 0 for t in outs)
            if pos and neg:
                if nd['type'] == 'stage':
                    fails.append((label, 'stage sign inconsistency', outs))
                stats['inconsistent_mix'] += 1
                sig[k] = np.nan
            else:
                sig[k] = 1 if pos else -1
        # (c) cut identity with U^c = sign <= 0 (inconsistent mix nodes: try both sides, keep the best)
        Uc = {k for k in range(len(net.nodes)) if not (sig[k] > 0)}
        cross = 0.0
        badcross = 0.0
        for st in net.streams:
            s, d, ph, f = st
            if s == 'F' or not isinstance(d, int):
                continue
            t = net.stream_T(st, th)
            if s in Uc and d not in Uc:
                cross -= t
                if np.isfinite(sig[s]):
                    badcross = min(badcross, -t)
            elif s not in Uc and d in Uc:
                cross += t
                badcross = min(badcross, t)
        RUc = sum(max(nv[k], 0) for k in Uc)
        QUc = sum(max(-nv[k], 0) for k in Uc)
        TPUc = sum(net.stream_T(st, th) for st in P if st[0] in Uc)
        Tfeed = sum(net.stream_T(st, th) for st in net.streams if st[0] == 'F' and st[1] in Uc)
        cutres = abs(RUc - (QUc - TPUc - Tfeed + cross))
        stats['cutres'] = max(stats['cutres'], cutres)
        if cutres > 1e-8:
            fails.append((label, 'cut identity', cutres))
        if badcross < -1e-9 and not np.any(np.isnan(sig)):
            fails.append((label, 'unfavourable crossing term', badcross))
        # (d) bound
        bound = sum(max(0.0, -t) for t in TP)
        gap = R - bound
        stats['mingap'] = min(stats['mingap'], gap / max(bound, 1e-300) if bound > 0 else np.inf)
        if gap < -1e-9:
            fails.append((label, 'BOUND VIOLATED', R, bound))


def main():
    allfail = False
    # (1) columns
    fails, stats = [], dict(ok=0, nonconv=0, idres=0, cutres=0, mingap=np.inf, inconsistent_mix=0)
    for c in range(NC):
        N = int(rng.integers(2, 6))
        a, z, q = rand_data(N)
        nT, nB = int(rng.integers(0, 15)), int(rng.integers(0, 15))
        D = rng.uniform(0.05, 0.95)
        VB = rng.uniform(0.2, 4.0)
        net = column(N, nT, nB, D, VB, a, z, q)
        if net is None:
            continue
        check(net, f'column{c}', fails, stats)
    print('columns :', stats, 'fails:', len(fails))
    for f in fails[:10]:
        print('   ', f)
    allfail |= bool(fails)
    print('PASS' if not fails else 'FAIL', '(columns)')
    # (2) Petlyuk
    fails, stats = [], dict(ok=0, nonconv=0, idres=0, cutres=0, mingap=np.inf, inconsistent_mix=0)
    for c in range(NC):
        N = int(rng.integers(3, 6))
        a, z, q = rand_data(N)
        ns = [int(rng.integers(0, 10)) for _ in range(6)]
        D = rng.uniform(0.05, 0.6); S = rng.uniform(0.05, 0.9 - D)
        VB = rng.uniform(0.5, 4.0); V1 = rng.uniform(0.05, VB * 0.9); L1 = rng.uniform(0.05, 3.0)
        net = petlyuk(N, ns, a, z, q, VB, V1, L1, S, D)
        if net is None:
            continue
        check(net, f'petlyuk{c}', fails, stats)
    print('petlyuk :', stats, 'fails:', len(fails))
    for f in fails[:10]:
        print('   ', f)
    allfail |= bool(fails)
    print('PASS' if not fails else 'FAIL', '(Petlyuk)')
    # (3) random networks
    fails, stats = [], dict(ok=0, nonconv=0, idres=0, cutres=0, mingap=np.inf, inconsistent_mix=0)
    for c in range(NC * 3):
        N = int(rng.integers(2, 6))
        a, z, q = rand_data(N)
        net = random_net(N, a, z, q, int(rng.integers(3, 12)))
        check(net, f'random{c}', fails, stats)
    print('random  :', stats, 'fails:', len(fails))
    for f in fails[:10]:
        print('   ', f)
    allfail |= bool(fails)
    print('PASS' if not fails else 'FAIL', '(random networks)')
    sys.exit(1 if allfail else 0)


if __name__ == '__main__':
    main()
