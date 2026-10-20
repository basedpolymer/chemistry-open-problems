#!/usr/bin/env python3
"""Numerical check of the stage-level cut theorem (Agent N, network.md, Theorem 1) on steady states of networks of
equilibrium stages: random section/junction networks (splitters, mixing of arbitrary streams, recycles, several feed
points, products in either phase, heaters/coolers anywhere, reverse net flows) and a ternary Petlyuk column.

For each steady state and each feed root theta_j:
  (i)   the cut identity  R(U^c) - Q(U^c) = H'_j + Phi(U^c->U) - Phi(U->U^c)  with U = {stages with u(x;theta_j) >= 0}
        holds (floating point residual);
  (ii)  every stream U^c -> U has flux >= 0 and every stream U -> U^c has flux <= 0;
  (iii) hence R >= H'_j (H'_j = flux of the products leaving U^c; = H_j for pure products).
Floating point: evidence for bookkeeping/sign conventions only; the proof is analytic (network.md).

Environment: Python 3, numpy.  Runtime: about 1-3 minutes.  Usage: python cut_bound_check.py [ntrials] [seed]
"""
import sys
import numpy as np
from stagenet import feed_roots, cut_check
from secnet import build


def peaks(alpha, z, q):
    th = feed_roots(alpha, z, q)
    H = [sum(alpha[p] * z[p] / (t - alpha[p]) for p in range(j + 1, len(alpha))) for j, t in enumerate(th)]
    return th, np.array(H)


def random_network(rng):
    N = int(rng.integers(2, 6))
    alpha = np.sort(rng.uniform(1.0, 12.0, N))[::-1]
    alpha = alpha / alpha[-1]
    z = rng.dirichlet(np.ones(N))
    q = float(rng.uniform(-0.3, 1.3))
    nj = int(rng.integers(2, 7))
    nsec = int(rng.integers(nj - 1, nj + 4))
    # random connected multigraph of sections between junctions
    edges = []
    for j in range(1, nj):
        i = int(rng.integers(0, j))
        edges.append((i, j) if rng.random() < 0.5 else (j, i))
    while len(edges) < nsec:
        i, j = rng.choice(nj, 2, replace=False)
        edges.append((int(i), int(j)))
    # net upward flows D_s: route the feed to product junctions along random tree paths, plus circulations
    feedj = int(rng.integers(0, nj))
    nprod = int(rng.integers(1, min(nj, 4) + 1))
    prodj = list(rng.choice(nj, nprod, replace=False))
    shares = rng.dirichlet(np.ones(nprod))
    D = np.zeros(len(edges))
    # undirected adjacency for path finding
    import collections
    adj = collections.defaultdict(list)
    for si, (b, t) in enumerate(edges):
        adj[b].append((t, si, +1))   # moving b -> t is upward along s
        adj[t].append((b, si, -1))

    def path(src, dst):
        prev = {src: None}
        dq = collections.deque([src])
        while dq:
            v = dq.popleft()
            if v == dst:
                break
            order = list(adj[v])
            rng.shuffle(order)
            for (w, si, sg) in order:
                if w not in prev:
                    prev[w] = (v, si, sg)
                    dq.append(w)
        out = []
        v = dst
        while prev[v] is not None:
            u_, si, sg = prev[v]
            out.append((si, sg))
            v = u_
        return out
    for pj, sh in zip(prodj, shares):
        for si, sg in path(feedj, pj):
            D[si] += sg * sh
    # circulations (recycles): random cycles
    for _ in range(int(rng.integers(0, 3))):
        a_, b_ = rng.choice(nj, 2, replace=False)
        amt = rng.uniform(-0.5, 0.5)
        p1 = path(int(a_), int(b_))
        p2 = path(int(b_), int(a_))
        for si, sg in p1 + p2:
            D[si] += sg * amt
    sections = []
    for si, (b, t) in enumerate(edges):
        V = max(D[si], 0) + rng.uniform(0.05, 2.0)
        L = V - D[si]
        sections.append((b, t, int(rng.integers(0, 6)), V, L))
    products = []
    for k, (pj, sh) in enumerate(zip(prodj, shares)):
        ph = 'V' if rng.random() < 0.25 else 'L'
        products.append((int(pj), ph, float(sh), f"P{k}"))
    feed_split = None
    if rng.random() < 0.3:
        # split the feed between two junctions: adjust by routing the second part along a path
        j2 = int(rng.integers(0, nj))
        fr = float(rng.uniform(0.2, 0.8))
        feed_split = [(feedj, 1 - fr), (j2, fr)]
        for si, sg in path(j2, feedj):
            # the part fr now enters at j2 instead of feedj: shift it along a path j2 -> feedj
            D_shift = sg * fr
            sections[si] = (sections[si][0], sections[si][1], sections[si][2],
                            sections[si][3] + max(0, D_shift) + 0.0, sections[si][4] + max(0, D_shift) - D_shift)
    net, J, _ = build(alpha, z, q, nj, sections, feedj, products, feed_split)
    return net, alpha, z, q


def petlyuk(alpha, z, q, nst, VT_pre, split_frac_liq, VB_main_extra):
    """Ternary Petlyuk: junctions 0=C (condenser), 1=J1, 2=J2 (side draw B), 3=J3, 4=Rb (reboiler), 5=F (prefrac feed)."""
    # products: A at C (liquid), B at J2 (liquid side draw), C at Rb (liquid)
    dA, dB, dC = z
    # prefractionator: top section F->J1 with V_pt, L_pt; bottom section J3->F with V_pb, L_pb
    V_pt = VT_pre
    L_pt = split_frac_liq
    D_pt = V_pt - L_pt
    D_pb = D_pt - 1.0          # net up through the bottom section: feed balance at F: D_pt - D_pb = 1 (all to sections)
    V_pb = V_pt - (1 - q)
    L_pb = V_pb - D_pb
    # main column: section 1: J1 -> C, D1 = dA ; section 2: J2 -> J1, D2 = dA - D_pt ; section 3: J3 -> J2, D3 = D2 + dB
    # section 4: Rb -> J3, D4 = D3 + D_pb ... check balance at Rb: -D4 - dC = 0 -> D4 = -dC
    V4 = VB_main_extra
    D4 = -dC
    L4 = V4 - D4
    V3 = V4 - V_pb
    D3 = D4 - D_pb
    L3 = V3 - D3
    V2 = V3
    D2 = D3 - dB
    L2 = V2 - D2
    V1 = V2 + V_pt
    D1 = D2 + D_pt
    L1 = V1 - D1
    sections = [(1, 0, nst, V1, L1), (2, 1, nst, V2, L2), (3, 2, nst, V3, L3), (4, 3, nst, V4, L4),
                (5, 1, nst, V_pt, L_pt), (3, 5, nst, V_pb, L_pb)]
    products = [(0, 'L', dA, 'A'), (2, 'L', dB, 'B'), (4, 'L', dC, 'C')]
    if min(L1, L2, L3, L4, L_pt, L_pb, V1, V2, V3, V4, V_pt, V_pb) <= 0:
        return None
    net, J, _ = build(alpha, z, q, 6, sections, 5, products)
    return net


def check_one(net, alpha, z, q, tag, out):
    X, res = net.solve()
    if not np.isfinite(res) or res > 1e-9 or np.min(X) < -1e-12:
        out['skipped'] += 1
        return None
    th, H = peaks(alpha, z, q)
    for j, t in enumerate(th):
        c = cut_check(net, X, t)
        out['n'] += 1
        out['max_id_res'] = max(out['max_id_res'], abs(c['identity_residual']))
        if c['min_in'] < -1e-9 or c['max_out'] > 1e-9:
            out['sign_fail'] += 1
        gap = c['R'] - c['Hprime']
        out['min_gap'] = min(out['min_gap'], gap)
        if gap < -1e-9:
            out['bound_fail'] += 1
            print("BOUND FAIL", tag, j, c['R'], c['Hprime'])
    return X


def main():
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    rng = np.random.default_rng(seed)
    out = dict(n=0, skipped=0, max_id_res=0.0, sign_fail=0, bound_fail=0, min_gap=np.inf)
    neg_flow_nets = 0
    for t in range(ntr):
        net, alpha, z, q = random_network(rng)
        Lout, Vout, Lin, Vin, Fl, Fv, delta, massres = net.flows()
        if np.max(np.abs(massres)) > 1e-9 or np.min(Lout) < 0 or np.min(Vout) < 0:
            out['skipped'] += 1
            continue
        check_one(net, alpha, z, q, f"random{t}", out)
    print(f"random networks: {out['n']} (network, theta_j) tests, {out['skipped']} networks skipped (solver or "
          f"inconsistent flows); max identity residual {out['max_id_res']:.2e}; min R - H'_j = {out['min_gap']:.4f}")
    ok = out['max_id_res'] < 1e-8 and out['sign_fail'] == 0 and out['bound_fail'] == 0 and out['n'] > 0
    print(("PASS" if out['max_id_res'] < 1e-8 else "FAIL") + "  cut identity (random networks)")
    print(("PASS" if out['sign_fail'] == 0 else "FAIL") + "  crossing streams have the right flux sign (random networks)")
    print(("PASS" if out['bound_fail'] == 0 else "FAIL") + "  R >= H'_j (random networks)")
    # Petlyuk near its minimum
    alpha = np.array([4.0, 2.0, 1.0])
    z = np.array([1 / 3, 1 / 3, 1 / 3])
    q = 1.0
    th, H = peaks(alpha, z, q)
    out2 = dict(n=0, skipped=0, max_id_res=0.0, sign_fail=0, bound_fail=0, min_gap=np.inf)
    best = None
    for nst in (8, 20, 40):
        for VTp in np.linspace(0.6, 1.4, 4):
            for Lp in np.linspace(0.1, 0.9, 4):
                for V4 in np.linspace(1.2, 2.2, 6):
                    net = petlyuk(alpha, z, q, nst, VTp, Lp, V4)
                    if net is None:
                        continue
                    X = check_one(net, alpha, z, q, "petlyuk", out2)
                    if X is None:
                        continue
                    # purity
                    pur = min(X[0][0], X[2][1], X[4][2])
                    if pur > 0.99 and (best is None or V4 < best[0]):
                        best = (V4, nst, pur)
    print(f"Petlyuk ternary (alpha=4,2,1, z=1/3, q=1, H = {H.round(5)}, max H = {H.max():.5f}): {out2['n']} tests, "
          f"max identity residual {out2['max_id_res']:.2e}, min R - H'_j = {out2['min_gap']:.4f}; "
          f"lowest reboil with all purities > 0.99 in the grid: {best}")
    ok2 = out2['max_id_res'] < 1e-8 and out2['sign_fail'] == 0 and out2['bound_fail'] == 0 and out2['n'] > 0
    print(("PASS" if ok2 else "FAIL") + "  cut identity, flux signs and R >= H'_j (Petlyuk grid)")
    sys.exit(0 if (ok and ok2) else 1)


if __name__ == "__main__":
    main()
