#!/usr/bin/env python3
"""Exact checks of the two algebraic ingredients of the stage-level cut theorem (Agent N, network.md, Sections 2-4).

(A) Cut identity (Proposition 2), exact rational arithmetic, on random unit graphs: random directed multigraphs of
    units with random phases and random nonnegative rational component flows on every stream, built so that every
    unit conserves every component (a feed -> products flow decomposition plus random circulations = recycles);
    the CMO energy balance then DEFINES each unit's duty Q_u = sum_out v - sum_in v - (feed vapor).  For random
    rational theta (not only feed roots) and random unit subsets W it checks
        sum_{u in W} Q_u = sum_{e: W -> outside or product} phi_e - sum_{e: outside -> W} phi_e - phi_F(W),
    with phi_e(theta) = v_e - sum_i alpha_i n_{e,i}/(alpha_i - theta) and phi_F(W) the flux of the feed portions
    entering W.  No equilibrium relation is needed for the identity.
(B) Sign uniformity (Lemma 3), symbolic (sympy) for N = 2..6: for a liquid of composition x and the vapor
    y = alpha x / S(x) in equilibrium with it, phi(vapor of flow V) = -(theta V / S(x)) G_x(theta) and
    phi(liquid of flow L) = -L G_x(theta), G_x(theta) = sum alpha_i x_i/(alpha_i - theta), identically in
    (alpha, x, theta).  Also: a vapor of the SAME composition x as the liquid (outlet of a total vaporizer,
    sum x = 1) has phi = V (1 - G_x(theta)), whose sign differs from that of the liquid when 0 < G_x < 1; this is
    why a total vaporizer (or condenser) is a unit of its own, with a single outflow composition.

Environment: Python 3, sympy.  Runtime: about 20 s.  Prints PASS/FAIL per claim, exits 1 on failure.
"""
import random
import sys
from fractions import Fraction as Fr
import sympy as sp


def phi(stream, alpha, th):
    n, vap = stream['n'], stream['vap']
    return (sum(n) if vap else 0) - sum(a * ni / (a - th) for a, ni in zip(alpha, n))


def random_graph(rng):
    N = rng.randint(2, 5)
    alpha = sorted({Fr(rng.randint(11, 99), 10) for _ in range(N)}, reverse=True)
    N = len(alpha)
    z = [Fr(rng.randint(1, 9)) for _ in range(N)]
    s = sum(z)
    z = [zi / s for zi in z]
    q = Fr(rng.randint(-5, 15), 10)
    nu = rng.randint(2, 9)
    streams = {}  # (src, dst) multigraph: list of dicts; src/dst: unit index, 'F' (feed) or ('P', k)

    def add(src, dst, n, vap):
        streams.setdefault((src, dst, len(streams)), {'n': list(n), 'vap': vap})

    # component-wise flow decomposition: feed portions -> random walks -> product outlets
    feed_portions = {}
    nprod = rng.randint(1, 4)
    for p in range(N):
        remaining = z[p]
        parts = rng.randint(1, 3)
        for t in range(parts):
            amt = remaining if t == parts - 1 else remaining * Fr(rng.randint(1, 9), 10)
            remaining -= amt
            if amt == 0:
                continue
            u0 = rng.randrange(nu)
            feed_portions.setdefault(u0, [Fr(0)] * N)
            feed_portions[u0][p] += amt
            cur = u0
            for _ in range(rng.randint(0, 5)):
                nxt = rng.randrange(nu)
                if nxt == cur:
                    continue
                n = [Fr(0)] * N
                n[p] = amt
                add(cur, nxt, n, rng.random() < 0.5)
                cur = nxt
            n = [Fr(0)] * N
            n[p] = amt
            add(cur, ('P', rng.randrange(nprod)), n, rng.random() < 0.3)
    # circulations (recycles) carrying several components
    for _ in range(rng.randint(0, 4)):
        cyc = rng.sample(range(nu), rng.randint(2, min(nu, 4)))
        n = [Fr(rng.randint(0, 5), 7) for _ in range(N)]
        for a_, b_ in zip(cyc, cyc[1:] + cyc[:1]):
            add(a_, b_, n, rng.random() < 0.5)
    # the feed: portions have the composition z only in total; here each portion carries the vector assigned;
    # to respect "feed portions have composition z", rescale: merge all portions into one feed unit per portion
    # proportional to z (a feed splitter unit 'FS' sends z-proportional portions).
    return alpha, z, q, nu, streams, feed_portions


def check_identity(rng, ntests=300):
    ok = True
    count = 0
    for _ in range(ntests):
        alpha, z, q, nu, streams, fp = random_graph(rng)
        N = len(alpha)
        # Build an explicit feed splitter: unit index nu ('FS') receives the whole feed (composition z, vapor 1-q)
        # and sends to each unit u0 a stream of liquid carrying fp[u0]; the feed splitter itself then has a duty
        # (it is just another unit for the identity).  This keeps the identity exact for an arbitrary assignment.
        FS = nu
        for u0, n in fp.items():
            streams[(FS, u0, len(streams))] = {'n': list(n), 'vap': False}
        units = list(range(nu + 1))
        # duties Q_u = sum_out v - sum_in v - feed vapor ; check component balances first
        for u in units:
            bal = [Fr(0)] * N
            for (s, d, _), st in streams.items():
                if s == u:
                    bal = [b + x for b, x in zip(bal, st['n'])]
                if d == u:
                    bal = [b - x for b, x in zip(bal, st['n'])]
            if u == FS:
                bal = [b - zi for b, zi in zip(bal, z)]
            if any(b != 0 for b in bal):
                print("component balance violated at unit", u)
                ok = False
        Q = {}
        for u in units:
            vout = sum(sum(st['n']) for (s, d, _), st in streams.items() if s == u and st['vap'])
            vin = sum(sum(st['n']) for (s, d, _), st in streams.items() if d == u and st['vap'])
            Q[u] = vout - vin - ((1 - q) if u == FS else 0)
        for _t in range(6):
            th = Fr(rng.randint(1, 200), 13)
            if th in alpha:
                continue
            phiF = (1 - q) - sum(a * zi / (a - th) for a, zi in zip(alpha, z))
            W = {u for u in units if rng.random() < 0.5}
            lhs = sum(Q[u] for u in W)
            rhs = Fr(0)
            for (s, d, _), st in streams.items():
                f = phi(st, alpha, th)
                if s in W and (d not in W):
                    rhs += f
                if (s not in W) and d in W:
                    rhs -= f
            if FS in W:
                rhs -= phiF
            count += 1
            if lhs != rhs:
                ok = False
    return ok, count


def check_sign_symbolic():
    ok = True
    th = sp.Symbol('theta', positive=True)
    for N in range(2, 7):
        a = sp.symbols(f'a1:{N + 1}', positive=True)
        x = sp.symbols(f'x1:{N + 1}', positive=True)
        V, L = sp.symbols('V L', positive=True)
        S = sum(ai * xi for ai, xi in zip(a, x))
        y = [ai * xi / S for ai, xi in zip(a, x)]
        G = sum(ai * xi / (ai - th) for ai, xi in zip(a, x))
        Xs = sum(x)  # liquid composition need not be normalized for the identity; use sum x = 1 below
        phi_vap = V * sum(y) - sum(ai * V * yi / (ai - th) for ai, yi in zip(a, y))
        phi_liq = -sum(ai * L * xi / (ai - th) for ai, xi in zip(a, x))
        e1 = sp.simplify(sp.together(phi_vap + th * V / S * G))
        e2 = sp.simplify(phi_liq + L * G)
        # vapor of the same composition as the liquid (normalized: sum x = 1): phi = V(1 - G) ... with sum x = 1
        phi_same = V * Xs - sum(ai * V * xi / (ai - th) for ai, xi in zip(a, x))
        e3 = sp.simplify(phi_same - V * (Xs - G))
        if e1 != 0 or e2 != 0 or e3 != 0:
            ok = False
            print("symbolic check failed for N =", N, e1, e2, e3)
    return ok


def main():
    rng = random.Random(20261020)
    ok1, cnt = check_identity(rng)
    print(f"cut identity: {cnt} (graph, theta, subset) tests in exact rationals")
    print(("PASS" if ok1 else "FAIL") + "  node balances sum to the cut identity for arbitrary unit subsets (exact)")
    ok2 = check_sign_symbolic()
    print(("PASS" if ok2 else "FAIL") + "  equilibrium outflows: phi_vap = -(theta V/S) G_x, phi_liq = -L G_x (symbolic, N = 2..6)")
    sys.exit(0 if (ok1 and ok2) else 1)


if __name__ == "__main__":
    main()
