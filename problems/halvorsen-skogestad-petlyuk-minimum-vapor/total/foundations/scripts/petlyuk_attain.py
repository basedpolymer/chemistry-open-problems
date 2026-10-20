"""Independent check (Agent F) of attainment for the ternary Petlyuk arrangement in the finite-stage MESH model.

Agent F, Claude Opus 5.5 (effort xhigh). Environment: Python 3 + numpy + scipy (WSL Fedora sci env).
Usage: python petlyuk_attain.py [n1,n2,...] [eps1,eps2,...] [x0 as 5 comma-separated flows]
       python petlyuk_attain.py fixed [delta] [n1,n2,...] [q]   (fixed flows near the Underwood attainment point)
Runtime: minutes to tens of minutes (finite-difference SLSQP). Deterministic.

Arrangement (all units of class A, built with unitnet.py): prefractionator C1 (top section T1, feed stage F1,
bottom section B1, n stages each) thermally coupled to the main column at junction stages J1 (top) and J2 (bottom);
main column: total condenser, section S1, J1, section S2, side-draw stage SD (liquid product B), section S3, J2,
section S4, partial reboiler (the only heat input, R = its vapor). Free flows: V1T, L1T (C1 top), V4 (= R), B_C,
S_B; the others follow from CMO. For each (n, eps): minimize R subject to every product purity >= 1 - eps (SLSQP
on the five flows, each evaluation a Newton steady-state solve warm-started from the previous one). Reports R, the
stage-level bound max_j sum_P max(0, phi_P(theta_j)) at the optimum, and max_j H_j. Numerical evidence only.
"""
import sys
import numpy as np
from scipy.optimize import minimize
from unitnet import Net

ALPHA, Z, Q = (4.0, 2.0, 1.0), (1 / 3, 1 / 3, 1 / 3), 1.0


def build_flows(flows):
    V1T, L1T, V4, BC, SB = flows
    q = Q
    V1B, L1B = V1T - (1 - q), L1T + q
    L4 = V4 + BC
    V3, L3 = V4 - V1B, L4 - L1B
    V2, L2 = V3, L3 + SB
    V1, L1 = V2 + V1T, L2 + L1T
    return dict(V1B=V1B, L1B=L1B, V3=V3, L3=L3, L2=L2, L1=L1, D=V1 - L1, L1T=L1T)


def build(flows, n):
    V1T, L1T, V4, BC, SB = flows
    q = Q
    V1B, L1B = V1T - (1 - q), L1T + q
    L4 = V4 + BC
    V3, L3 = V4 - V1B, L4 - L1B
    V2, L2 = V3, L3 + SB
    V1, L1 = V2 + V1T, L2 + L1T
    D = V1 - L1
    fl = dict(V1B=V1B, L1B=L1B, L4=L4, V3=V3, L3=L3, V2=V2, L2=L2, V1=V1, L1=L1, D=D, V1T=V1T, L1T=L1T, V4=V4)
    if min(fl.values()) <= 1e-6 or BC <= 1e-6 or SB <= 1e-6:
        return None, fl
    net = Net(ALPHA, Z, Q)
    U = net.unit
    cond = U('C')
    S1 = [U('E') for _ in range(n)]
    J1 = U('E')
    S2 = [U('E') for _ in range(n)]
    SD = U('E')
    S3 = [U('E') for _ in range(n)]
    J2 = U('E')
    S4 = [U('E') for _ in range(n)]
    RB = U('E')
    T1 = [U('E') for _ in range(n)]
    F1 = U('E')
    B1 = [U('E') for _ in range(n)]

    def chain(top, stages, bot, L, V):
        seq = [top] + stages + [bot]
        for a, b in zip(seq[:-1], seq[1:]):
            net.s(a, b, 'L', L)
            net.s(b, a, 'V', V)
    # main column (top to bottom); condenser: vapor V1 in, reflux L1 to the top stage, product D
    seq = S1 + [J1]
    net.s(seq[0], cond, 'V', V1)
    net.s(cond, seq[0], 'L', L1)
    net.s(cond, ('P', 0), 'L', D)
    for a, b in zip(seq[:-1], seq[1:]):
        net.s(a, b, 'L', L1)
        net.s(b, a, 'V', V1)
    chain(J1, S2, SD, L2, V2)
    net.s(SD, ('P', 1), 'L', SB)
    chain(SD, S3, J2, L3, V3)
    chain(J2, S4, RB, L4, V4)
    net.s(RB, ('P', 2), 'L', BC)
    # prefractionator: T1 top stage <-> J1, B1 bottom stage <-> J2
    net.s(T1[0], J1, 'V', V1T)
    net.s(J1, T1[0], 'L', L1T)
    # explicit flows in C1: above F1 (T1 and T1->F1) L1T/V1T, below F1 L1B/V1B
    seq_top = T1 + [F1]
    for a, b in zip(seq_top[:-1], seq_top[1:]):
        net.s(a, b, 'L', L1T)
        net.s(b, a, 'V', V1T)
    seq_bot = [F1] + B1
    for a, b in zip(seq_bot[:-1], seq_bot[1:]):
        net.s(a, b, 'L', L1B)
        net.s(b, a, 'V', V1B)
    net.s(seq_bot[-1], J2, 'L', L1B)
    net.s(J2, seq_bot[-1], 'V', V1B)
    net.s('F', F1, 'F', 1.0)
    return net, fl


class Model:
    def __init__(self, n):
        self.n = n
        self.X = None
        self.cache = {}

    def solve(self, flows):
        key = tuple(np.round(flows, 12))
        if key in self.cache:
            return self.cache[key]
        net, fl = build(flows, self.n)
        if net is None:
            return None
        try:
            if self.X is not None and self.X.shape[0] == len(net.types):
                net.X = self.X.copy()
            else:
                net.solve(iters=3000)
            net.newton(iters=60)
            if not np.isfinite(net.closure) or net.closure > 1e-9:
                net.solve(iters=3000)
                net.newton(iters=100)
        except np.linalg.LinAlgError:
            return None
        if not np.isfinite(net.closure) or net.closure > 1e-9 or net.X.min() < -1e-9:
            return None
        self.X = net.X.copy()
        pur = [net.comp(s[0], s[2])[s[1][1]] for s in net.streams if not isinstance(s[1], int)]
        res = (net, pur)
        self.cache[key] = res
        return res


X0 = None


def study(ns, epss):
    H = None
    rows = []
    x0 = np.array([float(v) for v in X0.split(',')]) if X0 else np.array([0.9, 0.4, 1.7, 1 / 3, 1 / 3])
    for eps in epss:
        for n in ns:
            m = Model(n)
            # continuation in n: start from the previous optimum
            def f(x):
                return x[2]

            def cons(x):
                r = m.solve(x)
                if r is None:
                    return -np.ones(3)
                return np.array(r[1]) - (1 - eps)
            def flowcons(x):
                fl = build_flows(x)
                return np.array(list(fl.values())) - 0.02
            res = minimize(f, x0, method='SLSQP', constraints=[{'type': 'ineq', 'fun': flowcons},
                                                              {'type': 'ineq', 'fun': cons}],
                           bounds=[(0.05, 5), (0.01, 5), (0.5, 5), (0.2, 0.45), (0.2, 0.45)],
                           options=dict(maxiter=200, ftol=1e-12, eps=1e-7))
            r = m.solve(res.x)
            if r is None:
                print(f"n={n} eps={eps}: no solution", flush=True)
                continue
            net, pur = r
            R, out, resid, viol = net.check(verbose=False)
            bound = max(b for (_, b, _) in out)
            H = max(h for (_, _, h) in out)
            ok = min(pur) >= 1 - eps - 1e-7
            rows.append((n, eps, R, bound, H, min(pur), ok, resid, viol))
            print(f"n={n:3d} eps={eps:g}: R={R:.6f}  bound={bound:.6f}  R-bound={R - bound:.2e}  max_j H_j={H:.6f}"
                  f"  R/maxH={R / H:.5f}  min purity={min(pur):.6f}  residual={resid:.1e}  sign-viol={viol}"
                  f"  {'feasible' if ok else 'INFEASIBLE'}  flows={np.round(res.x, 8).tolist()}", flush=True)
            if ok:
                x0 = res.x
    return rows


def ternary_point(q):
    """Underwood attainment point for alpha = (4,2,1), z = 1/3 and feed quality q: prefractionator preferred split
    (d_A = z_A, d_C = 0, both feed roots active) and max_j H_j."""
    from scipy.optimize import brentq
    a, z = np.array(ALPHA), np.array(Z)
    g = lambda t: np.sum(a * z / (a - t)) - (1 - q)
    th = [brentq(g, a[1] + 1e-12, a[0] - 1e-12), brentq(g, a[2] + 1e-12, a[1] - 1e-12)]
    # V = a1 z1/(a1 - th) + a2 dB/(a2 - th) at both roots
    A = np.array([[-a[1] / (a[1] - th[0]), 1.0], [-a[1] / (a[1] - th[1]), 1.0]])
    rhs = np.array([a[0] * z[0] / (a[0] - th[0]), a[0] * z[0] / (a[0] - th[1])])
    dB, V1T = np.linalg.solve(A, rhs)
    L1T = V1T - (z[0] + dB)
    H = [np.sum(a[j + 1:] * z[j + 1:] / (th[j] - a[j + 1:])) for j in range(2)]
    return V1T, L1T, max(H), dB


def fixed_flow_study(delta, ns):
    """Flows fixed at the Underwood attainment point (prefractionator at its preferred split, R = max_j H_j) plus
    delta on both prefractionator flows and 2 delta on the reboiler vapor; report product impurities vs n."""
    V1T, L1T, H2, dB = ternary_point(Q)
    ok = True
    for n in ns:
        net, fl = build(np.array([V1T + delta, L1T + delta, H2 + 2 * delta, 1 / 3, 1 / 3]), n)
        net.solve(iters=4000)
        net.newton(iters=100)
        pur = [net.comp(e[0], e[2])[e[1][1]] for e in net.streams if not isinstance(e[1], int)]
        R, out, resid, viol = net.check(verbose=False)
        ok &= viol == 0 and all(R >= b - 1e-9 for (_, b, _) in out)
        print(f"fixed flows, delta={delta:g}, n={n:3d}: R={R:.6f} (max_j H_j={H2:.6f}, q={Q:g}), closure={net.closure:.1e}, "
              f"impurities={[float('%.2e' % (1 - p)) for p in pur]}, sign-viol={viol}", flush=True)
    return ok


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'fixed':
        if len(sys.argv) > 4:
            Q = float(sys.argv[4])
        okf = fixed_flow_study(float(sys.argv[2]) if len(sys.argv) > 2 else 0.02,
                               [int(v) for v in (sys.argv[3] if len(sys.argv) > 3 else '10,20,40,80').split(',')])
        print("fixed-flow run: bound respected, no sign violation: " + ("PASS" if okf else "FAIL"))
        sys.exit(0 if okf else 1)
    ns = [int(v) for v in (sys.argv[1] if len(sys.argv) > 1 else '15,30,60').split(',')]
    epss = [float(v) for v in (sys.argv[2] if len(sys.argv) > 2 else '1e-2,1e-3').split(',')]
    X0 = sys.argv[3] if len(sys.argv) > 3 else None
    rows = study(ns, epss)
    bad = [r for r in rows if r[2] < r[3] - 1e-9 or r[8] > 0]
    print("bound respected at every optimum, no sign-uniformity violation: " + ("PASS" if not bad and rows else "FAIL"))
    sys.exit(0 if rows and not bad else 1)
