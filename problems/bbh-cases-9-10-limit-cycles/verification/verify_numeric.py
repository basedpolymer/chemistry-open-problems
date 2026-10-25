#!/usr/bin/env python3
"""Verifier's numerical search for counterexamples (floating point; it can refute, it cannot prove).

Written independently of scripts/numerical_sanity.py. Everything is done in the normalised coordinates of (F3)
(equilibrium at (1,1), one parameter rho = xbar/ybar > 0), which I re-derived in verify_symbolic.py (V1).

Section. On S = {x = 1, y > 1} we have xdot = alpha1 (1 - y) (case 9) or alpha2 (1 - y) (case 10), so S is
transversal, and every periodic orbit surrounding (1,1) crosses it (it must cross x = 1 above the equilibrium).
Return map: start at (1, 1 + r), integrate forward until the next crossing of S; delta(r) = y_return - 1 - r.
A periodic orbit through (1, 1 + r) means delta(r) = 0; since the equilibrium is attracting for the non-Hopf
families (tr J < 0 for every rho), a periodic orbit forces delta to reach 0 from below somewhere.

R1  random members of 9a, 9b, 10a (coefficients <= 6) plus the BBH example and members with k = 1 exactly and k just
    above 1; rho on a log grid 1e-4 .. 1e3; r on a log grid 1e-3 .. 1e3. Claim: delta(r) < 0 wherever the orbit
    returns (no sign change, no near-zero delta/r).
R2  positive controls: Hopf networks of cases 9 and 10, rho a few per cent beyond the Hopf value (tr J > 0):
    delta changes sign (a stable cycle exists), so the detector works.
R3  continuation across the Hopf boundary in network space: pairs (Hopf network with k < 1, non-Hopf network with
    k >= 1) differing in one coefficient; scan rho finely around the Hopf value of the Hopf member.
R4  hypothesis (H2) of the comparison lemma by direct quadrature of f and g (mpmath, 30 digits), not using the
    closed forms of the paper: F(a) > F(b) at G-paired points a < 1 < b.
Runtime: about 4.5 minutes on 12 processes. Environment: WSL Fedora, Python 3.12, numpy 2.5, scipy 1.18, mpmath.
"""
import math
import multiprocessing as mp
import random
import sys
import time

import mpmath as mpm
import numpy as np
from scipy.integrate import solve_ivp

SRC = {9: [(2, 0), (1, 1), (0, 1)], 10: [(2, 0), (1, 1), (0, 0)]}


def cross(c, d):
    return (c[1] * d[2] - c[2] * d[1], c[2] * d[0] - c[0] * d[2], c[0] * d[1] - c[1] * d[0])


def make(case, c, d, rho):
    u = cross(c, d)
    cu = [c[i] * u[i] for i in range(3)]
    du = [rho * d[i] * u[i] for i in range(3)]
    src = SRC[case]

    def rhs(s, z):
        x, y = z
        mm = [x ** a * y ** b for (a, b) in src]
        fx = cu[0] * mm[0] + cu[1] * mm[1] + cu[2] * mm[2]
        fy = du[0] * mm[0] + du[1] * mm[1] + du[2] * mm[2]
        sc = 1.0 + x * x + y * y          # orbit-preserving time change, avoids blow-up
        return [fx / sc, fy / sc]
    return rhs


def trace(case, c, d, rho):
    u = cross(c, d)
    A = c[0] * u[0] - c[2] * u[2]
    return A - rho * d[0] * u[0] if case == 9 else A + rho * d[1] * u[1]


def ret(args):
    """return (r, delta or None, status)"""
    case, c, d, rho, r = args
    rhs = make(case, c, d, rho)

    def ev_sec(s, z):
        return z[0] - 1.0
    ev_sec.direction = -1.0       # downward crossings of x = 1 happen only at y > 1 (xdot = const*(1-y) there)
    ev_sec.terminal = True

    def ev_esc(s, z):
        return 1e6 - (z[0] + z[1])
    ev_esc.terminal = True

    def ev_wall(s, z):
        return min(z[0], z[1])
    ev_wall.terminal = True

    eps = min(1e-6, 1e-3 * r)

    def ev_eq(s, z):            # entered a tiny neighbourhood of the attracting equilibrium: no return possible
        return math.hypot(z[0] - 1.0, z[1] - 1.0) - eps
    ev_eq.terminal = True
    ev_eq.direction = -1.0
    z0 = [1.0, 1.0 + r]
    # start slightly past the section so the initial point is not detected
    sol = solve_ivp(rhs, (0, 1e-6 * max(r, 1e-3)), z0, method='DOP853', rtol=1e-12, atol=1e-14)
    z1 = sol.y[:, -1]
    smax = 1e8
    try:
        sol = solve_ivp(rhs, (0, smax), z1, method='LSODA', rtol=1e-10, atol=1e-13,
                        events=(ev_sec, ev_esc, ev_wall, ev_eq))
    except Exception as e:  # noqa
        return (r, None, 'error')
    for ze in sol.y_events[0]:
        if ze[1] > 1.0:
            return (r, ze[1] - 1.0 - r, 'return')
    if len(sol.t_events[1]):
        return (r, None, 'escape')
    if len(sol.t_events[2]):
        return (r, None, 'axis')
    if len(sol.t_events[3]):
        return (r, None, 'to equilibrium')
    return (r, None, 'time limit')


def family_members():
    random.seed(4242)
    nets = []
    while len([n for n in nets if n[0] == '9a']) < 8:
        c1, c3 = random.randint(1, 6), random.randint(1, 6)
        d1, d3 = random.randint(c1 + 1, 12), random.randint(-1, c3 - 1)
        if d1 * c3 + d3 * c1 >= 2 * c1 * c3:
            nets.append(('9a', 9, (c1, -1, c3), (d1, -1, d3)))
    while len([n for n in nets if n[0] == '9b']) < 6:
        c1, c3, d1 = random.randint(1, 6), random.randint(1, 6), random.randint(1, 6)
        if c3 * d1 >= c1:
            nets.append(('9b', 9, (c1, -1, c3), (d1, 0, -1)))
    while len([n for n in nets if n[0] == '10a']) < 8:
        c1, c3 = random.randint(1, 6), random.randint(1, 6)
        d1, d3 = random.randint(c1 + 1, 12), random.randint(0, c3 - 1)
        if d1 * c3 + d3 * c1 >= 2 * c1 * c3:
            nets.append(('10a', 10, (c1, -1, c3), (d1, -1, d3)))
    # BBH example and boundary members (k = 1 exactly, k slightly above 1)
    nets += [('9a BBH example k=1', 9, (1, -1, 1), (2, -1, 0)),
             ('9a k=1 (d3=-1)', 9, (1, -1, 1), (3, -1, -1)),
             ('9a k=1.1', 9, (10, -1, 1), (21, -1, 0)),
             ('9b k=1', 9, (2, -1, 1), (2, 0, -1)),
             ('9b k=1 (c1=c3=1)', 9, (1, -1, 1), (1, 0, -1)),
             ('10a k=1', 10, (1, -1, 1), (2, -1, 0)),
             ('10a k=1.1', 10, (10, -1, 1), (21, -1, 0))]
    return nets


def kval(c, d):
    u = cross(c, d)
    return c[2] * u[2] / (c[0] * u[0])


def scan(case, c, d, rhos, rs, pool):
    jobs = [(case, c, d, rho, r) for rho in rhos for r in rs]
    out = pool.map(ret, jobs, chunksize=4)
    res = {}
    for (case_, c_, d_, rho, r), o in zip(jobs, out):
        res.setdefault(rho, []).append(o)
    return res


def main():
    t0 = time.time()
    fails = 0
    pool = mp.Pool(min(16, mp.cpu_count()))
    rhos = list(np.logspace(-4, 3, 15))
    rs = list(np.logspace(-3, 3, 37))
    print('R1 non-Hopf members: max over rho, r of delta(r)/r (must be < 0) and outcome counts')
    worst = -np.inf
    for name, case, c, d in family_members():
        u = cross(c, d)
        assert all(t > 0 for t in u), (name, u)
        k = kval(c, d)
        assert k >= 1 - 1e-12, (name, k)
        assert all(trace(case, c, d, rho) < 0 for rho in rhos)
        res = scan(case, c, d, rhos, rs, pool)
        mx = -np.inf
        stat = {}
        for rho, lst in res.items():
            for (r, dl, st) in lst:
                stat[st] = stat.get(st, 0) + 1
                if dl is not None:
                    mx = max(mx, dl / r)
        worst = max(worst, mx)
        print(f'   {name:20s} c={c} d={d} k={k:.4g}: max delta/r = {mx:.3e}; {stat}', flush=True)
        if not mx < 0:
            fails += 1
            print('FAIL   ^ possible periodic orbit')
    print(('PASS' if worst < 0 else 'FAIL') + f' R1 no sign change of the return-map displacement (worst delta/r = {worst:.3e})')

    print('R2 positive controls (Hopf networks, rho 3% past the Hopf value, tr J > 0)')
    controls = [(9, (2, -1, 1), (3, -1, 0)), (9, (3, -1, 2), (4, -1, 0)), (9, (1, -1, 0), (1, 0, -1)),
                (10, (2, -1, 1), (3, -1, 0)), (10, (3, -1, 2), (4, -1, 1))]
    okc = True
    for case, c, d in controls:
        u = cross(c, d)
        A = c[0] * u[0] - c[2] * u[2]
        assert A > 0
        rhoH = A / (d[0] * u[0]) if case == 9 else A / (-d[1] * u[1])
        rho = rhoH * (0.97 if case == 9 else 0.97)   # tr J = A - rho*(...) > 0 for rho < rhoH
        assert trace(case, c, d, rho) > 0
        res = scan(case, c, d, [rho], list(np.logspace(-3, 1, 41)), pool)[rho]
        sg = [np.sign(dl) for (r, dl, st) in res if dl is not None]
        change = any(sg[i] != sg[i + 1] for i in range(len(sg) - 1))
        okc &= change
        print(f'   case {case} c={c} d={d} rho={rho:.4g} (rho_H={rhoH:.4g}): sign change of delta: {change}')
    print(('PASS' if okc else 'FAIL') + ' R2 the detector finds the stable cycle in every control')
    fails += (not okc)

    print('R3 continuation across the Hopf boundary (Hopf member k<1 vs non-Hopf neighbour k>=1)')
    pairs = [(9, (10, -1, 1), (19, -1, 0), (20, -1, 0)), (9, (10, -1, 1), (19, -1, 0), (21, -1, 0)),
             (9, (3, -1, 1), (5, -1, 0), (6, -1, 0)), (9, (1, -1, 2), (2, -1, -1), (2, -1, 0)),
             (10, (10, -1, 1), (19, -1, 0), (20, -1, 0)), (10, (3, -1, 2), (4, -1, 1), (5, -1, 1))]
    okr = True
    for case, c, dH, dN in pairs:
        uH = cross(c, dH)
        AH = c[0] * uH[0] - c[2] * uH[2]
        rhoH = AH / (dH[0] * uH[0]) if case == 9 else AH / (-dH[1] * uH[1])
        kH, kN = kval(c, dH), kval(c, dN)
        assert kH < 1 <= kN and all(t > 0 for t in cross(c, dN))
        rr = list(rhoH * np.linspace(0.5, 2.0, 16)) + list(np.logspace(-4, -1, 4))
        res = scan(case, c, dN, rr, list(np.logspace(-3, 2, 26)), pool)
        mx = max((dl / r for lst in res.values() for (r, dl, st) in lst if dl is not None), default=-np.inf)
        okr &= mx < 0
        print(f'   case {case} c={c}: Hopf d={dH} (k={kH:.3g}, rho_H={rhoH:.4g}) -> non-Hopf d={dN} (k={kN:.3g}): '
              f'max delta/r = {mx:.3e}')
    print(('PASS' if okr else 'FAIL') + ' R3 no cycle in the non-Hopf neighbours near the Hopf value of their Hopf twin')
    fails += (not okr)
    pool.close()

    print('R4 (H2) by quadrature (mpmath, 30 digits)')
    mpm.mp.dps = 30
    random.seed(99)
    okh = True
    nchk = 0
    for name, case, c, d in family_members():
        u = cross(c, d)
        al1, al2, al3 = c[0] * u[0], u[1], c[2] * u[2]
        be1, be2, be3 = d[0] * u[0], -d[1] * u[1], d[2] * u[2]
        if case == 9:
            Af = lambda x: al1 * x ** 2; Bf = lambda x: al3 - al2 * x
            Cf = lambda x: be1 * x ** 2; Df = lambda x: be3 - be2 * x
            x0 = mpm.mpf(al3) / al2
        else:
            Af = lambda x: al1 * x ** 2 + al3; Bf = lambda x: -al2 * x
            Cf = lambda x: be1 * x ** 2 + be3; Df = lambda x: -be2 * x
            x0 = mpm.mpf(0)
        Nf = lambda x: -Af(x) / Bf(x)
        for rho in (mpm.mpf('0.001'), mpm.mpf('0.3'), mpm.mpf(5)):
            # F = int_1^x f with f = N' - rho D/B, i.e. F = N(x) - N(1) - rho int_1^x D/B
            g = lambda x: -rho * (Bf(x) * Cf(x) - Af(x) * Df(x)) / Bf(x) ** 2
            F = lambda x: Nf(x) - Nf(mpm.mpf(1)) - rho * mpm.quad(lambda s: Df(s) / Bf(s), [1, x])
            G = lambda x: mpm.quad(g, [1, x])
            for frac in (mpm.mpf('0.02'), mpm.mpf('0.3'), mpm.mpf('0.7'), mpm.mpf('0.97')):
                a = x0 + (1 - x0) * frac
                Ga = G(a)
                hi = mpm.mpf(2)
                while G(hi) < Ga:
                    hi *= 2
                b = mpm.findroot(lambda t: G(t) - Ga, (mpm.mpf(1) + mpm.mpf('1e-20'), hi), solver='illinois')
                diff = F(a) - F(b)
                nchk += 1
                if not diff > 0:
                    okh = False
                    print(f'   FAIL {name} rho={rho} a={mpm.nstr(a, 8)} b={mpm.nstr(b, 8)} F(a)-F(b)={mpm.nstr(diff, 8)}')
    print(('PASS' if okh else 'FAIL') + f' R4 F(a) > F(b) at {nchk} G-paired points (3 rho values, 4 a values, all members)')
    fails += (not okh)
    print(f'runtime {time.time() - t0:.0f} s')
    print('ALL PASS' if fails == 0 else f'{fails} FAIL')
    return fails


if __name__ == '__main__':
    sys.exit(1 if main() else 0)
