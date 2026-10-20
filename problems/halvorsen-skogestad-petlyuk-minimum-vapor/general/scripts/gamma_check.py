"""Numerical stress test of the general-N lemmas of general/general-N.md (Agent G).  Floating point, not a
certificate: it tests the inequalities on optimal and random feasible points of the PERMISSIVE (TG literal)
model and of the STRICT model (computation/scripts/nlp.py, used read-only).

Checked at every point, for every split j and every mixed column c of M_j (theta in I_j = (a_{j+1}, a_j)):
  (G1)  recovery ordering d_p/f_p >= d_q/f_q for distributing p < q (every column, not only mixed ones);
  (SR)  sign rule: not (T_c(theta_j) < 0 and B_c(theta_j) > 0);
  (LOC) local Gamma inequality (Lemma G4 / proof of Theorem G3):
            l_c + kappa_out + sum_{mixed out-links e} Gamma(y_e) >= Gamma(a_c),
        where Gamma(y) = G_y(zeta_y) - G_y(theta_j), zeta_y = zero of y in I_j, G_y = light part;
  (MIN) local min inequality (proof of Theorem G2): l_c + kappa_out + sum_e min(0, y_e(theta_j)) >= min(0, a_c(theta_j))
        reported separately for type-(d) columns with excess (expected to fail there in PERMISSIVE);
  (TREE) R - H_j >= 0 when M_j has no merge node (Theorem G3), and R - H_j >= 0 in general (conjecture);
  merge defects  Gamma(y') + Gamma(y'') - Gamma(y'+y'') at merge nodes are recorded (statistics only).
Environment: WSL Fedora, Python 3, pyomo + ipopt.  Runtime: ~1-3 min for 100 trials.
Usage: python gamma_check.py N ntrials seed permissive(0/1) [config ids...]
Prints PASS if (G1), (SR), (LOC), (TREE) hold within tolerance at all tested points, FAIL otherwise.
"""
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'computation', 'scripts'))

import pyomo.environ as pyo  # noqa: E402
from pyomo.opt import TerminationCondition  # noqa: E402

from hsmodel import all_configurations, config_name, feed_roots, petlyuk_peaks  # noqa: E402
from nlp_g import build_model, random_start, solver, max_violation, snapshot, solve_config, structure  # noqa: E402


TOL = 1e-6


def random_feasible(cfg, alpha, F, q, rng, permissive):
    m = build_model(cfg, alpha, F, q, permissive=permissive)
    cols = list(m.VT.keys())
    m.obj.deactivate()
    w = {c: rng.uniform(-0.2, 1.0) for c in cols}
    w2 = {c: rng.random() for c in cols}
    wd = {k: rng.uniform(-1, 1) for k in m.d.keys()}
    m.obj2 = pyo.Objective(expr=sum(w[c] * m.VB[c] + w2[c] * m.VT[c] for c in cols)
                           + sum(wd[k] * m.d[k] for k in m.d.keys()))
    random_start(m, rng, 2.0)
    try:
        res = solver().solve(m, load_solutions=False)
    except Exception:  # noqa: BLE001
        return None
    if res.solver.termination_condition != TerminationCondition.optimal:
        return None
    m.solutions.load_from(res)
    if max_violation(m) > 1e-8:
        return None
    return snapshot(m)


class Point:
    def __init__(self, cfg, sol, alpha, F, q):
        self.cfg, self.alpha, self.F, self.q = cfg, alpha, F, q
        self.N = cfg['N']
        self.cols, self.info, self.reb, self.side = structure(cfg)
        self.al = {p: alpha[p - 1] for p in range(1, self.N + 1)}
        self.d = {eval(k): v for k, v in sol['d'].items()}
        self.b = {eval(k): v for k, v in sol['b'].items()}
        self.VT = {eval(k): v for k, v in sol['VT'].items()}
        self.VB = {eval(k): v for k, v in sol['VB'].items()}
        self.R = sum(self.VB[u] for u in self.reb)

    def dd(self, s, p):
        return max(self.d.get((s[0], s[1], p), 0.0), 0.0)

    def bb(self, s, p):
        return max(self.b.get((s[0], s[1], p), 0.0), 0.0)

    def feed(self, s):
        """(U, {p: f_p}) of column s."""
        N = self.N
        I = self.info[s]
        if s == (1, N):
            return (1 - self.q) * sum(self.F), {p: self.F[p - 1] for p in range(1, N + 1)}
        f = {p: 0.0 for p in range(I['i'], I['j'] + 1)}
        U = 0.0
        if I['top'] is not None:
            t = I['top']
            for p in f:
                f[p] += self.dd(t, p)
            if I['kind'] == 'HX':
                U += sum(self.dd(t, p) for p in range(self.info[t]['i'], self.info[t]['k'] + 1))
            else:
                U += self.VT[t]
        if I['bot'] is not None:
            u = I['bot']
            for p in f:
                f[p] += self.bb(u, p)
            if I['kind'] != 'HX':
                U -= self.VB[u]
        return U, f


def fun(U, w, al):
    """y(theta) = U - sum_p a_p w_p/(a_p - theta) as a closure, with its light-part G(theta; j)."""
    return lambda th: U - sum(al[p] * x / (al[p] - th) for p, x in w.items())


def zero_in(y, lo, hi):
    """Zero of a non-increasing function y on (lo, hi); returns lo/hi clamp if no sign change."""
    a, b = lo + 1e-14 * (hi - lo), hi - 1e-14 * (hi - lo)
    ya, yb = y(a), y(b)
    if ya <= 0:
        return a
    if yb >= 0:
        return b
    for _ in range(200):
        m = 0.5 * (a + b)
        if y(m) > 0:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def Glight(w, al, j):
    return lambda th: sum(al[p] * x / (al[p] - th) for p, x in w.items() if p <= j)


def analyse(P, j, thj):
    """Return dict of checks for split j."""
    N, al, info = P.N, P.al, P.info
    lo, hi = al[j + 1], al[j]
    M = [s for s in P.cols if s[0] <= j < s[1]]
    out = {'loc_min': 1e9, 'min_nonD': 1e9, 'min_D': 1e9, 'sr_viol': 0, 'merge_def': [], 'tree': True}
    # per-column functions
    data = {}
    for c in M:
        I = info[c]
        U, f = P.feed(c)
        top = {p: P.dd(c, p) for p in range(I['i'], I['k'] + 1)}
        bot = {p: P.bb(c, p) for p in range(I['l'], I['j'] + 1)}
        T = fun(P.VT[c], top, al)
        B = lambda th, c=c, bot=bot: -P.VB[c] - sum(al[p] * x / (al[p] - th) for p, x in bot.items())
        a = fun(U, f, al)
        data[c] = dict(U=U, f=f, top=top, bot=bot, T=T, B=B, a=a)
    def Gamma_of(y, w):
        z = zero_in(y, lo, hi)
        G = Glight(w, al, j)
        return G(z) - G(thj), z
    # links into mixed children
    links = {}  # child -> list of (parent, kind, y, w, kappa)
    for c in M:
        I = info[c]
        dist, res = (I['i'], I['k']), (I['l'], I['j'])
        D = data[c]
        if dist in data:
            kap = 0.0
            if info[dist]['kind'] == 'HX':
                kap = P.VT[c] - sum(D['top'].values())
            y = (lambda th, T=D['T'], kap=kap: T(th) - kap)
            links.setdefault(dist, []).append((c, 'dist', y, D['top'], kap))
        if res in data:
            rho = P.VB[c] if info[res]['kind'] == 'HX' else 0.0
            y = (lambda th, B=D['B'], rho=rho: B(th) + rho)
            links.setdefault(res, []).append((c, 'res', y, D['bot'], 0.0))
    for c in M:
        I = info[c]
        D = data[c]
        dist, res = (I['i'], I['k']), (I['l'], I['j'])
        t, bval, aval = D['T'](thj), D['B'](thj), D['a'](thj)
        if t < -TOL and bval > TOL:
            out['sr_viol'] += 1
        ell = t if I['k'] <= j else 0.0
        gam_a, _ = Gamma_of(D['a'], D['f'])
        lhs_g = ell
        lhs_m = ell
        for child in (dist, res):
            for (par, kind, y, w, kap) in links.get(child, []):
                if par != c:
                    continue
                g, _ = Gamma_of(y, w)
                lhs_g += g + kap
                lhs_m += min(0.0, y(thj)) + kap
        out['loc_min'] = min(out['loc_min'], lhs_g - gam_a)
        typ_d = I['k'] >= j + 1 and I['l'] <= j
        mm = lhs_m - min(0.0, aval)
        if typ_d:
            out['min_D'] = min(out['min_D'], mm)
        else:
            out['min_nonD'] = min(out['min_nonD'], mm)
    for s, L in links.items():
        if len(L) == 2:
            out['tree'] = False
            gsum = 0.0
            zs = []
            for (par, kind, y, w, kap) in L:
                g, z = Gamma_of(y, w)
                gsum += g
                zs.append((kind, z))
            g_s, z_s = Gamma_of(data[s]['a'], data[s]['f'])
            out['merge_def'].append((gsum - g_s, dict(zs), z_s))
    return out


def recovery_ok(P):
    worst = 1e9
    for c in P.cols:
        I = P.info[c]
        U, f = P.feed(c)
        r = [P.dd(c, p) / f[p] for p in range(I['l'], I['k'] + 1) if f[p] > 1e-9]
        for x, y in zip(r, r[1:]):
            worst = min(worst, x - y)
    return worst


def main():
    N = int(sys.argv[1])
    ntr = int(sys.argv[2])
    rng = random.Random(int(sys.argv[3]))
    permissive = bool(int(sys.argv[4]))
    ids = [int(x) for x in sys.argv[5:]]
    C = all_configurations(N)
    stats = dict(points=0, g1=1e9, sr=0, loc=1e9, minD=1e9, minnonD=1e9, tree=1e9, gap=1e9,
                 mdef_pos=0, mdef_n=0, mdef_max=-1e9)
    for tr in range(ntr):
        c = C[rng.choice(ids)] if ids else rng.choice(C)
        rat = [rng.uniform(1.05, 4) for _ in range(N - 1)]
        al = [1.0]
        for r in reversed(rat):
            al.insert(0, al[0] * r)
        wz = [rng.random() + 0.05 for _ in range(N)]
        z = [x / sum(wz) for x in wz]
        q = rng.choice([1.0, 0.0, rng.uniform(-0.3, 1.3)])
        ths = [float(x) for x in feed_roots(al, z, q)]
        Pk = [float(x) for x in petlyuk_peaks(al, z, q)]
        for mode in ('opt', 'rand'):
            if mode == 'opt':
                best, _ = solve_config(c, al, z, q, starts=3, seed=tr, permissive=permissive)
                sol = best[1] if best else None
            else:
                sol = random_feasible(c, al, z, q, rng, permissive)
            if sol is None:
                continue
            P = Point(c, sol, al, z, q)
            stats['points'] += 1
            stats['g1'] = min(stats['g1'], recovery_ok(P))
            for j in range(1, N):
                o = analyse(P, j, ths[j - 1])
                gap = P.R - (Pk[j - 1] - (1 - q))
                stats['gap'] = min(stats['gap'], gap)
                stats['sr'] += o['sr_viol']
                stats['loc'] = min(stats['loc'], o['loc_min'])
                stats['minD'] = min(stats['minD'], o['min_D'])
                stats['minnonD'] = min(stats['minnonD'], o['min_nonD'])
                if o['tree']:
                    stats['tree'] = min(stats['tree'], gap)
                for (dv, zs, zsn) in o['merge_def']:
                    stats['mdef_n'] += 1
                    stats['mdef_max'] = max(stats['mdef_max'], dv)
                    if dv > 1e-7:
                        stats['mdef_pos'] += 1
                if o['loc_min'] < -1e-5 or o['min_nonD'] < -1e-5 or gap < -1e-6:
                    print('ALERT', mode, c['id'], config_name(c), 'j', j, 'alpha', al, 'z', z, 'q', q,
                          'loc', o['loc_min'], 'minnonD', o['min_nonD'], 'gap', gap, flush=True)
    for k, v in stats.items():
        print(f'  {k}: {v}')
    ok = (stats['g1'] > -1e-6 and stats['sr'] == 0 and stats['loc'] > -1e-5 and stats['minnonD'] > -1e-5
          and stats['tree'] > -1e-6)
    print('PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
