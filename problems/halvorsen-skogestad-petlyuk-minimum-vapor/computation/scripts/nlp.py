"""Minimum total reboiler vapor of a given configuration (TG / Nallasivam model), solved with ipopt via pyomo.

Environment: WSL Fedora, Python 3, pyomo 6.10, ipopt 3.14 (/usr/sbin/ipopt).  Floating point: exploration only.
Runtime: 0.03-0.2 s per local solve (N = 3..5).

Model for configuration c (see hsmodel.py for the configuration space), feed flows F_p, liquid fraction q:
  For each mixture stream s = [i,j] with split [i,k]/[l,j]:
     d[s,p] >= 0 (p = i..k), b[s,p] >= 0 (p = l..j): net component flows up (top section) / down (bottom)
     VT[s], VB[s] >= 0: vapor in top / bottom section
     th[s,r], r = l-1..k, alpha_{r+1} < th < alpha_r: common Underwood roots of the column's feed
  feed of s:  f_p = F_p for s = [1,N]; else sum of d[t,p] (top producer t) and b[u,p] (bottom producer u)
  vapor feed: vf = (1-q) sum F for [1,N]; else  top part - bottom part, where
     top part    = VT[t] (thermal coupling or side draw)  or  sum_p d[t,p] (condenser: saturated-vapor transfer)
     bottom part = VB[u] (thermal coupling or side draw)  or  0          (reboiler: saturated-liquid transfer)
  balances: f_p = d_p [p<=k] + b_p [p>=l];   VT - VB = vf
  Underwood: sum_p alpha_p f_p/(alpha_p - th_r) = vf  (all r)
             VT >= sum_{p<=k} alpha_p d_p/(alpha_p - th_r)  (r = l-1, k), with equality for l <= r <= k-1
  pure product p produced both as distillate of t and as residue of u (liquid side draw): VT[t] = VB[u]
  objective: sum of VB[u] over reboilers (submixture HX residues and pure products produced only as residues)
"""
from __future__ import annotations

import math
import random

import pyomo.environ as pyo
from pyomo.opt import SolverFactory, TerminationCondition

from hsmodel import split_of  # noqa: F401  (re-export convenience)
from nlp_structure import structure  # noqa: F401


def build_model(cfg, alpha, F, q, eps=1e-9, permissive=False):
    """permissive=False: strict infinite-stage model (VT = Underwood value at interior roots, Halvorsen 3.2.4);
    permissive=True: TG/Nallasivam literal model (Upsilon = value at interior roots >= value at outer roots,
    VT >= Upsilon)."""
    N = cfg['N']
    al = {p: float(alpha[p - 1]) for p in range(1, N + 1)}
    Fd = {p: float(F[p - 1]) for p in range(1, N + 1)}
    cols, info, reb, side = structure(cfg)
    m = pyo.ConcreteModel()
    Ftot = sum(Fd.values())
    m.d = pyo.Var([(s, p) for s in cols for p in range(info[s]['i'], info[s]['k'] + 1)], bounds=(0, Ftot))
    m.b = pyo.Var([(s, p) for s in cols for p in range(info[s]['l'], info[s]['j'] + 1)], bounds=(0, Ftot))
    m.VT = pyo.Var(cols, bounds=(0, 1e3 * Ftot))
    m.VB = pyo.Var(cols, bounds=(0, 1e3 * Ftot))
    th_idx = [(s, r) for s in cols for r in range(info[s]['l'] - 1, info[s]['k'] + 1)]

    def thb(mm, i, j, r):
        w = al[r] - al[r + 1]
        return (al[r + 1] + eps * w, al[r] - eps * w)
    m.th = pyo.Var(th_idx, bounds=thb)

    def feed(s, p):
        if s == (1, N):
            return Fd[p]
        d = info[s]
        e = 0
        if d['top'] is not None:
            e = e + m.d[d['top'], p]
        if d['bot'] is not None:
            e = e + m.b[d['bot'], p]
        return e

    def vfeed(s):
        if s == (1, N):
            return (1 - q) * Ftot
        d = info[s]
        e = 0
        if d['top'] is not None:
            t = d['top']
            if d['kind'] == 'HX':
                e = e + sum(m.d[t, p] for p in range(info[t]['i'], info[t]['k'] + 1))
            else:
                e = e + m.VT[t]
        if d['bot'] is not None:
            u = d['bot']
            if d['kind'] != 'HX':
                e = e - m.VB[u]
        return e

    m.cons = pyo.ConstraintList()
    for s in cols:
        d = info[s]
        for p in range(d['i'], d['j'] + 1):
            rhs = 0
            if p <= d['k']:
                rhs = rhs + m.d[s, p]
            if p >= d['l']:
                rhs = rhs + m.b[s, p]
            m.cons.add(feed(s, p) == rhs)
        m.cons.add(m.VT[s] - m.VB[s] == vfeed(s))
        has_int = d['l'] <= d['k'] - 1
        if permissive and has_int:
            ups = pyo.Var(bounds=(0, 1e3 * Ftot))
            m.add_component(f"ups_{s[0]}_{s[1]}", ups)
            m.cons.add(m.VT[s] >= ups)
        for r in range(d['l'] - 1, d['k'] + 1):
            m.cons.add(sum(al[p] * feed(s, p) / (al[p] - m.th[s, r]) for p in range(d['i'], d['j'] + 1))
                       == vfeed(s))
            g = sum(al[p] * m.d[s, p] / (al[p] - m.th[s, r]) for p in range(d['i'], d['k'] + 1))
            tgt = ups if (permissive and has_int) else m.VT[s]
            if d['l'] <= r <= d['k'] - 1:
                m.cons.add(tgt == g)
            else:
                m.cons.add(tgt >= g)
    for (t, u) in side:
        m.cons.add(m.VT[t] == m.VB[u])
    m.obj = pyo.Objective(expr=sum(m.VB[u] for u in reb), sense=pyo.minimize)
    m._info = (cols, info, reb, side, al, Fd, q)
    return m


def _roots_of(al_list, f_list, rhs, idx_pairs):
    """Roots of sum al f/(al - th) = rhs in given intervals (float bisection); returns dict r->root or None."""
    out = {}
    for (r, lo, hi) in idx_pairs:
        a, b = lo, hi
        fa = None
        for _ in range(200):
            mid = 0.5 * (a + b)
            v = sum(x * y / (x - mid) for x, y in zip(al_list, f_list)) - rhs
            if v > 0:
                b = mid
            else:
                a = mid
        out[r] = 0.5 * (a + b)
    return out


def random_start(m, rng, vscale=1.0):
    """Heuristic initial point: random component splits, roots from feed eq., vapors from Underwood + margin.
    Processes columns from the longest stream down (producers before consumers)."""
    cols, info, reb, side, al, Fd, q = m._info
    N = max(al)
    order = sorted(cols, key=lambda s: -(s[1] - s[0]))
    val = {}
    Ntot = sum(Fd.values())
    for s in order:
        d = info[s]
        if s == (1, max(p for p in al)):
            f = {p: Fd[p] for p in range(d['i'], d['j'] + 1)}
            vf = (1 - q) * Ntot
        else:
            f = {}
            for p in range(d['i'], d['j'] + 1):
                v = 0.0
                if d['top'] is not None:
                    v += val[('d', d['top'], p)] if ('d', d['top'], p) in val else 0.0
                if d['bot'] is not None:
                    v += val[('b', d['bot'], p)] if ('b', d['bot'], p) in val else 0.0
                f[p] = v
            vf = 0.0
            if d['top'] is not None:
                t = d['top']
                if d['kind'] == 'HX':
                    vf += sum(val[('d', t, p)] for p in range(info[t]['i'], info[t]['k'] + 1))
                else:
                    vf += val[('VT', t)]
            if d['bot'] is not None and d['kind'] != 'HX':
                vf -= val[('VB', d['bot'])]
        # split
        for p in range(d['i'], d['j'] + 1):
            if p < d['l']:
                val[('d', s, p)] = f[p]
            elif p > d['k']:
                val[('b', s, p)] = f[p]
            else:
                x = rng.random()
                val[('d', s, p)] = x * f[p]
                val[('b', s, p)] = (1 - x) * f[p]
        alist = [al[p] for p in range(d['i'], d['j'] + 1)]
        flist = [max(f[p], 1e-12) for p in range(d['i'], d['j'] + 1)]
        pairs = [(r, al[r + 1], al[r]) for r in range(d['l'] - 1, d['k'] + 1)]
        rts = _roots_of(alist, flist, vf, pairs)
        req = 0.0
        for r, th in rts.items():
            val[('th', s, r)] = th
            g = sum(al[p] * val[('d', s, p)] / (al[p] - th) for p in range(d['i'], d['k'] + 1))
            req = max(req, g)
        VT = max(req, 0.0) * (1 + vscale * rng.random()) + 1e-6
        VB = VT - vf
        if VB < 0:
            VT -= VB - 1e-6
            VB = 1e-6
        val[('VT', s)] = VT
        val[('VB', s)] = VB
    for key, v in val.items():
        if key[0] == 'd':
            m.d[key[1], key[2]].value = v
        elif key[0] == 'b':
            m.b[key[1], key[2]].value = v
        elif key[0] == 'th':
            lo, hi = m.th[key[1], key[2]].bounds
            m.th[key[1], key[2]].value = min(max(v, lo), hi)
        elif key[0] == 'VT':
            m.VT[key[1]].value = v
        elif key[0] == 'VB':
            m.VB[key[1]].value = v


_SOLVER = None


def solver():
    global _SOLVER
    if _SOLVER is None:
        _SOLVER = SolverFactory('ipopt', executable='/usr/sbin/ipopt')
        _SOLVER.options['tol'] = 1e-10
        _SOLVER.options['acceptable_tol'] = 1e-8
        _SOLVER.options['max_iter'] = 3000
        _SOLVER.options['print_level'] = 0
        _SOLVER.options['mu_strategy'] = 'adaptive'
        _SOLVER.options['constr_viol_tol'] = 1e-10
        _SOLVER.options['bound_relax_factor'] = 1e-10
    return _SOLVER


def max_violation(m):
    v = 0.0
    for c in m.component_data_objects(pyo.Constraint, active=True):
        body = pyo.value(c.body)
        lo = pyo.value(c.lower) if c.has_lb() else None
        hi = pyo.value(c.upper) if c.has_ub() else None
        if lo is not None:
            v = max(v, lo - body)
        if hi is not None:
            v = max(v, body - hi)
    return v


def solve_config(cfg, alpha, F, q, starts=10, seed=0, vscale=1.0, eps=1e-9, tee=False, permissive=False):
    """Multi-start local minimization; returns (best objective, best point dict, list of all results)."""
    m = build_model(cfg, alpha, F, q, eps=eps, permissive=permissive)
    rng = random.Random(seed)
    best = None
    results = []
    for st in range(starts):
        random_start(m, rng, vscale)
        try:
            res = solver().solve(m, tee=tee, load_solutions=False)
        except Exception as ex:  # noqa: BLE001
            results.append(('error', str(ex)))
            continue
        tc = res.solver.termination_condition
        if tc in (TerminationCondition.optimal, TerminationCondition.locallyOptimal,
                  TerminationCondition.feasible) or str(tc) == 'acceptable':
            m.solutions.load_from(res)
            viol = max_violation(m)
            obj = pyo.value(m.obj)
            results.append((obj, viol, str(tc)))
            if viol < 1e-6 and (best is None or obj < best[0]):
                snap = snapshot(m)
                vv, obj2 = verify_solution(cfg, alpha, F, q, snap, permissive=permissive)
                if vv < 1e-7:
                    best = (obj, snap)
                else:
                    results[-1] = (obj, viol, str(tc) + f' REJECTED(verify {vv:.1e})')
        else:
            results.append(('fail', str(tc)))
    return best, results


def verify_solution(cfg, alpha, F, q, sol, permissive=False):
    """Independent a-posteriori check of a solution snapshot: recompute every column's feed (from the producers'
    flows and vapors), its Underwood roots by bisection in each interval, and all constraints, with RELATIVE
    tolerances.  Returns (max relative violation, recomputed objective)."""
    N = cfg['N']
    cols, info, reb, side = structure(cfg)
    al = {p: float(alpha[p - 1]) for p in range(1, N + 1)}
    Fd = {p: float(F[p - 1]) for p in range(1, N + 1)}
    get = lambda dct, *k: dct.get(str(k), 0.0) if len(k) > 1 else dct[str(k[0])]
    d = {eval(k): v for k, v in sol['d'].items()}
    b = {eval(k): v for k, v in sol['b'].items()}
    VT = {eval(k): v for k, v in sol['VT'].items()}
    VB = {eval(k): v for k, v in sol['VB'].items()}
    dd = lambda s, p: d.get((s[0], s[1], p), 0.0)
    bb = lambda s, p: b.get((s[0], s[1], p), 0.0)
    scale = max(1.0, max(VT.values()))
    viol = 0.0
    for s in cols:
        I = info[s]
        if s == (1, N):
            f = {p: Fd[p] for p in range(I['i'], I['j'] + 1)}
            vf = (1 - q) * sum(Fd.values())
        else:
            f = {p: 0.0 for p in range(I['i'], I['j'] + 1)}
            vf = 0.0
            if I['top'] is not None:
                t = I['top']
                for p in range(info[t]['i'], info[t]['k'] + 1):
                    f[p] += dd(t, p)
                vf += sum(dd(t, p) for p in range(info[t]['i'], info[t]['k'] + 1)) if I['kind'] == 'HX' else VT[t]
            if I['bot'] is not None:
                u = I['bot']
                for p in range(info[u]['l'], info[u]['j'] + 1):
                    f[p] += bb(u, p)
                if I['kind'] != 'HX':
                    vf -= VB[u]
        for p in range(I['i'], I['j'] + 1):
            rhs = (dd(s, p) if p <= I['k'] else 0.0) + (bb(s, p) if p >= I['l'] else 0.0)
            viol = max(viol, abs(f[p] - rhs) / max(1e-12, Fd[p]))
            if p <= I['k']:
                viol = max(viol, -dd(s, p) / Fd[p])
            if p >= I['l']:
                viol = max(viol, -bb(s, p) / Fd[p])
        viol = max(viol, abs(VT[s] - VB[s] - vf) / scale, -VT[s] / scale, -VB[s] / scale)
        gvals = {}
        for r in range(I['l'] - 1, I['k'] + 1):
            lo, hi = al[r + 1], al[r]
            fpos = {p: f[p] for p in f if f[p] > 0}
            fun = lambda th: sum(al[p] * fpos[p] / (al[p] - th) for p in fpos) - vf
            a, c_ = lo + 1e-10 * (hi - lo), hi - 1e-10 * (hi - lo)
            try:
                fa, fc = fun(a), fun(c_)
            except ZeroDivisionError:
                return float('inf'), None
            if not (fa < 0 < fc):
                return float('inf'), None   # no root in the interval (or root within 1e-10 of an end: rejected)
            for _ in range(200):
                mid = 0.5 * (a + c_)
                if fun(mid) > 0:
                    c_ = mid
                else:
                    a = mid
            th = 0.5 * (a + c_)
            gvals[r] = sum(al[p] * dd(s, p) / (al[p] - th) for p in range(I['i'], I['k'] + 1) if dd(s, p) > 0)
        interior = [r for r in gvals if I['l'] <= r <= I['k'] - 1]
        outer = [r for r in gvals if r not in interior]
        if interior and permissive:
            ups = max(gvals[r] for r in interior)
            for r in interior:
                viol = max(viol, abs(gvals[r] - ups) / scale)
            for r in outer:
                viol = max(viol, (gvals[r] - ups) / scale)
            viol = max(viol, (ups - VT[s]) / scale)
        else:
            for r in interior:
                viol = max(viol, abs(VT[s] - gvals[r]) / scale)
            for r in outer:
                viol = max(viol, (gvals[r] - VT[s]) / scale)
    for (t, u) in side:
        viol = max(viol, abs(VT[t] - VB[u]) / scale)
    return viol, sum(VB[u] for u in reb)


def snapshot(m):
    return {'d': {str(k): pyo.value(v) for k, v in m.d.items()},
            'b': {str(k): pyo.value(v) for k, v in m.b.items()},
            'VT': {str(k): pyo.value(v) for k, v in m.VT.items()},
            'VB': {str(k): pyo.value(v) for k, v in m.VB.items()},
            'th': {str(k): pyo.value(v) for k, v in m.th.items()}}
