#!/usr/bin/env python3
"""paper_tables.py -- Agent C.  Reproducible high-precision demonstration for the numerical section of the paper.

Table 1  N = 3 (alpha = (4,2,1), z = 1/3, q = 1, delta = 0.05, prefractionator slack 0.025, pattern BT with
         eta_1 = eta_2 = 0.003; the data of network.md 11.3): for n = 16 .. 256 stages in every section,
         product impurities, junction compositions versus the limit configuration, |v_n - v_inf|, the trace net flows of
         the prefractionator, sigma_min(DG_n), and the stage-level verification (max relative residual of all stage
         and unit balances, min mole fraction, min S).
Table 2  N = 4 (alpha = (8,4,2,1), z = 1/4, q = 1, delta = 0.05, eta_nb = 0.005, pattern BBB, eta_p = 0.003; D's
         preset n4_1): same quantities.
Table 3  approach to the Petlyuk minimum (N = 3 data of Table 1): (delta, eta) -> 0 jointly; R - H* = delta exactly,
         at a stage number large enough for the limit regime: impurities of all products.
Writes runs/paper_tables.md (Markdown tables) and runs/paper_tables.json; prints PASS/FAIL (every row must be a
converged, physical steady state).  Numerical evidence only (high-precision floating point, not interval arithmetic).
Usage: python paper_tables.py
Environment: Python 3.12, mpmath 1.3, numpy, python-flint (WSL Fedora sci env).  Runtime: about 5 minutes.
"""
import os
import sys
import json
import time
import math
import mpmath as mp
from petnet import Problem, Design, limit_config, System, svals, impurities, solve_audited

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, '..', 'runs')
t0 = time.time()
fails = []
md = []
js = {}


def e_(x, d=2):
    x = float(x)
    if x == 0:
        return "0"
    return f"{x:.{d}e}".replace("e-0", "e-").replace("e+0", "e+")


def run_table(title, build, nlist, junctions, pre_traces=None):
    md.append(f"### {title}\n")
    mp.mp.dps = 60
    pb, des, lim = build()
    N = pb.N
    md.append(f"R = H* + δ = {mp.nstr(des.R, 10)} (H* = {mp.nstr(pb.Hmax, 10)}); product flows "
              f"{[mp.nstr(x, 6) for x in des.P]}; limit product impurities "
              f"{[e_(1 - des.xprod[p][p], 4) for p in range(N)]} (designed: O(η)).\n")
    sl = System(des, lim, None)
    sv = svals(sl.jacobian(sl.start()))[-1]
    md.append(f"σ_min(DG_∞(v_∞)) = {mp.nstr(sv, 5)}; limit junction compositions: " +
              "; ".join(f"J[{u[0] + 1},{u[1] + 1}] = ({', '.join(mp.nstr(x, 8) for x in lim['X'][u])})" for u in junctions) + "\n")
    head = ("| n | impurity of products 1..N | " + " | ".join(f"‖x_J[{u[0]+1},{u[1]+1}] − limit‖∞" for u in junctions)
            + " | ‖v_n − v_∞‖ (majors) | " + ("prefractionator trace flows | " if pre_traces else "")
            + "σ_min(DG_n) | max rel. stage residual | min mole fraction | min S | digits |")
    md.append(head)
    md.append("|" + "---|" * (head.count("|") - 1))
    rows = []
    vprev = None
    for n in sorted(nlist, reverse=True):
        dps = int(40 + 1.25 * (n + 1) * math.log10(float(max(max(r) / min(r) for r in lim['roots'].values()))))
        res = solve_audited(build, n, max(60, dps), start=None)
        if not res['phys'] and vprev is not None:
            res = solve_audited(build, n, max(60, dps), start=vprev)
        sn, v = res['sn'], res['v']
        X, W = sn.unpack(v)
        Xl = res['lim']['X']
        vl = sn.start()
        dmaj = max(abs(v[j] - vl[j]) for j in range(len(v)) if vl[j] != 0)
        dj = [max(abs(a - b) for a, b in zip(X[u], Xl[u])) for u in junctions]
        smin = svals(sn.jacobian(v))[-1]
        imp = impurities(sn, v)
        chk = res['chk']
        worst = max(chk['worst_rel'], chk['unit_rel'])
        pt = ""
        if pre_traces:
            pt = ", ".join(f"{lab} {e_(f(W))}" for lab, f in pre_traces) + " | "
        row = dict(n=n, impurities=[str(x) for x in imp], junction_dist=[float(x) for x in dj], dist=float(dmaj),
                   sigma_min=float(smin), worst_rel=float(worst), minx=str(chk['minx']), minS=float(chk['minS']),
                   dps=res['dps'], physical=res['phys'])
        rows.append((n, f"| {n} | {', '.join(e_(x) for x in imp)} | " + " | ".join(e_(x) for x in dj) +
                     f" | {e_(dmaj)} | {pt}{mp.nstr(smin, 5)} | {e_(worst, 1)} | {e_(chk['minx'], 1)} | "
                     f"{mp.nstr(chk['minS'], 4)} | {res['dps']} |", row))
        print(f"  {title[:20]} n = {n}: physical {res['phys']}, impurities {[e_(x) for x in imp]}, dist {e_(dmaj)}, "
              f"sigma_min {mp.nstr(smin, 5)}, worst rel {e_(worst, 1)}", flush=True)
        if not res['phys']:
            fails.append(f"{title} n={n}")
        else:
            vprev = [str(x) for x in v]
        mp.mp.dps = 60
    rows.sort(key=lambda t: t[0])
    md.extend(r[1] for r in rows)
    md.append("")
    # rate
    a, b = rows[-2][2], rows[-1][2]
    rate = math.exp((math.log(b['dist']) - math.log(a['dist'])) / (b['n'] - a['n']))
    rho = max(max(float(r[k] / r[j]) for j in range(k)) if k > 0 else 0 for r, k in
              [(lim['roots'][s], lim['kstar'][s]) for s in pb.secs])
    rho2 = max(max([float(r[j] / r[k]) for j in range(k + 1, N)] + [0]) for r, k in
               [(lim['roots'][s], lim['kstar'][s]) for s in pb.secs])
    md.append(f"Observed rate ‖v_n − v_∞‖ ∝ {rate:.4f}^n; predicted ρ_∞ = max adjacent mode ratio at the pinch = "
              f"{max(rho, rho2):.4f}.\n")
    return [r[2] for r in rows]


case3 = dict(alpha=[4, 2, 1], z=[1, 1, 1], q=1)


def build3(delta='0.05', eta_pre='0.025', eta='0.003'):
    def b():
        pb = Problem([4, 2, 1], [1, 1, 1], 1)
        des = Design(pb, mp.mpf(delta), eta_nb_default=mp.mpf(eta_pre),
                     eta_bin={(0, 1): ('B', mp.mpf(eta)), (1, 2): ('T', mp.mpf(eta))})
        return pb, des, limit_config(des)
    return b


def build4():
    pb = Problem([8, 4, 2, 1], [1, 1, 1, 1], 1)
    des = Design(pb, mp.mpf('0.05'), eta_nb_default=mp.mpf('0.005'),
                 eta_bin={(0, 1): ('B', mp.mpf('0.003')), (1, 2): ('B', mp.mpf('0.003')), (2, 3): ('B', mp.mpf('0.003'))})
    return pb, des, limit_config(des)


md.append("## Numerical demonstration (Agent C): finite-stage steady states of the forced generalized Petlyuk "
          "arrangement at fixed flows\n")
md.append("All numbers below are numerical evidence (high-precision floating point, mpmath), not a proof. Every row is a "
          "steady state of the full stage-by-stage model (every stage and unit component balance re-checked with the "
          "actual stream compositions; all mole fractions > 0; S > 0 everywhere).\n")
js['N3'] = run_table("Table 1. N = 3, α = (4, 2, 1), z = (1/3, 1/3, 1/3), q = 1, δ = 0.05, η_[1,3] = 0.025, "
                     "pattern (B, T), η_1 = η_2 = 0.003", build3(), [16, 24, 32, 48, 64, 96, 128, 192, 256],
                     [(0, 1), (0, 2), (1, 2)],
                     pre_traces=[("ε_C = w_C(T[1,3])", lambda W: W[('T', 0, 2)][2]),
                                 ("ε_A = −w_A(B[1,3])", lambda W: -W[('B', 0, 2)][0])])
js['N4'] = run_table("Table 2. N = 4, α = (8, 4, 2, 1), z = 1/4, q = 1, δ = 0.05, η_nb = 0.005, pattern (B, B, B), "
                     "η_p = 0.003", build4, [16, 32, 64, 96, 128, 192],
                     [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3), (0, 3)])

md.append("### Table 3. Approach to the Petlyuk minimum H* (N = 3 data of Table 1, pattern (B, T)): "
          "(δ, η_[1,3], η_p) → 0 jointly, R − H* = δ exactly\n")
md.append("| δ | η_[1,3] | η_1 = η_2 | n | R − H* | impurity of products 1, 2, 3 | designed limit impurity of product 2 | σ_min(DG_∞) |")
md.append("|---|---|---|---|---|---|---|---|")
js['approach'] = []
for delta, epre, eta, n in [('0.05', '0.025', '0.003', 128), ('0.02', '0.01', '1e-3', 192),
                            ('0.005', '0.0025', '2e-4', 256), ('0.001', '5e-4', '3e-5', 384)]:
    b = build3(delta, epre, eta)
    mp.mp.dps = 60
    pb, des, lim = b()
    sl = System(des, lim, None)
    sv = svals(sl.jacobian(sl.start()))[-1]
    res = solve_audited(b, n, int(40 + 0.8 * n))
    imp = impurities(res['sn'], res['v'])
    md.append(f"| {delta} | {epre} | {eta} | {n} | {mp.nstr(res['des'].R - res['pb'].Hmax, 6)} | "
              f"{', '.join(e_(x) for x in imp)} | {e_(1 - res['des'].xprod[1][1], 3)} | {mp.nstr(sv, 4)} |")
    js['approach'].append(dict(delta=delta, eta=eta, n=n, physical=res['phys'], impurities=[str(x) for x in imp],
                               sigma_min_limit=float(sv)))
    print(f"  approach delta {delta}: physical {res['phys']}, impurities {[e_(x) for x in imp]}", flush=True)
    if not res['phys']:
        fails.append(f"approach delta={delta}")
md.append("")
md.append("Reading: at fixed (δ, η) the products 1 and N purify geometrically in n while the side product tends to its "
          "designed impurity (η_1 + η_2)/P_2; letting n → ∞, then η → 0, then δ → 0 gives purity → 1 with "
          "R − H* = δ → 0 (the order of limits of the theorem). The stage numbers in Table 3 grow because "
          "σ_min(DG_∞) ∝ η in the product-form normalization (n_0(η) ~ log(1/η)/log(1/ρ)).")
with open(os.path.join(RUNS, 'paper_tables.md'), 'w', encoding='utf-8') as fh:
    fh.write("\n".join(md) + "\n")
with open(os.path.join(RUNS, 'paper_tables.json'), 'w') as fh:
    json.dump(js, fh, indent=1, default=str)
print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
sys.exit(1 if fails else 0)
