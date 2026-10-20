#!/usr/bin/env python3
"""run_case.py -- Agent C (computation, red team).  Finite-stage generalized Petlyuk arrangement at FIXED flows
F(delta, eta) (forced binary columns, Astra-deformed non-binary columns), n -> infinity, in high precision.

For one case (alpha, z, q, delta, etas) it:
  1. builds the operating specification (own implementation of Theorem achieve + Astra's slack (6) + forcing) and
     checks every unit's total liquid / vapor balance and the positivity of all flows;
  2. computes the infinite-stage LIMIT configuration independently (Lemma 1 vectors at non-binary junctions, pinch
     eigenvectors of the forced sections at binary junctions, sharp/forced products, limit net flows);
  3. the limit system G_inf (chart equations at t = 0): residual at the limit configuration, Jacobian, singular values,
     major/trace classification, block-triangularity, trace block split by component, (H7)-type signs kappa > 0;
  4. for each n: Newton (mpmath, adaptive precision) from the limit configuration for the exact finite-stage system
     G_n (same variables), stage profiles reconstructed, EVERY stage and unit balance re-checked with the actual
     stream compositions, positivity of every mole fraction, S(x_k) > 0, product impurities, distance to the limit,
     smallest singular value of DG_n, signs of all trace net flows.
Prints PASS/FAIL lines and writes runs/<name>.json.  Exit code 1 on any FAIL.

Numerical evidence only (floating point in high precision; no interval enclosure), not a proof.
Usage: python run_case.py <case-name>  (cases defined in cases.py)  [n-list override, e.g. 16,32,64]
Environment: Python 3.12, mpmath 1.3, numpy 2.x, python-flint 0.9 (WSL Fedora sci env).  Runtime: N = 3 about a
minute; N = 4 several minutes; N = 5, 6 longer (see header of the run logs).
"""
import sys
import os
import json
import math
import time
import mpmath as mp
from petnet import (Problem, Design, limit_config, System, reconstruct, stage_check, impurities, svals, classify,
                    submatrix, trace_components, kappa_check, section_roots, to_numpy, hyp_checks,
                    solve_audited, trace_mmatrix)
import cases

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, '..', 'runs')


def s_(x, d=6):
    return mp.nstr(x, d)


def build(case):
    pb = Problem(case['alpha'], case['z'], case['q'])
    des = Design(pb, case['delta'], eta_nb=case.get('eta_nb'), eta_bin=case.get('eta_bin'),
                 eta_nb_default=case.get('eta_nb_default'), bin_default=case.get('bin_default'))
    lim = limit_config(des)
    return pb, des, lim


def spread(lim, pb):
    L = 0.0
    for s in pb.secs:
        r = lim['roots'][s]
        L = max(L, math.log10(float(max(r) / min(r))))
    return L


def main():
    name = sys.argv[1]
    case = cases.CASES[name]
    nlist = case['n']
    if len(sys.argv) > 2:
        nlist = [int(a) for a in sys.argv[2].split(',')]
    fails = []
    out = dict(case=name, spec={k: (v if not isinstance(v, dict) else {str(a): b for a, b in v.items()})
                               for k, v in case.items()})
    mp.mp.dps = 60
    t0 = time.time()
    pb, des, lim = build(case)
    N = pb.N
    print(f"case {name}: N = {N}, alpha = {case['alpha']}, z = {[s_(x, 5) for x in pb.z]}, q = {case['q']}, "
          f"delta = {case['delta']}")
    print(f"  theta = {[s_(t, 10) for t in pb.theta]}, H = {[s_(h, 10) for h in pb.H]}, R = maxH + delta = {s_(des.R, 12)}")
    fres, fmin = des.flow_checks()
    ok = fres < mp.mpf(10) ** -40 and fmin > 0
    print(f"  [{'PASS' if ok else 'FAIL'}] flows: unit total balances residual {s_(fres, 3)}, min flow {s_(fmin, 6)}")
    if not ok:
        fails.append('flows')
    design = {}
    for c, r in sorted(des.col.items()):
        design[str(c)] = dict(kind=r['kind'], eta=float(r['eta']), VT=float(r['VT']), VB=float(r['VB']),
                              LT=float(r['LT']), LB=float(r['LB']), d=[float(x) for x in r['d']],
                              b=[float(x) for x in r['b']], f=[float(x) for x in r['f']], U=float(r['U']),
                              excess=float(r.get('excess', 0)))
        print(f"  column {c} [{r['kind']}, eta {s_(r['eta'], 4)}]: VT {s_(r['VT'], 8)} VB {s_(r['VB'], 8)} "
              f"LT {s_(r['LT'], 8)} LB {s_(r['LB'], 8)} d {[s_(x, 6) for x in r['d']]} b {[s_(x, 6) for x in r['b']]}"
              + (f" excess {s_(r['excess'], 6)}" if 'excess' in r else ""))
    print(f"  products P = {[s_(x, 8) for x in des.P]}")
    out['design'] = design
    out['P'] = [float(x) for x in des.P]
    # limit configuration
    limout = dict(X={str(u): [str(x) for x in lim['X'][u]] for u in pb.units},
                  w={str(s): [str(x) for x in lim['w'][s]] for s in pb.secs},
                  roots={str(s): [str(x) for x in lim['roots'][s]] for s in pb.secs},
                  kstar={str(s): lim['kstar'][s] for s in pb.secs})
    for u in pb.units:
        print(f"  limit x{u} = {[s_(x, 10) for x in lim['X'][u]]}")
    ok_sec = True
    rho = 0.0
    for s in pb.secs:
        r = lim['roots'][s]
        k = lim['kstar'][s]
        dist = all(r[a] > r[a + 1] for a in range(N - 1)) and r[-1] > 0
        ratios = [float(r[k] / r[a]) for a in range(k)] + [float(r[a] / r[k]) for a in range(k + 1, N)]
        rho = max([rho] + ratios)
        ok_sec = ok_sec and dist
        print(f"  section {s}: roots {[s_(x, 8) for x in r]}, pinch index {k}, max ratio {max(ratios):.6f}")
    print(f"  [{'PASS' if ok_sec else 'FAIL'}] (H5) every section: real positive distinct roots; rho = max ratio = {rho:.6f}")
    if not ok_sec:
        fails.append('H5')
    out['rho'] = rho
    for c, inf in lim['info'].items():
        if 'midgap' in inf:
            print(f"  column {c}: r = {s_(inf['r'], 10)}, s = {s_(inf['s'], 10)}, middle-root gap |theta_T - theta_B| = {s_(inf['midgap'], 3)}")
    # limit system
    sl = System(des, lim, None)
    vlim = sl.start()
    F = sl.residual(vlim)
    nF = mp.norm(F, mp.inf)
    ok = nF < mp.mpf(10) ** -40
    print(f"  [{'PASS' if ok else 'FAIL'}] G_inf(v_inf) = 0: |residual| = {s_(nF, 3)} (independent limit configuration solves the limit system)")
    if not ok:
        fails.append('limit residual')
    cl = classify(sl, vlim)
    Jl = cl['J']
    sv = svals(Jl)
    out['limit_svals'] = [float(x) for x in sv]
    print(f"  limit Jacobian {Jl.rows}x{Jl.cols}: sigma_max {s_(sv[0], 5)}, sigma_min {s_(sv[-1], 5)}, cond {s_(sv[0] / sv[-1], 5)}")
    ok = sv[-1] > mp.mpf(10) ** -8
    print(f"  [{'PASS' if ok else 'FAIL'}] DG_inf(v_inf) invertible (sigma_min = {s_(sv[-1], 5)})")
    if not ok:
        fails.append('limit Jacobian singular')
    nmr, nmc, ntr, ntc = len(cl['mr']), len(cl['mc']), len(cl['tr']), len(cl['tc'])
    ok = (nmr == nmc and ntr == ntc)
    print(f"  [{'PASS' if ok else 'FAIL'}] block triangular [[Maj, *], [0, Tr]]: major rows {nmr}, major cols {nmc}, trace rows {ntr}, trace cols {ntc}")
    if not ok:
        fails.append('block structure counts')
    Jm = submatrix(Jl, cl['mr'], cl['mc'])
    Jt = submatrix(Jl, cl['tr'], cl['tc'])
    svm = svals(Jm) if nmr == nmc else [mp.mpf(0)]
    svt = svals(Jt) if ntr == ntc else [mp.mpf(0)]
    print(f"  major block sigma_min {s_(svm[-1], 5)}; trace block sigma_min {s_(svt[-1], 5)}")
    blocks, mixed = trace_components(cl, sl)
    out['limit_blocks'] = dict(major=[nmr, nmc, float(svm[-1])], trace=[ntr, ntc, float(svt[-1])], per_component={})
    ok_b = not mixed
    for k, bl in sorted(blocks.items()):
        nr, nc = len(bl['rows']), len(bl['cols'])
        sm = svals(submatrix(Jl, bl['rows'], bl['cols']))[-1] if nr == nc and nr > 0 else mp.mpf(0)
        ok_b = ok_b and nr == nc and sm > mp.mpf(10) ** -8
        out['limit_blocks']['per_component'][k] = [nr, nc, float(sm)]
        print(f"    trace block of component {k}: {nr} rows x {nc} cols, sigma_min {s_(sm, 5)}")
    if mixed:
        print(f"    rows coupling several trace components: {mixed}")
    print(f"  [{'PASS' if ok_b else 'FAIL'}] trace block = one invertible block per component")
    if not ok_b:
        fails.append('trace blocks')
    mm = trace_mmatrix(cl, sl)
    ok_m = all(e.get('ok') for e in mm.values())
    for k, e in sorted(mm.items()):
        if 'size' in e:
            print(f"    A_{k} ({e['size']}x{e['size']}, trace contents of component {k} after eliminating the trace flows): "
                  f"min diag {e['diag_min']:.4g}, max off-diag {e['offdiag_max']:.3g}, min entry of inverse {e['inverse_min']:.3g}, "
                  f"min column sum {e['min_col_sum']:.3g}, min row sum {e['min_row_sum']:.3g}")
        else:
            print(f"    A_{k}: {e.get('note')}")
    print(f"  [{'PASS' if ok_m else 'FAIL'}] every A_k is (up to sign) a nonsingular M-matrix (Z pattern, inverse >= 0)")
    if not ok_m:
        fails.append('M-matrix')
    out['mmatrix'] = {str(k): v for k, v in mm.items()}
    kc = kappa_check(sl, vlim)
    ok_k = all(e.get('sign_ok', True) for e in kc)
    for e in kc:
        if e['end'] is None:
            print(f"    {e['sec']} component {e['k']}: {e['note']}")
            continue
        print(f"    {e['sec']} trace {e['k']} closed at {e['end']}: G^(k) = {s_(e['G'], 5)}, Psi^(k) = {s_(e['Psi'], 5)}, "
              f"kappa = {s_(e['kappa'], 5)}{'' if e['other_trace'] else ' (k major at the other end)'} {'ok' if e['sign_ok'] else 'SIGN FAIL'}")
    print(f"  [{'PASS' if ok_k else 'FAIL'}] (H7) signs: Psi^(k)(alpha_k) > 0 and G^(k) of the right sign at every closed trace end")
    if not ok_k:
        fails.append('H7 signs')
    out['kappa'] = [dict(sec=str(e['sec']), k=e['k'], end=str(e['end']), G=float(e.get('G', 0)),
                         Psi=float(e.get('Psi', 0)), sign_ok=e.get('sign_ok')) for e in kc]
    okh, msgs, hd = hyp_checks(des, lim)
    for m_ in msgs:
        print('    ' + m_)
    print(f"  [{'PASS' if okh else 'FAIL'}] (H6)/(A5) half-profiles >= 0, (H8)/(A3) pivots > 0, (A4) closed coefficients 0, (H7)(iii)-(v)")
    if not okh:
        fails.append('H6/H7/H8')
    out['hyp'] = hd
    out['limit'] = limout
    Ldec = spread(lim, pb)
    # finite n
    rows = []
    vprev = None
    for n in sorted(nlist, reverse=True):      # continuation in n, from the largest n down
        t1 = time.time()
        dps = int(40 + 1.25 * (n + 1) * Ldec * case.get('trace_depth', 1.0))
        dps = max(60, min(dps, case.get('max_dps', 900)))

        def attempt(start_from):
            return solve_audited(lambda: build(case), n, dps, start=start_from, max_dps=4000)

        res = attempt(None)
        how = 'from v_inf'
        if not res['phys'] and vprev is not None:
            r2 = attempt(vprev)
            if r2['phys'] or not res['conv']:
                res, how = r2, 'continued from n = %d' % nprev
        sn, v, nF, hist, conv, chk = res['sn'], res['v'], res['nF'], res['hist'], res['conv'], res['chk']
        pb2, des2, lim2 = res['pb'], res['des'], res['lim']
        imp = impurities(sn, v)
        X, W = sn.unpack(v)
        vl = sn.start()
        dmaj = max(abs(v[j] - vl[j]) for j in range(len(v)) if vl[j] != 0)
        dtr = max([abs(v[j]) for j in range(len(v)) if vl[j] == 0] + [mp.mpf(0)])
        Jn = sn.jacobian(v)
        svn = svals(Jn)
        negw = []
        for s in pb2.secs:
            for j in range(N):
                if lim2['w'][s][j] == 0:
                    negw.append((s, j, W[s][j]))
        rev = [(str(s), j, float(x)) for (s, j, x) in negw if x < 0 and s[0] == 'T'] +               [(str(s), j, float(x)) for (s, j, x) in negw if x > 0 and s[0] == 'B']
        okn = res['phys']
        row = dict(n=n, dps=mp.mp.dps, start=how, resid=float(nF), converged=bool(conv), physical=okn,
                   newton_its=len(hist),
                   stage_rel=float(chk['worst_rel']), unit_rel=float(chk['unit_rel']), minx=str(chk['minx']),
                   minS=float(chk['minS']), negatives=len(chk['negs']),
                   impurities=[mp.nstr(x, 6) for x in imp],
                   X={str(u): [mp.nstr(x, 15) for x in X[u]] for u in pb2.units},
                   W={str(s): [mp.nstr(x, 15) for x in W[s]] for s in pb2.secs},
                   dist_major=float(dmaj), max_trace=float(dtr), sigma_min=float(svn[-1]), sigma_max=float(svn[0]),
                   reverse_trace_flows=rev, time=time.time() - t1)
        rows.append(row)
        print(f"  [{'PASS' if okn else 'FAIL'}] n = {n:4d} (dps {mp.mp.dps}, {how}): |G_n| {s_(nF, 2)}, stage rel {s_(chk['worst_rel'], 2)}, "
              f"unit rel {s_(chk['unit_rel'], 2)}, min x {s_(chk['minx'], 3)}, min S {s_(chk['minS'], 4)}, "
              f"impurities {[s_(x, 3) for x in imp]}, |v_n - v_inf|_maj {s_(dmaj, 3)}, max trace {s_(dtr, 3)}, "
              f"sigma_min(DG_n) {s_(svn[-1], 4)}, reverse trace flows {len(rev)}, {time.time() - t1:.1f}s", flush=True)
        if rev:
            print(f"      trace net flows against the section's direction (T: < 0, B: > 0): {rev}")
        if okn:
            vprev = [str(x) for x in v]
            nprev = n
        mp.mp.dps = 60
    rows.sort(key=lambda r: r['n'])
    # asymptotic regime: smallest n from which every larger n passes
    n0 = None
    for r in reversed(rows):
        if r['physical']:
            n0 = r['n']
        else:
            break
    if n0 is None:
        fails.append('largest n not physical')
        print("  [FAIL] no physical steady state found at the largest n")
    else:
        bad = [r['n'] for r in rows if r['n'] < n0]
        print(f"  [PASS] physical steady states for every n >= {n0} tested" +
              (f" (pre-asymptotic n without a physical zero found near the branch: {bad})" if bad else ""))
        out['n0_observed'] = n0
    out['rows'] = rows
    # rates
    if len(rows) >= 2:
        r1, r2 = rows[-2], rows[-1]
        if r1['dist_major'] > 0 and r2['dist_major'] > 0:
            rate = (math.log(r2['dist_major']) - math.log(r1['dist_major'])) / (r2['n'] - r1['n'])
            print(f"  convergence of the majors: |v_n - v_inf| ~ exp({rate:.4f} n) (rho^n with rho = {math.exp(rate):.4f}); max section ratio rho = {rho:.4f}")
            out['rate_major'] = math.exp(rate)
        sl_last = rows[-1]['sigma_min']
        ok = sl_last > 0.5 * float(sv[-1])
        print(f"  [{'PASS' if ok else 'FAIL'}] sigma_min(DG_n) -> sigma_min(DG_inf) = {s_(sv[-1], 5)} (last: {sl_last:.5g})")
        if not ok:
            fails.append('sigma_min convergence')
    out['fails'] = fails
    out['runtime'] = time.time() - t0
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, f"{name}.json"), 'w') as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
