#!/usr/bin/env python3
"""compare_D.py -- Agent C.  Independent recomputation of Agent D's exported specifications and limit
configurations ("hs-petlyuk-limit-config v1" JSON, produced by D's design/scripts/export_config.py; copies in
../runs/D_exports/ made with run_D_export.py).  From the JSON's data (alpha, z, q, delta) and parameters (pattern,
eta's) only, my own code (petnet.py: Theorem achieve + Astra's slack + forcing; Lemma 1 vectors and pinch
eigenvectors; my own eigen-solver) recomputes everything and compares:
  flows V, L of every section; product flows; limit net flows; limit node compositions; section modes; pinch index;
  pivots (coefficients of the end compositions on the normalized pinch vector); kappa at closed trace ends; Astra's
  cone weights A, B, C, D of every non-binary column; then solves my limit system at D's limit point
  (residual of G_inf at D's configuration) and checks invertibility of DG_inf there.
Prints the max relative discrepancy per quantity and PASS/FAIL (tolerance 1e-40).  Numerical evidence only.
Usage: python compare_D.py [json files ...]   (default: every ../runs/D_exports/*.json)
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env).  Runtime: about a minute.
"""
import os
import sys
import glob
import json
import time
import mpmath as mp
from petnet import (Problem, Design, limit_config, System, lhat, dot, right_vec_lim, svals)

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 60
TOL = mp.mpf(10) ** -40


def parse_col(s):            # "[i,m]" -> (i-1, m-1)
    a, b = s.strip('[]').split(',')
    return int(a) - 1, int(b) - 1


def node_name(u):
    i, m = u
    return f"P[{i + 1}]" if i == m else f"J[{i + 1},{m + 1}]"


def sec_name(s):
    return f"{s[0]}[{s[1] + 1},{s[2] + 1}]"


def rel(a, b):
    a, b = mp.mpf(a), mp.mpf(b)
    return abs(a - b) / max(mp.mpf(1), abs(a), abs(b))


def compare(path):
    D = json.load(open(path))
    data, par = D['data'], D['params']
    N = int(data['N'])
    pb = Problem([mp.mpf(x) for x in data['alpha']], [mp.mpf(x) for x in data['z']], mp.mpf(data['q']))
    pat = par['pattern']
    eta_nb = {parse_col(k): mp.mpf(v) for k, v in par['eta_nb'].items()}
    eta_bin = {(int(p) - 1, int(p)): (pat[int(p) - 1], mp.mpf(v)) for p, v in par['eta_bin'].items()}
    des = Design(pb, mp.mpf(data['delta']), eta_nb=eta_nb, eta_bin=eta_bin)
    lim = limit_config(des)
    err = {}

    def upd(key, e):
        err[key] = max(err.get(key, mp.mpf(0)), e)
    upd('theta', max(rel(a, b) for a, b in zip(pb.theta, data['theta'])))
    upd('H', max(rel(a, b) for a, b in zip(pb.H, data['H'])))
    upd('R', rel(des.R, data['R']))
    for s in pb.secs:
        Ds = D['sections'][sec_name(s)]
        upd('V, L', max(rel(des.V[s], Ds['V']), rel(des.L[s], Ds['L'])))
        upd('limit net flows w', max(rel(a, b) for a, b in zip(lim['w'][s], Ds['w'])))
        modes = sorted([mp.mpf(m_['value']) for m_ in Ds['modes']], reverse=True)
        upd('modes', max(rel(a, b) for a, b in zip(lim['roots'][s], modes)))
        if Ds['pinch'] != lim['kstar'][s]:
            upd('pinch index', mp.mpf(1))
        else:
            upd('pinch index', mp.mpf(0))
        # pivots: coefficients of x_bot, x_top on e_* (1.e_* = 1, l_*.e_* = 1)
        bu, tu = pb.sec_ends(s)
        k = lim['kstar'][s]
        ph = lim['roots'][s][k]
        e = right_vec_lim(pb.al, des.V[s], lim['w'][s], ph)
        se = mp.fsum(e)
        e = [x / se for x in e]
        l = lhat(pb.al, ph)
        le = dot(l, e)
        piv = [dot(l, lim['X'][bu]) / le, dot(l, lim['X'][tu]) / le]
        upd('pivots', max(rel(a, b) for a, b in zip(piv, Ds['pivots'])))
        for tr in Ds['traces']:
            kk = tr['k'] - 1
            end = tu if tr['end'] == 'top' else bu
            xJ = lim['X'][end]
            G = mp.fsum(pb.al[i] * xJ[i] / (pb.al[i] - pb.al[kk]) for i in range(N) if i != kk)
            Psi = des.V[s] - mp.fsum(pb.al[i] * lim['w'][s][i] / (pb.al[i] - pb.al[kk]) for i in range(N) if i != kk)
            upd('kappa, G_k, Psi_k', max(rel(G, tr['G_k']), rel(Psi, tr['Psi_k']), rel(abs(G) / Psi, tr['kappa'])))
    for u in pb.units:
        upd('limit node compositions', max(rel(a, b) for a, b in zip(lim['X'][u], D['nodes'][node_name(u)])))
    for p in range(N):
        upd('product flows', rel(des.P[p], D['products'][str(p + 1)]['P']))
    # Astra cone weights: X = A e_c + B e_u (top section), = C e_v + D e_a (bottom section)
    for c, rec in des.col.items():
        if rec['kind'] != 'nb':
            continue
        i, m = c
        Dc = D['columns'][f"[{i + 1},{m + 1}]"]
        X = lim['X'][c]
        out = []
        for s, mode_idx in ((('T', i, m), 'T'), (('B', i, m), 'B')):
            rts = lim['roots'][s]
            k = lim['kstar'][s]
            # the trace mode used: alpha_m (top) or alpha_i (bottom)
            t_al = pb.al[m] if mode_idx == 'T' else pb.al[i]
            kt = min(range(N), key=lambda r: abs(rts[r] - t_al))
            coeffs = {}
            for r in (k, kt):
                e = right_vec_lim(pb.al, des.V[s], lim['w'][s], rts[r])
                se = mp.fsum(e)
                e = [x / se for x in e]
                l = lhat(pb.al, rts[r])
                coeffs[r] = dot(l, X) / dot(l, e)
            if mode_idx == 'T':
                out += [coeffs[kt], coeffs[k]]      # A (e_c, trace mode alpha_m), B (e_u, pinch)
            else:
                out += [coeffs[k], coeffs[kt]]      # C (e_v, pinch), D (e_a, trace mode alpha_i)
        upd('cone weights A, B, C, D', max(rel(a, b) for a, b in zip(out, Dc['ABCD'])))
        if min(out) <= 0:
            upd('cone weights positive', mp.mpf(1))
    # my limit system at D's limit point
    sl = System(des, lim, None)
    Xd = {u: [mp.mpf(x) for x in D['nodes'][node_name(u)]] for u in pb.units}
    Wd = {s: [mp.mpf(x) for x in D['sections'][sec_name(s)]['w']] for s in pb.secs}
    vD = sl.pack(Xd, Wd)
    FD = mp.norm(sl.residual(vD), mp.inf)
    sv = svals(sl.jacobian(vD))
    ok = all(v_ < TOL for v_ in err.values()) and FD < TOL and sv[-1] > mp.mpf(10) ** -10
    return ok, err, FD, sv[-1], N, pat


def main():
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(HERE, '..', 'runs', 'D_exports', '*.json')))
    t0 = time.time()
    fails = []
    for f in files:
        ok, err, FD, smin, N, pat = compare(f)
        name = os.path.basename(f)
        worst = max(err.items(), key=lambda kv: kv[1])
        print(f"[{'PASS' if ok else 'FAIL'}] {name} (N = {N}, pattern {pat}): max rel. discrepancy {mp.nstr(worst[1], 3)} "
              f"({worst[0]}); |G_inf(D's limit point)| = {mp.nstr(FD, 3)}; sigma_min(DG_inf) = {mp.nstr(smin, 4)}")
        for k_, v_ in sorted(err.items()):
            if v_ >= TOL:
                print(f"      {k_}: {mp.nstr(v_, 4)}")
        if not ok:
            fails.append(name)
    print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
