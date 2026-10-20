#!/usr/bin/env python3
"""redteam_design.py -- Agent C (red team) on Agent D's design claims (design/hypotheses.md v1, section 4).

Part 1 (margin identities, D's finding (b)): for every binary column, kappa_p(eta_nb, delta) := excess of the binary
  column over its own minimum vapor (= Psi_T[p,p+1](binary feed root)) at eta_p = 0 satisfies
  kappa_p(eta_nb, delta) - kappa_p(eta_nb, 0) = delta and kappa_p(0, 0) = e_p = H* - H_p.  Also reports
  d kappa_p / d eta_nb at 0 (common eta on all non-binary columns) and the largest common eta_nb for which
  min_p kappa_p >= delta/2 (an explicit eta_bar(delta) for these data).
Part 2 (forcing only increases the margin, D's (c)): the (H4) quantity Psi_B(phi^T) (B-pinned) / Psi_T(phi^B)
  (T-pinned) as a function of eta_p is increasing and positive.
Part 3 (D's explicit bounds on eta_p): at 0.97 x bound the limit configuration passes every check; at 1.03 x bound
  the forced root crosses the neighbouring trace mode alpha_{p-1} (resp. alpha_{p+2}) and the pinch rule of the mode
  table fails ((H5)); a finite-n solve is then attempted with the old mode order to see what breaks.
Prints PASS/FAIL; numerical evidence only.  Usage: python redteam_design.py
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env).  Runtime: a few minutes.
"""
import sys
import time
import mpmath as mp
from petnet import (Problem, Design, limit_config, System, hyp_checks, kappa_check, svals, solve_audited,
                    section_roots, impurities)

t0 = time.time()
fails = []
mp.mp.dps = 50
DATA = [([4, 2, 1], [1, 1, 1], 1), ([4, 2, 1], [1, 1, 1], -0.5), ([8, 4, 2, 1], [1, 1, 1, 1], 1),
        ([6, 3.5, 1.8, 1], [0.3, 0.2, 0.15, 0.35], 0), ([16, 8, 4, 2, 1], [1] * 5, 0.5),
        ([7, 5, 3.4, 2.3, 1.5, 1], [0.2, 0.1, 0.25, 0.15, 0.1, 0.2], 0.8)]


def pattern_bins(N, pat, eta):
    return {(p, p + 1): (pat[p], eta) for p in range(N - 1)}


print("Part 1: margin identities kappa_p(eta_nb, delta) = kappa_p(eta_nb, 0) + delta, kappa_p(0, 0) = e_p")
ok1 = True
for al, z, q in DATA:
    pb = Problem(al, z, q)
    N = pb.N
    e = [pb.Hmax - h for h in pb.H]
    unf = pattern_bins(N, 'N' * (N - 1), 0)
    k00 = Design(pb, 0, eta_nb_default=0, eta_bin=unf)
    err0 = max(abs(k00.col[(p, p + 1)]['excess'] - e[p]) for p in range(N - 1))
    errd = mp.mpf(0)
    for eta_nb in ['0', '0.001', '0.01']:
        for delta in ['0.003', '0.05']:
            dA = Design(pb, mp.mpf(delta), eta_nb_default=mp.mpf(eta_nb), eta_bin=unf)
            dB = Design(pb, 0, eta_nb_default=mp.mpf(eta_nb), eta_bin=unf)
            errd = max(errd, max(abs(dA.col[(p, p + 1)]['excess'] - dB.col[(p, p + 1)]['excess'] - mp.mpf(delta))
                                 for p in range(N - 1)))
    # slope in eta_nb at 0 (delta = 0)
    h = mp.mpf('1e-12')
    dH = Design(pb, 0, eta_nb_default=h, eta_bin=unf)
    slope = [(dH.col[(p, p + 1)]['excess'] - k00.col[(p, p + 1)]['excess']) / h for p in range(N - 1)]
    # eta_bar(delta = 0.05): largest common eta_nb with min_p kappa_p >= delta/2 (bisection on [0, 0.2])
    delta = mp.mpf('0.05')

    def minkap(x):
        try:
            d = Design(pb, delta, eta_nb_default=x, eta_bin=unf)
            return min(d.col[(p, p + 1)]['excess'] - e[p] * 0 for p in range(N - 1)) - delta / 2 \
                if all(min(d.col[c]['d'][j] for j in range(c[0], c[1])) > 0 and min(d.col[c]['b'][j] for j in range(c[0] + 1, c[1] + 1)) > 0 for c in d.col if d.col[c]['kind'] == 'nb') else mp.mpf(-1)
        except Exception:
            return mp.mpf(-1)
    lo, hi = mp.mpf(0), mp.mpf('0.2')
    if minkap(hi) >= 0:
        ebar = hi
    else:
        for _ in range(40):
            m = (lo + hi) / 2
            if minkap(m) >= 0:
                lo = m
            else:
                hi = m
        ebar = lo
    ok = err0 < mp.mpf(10) ** -40 and errd < mp.mpf(10) ** -40
    ok1 = ok1 and ok
    print(f"  N={N} alpha={al} q={q}: |kappa(0,0) - e| = {mp.nstr(err0, 3)}, |kappa(eta,delta) - kappa(eta,0) - delta| = "
          f"{mp.nstr(errd, 3)}; d kappa_p/d eta_nb = {[mp.nstr(x, 4) for x in slope]}; e_p = {[mp.nstr(x, 4) for x in e]}; "
          f"eta_bar(0.05) (min_p kappa_p >= delta/2, splits > 0) = {mp.nstr(ebar, 4)}")
print(f"[{'PASS' if ok1 else 'FAIL'}] margin identities")
if not ok1:
    fails.append('Part 1')

print("Part 2: (H4) quantity versus the forcing eta_p (N = 3, 4; all patterns)")
ok2 = True


def h4_quantity(des, p):
    """B-pinned: Psi_B(phi^T) with phi^T = alpha_p L_T/V_T; T-pinned: Psi_T(phi^B), phi^B = alpha_{p+1} L_B/V_B."""
    pb = des.pb
    al = pb.al
    c = des.col[(p, p + 1)]
    if c['kind'] == 'B':
        phiT = al[p] * c['LT'] / c['VT']
        w = des.wlim[('B', p, p + 1)]
        return des.V[('B', p, p + 1)] - mp.fsum(al[i] * w[i] / (al[i] - phiT) for i in range(pb.N) if w[i] != 0)
    phiB = al[p + 1] * c['LB'] / c['VB']
    w = des.wlim[('T', p, p + 1)]
    return des.V[('T', p, p + 1)] - mp.fsum(al[i] * w[i] / (al[i] - phiB) for i in range(pb.N) if w[i] != 0)


for al, z, q, pats in [([4, 2, 1], [1, 1, 1], 1, ['BB', 'BT', 'TB', 'TT']),
                       ([8, 4, 2, 1], [1, 1, 1, 1], 1, ['BBB', 'TBT', 'BTB', 'TTT'])]:
    pb = Problem(al, z, q)
    N = pb.N
    for pat in pats:
        vals = []
        for eta in ['0', '0.001', '0.003', '0.01', '0.03', '0.06']:
            des = Design(pb, 0.05, eta_nb_default=0.005, eta_bin=pattern_bins(N, pat, mp.mpf(eta)))
            vals.append([h4_quantity(des, p) for p in range(N - 1)])
        inc = all(vals[k + 1][p] >= vals[k][p] for k in range(len(vals) - 1) for p in range(N - 1))
        pos = all(v > 0 for row in vals for v in row)
        ok2 = ok2 and inc and pos
        print(f"  N={N} pattern {pat}: (H4) quantity at eta = 0, 1e-3, 3e-3, 1e-2, 3e-2, 6e-2: "
              + "; ".join("[" + ", ".join(mp.nstr(vals[k][p], 5) for k in range(len(vals))) + "]" for p in range(N - 1))
              + f"  increasing: {inc}, positive: {pos}")
print(f"[{'PASS' if ok2 else 'FAIL'}] forcing increases the (H4) margin")
if not ok2:
    fails.append('Part 2')

print("Part 3: D's explicit bounds on eta_p (N = 4, z = 1/4, q = 1, eta_nb = 0.005, delta = 0.05; volatilities chosen so that")
print("        the root-crossing bound is below the cap f_p: alpha = (4.4, 4, 2, 1) for B-pinned [2,3], (8, 4, 2, 1.8) for T-pinned [2,3])")
ok3 = True
for alx, pat, p in [([4.4, 4, 2, 1], 'BBB', 1), ([4.4, 4, 2, 1], 'TBB', 1), ([8, 4, 2, 1.8], 'BTB', 1), ([8, 4, 2, 1.8], 'TTB', 1)]:
    pb = Problem(alx, [1, 1, 1, 1], 1)
    al = pb.al
    base = Design(pb, 0.05, eta_nb_default=0.005, eta_bin=pattern_bins(4, pat, mp.mpf('0.003')))
    c = base.col[(p, p + 1)]
    f = c['f']
    if pat[p] == 'B':
        V = c['VB']
        Psi = V - al[p + 1] * f[p + 1] / (al[p - 1] - al[p + 1])
        bound = (al[p - 1] - al[p]) * Psi / al[p]
        cap = f[p]
    else:
        V = c['VT']
        Psi = V - al[p] * f[p] / (al[p] - al[p + 2])
        bound = (al[p + 1] - al[p + 2]) * Psi / al[p + 1]
        cap = f[p + 1]
    print(f"  column [{p + 1},{p + 2}] {pat[p]}-pinned: D's bound {mp.nstr(bound, 6)}, cap f = {mp.nstr(cap, 6)}")
    tests = [('0.97', mp.mpf('0.97') * bound), ('1.03', mp.mpf('1.03') * bound)]
    if bound >= cap:
        print(f"    D's bound is not binding (>= cap f); stress test of the forcing near the cap instead")
        tests = [('0.97', mp.mpf('0.5') * cap), ('0.97', mp.mpf('0.95') * cap)]
    for fac, eta in tests:
        if eta >= cap:
            print(f"    eta = {mp.nstr(eta, 5)} >= cap: skipped")
            continue
        eb = pattern_bins(4, pat, mp.mpf('0.003'))
        eb[(p, p + 1)] = (pat[p], eta)
        des = Design(pb, 0.05, eta_nb_default=0.005, eta_bin=eb)
        sec = ('B', p, p + 1) if pat[p] == 'B' else ('T', p, p + 1)
        rts = section_roots(al, des.V[sec], des.wlim[sec])
        kst = p if pat[p] == 'B' else p + 1          # pinch index of the forcing section in the mode table
        crossed = min(abs(rts[kst] - a_) for a_ in al) < mp.mpf(10) ** -30   # a trace mode sits at the pinch slot
        try:
            lim = limit_config(des)
            okh, msgs, hd = hyp_checks(des, lim, K=200)
            kc = kappa_check(System(des, lim, None), System(des, lim, None).start())
            ksig = all(e_.get('sign_ok', True) for e_ in kc)
            sl = System(des, lim, None)
            smin = svals(sl.jacobian(sl.start()))[-1]
            bad = [m_ for m_ in msgs if 'FAIL' in m_]
            print(f"    eta = {mp.nstr(eta, 5)}: forced root {'CROSSED' if crossed else 'below'} the neighbouring alpha; roots "
                  f"{[mp.nstr(x, 6) for x in rts]}; hypothesis checks {'pass' if okh else 'FAIL'}"
                  f"{(': ' + bad[0][:110]) if bad else ''}; kappa signs {'ok' if ksig else 'WRONG'}; sigma_min(DG_inf) {mp.nstr(smin, 4)}")
            expect = (fac == '0.97')
            if expect and not (okh and ksig and not crossed):
                ok3 = False
            if True:
                # beyond the bound: try a finite-n solve with the (now inconsistent) mode table
                def build_(alx=alx):
                    pb_ = Problem(alx, [1, 1, 1, 1], 1)
                    d_ = Design(pb_, 0.05, eta_nb_default=0.005, eta_bin=eb)
                    return pb_, d_, limit_config(d_)
                res = solve_audited(build_, 64, 120)
                if expect and not res['phys']:
                    ok3 = False
                print(f"      finite n = 64: converged {res['conv']}, physical {res['phys']}, min x "
                      f"{mp.nstr(res['chk']['minx'], 4)}, impurities {[mp.nstr(x, 3) for x in impurities(res['sn'], res['v'])]}")
                mp.mp.dps = 50
        except Exception as ex:
            print(f"    eta = {fac} x bound: limit configuration not constructible ({type(ex).__name__}: {str(ex)[:80]}); "
                  f"forced root {'CROSSED' if crossed else 'below'}")
            if fac == '0.97':
                ok3 = False
print(f"[{'PASS' if ok3 else 'FAIL'}] below D's bounds every check passes")
if not ok3:
    fails.append('Part 3')
print(f"{'ALL PASS' if not fails else 'FAIL: ' + ', '.join(fails)}  ({time.time() - t0:.1f}s)")
sys.exit(1 if fails else 0)
