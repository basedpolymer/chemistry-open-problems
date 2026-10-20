"""
ftc_fixed_flows.py -- attainment at FIXED flows: the finite-stage FTC (generalized Petlyuk) operated at the
Underwood attainment flows of result 12 (Theorem achieve: every column with >= 3 components at its preferred split,
binary columns with excess e_r = max H - H_r) plus a vapor margin delta circulated through every column
(both sections of every column get +delta in V and L; net component flows and products unchanged), so that
R = max_j H_j + (number of columns) * delta exactly.  The stage-by-stage steady state (stagenet.py, independent
of Underwood) is computed for n = n1, n2, ... stages per section, and the product impurities are recorded:
attainment predicts impurities -> 0 (geometrically) as n -> infinity for every fixed delta > 0.

The Underwood construction is used ONLY to choose the flows; purities are physical (MESH) results.
Agent P (Claude Opus 5.5, effort xhigh).
Usage: python ftc_fixed_flows.py DATASET "delta1,delta2" "n1,n2,..."   ->  ../runs/fixed_DATASET.json
"""
import sys, json, time
import mpmath as mp
import numpy as np
from arrangements import regular_network, ftc_columns, lab, LET
from refvals import feed_roots, peaks
from attain import DATA, cex5

mp.mp.dps = 50


def attainment_flows(alpha, z, q):
    """Result 12, Theorem achieve (FTC): per-column d, b, V_T, V_B (mpmath)."""
    N = len(alpha)
    a = [mp.mpf(v) for v in alpha]
    th = feed_roots(alpha, z, q)
    H, VF, _ = peaks(alpha, z, q)
    Hm = [mp.mpf(h) for h in H]
    Hmax = max(Hm)
    cols = sorted(ftc_columns(N), key=lambda c: -(c[1] - c[0]))
    feed = {(0, N - 1): [mp.mpf(v) for v in z]}
    U = {(0, N - 1): 1 - mp.mpf(q)}
    res = {}
    for (i, m) in cols:
        f = feed[(i, m)]
        Uc = U[(i, m)]
        d = [mp.mpf(0)] * N; b = [mp.mpf(0)] * N
        if m - i >= 2:
            def coef(p):
                num = mp.mpf(1)
                for r in range(i, m):
                    num *= (a[p] - th[r])
                den = mp.mpf(1)
                for pp in range(i, m):
                    if pp != p:
                        den *= (a[p] - a[pp])
                return num / den / a[p]
            VT = f[i] / coef(i)
            for p in range(i, m):
                d[p] = VT * coef(p)
            for p in range(i + 1, m + 1):
                b[p] = f[p] - d[p]
            VB = VT - Uc
        else:
            r = i
            e = Hmax - Hm[r]
            d[r] = f[r]; b[r + 1] = f[r + 1]
            VT = a[r] * f[r] / (a[r] - th[r]) + e
            VB = a[r + 1] * f[r + 1] / (th[r] - a[r + 1]) + e
        res[(i, m)] = dict(d=d, b=b, VT=VT, VB=VB)
        # children: distillate [i, m-1] and residue [i+1, m]
        for child, vec, sign in (((i, m - 1), d, +1), ((i + 1, m), b, -1)):
            if child[1] - child[0] >= 1:
                feed.setdefault(child, [mp.mpf(0)] * N)
                feed[child] = [feed[child][p] + vec[p] for p in range(N)]
                U.setdefault(child, mp.mpf(0))
                U[child] += VT if sign > 0 else -VB
    return res, Hmax


def column_path(N, c):
    """Sections (names) on a vapor path reboiler -> ... -> column c -> ... -> condenser through both sections of c."""
    i, m = c
    secs = []
    # upward from C[i,m]: T[i,m], T[i,m-1], ..., T[i,i+1] reaches P_i; then B[i-1,i], T[i-1,i], ... to P_1
    k = m
    while k > i:
        secs.append('T' + lab(i, k)); k -= 1
    p = i
    while p > 0:
        secs.append('B' + lab(p - 1, p)); secs.append('T' + lab(p - 1, p)); p -= 1
    # downward from C[i,m]: B[i,m], B[i+1,m], ..., B[m-1,m] reaches P_m; then T[m,m+1], B[m,m+1], ... to P_N
    k = i
    while k < m:
        secs.append('B' + lab(k, m)); k += 1
    p = m
    while p < N - 1:
        secs.append('T' + lab(p, p + 1)); secs.append('B' + lab(p, p + 1)); p += 1
    return secs


def flows_for(net, alpha, z, q, delta):
    N = len(alpha)
    res, Hmax = attainment_flows(alpha, z, q)
    f = np.zeros(net.nf)
    for (i, m), r in res.items():
        c = lab(i, m)
        D = sum(r['d']); B = sum(r['b'])
        f[net.fidx['V:T' + c]] = float(r['VT']); f[net.fidx['L:T' + c]] = float(r['VT'] - D)
        f[net.fidx['V:B' + c]] = float(r['VB']); f[net.fidx['L:B' + c]] = float(r['VB'] + B)
    for p in range(N):
        f[net.fidx['P:' + LET[p]]] = float(z[p])
    for c in res:
        for s in column_path(N, c):
            f[net.fidx['V:' + s]] += delta; f[net.fidx['L:' + s]] += delta
    f[net.fidx['Q:P' + LET[N - 1]]] = f[net.fidx['V:B' + lab(N - 2, N - 1)]]
    f[net.fidx['Q:P' + LET[0]]] = -f[net.fidx['V:T' + lab(0, 1)]]
    return f, float(Hmax), len(res)


if __name__ == '__main__':
    ds = sys.argv[1]
    deltas = [float(v) for v in sys.argv[2].split(',')]
    ns = [int(v) for v in sys.argv[3].split(',')]
    alpha, z, q = cex5() if ds == 'CEX5' else DATA[ds]
    N = len(alpha)
    out = dict(dataset=ds, alpha=alpha, z=z, q=q, runs=[])
    for delta in deltas:
        x = None; prev = None
        for n in ns:
            net = regular_network(alpha, z, q, ftc_columns(N), {}, n=n)
            f, Hmax, ncol = flows_for(net, alpha, z, q, delta)
            bal = net.check_flows(f)
            t = time.time()
            x0 = net.interpolate_from(prev[0], prev[1]) if prev is not None else None
            x, info = net.solve(f, x0=x0)
            if not info['ok'] and x0 is not None:
                x, info = net.solve(f)
            if info['ok']:
                prev = (net, x)
            imp = net.impurities(x)
            row = dict(delta=delta, n=n, R=net.R(f), V_FTC=Hmax, R_minus_VFTC=net.R(f) - Hmax, balance=bal,
                       solve_ok=bool(info['ok']), max_imp=max(imp.values()), impurities=imp,
                       min_flow=float(min(f[i] for i, nm in enumerate(net.fnames) if nm[0] in 'VL')),
                       sec=round(time.time() - t, 1))
            out['runs'].append(row)
            print(json.dumps(row), flush=True)
            json.dump(out, open('../runs/fixed_%s.json' % ds, 'w'), indent=1)
            if info['ok'] and row['max_imp'] < 1e-12:
                break       # beyond double-precision resolution of the global balances
