"""
permissive_physical.py -- physical (stage-by-stage) test of the permissive-model counterexample of result 12.

Agent P (Claude Opus 5.5, effort xhigh).  The certified permissive point (verification/permissive_counterexample_N5.json:
basic configuration 161, all thermal couplings, N = 5, q = 0, R = 0.4696851 < V_FTC = 0.4699411) prescribes every
section's vapor and liquid flow (V_T, V_B and the net flows d, b of every column).  We build the same arrangement
with n equilibrium stages per section (arrangements.regular_network, result-12 conventions), impose exactly these
flows (so R = 0.4696851), compute the physical steady state for n = n1, n2, ... and report the product impurities
and the stage-level bound max_j sum_P max(0, phi_P(theta_j)) (which must be <= R at every steady state).
If the point were physically realizable with pure products, impurities would -> 0 as n -> infinity; the
stage-level theorem says they cannot (R < max_j H_j).
Usage: python permissive_physical.py "n1,n2,..."   ->  ../runs/permissive_physical.json
"""
import sys, json, time
from fractions import Fraction
import numpy as np
from arrangements import regular_network, lab, LET
from attain import CONF, cex5
from study import stage_bound
from refvals import peaks

d = json.load(open('../../../verification/permissive_counterexample_N5.json'))
c = {k: float(v) for k, v in d['center'].items()}
alpha, z, q = cex5()
N = 5
H, VF, th = peaks(alpha, z, q)
cols, tr = CONF['C161']


def flows(net):
    f = np.zeros(net.nf)
    for (i, m), (k, l) in cols.items():
        key = '(%d, %d)' % (i + 1, m + 1)
        D = sum(v for kk, v in c.items() if kk.startswith('d' + key))
        B = sum(v for kk, v in c.items() if kk.startswith('b' + key))
        VT, VB = c['VT' + key], c['VB' + key]
        s = lab(i, m)
        f[net.fidx['V:T' + s]] = VT; f[net.fidx['L:T' + s]] = VT - D
        f[net.fidx['V:B' + s]] = VB; f[net.fidx['L:B' + s]] = VB + B
    for p in range(N):
        f[net.fidx['P:' + LET[p]]] = z[p]
    f[net.fidx['Q:P' + LET[N - 1]]] = f[net.fidx['V:B' + lab(3, 4)]]
    f[net.fidx['Q:P' + LET[0]]] = -f[net.fidx['V:T' + lab(0, 1)]]
    return f


ns = [int(v) for v in sys.argv[1].split(',')] if len(sys.argv) > 1 else [20, 40, 80, 160]
out = dict(V_FTC=VF, H=H, R_permissive=None, runs=[])
x = None; prev = None
for n in ns:
    net = regular_network(alpha, z, q, cols, tr, n=n)
    f = flows(net)
    bal = net.check_flows(f)
    t = time.time()
    x0 = net.interpolate_from(prev[0], prev[1]) if prev is not None else None
    x, info = net.solve(f, x0=x0)
    if not info['ok'] and x0 is not None:
        x, info = net.solve(f)
    if info['ok']:
        prev = (net, x)
    imp = net.impurities(x)
    b, per = stage_bound(net, f, x, alpha, z, q)
    out['R_permissive'] = net.R(f)
    row = dict(n=n, R=net.R(f), R_over_VFTC=net.R(f) / VF, balance=bal, solve_ok=bool(info['ok']),
               impurities=imp, max_imp=max(imp.values()), bound=b, bound_per_j=per, R_minus_bound=net.R(f) - b,
               sec=round(time.time() - t, 1))
    out['runs'].append(row)
    print(json.dumps(row), flush=True)
    json.dump(out, open('../runs/permissive_physical.json', 'w'), indent=1)
