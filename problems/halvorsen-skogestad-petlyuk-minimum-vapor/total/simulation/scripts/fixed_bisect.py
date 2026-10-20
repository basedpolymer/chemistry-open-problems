"""
fixed_bisect.py -- upper bound on the finite-stage FTC minimum vaporization at purity 1-eps with n stages per
section, along the one-parameter family of flows of ftc_fixed_flows.py (result 12's Underwood attainment flows +
a circulation delta through every column; R = max_j H_j + K delta, K = number of columns):
    R_up(n, eps) = max_j H_j + K delta*(n, eps),  delta*(n, eps) = min{delta : max product impurity <= eps}
(bisection on log delta; every evaluation is a converged MESH steady state, warm-started).  Together with the
stage-level lower bound R >= max_j H_j - c eps (class.md, valid for every arrangement), this brackets the minimum.

Agent P (Claude Opus 5.5, effort xhigh).
Usage: python fixed_bisect.py DATASET "n1,n2,..." "eps1,eps2,..."  ->  ../runs/bisect_DATASET.json
"""
import sys, json, time
import numpy as np
from arrangements import regular_network, ftc_columns
from ftc_fixed_flows import flows_for
from refvals import peaks
from attain import DATA, cex5

ds = sys.argv[1]
ns = [int(v) for v in sys.argv[2].split(',')]
epss = [float(v) for v in sys.argv[3].split(',')]
alpha, z, q = cex5() if ds == 'CEX5' else DATA[ds]
N = len(alpha)
H, VF, th = peaks(alpha, z, q)
out = dict(dataset=ds, alpha=alpha, z=z, q=q, H=H, V_FTC=VF, runs=[])
for n in ns:
    net = regular_network(alpha, z, q, ftc_columns(N), {}, n=n)
    cache = {}
    xw = [None]

    def imp(delta):
        if delta in cache:
            return cache[delta]
        f, Hmax, K = flows_for(net, alpha, z, q, delta)
        x, info = net.solve(f, x0=xw[0])
        if not info['ok']:
            x, info = net.solve(f)
        if info['ok']:
            xw[0] = x
            v = max(net.impurities(x).values())
        else:
            v = np.nan
        cache[delta] = (v, info['ok'], K, net.R(f))
        return cache[delta]

    for eps in epss:
        t = time.time()
        lo, hi = 1e-6, 1.0                      # delta range (log bisection)
        vhi = imp(hi)[0]
        if not (vhi <= eps):
            print(json.dumps(dict(n=n, eps=eps, error='not reachable with delta <= 1', imp_at_1=vhi)), flush=True)
            continue
        vlo = imp(lo)[0]
        if vlo <= eps:
            hi = lo
        else:
            for _ in range(40):
                mid = np.sqrt(lo * hi)
                v, ok, K, _ = imp(mid)
                if ok and v <= eps:
                    hi = mid
                else:
                    lo = mid
                if hi / lo < 1 + 1e-4:
                    break
        v, ok, K, Rh = imp(hi)
        row = dict(n=n, eps=eps, delta_star=hi, K=K, R_up=Rh, ratio_up=Rh / VF, max_imp=v,
                   sec=round(time.time() - t, 1))
        out['runs'].append(row)
        print(json.dumps(row), flush=True)
        json.dump(out, open('../runs/bisect_%s.json' % ds, 'w'), indent=1)
