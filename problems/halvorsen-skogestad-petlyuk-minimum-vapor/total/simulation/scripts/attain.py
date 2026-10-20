"""
attain.py -- purity/stage study of the minimum total vaporization of an arrangement in the MESH model.

Agent P (Claude Opus 5.5, effort xhigh).
Usage: python attain.py DATASET CONFIG "n1,n2,..." "eps1,eps2,..." [tag]
  DATASET: T61 T61q0 T61q05 T61q13 P2 E8 M05 Q04 Q13 QN N13 N03 CEX5
  CONFIG : FTC | K (N=4 family K, all TC) | C161 (N=5 basic configuration 161, all TC) | direct | indirect
For every (eps, n): IPOPT minimum of R at purity >= 1-eps (all products), re-solved by stagenet's Newton solver;
records R, R/V_FTC (V_FTC = max_j H_j from refvals), the stage-level bound at the optimum, the true impurities.
Output: ../runs/attain_DATASET_CONFIG[_tag].json
"""
import sys, json, time
from fractions import Fraction
import numpy as np
from study import run_config
from arrangements import ftc_columns
from refvals import peaks

DATA = {
    'T61': ((4.0, 2.0, 1.0), (1 / 3, 1 / 3, 1 / 3), 1.0),
    'T61q0': ((4.0, 2.0, 1.0), (1 / 3, 1 / 3, 1 / 3), 0.0),
    'T61q05': ((4.0, 2.0, 1.0), (1 / 3, 1 / 3, 1 / 3), 0.5),
    'T61q13': ((4.0, 2.0, 1.0), (1 / 3, 1 / 3, 1 / 3), 1.3),
    'P2': ((47.5774, 25.5263, 15.1071, 1.0), (0.1358, 0.5169, 0.0923, 0.255), 1.0),
    'E8': ((8.0, 4.0, 2.0, 1.0), (0.25, 0.25, 0.25, 0.25), 1.0),
    'Q04': ((6.0, 3.0, 1.6, 1.0), (0.15, 0.35, 0.35, 0.15), 0.4),
    'Q13': ((5.0, 2.2, 1.5, 1.0), (0.2, 0.3, 0.3, 0.2), 1.3),
    'QN': ((7.0, 4.0, 1.4, 1.0), (0.1, 0.4, 0.4, 0.1), -0.3),
    'M05': ((10.0, 5.0, 4.2, 1.0), (0.15, 0.35, 0.35, 0.15), 0.5),
    'N13': ((14.5262, 5.957, 3.327, 1.0), (0.1368, 0.2208, 0.3374, 0.305), 1.3),
    'N03': ((5.9208, 4.3871, 3.4747, 1.0), (0.4139, 0.2266, 0.1566, 0.2029), 0.3),
}


def cex5():
    d = json.load(open('../../../verification/permissive_counterexample_N5.json'))
    a = tuple(float(Fraction(s)) for s in d['alpha'])
    z = tuple(float(Fraction(s)) for s in d['z'])
    return a, z, float(Fraction(d['q']))


CONF = {
    'K': ({(0, 3): (2, 1), (0, 2): (1, 2), (1, 3): (1, 2), (0, 1): (0, 1), (2, 3): (2, 3)}, {}),
    'C161': ({(0, 4): (3, 1), (0, 3): (2, 1), (1, 4): (3, 2), (0, 2): (1, 1), (1, 3): (2, 3), (2, 4): (2, 3),
              (0, 1): (0, 1), (1, 2): (1, 2), (3, 4): (3, 4)}, {}),
}


def config(name, N):
    if name == 'FTC':
        return ftc_columns(N), {}
    if name == 'direct':
        return {(i, N - 1): (i, i + 1) for i in range(N - 1)}, {LETS(i + 1, N - 1): 'HX' for i in range(N - 2)}
    if name == 'indirect':
        return {(0, m): (m - 1, m) for m in range(1, N)}, {LETS(0, m - 1): 'HX' for m in range(2, N)}
    return CONF[name]


def LETS(i, m):
    return 'ABCDEFGH'[i:m + 1]


if __name__ == '__main__':
    ds, cf = sys.argv[1], sys.argv[2]
    ns = [int(v) for v in sys.argv[3].split(',')]
    epss = [float(v) for v in sys.argv[4].split(',')]
    tag = ('_' + sys.argv[5]) if len(sys.argv) > 5 else ''
    alpha, z, q = cex5() if ds == 'CEX5' else DATA[ds]
    N = len(alpha)
    H, VF, th = peaks(alpha, z, q)
    cols, tr = config(cf, N)
    out = dict(dataset=ds, config=cf, alpha=alpha, z=z, q=q, H=H, V_FTC=VF, runs=[])
    print('V_FTC', VF, 'H', H); sys.stdout.flush()
    for eps in epss:
        for n in ns:
            try:
                r = run_config(alpha, z, q, cols, tr, n, eps)
            except Exception as ex:
                print(json.dumps(dict(n=n, eps=eps, error=str(ex)))); sys.stdout.flush()
                continue
            row = dict(n=n, eps=eps, R=r['R'], ratio=r['R'] / VF, bound=r['bound'], bound_per_j=r['bound_per_j'],
                       R_minus_bound=r['R'] - r['bound'], status=r['status'], max_imp=r['max_imp'],
                       resolve_ok=bool(r['resolve_ok']), sec=round(r['sec'], 1))
            out['runs'].append(row); print(json.dumps(row)); sys.stdout.flush()
            json.dump(out, open('../runs/attain_%s_%s%s.json' % (ds, cf, tag), 'w'), indent=1)
