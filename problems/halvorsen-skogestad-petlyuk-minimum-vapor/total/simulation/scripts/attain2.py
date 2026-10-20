"""
attain2.py -- purity/stage study of the minimum total vaporization of an arrangement in the MESH model, with
continuation (reduced-space SQP of opt.py: every iterate is a converged steady state, exact adjoint gradients).

Agent P (Claude Opus 5.5, effort xhigh).
Usage: python attain2.py DATASET CONFIG "n1,n2,..." "eps1,eps2,..." [tag]
  DATASET, CONFIG: see attain.py (T61..., P2, E8, ..., CEX5; FTC, K, C161, direct, indirect).
For each eps (in the given order) and each n (ascending), minimizes R at purity >= 1-eps for every product,
starting from the optimum at the previous n (or previous eps).  Records R, R/V_FTC, the stage-level bound
max_j sum_P max(0, phi_P(theta_j)) at the optimum, the impurities of the converged steady state.
Output: ../runs/attain2_DATASET_CONFIG[_tag].json (+ the optimal flows of every run).
"""
import sys, json, time
import numpy as np
from arrangements import regular_network, init_extra
from opt import MinR, lp_init
from study import stage_bound
from refvals import peaks
from attain import DATA, cex5, config

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
    print('V_FTC', VF, 'H', H, flush=True)
    f_eps = None
    for eps in epss:
        f_prev = f_eps
        first = True
        for n in ns:
            net = regular_network(alpha, z, q, cols, tr, n=n)
            if f_prev is None:
                f_prev = lp_init(net, extra=init_extra(net, N, z, cols), vmin=0.8 + max(0.0, 1 - q), lmin=0.5)
            t = time.time()
            opt = MinR(net, eps)
            try:
                r = opt.run(f_prev, maxiter=400)
            except Exception as ex:
                print(json.dumps(dict(n=n, eps=eps, error=repr(ex))), flush=True)
                continue
            if r['x'] is None:
                print(json.dumps(dict(n=n, eps=eps, error='no feasible point', status=r['status'])), flush=True)
                continue
            b, per = stage_bound(net, r['f'], r['x'], alpha, z, q)
            imp = net.impurities(r['x'])
            row = dict(n=n, eps=eps, R=r['R'], ratio=r['R'] / VF, bound=b, bound_per_j=per, R_minus_bound=r['R'] - b,
                       min_g=r['min_g'], max_imp=max(imp.values()), status=int(r['status']), nit=int(r['nit']),
                       nfev=int(r['nfev']), sec=round(time.time() - t, 1), flows=r['f'].tolist())
            out['runs'].append(row)
            print(json.dumps({k: v for k, v in row.items() if k != 'flows'}), flush=True)
            if r['min_g'] > -1e-6:
                f_prev = r['f']
                if first:
                    f_eps = r['f']; first = False
            json.dump(out, open('../runs/attain2_%s_%s%s.json' % (ds, cf, tag), 'w'), indent=1)
    out['fnames'] = net.fnames
    json.dump(out, open('../runs/attain2_%s_%s%s.json' % (ds, cf, tag), 'w'), indent=1)
