"""val_binary.py -- validation 1: binary column, minimum reboiler vapor at purity 1-eps with n stages per section,
compared with the exact infinite-stage (Underwood) value at the same purity.
Run: python val_binary.py   (~1 min, 1 thread)."""
import numpy as np, json, time, sys
from stagenet import Network
from opt import MinR, lp_init
from refvals import binary_vmin


def binary(n, alpha, z, q):
    net = Network(alpha)
    net.node('cond', 'cooler'); net.node('feed'); net.node('reb', 'heater')
    net.section('T', n, 'cond', 'feed'); net.section('B', n, 'feed', 'reb')
    net.feed('feed', 1.0, q, z)
    net.product('A', 'cond', 'L', 0); net.product('B', 'reb', 'L', 1)
    net.build()
    return net


out = []
for (alpha, z, q) in [((2.0, 1.0), (0.5, 0.5), 1.0), ((3.0, 1.0), (0.3, 0.7), 0.5), ((1.5, 1.0), (0.6, 0.4), 0.0)]:
    for eps in (1e-2, 1e-4, 1e-6):
        ref = binary_vmin(alpha, z, q, eps)
        fprev = None
        for n in (10, 20, 40, 80, 160):
            net = binary(n, alpha, z, q)
            # continuation in n: start from the optimum at the previous n (same flow layout)
            f0 = lp_init(net, extra=[({'P:A': 1.0}, z[0])], vmin=3.0, lmin=3.0) if fprev is None else fprev
            t = time.time()
            opt = MinR(net, eps)
            r = opt.run(f0)
            row = dict(alpha=alpha, z=z, q=q, eps=eps, n=n, R=r['R'], ref_inf=ref, ratio=r['R'] / ref,
                       min_g=r['min_g'], ok=bool(r['ok']), status=r['status'], nit=r['nit'], sec=time.time() - t)
            out.append(row)
            if r['ok'] and r['min_g'] > -1e-6:
                fprev = r['f']
            print(json.dumps(row)); sys.stdout.flush()
json.dump(out, open('../runs/val_binary.json', 'w'), indent=1)
