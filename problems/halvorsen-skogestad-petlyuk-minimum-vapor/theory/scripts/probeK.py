"""Exploration (floating point, not a certificate): config K variants (N=4) vs FTC value.
Environment: WSL Fedora, pyomo + ipopt (Agent C's nlp.py).  Runtime: about 10 minutes for 15 trials."""
import sys, random, math
sys.path.insert(0, '../../computation/scripts')
from hsmodel import all_configurations, config_name, petlyuk_value, petlyuk_peaks
from nlp import solve_config
A = all_configurations(4)
ids = list(range(88, 104))
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
worst = {}
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 20):
    # random alpha, z, q with P_2 the highest peak preferred
    r = [math.exp(rng.uniform(math.log(1.05), math.log(4))) for _ in range(3)]
    alpha = [r[0]*r[1]*r[2], r[1]*r[2], r[2], 1.0]
    z = [rng.uniform(0.02, 1) for _ in range(4)]; s = sum(z); z = [x/s for x in z]
    q = rng.uniform(0, 1)
    P = [float(x) for x in petlyuk_peaks(alpha, z, q)]
    vt, vr = petlyuk_value(alpha, z, q); vr = float(vr)
    for cid in ids:
        best, res = solve_config(A[cid], alpha, z, q, starts=8, seed=trial)
        if best is None: continue
        gap = (best[0] - vr) / vr
        if cid not in worst or gap < worst[cid][0]:
            worst[cid] = (gap, alpha, z, q, P)
for cid in ids:
    if cid in worst:
        g = worst[cid]
        print(cid, f"{g[0]:+.5f}", config_name(A[cid]).split('|')[1], 'peaks', [round(x,4) for x in g[4]])
