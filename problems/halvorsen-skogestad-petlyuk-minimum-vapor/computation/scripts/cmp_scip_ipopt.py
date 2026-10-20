"""Compare SCIP global optimum and ipopt multistart for all configurations at one parameter point.
Usage: python cmp_scip_ipopt.py N alpha(comma) z(comma) q outfile.json [timelimit]
Environment: WSL Fedora, Python 3, pyomo+ipopt, pyscipopt; 1 thread.
Runtime: about 30-60 min for the 152 configurations of N = 4, single thread."""
import json
import sys
import time
from hsmodel import all_configurations, config_name, petlyuk_value
from nlp import solve_config
from global_scip import scip_min

N = int(sys.argv[1])
al = [float(x) for x in sys.argv[2].split(',')]
z = [float(x) for x in sys.argv[3].split(',')]
q = float(sys.argv[4])
out = sys.argv[5]
tl = float(sys.argv[6]) if len(sys.argv) > 6 else 300
perm = len(sys.argv) > 7 and sys.argv[7] == 'permissive'
vt, vp = petlyuk_value(al, z, q)
vt, vp = float(vt), float(vp)
rows = []
for c in all_configurations(N):
    t0 = time.time()
    best, _ = solve_config(c, al, z, q, starts=6, seed=c['id'], permissive=perm)
    ti = time.time() - t0
    r = scip_min(c, al, z, q, vmax=20 * vt, timelimit=tl, permissive=perm)
    rows.append({'id': c['id'], 'name': config_name(c), 'ipopt': best[0] if best else None, 'ipopt_time': ti,
                 'scip': r, 'VP': vp})
    print(c['id'], best[0] if best else None, r['primal'], r['dual'], r['status'], round(r['time'], 1),
          flush=True)
    json.dump(rows, open(out, 'w'), indent=1)
