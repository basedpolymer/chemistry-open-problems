"""Summarise cmp_scip_ipopt.py outputs: SCIP status counts, agreement between the ipopt multistart value and the
SCIP primal/dual values, and the smallest SCIP dual bound relative to the Petlyuk value V_P (a dual bound >= V_P
certifies, up to SCIP's floating-point tolerances, that the configuration cannot beat the FTC at that point).
Usage: python summarize_scip.py file.json [...]   Environment: Python 3.  Runtime: seconds."""
import json
import sys
from collections import Counter

for fn in sys.argv[1:]:
    rows = json.load(open(fn))
    st = Counter(r['scip']['status'] for r in rows)
    vp = rows[0]['VP']
    diffs, missed, below, mind = [], [], [], float('inf')
    for r in rows:
        s = r['scip']
        if s['dual'] is not None:
            mind = min(mind, s['dual'] / vp)
        if r['ipopt'] is not None and s['primal'] is not None:
            diffs.append(abs(r['ipopt'] - s['primal']) / s['primal'])
            if r['ipopt'] > s['primal'] * (1 + 1e-6):
                missed.append((r['id'], r['ipopt'], s['primal']))
        if r['ipopt'] is not None and s['dual'] is not None and r['ipopt'] < s['dual'] * (1 - 1e-6):
            below.append((r['id'], r['ipopt'], s['dual']))
    nfail = sum(1 for r in rows if r['ipopt'] is None)
    print(f"{fn}: {len(rows)} configurations; SCIP status {dict(st)}; ipopt failures {nfail}; "
          f"max rel |ipopt - SCIP primal| = {max(diffs):.2e}; ipopt above SCIP primal by >1e-6: {missed}; "
          f"ipopt below SCIP dual by >1e-6: {below}; min SCIP dual / V_P = {mind:.9f}")
