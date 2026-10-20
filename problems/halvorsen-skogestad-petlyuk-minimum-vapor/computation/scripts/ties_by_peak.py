"""List, per dominant peak j*, the configurations observed to tie with the Petlyuk value (ratio < 1 + 1e-9) in
search_random.py outputs and satisfying the necessary tie condition of tie_condition.py.
Usage: python ties_by_peak.py out.json file.jsonl [...]   Environment: Python 3 (WSL).
Runtime: seconds."""
import json
import sys
from collections import defaultdict
from hsmodel import all_configurations, config_name
from tie_condition import can_tie

rows = []
for fn in sys.argv[2:]:
    rows += [json.loads(l) for l in open(fn)]
N = len(rows[0]['alpha'])
C = all_configurations(N)
ties = defaultdict(lambda: defaultdict(int))
npts = defaultdict(int)
for r in rows:
    P = r['peaks']
    s = sorted(P)
    if s[-1] - s[-2] < 1e-3 * s[-1]:
        continue
    js = P.index(max(P)) + 1
    npts[js] += 1
    for k, v in r['ratios'].items():
        if v is not None and v < 1 + 1e-9 and can_tie(C[int(k)], js):
            ties[js][int(k)] += 1
out = {}
for js in sorted(ties):
    lst = sorted(ties[js].items(), key=lambda x: -x[1])
    out[js] = {'points_with_this_dominant_peak': npts[js],
               'configs': [{'id': k, 'name': config_name(C[k]), 'tie_count': n} for k, n in lst]}
    print(f'j*={js}: {npts[js]} points, {len(lst)} configurations tie at least once')
    for k, n in lst:
        print(f'   {n:4d}  {k:4d}  {config_name(C[k])}')
json.dump(out, open(sys.argv[1], 'w'), indent=1)
