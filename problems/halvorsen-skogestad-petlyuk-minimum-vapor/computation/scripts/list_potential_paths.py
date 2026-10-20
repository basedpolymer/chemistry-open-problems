"""List, for N = 3..6 (full lists for N <= 5, counts for N = 6), the basic configurations admitting a sign-consistent light->heavy path in the sign graph
(see badpath_search.consistent), in the strict and in the permissive (TG literal) model.
Environment: Python 3 (WSL).  Pure combinatorics, deterministic.  Runtime ~1 min for N = 6.
Usage: python list_potential_paths.py outfile.json"""
import json
import sys
import badpath_search as bp
from hsmodel import basic_configurations, config_name

out = {}
for strict in (True, False):
    bp.STRICT = strict
    key = 'strict' if strict else 'permissive'
    out[key] = {}
    for N in (3, 4, 5, 6):
        rows = []
        for b, c in enumerate(basic_configurations(N)):
            c['hx'] = {s: 'TC' for s in c['streams'] if s != (1, N) and len(c['producers'][s]) == 1}
            c['id'] = b
            for j in range(1, N):
                for p in bp.potential_paths(c, j):
                    rows.append({'basic_id': b, 'name': config_name(c).split(' | ')[0], 'j': j,
                                 'reqs': [[list(x[0]), x[1], x[2]] for x in p]})
        if N <= 5:
            out[key][N] = rows
        else:  # N = 6: store only the number of paths per (basic configuration, j) to keep the file small
            cnt = {}
            for r in rows:
                cnt[f"{r['basic_id']},{r['j']}"] = cnt.get(f"{r['basic_id']},{r['j']}", 0) + 1
            out[key][N] = {'n_paths': len(rows), 'per_basic_j': cnt}
        print(key, N, 'paths', len(rows), 'configs', len(set(r['basic_id'] for r in rows)), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
