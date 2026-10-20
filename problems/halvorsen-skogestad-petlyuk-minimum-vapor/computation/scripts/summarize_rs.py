"""Summarise search_random.py output files: number of points, smallest V_c/V_P over all (point, configuration),
and per configuration the smallest ratio and the number of points where it ties with the Petlyuk value.
Usage: python summarize_rs.py out.json file1.jsonl [file2.jsonl ...]
Environment: Python 3 (WSL).  Floating point summary of numerical evidence.
Runtime: seconds.
"""
import json
import sys
from collections import defaultdict

from hsmodel import all_configurations, config_name

out = sys.argv[1]
rows = []
for fn in sys.argv[2:]:
    rows += [json.loads(l) for l in open(fn)]
N = len(rows[0]['alpha'])
C = all_configurations(N)
minr = defaultdict(lambda: float('inf'))
ties = defaultdict(int)
ties_by_peak = defaultdict(lambda: defaultdict(int))
global_min = float('inf')
nfail = 0
for r in rows:
    P = r['peaks']
    jstar = P.index(max(P)) + 1
    for k, v in r['ratios'].items():
        if v is None:
            nfail += 1
            continue
        k = int(k)
        minr[k] = min(minr[k], v)
        global_min = min(global_min, v)
        if v < 1 + 1e-6:
            ties[k] += 1
            ties_by_peak[k][jstar] += 1
summary = {'files': sys.argv[2:], 'N': N, 'points': len(rows), 'global_min_ratio': global_min,
           'failed_solves': nfail,
           'q_range': [min(r['q'] for r in rows), max(r['q'] for r in rows)],
           'configs': [{'id': k, 'name': config_name(C[k]), 'min_ratio': minr[k], 'ties': ties[k],
                        'ties_by_dominant_peak': dict(ties_by_peak[k])} for k in sorted(minr, key=lambda x: minr[x])]}
json.dump(summary, open(out, 'w'), indent=1)
print('points', len(rows), 'global min ratio', global_min, 'failed solves', nfail)
print('configs that tie at least once:', sum(1 for k in ties if ties[k] > 0), 'of', len(C))
never = [k for k in minr if ties[k] == 0]
print('closest non-tying configs:')
for k in sorted(never, key=lambda x: minr[x])[:10]:
    print('  ', f'{minr[k]:.6f}', k, config_name(C[k]))
