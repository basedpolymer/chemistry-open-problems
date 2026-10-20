"""Print the configurations closest to the Petlyuk value for each point of a search_random output file.
Usage: python inspect_rs.py file.jsonl [ntop]
Environment: Python 3 (WSL). Runtime: seconds."""
import json
import sys
from hsmodel import all_configurations, config_name

fn = sys.argv[1]
ntop = int(sys.argv[2]) if len(sys.argv) > 2 else 12
rows = [json.loads(l) for l in open(fn)]
N = len(rows[0]['alpha'])
C = all_configurations(N)
for r in rows:
    P = r['peaks']
    print('alpha', [round(a, 3) for a in r['alpha']], 'z', [round(a, 4) for a in r['z']], 'q', round(r['q'], 3),
          'peaks-(1-q)', [round(p - (1 - r['q']), 4) for p in P], 'argmax', P.index(max(P)) + 1)
    srt = sorted((v, int(k)) for k, v in r['ratios'].items() if v)
    for v, k in srt[:ntop]:
        print('   ', f'{v:.8f}', k, config_name(C[k]))
