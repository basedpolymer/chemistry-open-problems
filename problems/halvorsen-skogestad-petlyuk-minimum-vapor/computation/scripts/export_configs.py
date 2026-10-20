"""Export the enumerated configurations (Shah-Agrawal basic configurations x TC/HX choices) as JSON:
certificates/configurations_N{3,4,5}.json.  Each entry: id, basic_id, name, streams, splits [k,l] (split of
[i,j] is [i,k]/[l,j]), producers, transfer type of each single-producer submixture ('TC' or 'HX'), is_ftc.
Also prints the counts (expected 3+5=8, 18+134=152, 203+5925=6128, Shah & Agrawal 2010 Table 1) and exits
non-zero if they differ.  Environment: Python 3 (+ mpmath for hsmodel import).  Runtime < 5 s.
"""
import json
import sys
from hsmodel import all_configurations, basic_configurations, to_json

expected = {3: (3, 8), 4: (18, 152), 5: (203, 6128)}
ok = True
for N in (3, 4, 5):
    B = basic_configurations(N)
    A = all_configurations(N)
    json.dump([to_json(c) for c in A], open(f'../certificates/configurations_N{N}.json', 'w'))
    good = (len(B), len(A)) == expected[N]
    ok &= good
    print(f"{'PASS' if good else 'FAIL'}: N={N}: {len(B)} basic configurations, {len(A) - len(B)} additional "
          f"thermally coupled, {len(A)} in total (Shah-Agrawal 2010: {expected[N][0]} / {expected[N][1]})")
sys.exit(0 if ok else 1)
