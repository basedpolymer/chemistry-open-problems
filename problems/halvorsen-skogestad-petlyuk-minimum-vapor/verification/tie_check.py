"""Checks, on the verifier's N = 4 runs (logs/nlp_4_*.json), that every observed tie (R/V_FTC < 1 + 1e-6)
satisfies the necessary tie condition of Corollary 'tie': no reboiler with positive duty at a node inside [1, j*]
and no condenser with positive duty at a node inside [j*+1, N] (j* = index of the strict highest peak).
Also checks that no configuration is below V_FTC by more than 1e-7 (the feasibility tolerance).
Environment: Python 3, pyomo (only for importing indep_nlp).  Runtime: seconds."""
import json, glob, sys
from indep_nlp import configurations, name
cfgs = {name(c): c for c in configurations(4)}
ok = True; nt = 0; nrun = 0
for f in sorted(glob.glob('logs/nlp_4_*.json')):
    d = json.load(open(f)); nrun += 1
    H = d['H']; js = max(range(3), key=lambda j: H[j]) + 1
    for r in d['rows']:
        if r['R'] is None: continue
        if r['ratio'] < 1 - 1e-7: ok = False; print('BELOW', f, r)
        if r['ratio'] >= 1 + 1e-6: continue
        nt += 1
        c = cfgs[r['name']]
        bad = []
        for s, t in c['hx'].items():
            kind, P = c['prod'][s][0]
            if t == 'HX' and kind == 'bot' and s[1] <= js: bad.append(('reboiler on light', s))
            if t == 'HX' and kind == 'top' and s[0] >= js + 1: bad.append(('condenser on heavy', s))
        for p in range(1, 5):
            pr = c['prod'][(p, p)]
            if len(pr) == 1 and pr[0][0] == 'bot' and p <= js: bad.append(('reboiler on light product', p))
            if len(pr) == 1 and pr[0][0] == 'top' and p >= js + 1: bad.append(('condenser on heavy product', p))
        if bad:
            ok = False; print('TIE VIOLATES CONDITION', f, r['name'], bad)
print(f'{nrun} runs, {nt} ties checked against the tie condition')
print('PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
