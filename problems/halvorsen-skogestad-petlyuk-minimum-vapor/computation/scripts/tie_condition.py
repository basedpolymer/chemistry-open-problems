"""Necessary condition for a tie R = V_P (from the cut identity): if peak j* is the (strict) maximum, then
R - V_P >= FR(light-only nodes) + FC(heavy-only nodes) at every feasible point, and these duties are strictly
positive whenever present (a reboiler duty is the boil-up of a section with V_B > 0; a condenser duty on a
transfer stream is the reflux L_T > 0, on a product the full top vapor).  Hence a configuration can tie with the
Petlyuk value only if, for j*:
  (a) no light product p <= j* and no light-only submixture (stream within [1, j*]) carries a reboiler, and
  (b) no heavy product p > j* and no heavy-only submixture (stream within [j*+1, N]) carries a condenser.
This script checks (a)-(b) against the ties observed in search_random.py outputs (a tie = ratio < 1 + 1e-9; trace components can make strict gaps smaller than 1e-6).
Usage: python tie_condition.py file.jsonl [...]
Environment: Python 3 (WSL).  Deterministic analysis of stored numerical results.
Runtime: seconds.
"""
import json
import sys

from hsmodel import all_configurations, config_name


def heat_locations(c):
    """Return lists of streams carrying a reboiler and a condenser."""
    N = c['N']
    reb, cond = [], []
    for s in c['streams']:
        if s == (1, N):
            continue
        if c['hx'].get(s) == 'HX':
            kind = c['producers'][s][0][0]
            (cond if kind == 'top' else reb).append(s)
    for p in range(1, N + 1):
        kinds = {t for (t, _) in c['producers'][(p, p)]}
        if kinds == {'top'}:
            cond.append((p, p))
        elif kinds == {'bot'}:
            reb.append((p, p))
    return reb, cond


def can_tie(c, js):
    reb, cond = heat_locations(c)
    return not any(s[1] <= js for s in reb) and not any(s[0] >= js + 1 for s in cond)


TOL = 1e-9

if __name__ == "__main__":
    rows = []
    for fn in sys.argv[1:]:
        rows += [json.loads(l) for l in open(fn)]
    N = len(rows[0]['alpha'])
    C = all_configurations(N)
    viol, nt, nallow, nallow_tie = 0, 0, 0, 0
    for r in rows:
        P = r['peaks']
        srt = sorted(P)
        if srt[-1] - srt[-2] < 1e-4 * srt[-1]:
            continue  # skip near-ties between peaks
        js = P.index(max(P)) + 1
        for k, v in r['ratios'].items():
            c = C[int(k)]
            allowed = can_tie(c, js)
            nallow += allowed
            if v is not None and v < 1 + TOL:
                nt += 1
                nallow_tie += allowed
                if not allowed:
                    viol += 1
                    print('VIOLATION', k, config_name(c), 'j*', js, v)
    print(f'points {len(rows)}; observed ties {nt}; ties violating the necessary condition: {viol}; '
          f'(config, point) pairs satisfying the condition: {nallow}, of which tie: {nallow_tie}')
    print('PASS' if viol == 0 else 'FAIL')
