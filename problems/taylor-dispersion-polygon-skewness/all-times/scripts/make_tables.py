#!/usr/bin/env python3
"""make_tables.py: LaTeX table of the finite-n certificates from certificates/finite_4_200.json; prints summary statistics."""
import json, os, math

here = os.path.dirname(os.path.abspath(__file__))
res = json.load(open(os.path.join(here, '..', 'certificates', 'finite_4_200.json')))
by = {r['n']: r for r in res}
assert all(r['ok'] for r in res)
print('n range', min(by), max(by), 'count', len(res))
print('min margin', min(r['worst'][0] for r in res), 'at n =', min(res, key=lambda r: r['worst'][0])['n'])
print('min S', min(r['short']['S'] for r in res), 'at n =', min(res, key=lambda r: r['short']['S'])['n'])
print('worst point always tau0:', all(abs(r['worst'][1] - 0.015) < 1e-9 for r in res))
print('min slope - e1', min(r['slope'] - r['e1'] for r in res))
print('max Delta0', max(r['Delta0'] for r in res))


def sci(x, d=1):
    mant, exp = f"{x:.{d}e}".split('e')
    return f"{mant}\cdot10^{{{int(exp)}}}"


rows = [4, 5, 6, 8, 12, 16, 24, 32, 48, 64, 100, 150, 199]
out = [r"\begin{tabular}{rrrrrrrrr}", r"\toprule",
       r"\(n\) & \(M\) & \(\omega\) & \(\rho\) & \(L_1\) & \(\Delta_0\) & \(S(\tau_0)\) & \(\tilde K(\tau_0)\) & \(\mathcal E\)\\",
       r"\midrule"]
for n in rows:
    r = by[n]
    out.append(f"{n} & {r['M']} & \\({sci(r['omega_L'])}\\) & \\({sci(r['rho'])}\\) & \\({sci(r['L1'])}\\) & \\({sci(r['Delta0'])}\\) & "
               f"\\({sci(r['short']['S'],2)}\\) & \\({sci(r['worst'][2],2)}\\) & \\({sci(r['worst'][3],2)}\\)\\\\")
out += [r"\bottomrule", r"\end{tabular}"]
open(os.path.join(here, '..', 'table-finite.tex'), 'w').write("\n".join(out) + "\n")
