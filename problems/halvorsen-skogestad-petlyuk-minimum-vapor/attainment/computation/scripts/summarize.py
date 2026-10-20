#!/usr/bin/env python3
"""summarize.py -- Agent C.  Collects the JSON outputs of run_case.py (../runs/<case>.json) into one Markdown table
(../runs/summary.md) and prints it.  No computation.  Usage: python summarize.py
Environment: Python 3.12 (standard library only).  Runtime: < 1 s."""
import os
import json
import glob

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, '..', 'runs')
rows = []
for f in sorted(glob.glob(os.path.join(RUNS, '*.json'))):
    try:
        J = json.load(open(f))
    except Exception:
        continue
    if 'rows' not in J or 'case' not in J:
        continue
    spec = J['spec']
    N = len(spec['alpha'])
    pat = ''.join(v[0] if isinstance(v, (list, tuple)) else '?' for k, v in sorted(spec.get('eta_bin', {}).items()))
    rws = sorted(J['rows'], key=lambda r: r['n'])
    last = rws[-1]
    phys = [r['n'] for r in rws if r.get('physical', True)]
    blocks = J.get('limit_blocks', {})
    mm = J.get('mmatrix', {})
    mm_ok = all(v.get('ok') for v in mm.values()) if mm else None
    imp = last['impurities']
    rows.append(dict(case=J['case'], N=N, pattern=pat, q=f"{float(spec['q']):.4g}", delta=f"{float(spec['delta']):.4g}", rho=J.get('rho'),
                     smin=J['limit_svals'][-1] if 'limit_svals' in J else None,
                     maj=blocks.get('major', [None] * 3), tr=blocks.get('trace', [None] * 3),
                     mm=mm_ok, n0=J.get('n0_observed'), nmax=last['n'], dps=last['dps'], imp=imp,
                     dist=last['dist_major'], smin_n=last['sigma_min'], fails=J.get('fails', []),
                     rate=J.get('rate_major'), minS=min(r['minS'] for r in rws if r.get('physical', True)) if phys else None))
lines = ["| case | N | pattern | q | δ | ρ | σ_min(DG_∞) | majors / traces | A_k M-matrices | n0 obs. | largest n (digits) | impurities at largest n | ‖v_n − v_∞‖ | rate | verdict |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    imp = ", ".join(f"{float(x):.1e}" for x in r['imp'])
    lines.append(f"| {r['case']} | {r['N']} | {r['pattern']} | {r['q']} | {r['delta']} | "
                 f"{r['rho']:.4f} | {r['smin']:.3g} | {r['maj'][0]} / {r['tr'][0]} | "
                 f"{'yes' if r['mm'] else ('n/a' if r['mm'] is None else 'NO')} | {r['n0']} | {r['nmax']} ({r['dps']}) | "
                 f"{imp} | {r['dist']:.1e} | {(r['rate'] or float('nan')):.4f} | "
                 f"{'PASS' if not r['fails'] else 'FAIL: ' + ', '.join(r['fails'])} |")
txt = "\n".join(lines) + "\n"
with open(os.path.join(RUNS, 'summary.md'), 'w', encoding='utf-8') as fh:
    fh.write(txt)
print(txt)
print(f"{len(rows)} cases, {sum(1 for r in rows if not r['fails'])} PASS")
