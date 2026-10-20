#!/usr/bin/env python3
"""check_D_configs.py -- Theorem L (../limit.md) on Agent D's exported limit configurations.

Reads every JSON file "hs-petlyuk-limit-config v1" in ../../design/runs/ (format documented in
../../design/scripts/export_config.py), rebuilds from it ONLY the operating specification and the limit net flows
(alpha, z, q, pattern, per column d, b, V_T, L_T, V_B, L_B, feed roots; product flows), and then, with Agent L's own
code (limitsys.System): recomputes every section's modes, pinch and closed ends, the junction compositions of the
limit configuration (Definition 3.1 of limit.md), and compares them with D's exported values (modes, pinch index,
node compositions, kappa's); then runs the Theorem-L battery C1-C6 of check_limit_jacobian.py on it.
PASS/FAIL per file; exit code non-zero on failure.  Evidence only.

Usage: python check_D_configs.py [glob]
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env), 1 thread.  Runtime: about 1-3 minutes.
"""
import sys
import glob
import json
import os
import mpmath as mp
import limitsys as ls
from check_limit_jacobian import theorem_L_checks

mp.mp.dps = 60


def load(path):
    with open(path) as fh:
        J = json.load(fh)
    assert J['format'] == 'hs-petlyuk-limit-config v1'
    M = mp.mpf
    N = int(J['data']['N'])
    D = ls.Design()
    D.N = N
    D.al = [M(a) for a in J['data']['alpha']]
    D.z = [M(a) for a in J['data']['z']]
    D.q = M(J['data']['q'])
    D.H = [M(a) for a in J['data']['H']]
    D.R = M(J['data']['R'])
    pat = J['params']['pattern']
    D.pattern = {p: pat[p] for p in range(N - 1)}
    D.cols = sorted([(i, m) for i in range(N) for m in range(i + 1, N)], key=lambda c: (-(c[1] - c[0]), c[0]))
    D.d, D.b, D.f, D.VT, D.VB, D.LT, D.LB, D.chi = {}, {}, {}, {}, {}, {}, {}, {}
    for (i, m) in D.cols:
        e = J['columns'][f'[{i + 1},{m + 1}]']
        D.d[(i, m)] = [M(x) for x in e['d']]
        D.b[(i, m)] = [M(x) for x in e['b']]
        D.f[(i, m)] = [M(x) for x in e['feed']]
        D.VT[(i, m)], D.LT[(i, m)] = M(e['VT']), M(e['LT'])
        D.VB[(i, m)], D.LB[(i, m)] = M(e['VB']), M(e['LB'])
        if m - i >= 2:
            D.chi[(i, m)] = [M(x) for x in e['feed_roots']]
    D.P = [M(J['products'][str(p + 1)]['P']) for p in range(N)]
    return D, J


def compare(D, J):
    """Compare L's recomputed limit configuration with D's export.  Returns (ok, worst relative difference)."""
    S = ls.System(D, None)
    worst = mp.mpf(0)
    ok = True
    for c in D.cols:
        i, m = c
        Xd = [mp.mpf(x) for x in J['nodes'][f'J[{i + 1},{m + 1}]']]
        worst = max(worst, max(abs(a - b) for a, b in zip(S.X_inf[c], Xd)))
        for kind in ('T', 'B'):
            s = (kind, i, m)
            Js = J['sections'][f'{kind}[{i + 1},{m + 1}]']
            mine = S.sec_info[s]['modes']
            theirs = Js['modes']
            if len(mine) != len(theirs) or S.sec_info[s]['ipin'] != int(Js['pinch']):
                ok = False
                continue
            for (typ, val, k), t in zip(mine, theirs):
                tt = 'maj' if t['type'] == 'major' else 'tr'
                if tt != typ or (typ == 'tr' and int(t['k']) - 1 != k):
                    ok = False
                worst = max(worst, abs(val - mp.mpf(t['value'])) / abs(val))
    # flows consistency: V_T - L_T = sum d
    for c in D.cols:
        worst = max(worst, abs(D.VT[c] - D.LT[c] - sum(D.d[c])))
    return ok and worst < mp.mpf('1e-40'), worst


def main():
    pat = sys.argv[1] if len(sys.argv) > 1 else os.path.join('..', '..', 'design', 'runs', 'config_*.json')
    files = sorted(glob.glob(pat))
    if not files:
        print('no D export found:', pat)
        sys.exit(1)
    allok = True
    for f in files:
        D, J = load(f)
        okc, worst = compare(D, J)
        head = (f'[{os.path.basename(f)}] N={D.N} q={mp.nstr(D.q, 3)} pattern={"".join(D.pattern[p] for p in range(D.N - 1))} '
                f'R-H*={mp.nstr(D.R - max(D.H), 3)}  (L vs D: modes/pinch/X agree: {"yes" if okc else "NO"}, '
                f'max diff {mp.nstr(worst, 2)})')
        okL = theorem_L_checks(D, head)
        allok &= okc and okL
        sys.stdout.flush()
    print('ALL PASS' if allok else 'SOME FAIL')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
