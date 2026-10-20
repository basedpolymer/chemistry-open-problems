#!/usr/bin/env python3
"""Agent D: export the specification F(delta, eta, pi) and its limit configuration as JSON.

Usage:
  python export_config.py CASE [out.json] [dps]
CASE is one of the presets below, or "rand:N:seed[:pattern]" (random data as in check_design.py, kind generic,
delta = 0.05 H*, eta_nb found by halving from 0.1 delta, eta_bin = 0.25 x explicit bound).
Presets:
  n3_N     alpha=(4,2,1), z=1/3, q=1, delta=0.05, eta_[1,3]=0.025, eta_1=eta_2=0.003, pi=BT (network.md §11.3)
  n3_B     same data, pi=BB
  n4_1     alpha=(8,4,2,1), z=1/4, q=1, delta=0.05, eta_nb=0.005, eta_bin=0.003, pi=BBB
  n5_1     alpha=(16,8,4,2,1), z=1/5, q=0.5, delta=0.05, eta_nb=0.002, eta_bin=0.002, pi=BBBB
Environment: Python 3 + mpmath.  Runtime: seconds.

JSON format "hs-petlyuk-limit-config v1" (all numbers are decimal strings with dps significant digits;
ALL INDICES 1-BASED: components 1..N, columns [i,m], products p):
{
 "format": "hs-petlyuk-limit-config v1",
 "data": {"N", "alpha": [..], "z": [..], "q", "theta": [N-1 global feed roots], "H": [peaks], "Hstar", "delta", "R"},
 "params": {"pattern": "B..T", "eta_nb": {"[i,m]": eta}, "eta_bin": {"p": eta_p}},
 "columns": {"[i,m]": {"kind": "nb"|"bin", "feed": [N], "U", "VT", "LT", "VB", "LB", "D", "B",
             "d": [N] (limit net flow up the top section), "b": [N] (limit net flow down the bottom section),
             "eta", "feed_roots": [..], nb only: "VT0", "d0": [N], "h": [N], "ABCD": [A,B,C,D] (Astra's cone weights),
             bin only: "pi", "kappa" (margin at eta_p = 0)}},
 "sections": {"T[i,m]" | "B[i,m]": {"V", "L", "w": [N] (upward net flow), "bottom": node, "top": node,
              "modes": [{"type": "major"|"trace", "value" (eigenvalue of M/V), "k" (trace component, 1-based),
                         "right": [N] normalized pinch composition, "left": [N] left eigenvector
                         (alpha_i/(alpha_i - phi) for a major mode, e_k for a trace mode)}],   (decreasing value)
              "pinch": index into modes (0-based), "rho": max ratio of consecutive modes,
              "closed_bottom": [mode indices above the pinch], "closed_top": [mode indices below],
              "pivots": [coefficient of x_bottom, of x_top on the pinch vector] (normalized basis; > 0),
              "traces": [{"k", "end": "bottom"|"top", "node", "Psi_k" = Psi^(k)(alpha_k), "G_k" = G^(k)_{x_J}(alpha_k),
                          "kappa" = |G_k|/Psi_k}]}},
 "nodes": {"J[i,m]" | "P[p]": [N] limit liquid composition},
 "products": {"p": {"P": flow, "x": [N]}}
}
Node naming: J[i,m] = junction (feed) stage of column [i,m]; P[p] = product stage p (P[1] total condenser,
P[N] partial reboiler). Section T[i,m] runs from J[i,m] (bottom) up to J[i,m-1] or P[i]; B[i,m] from J[i+1,m]
or P[m] (bottom) up to J[i,m] (top).
"""
import sys
import json
import random
import mpmath as mp
from design_lib import build_spec, limit_config, coefficients
import check_design as cd


def s_(x, dps):
    return mp.nstr(x, dps, strip_zeros=False) if x != 0 else "0"


def nn(n):
    return f"J[{n[1]+1},{n[2]+1}]" if n[0] == 'J' else f"P[{n[1]+1}]"


def export(S, dps):
    N = S['N']
    v = lambda x: s_(x, dps)
    vec = lambda xs: [v(x) for x in xs]
    cd.check_instance(S, 'export')                # fills cone weights, pivots, traces (and checks)
    out = {'format': 'hs-petlyuk-limit-config v1',
           'data': {'N': N, 'alpha': vec(S['al']), 'z': vec(S['z']), 'q': v(S['q']), 'theta': vec(S['th']),
                    'H': vec(S['H']), 'Hstar': v(S['Hstar']), 'delta': v(S['delta']), 'R': v(S['R'])},
           'params': {'pattern': S['pattern'],
                      'eta_nb': {f"[{i+1},{m+1}]": v(c['eta']) for (i, m), c in S['cols'].items() if c['kind'] == 'nb'},
                      'eta_bin': {str(p + 1): v(S['cols'][(p, p + 1)]['eta']) for p in range(N - 1)}},
           'columns': {}, 'sections': {}, 'nodes': {}, 'products': {}}
    for (i, m), c in S['cols'].items():
        e = {'kind': c['kind'], 'feed': vec(c['f']), 'U': v(c['U']), 'VT': v(c['VT']), 'LT': v(c['LT']),
             'VB': v(c['VB']), 'LB': v(c['LB']), 'D': v(c['D']), 'B': v(c['B']), 'd': vec(c['d']), 'b': vec(c['b']),
             'eta': v(c['eta']), 'feed_roots': vec(c['roots_feed'])}
        if c['kind'] == 'nb':
            e['VT0'] = v(c['VT0'])
            e['d0'] = vec([c['d0'].get(p, 0) for p in range(N)])
            e['h'] = vec([c['h'].get(p, 0) for p in range(N)])
            e['ABCD'] = vec(c['ABCD'])
        else:
            e['pi'] = c['pi']
            e['kappa'] = v(c['kappa'])
        out['columns'][f"[{i+1},{m+1}]"] = e
    for key, s in S['secs'].items():
        ip = s['modes'].index(s['pinch'])
        modes = []
        for mo in s['modes']:
            d = {'type': mo['type'], 'value': v(mo['val']), 'right': vec(mo['right']), 'left': vec(mo['left'])}
            if mo['type'] == 'trace':
                d['k'] = mo['k'] + 1
            modes.append(d)
        out['sections'][f"{key[0]}[{key[1]+1},{key[2]+1}]"] = {
            'V': v(s['V']), 'L': v(s['L']), 'w': vec(s['w']), 'bottom': nn(s['bottom']), 'top': nn(s['top']),
            'modes': modes, 'pinch': ip, 'rho': v(s['rho']), 'closed_bottom': list(range(ip)),
            'closed_top': list(range(ip + 1, len(modes))), 'pivots': vec(s['pivots']),
            'traces': [{'k': t['k'] + 1, 'end': t['end'], 'node': nn(t['node']), 'Psi_k': v(t['Psi_k']),
                        'G_k': v(t['G_k']), 'kappa': v(t['kappa'])} for t in s['traces']]}
    for node, x in S['X'].items():
        out['nodes'][nn(node)] = vec(x)
    for p, pr in S['prods'].items():
        out['products'][str(p + 1)] = {'P': v(pr['P']), 'x': vec(pr['x'])}
    return out


PRESETS = {
    'n3_N': ([4, 2, 1], [mp.mpf(1) / 3] * 3, 1, '0.05', {(0, 2): '0.025'}, '0.003', 'BT'),
    'n3_B': ([4, 2, 1], [mp.mpf(1) / 3] * 3, 1, '0.05', {(0, 2): '0.025'}, '0.003', 'BB'),
    'n4_1': ([8, 4, 2, 1], [mp.mpf(1) / 4] * 4, 1, '0.05', '0.005', '0.003', 'BBB'),
    'n5_1': ([16, 8, 4, 2, 1], [mp.mpf(1) / 5] * 5, '0.5', '0.05', '0.002', '0.002', 'BBBB'),
}


def main():
    case = sys.argv[1]
    outp = sys.argv[2] if len(sys.argv) > 2 else f"../runs/config_{case.replace(':', '_')}.json"
    mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    if case in PRESETS:
        al, z, q, delta, enb, eb, pat = PRESETS[case]
        z = [mp.mpf(x) for x in z]
        z[-1] = 1 - mp.fsum(z[:-1])
        enb = {k: mp.mpf(x) for k, x in enb.items()} if isinstance(enb, dict) else mp.mpf(enb)
        S = limit_config(build_spec(al, z, mp.mpf(q), mp.mpf(delta), enb, mp.mpf(eb), pat))
    else:
        _, N, seed = case.split(':')[:3]
        N, seed = int(N), int(seed)
        rng = random.Random(seed)
        al, z, q = cd.rand_data(rng, N, 'generic')
        pat = case.split(':')[3] if len(case.split(':')) > 3 else 'B' * (N - 1)
        S, eta = cd.make_instance(al, z, q, mp.mpf('0.05'), pat)
    out = export(S, mp.mp.dps - 5)
    nf = sum(len(v) for v in cd.fails.values())
    with open(outp, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {outp}; checks: {'ALL PASS' if nf == 0 else str(nf) + ' FAIL'}")
    for k, vv in cd.fails.items():
        print('FAIL', k, vv[:3])
    sys.exit(0 if nf == 0 else 1)


if __name__ == '__main__':
    main()
