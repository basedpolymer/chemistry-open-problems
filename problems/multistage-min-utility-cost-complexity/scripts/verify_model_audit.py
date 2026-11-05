#!/usr/bin/env python3
"""Exact checks for the PARTIAL model audit, not a solution of the open problem.

Environment: Python 3.12+ standard library (tested with Windows Python 3.14).
Runtime: <1 second. Run: python scripts/verify_model_audit.py
All arithmetic is fractions.Fraction; no floating-point optimizer is used.
Each assertion is reported PASS/FAIL; any failure causes a nonzero exit.
General arguments are in PROGRESS.md. These finite checks corroborate them.

The source's (5m) index j in H is explicitly normalized to j in C and
F_i^H,F_j^C are identified with the given F_i,F_j. 'printed' retains the
printed <= signs in (5j,k); 'repaired' reverses them. 'physical' additionally
checks both zero-approach endpoint inequalities for q>0. The latter is an
explicit extension, NOT asserted to be the complete source model.
"""
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
REPORT = []


def check(label, condition):
    print(('PASS ' if condition else 'FAIL ') + label)
    REPORT.append({'claim': label, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(label)


def instance(hot, cold, stages=1):
    # Each stream is (capacity, inlet temperature, target temperature).
    return {'hot': [tuple(map(F, x)) for x in hot],
            'cold': [tuple(map(F, x)) for x in cold], 'stages': stages}


def blank(inst):
    n, m, ell = len(inst['hot']), len(inst['cold']), inst['stages']
    s = {'h': [[v[1]] * (ell + 1) for v in inst['hot']],
         'c': [[v[1]] * (ell + 1) for v in inst['cold']]}
    for key in ['fh', 'fc', 'a', 'b', 'q']:
        s[key] = [[[F(0) for _ in range(m)] for _ in range(n)] for _ in range(ell)]
    return s


def fill_idle(inst, s, stage):
    n, m = len(inst['hot']), len(inst['cold'])
    for i, (fi, _, _) in enumerate(inst['hot']):
        for j, (fj, _, _) in enumerate(inst['cold']):
            s['fh'][stage][i][j] = fi / m
            s['fc'][stage][i][j] = fj / n
            s['a'][stage][i][j] = s['h'][i][stage]
            s['b'][stage][i][j] = s['c'][j][stage + 1]


def finish(inst, s):
    s['cu'] = [fi * (s['h'][i][-1] - tout) for i, (fi, _, tout) in enumerate(inst['hot'])]
    s['hu'] = [fj * (tout - s['c'][j][0]) for j, (fj, _, tout) in enumerate(inst['cold'])]
    s['cost_unit_prices'] = sum(s['cu']) + sum(s['hu'])
    return s


def zero_heat(inst):
    s = blank(inst)
    for stage in range(inst['stages']):
        fill_idle(inst, s, stage)
    return finish(inst, s)


def violations(inst, s, mode='repaired'):
    errors = []

    def require(ok, label):
        if not ok:
            errors.append(label)

    n, m, ell = len(inst['hot']), len(inst['cold']), inst['stages']
    for i, (fi, tin, tout) in enumerate(inst['hot']):
        require(s['h'][i][0] == tin and tin >= s['h'][i][-1] >= tout, '5l')
        require(s['cu'][i] == fi * (s['h'][i][-1] - tout) and s['cu'][i] >= 0, '5b/5o')
    for j, (fj, tin, tout) in enumerate(inst['cold']):
        require(s['c'][j][-1] == tin and tout >= s['c'][j][0] >= tin, '5m (j in C)')
        require(s['hu'][j] == fj * (tout - s['c'][j][0]) and s['hu'][j] >= 0, '5c/5o')
    for t in s['h'] + s['c']:
        require(all(x >= 0 for x in t), '5n boundaries')
    for k in range(ell):
        for i, (fi, _, _) in enumerate(inst['hot']):
            require(sum(s['fh'][k][i]) == fi, '5d')
            require(fi * s['h'][i][k + 1] == sum(s['fh'][k][i][j] * s['a'][k][i][j] for j in range(m)), '5h')
            h0, h1 = s['h'][i][k:k + 2]
            require(h0 <= h1 if mode == 'printed' else h0 >= h1, '5j ' + mode)
            # The conservation identity used in the universal degeneracy proof.
            require(sum(s['q'][k][i]) == fi * (h0 - h1), 'derived hot conservation')
        for j, (fj, _, _) in enumerate(inst['cold']):
            require(sum(s['fc'][k][i][j] for i in range(n)) == fj, '5e')
            require(fj * s['c'][j][k] == sum(s['fc'][k][i][j] * s['b'][k][i][j] for i in range(n)), '5i')
            c0, c1 = s['c'][j][k:k + 2]
            require(c0 <= c1 if mode == 'printed' else c0 >= c1, '5k ' + mode)
            require(sum(s['q'][k][i][j] for i in range(n)) == fj * (c0 - c1), 'derived cold conservation')
        for i in range(n):
            for j in range(m):
                fh, fc, a, b, q = [s[key][k][i][j] for key in ['fh', 'fc', 'a', 'b', 'q']]
                require(all(v >= 0 for v in (fh, fc, a, b, q)), '5n/5o')
                require(q == fh * (s['h'][i][k] - a), '5f')
                require(q == fc * (b - s['c'][j][k + 1]), '5g')
                if mode == 'physical' and q > 0:
                    require(s['h'][i][k] >= b, 'hot inlet >= cold outlet')
                    require(a >= s['c'][j][k + 1], 'hot outlet >= cold inlet')
    return errors


def no_approach_optimum(inst):
    """Construct max-recovery solution of the sign-only repair in one stage."""
    n, m = len(inst['hot']), len(inst['cold'])
    remh = [fi * (tin - tout) for fi, tin, tout in inst['hot']]
    remc = [fj * (tout - tin) for fj, tin, tout in inst['cold']]
    s = blank(inst)
    q = s['q'][0]
    for i in range(n):
        for j in range(m):
            q[i][j] = min(remh[i], remc[j])
            remh[i] -= q[i][j]
            remc[j] -= q[i][j]
    rows = [sum(row) for row in q]
    cols = [sum(q[i][j] for i in range(n)) for j in range(m)]
    for i, (fi, tin, _) in enumerate(inst['hot']):
        s['h'][i][1:] = [tin - rows[i] / fi] * inst['stages']
    for j, (fj, tin, _) in enumerate(inst['cold']):
        s['c'][j][0] = tin + cols[j] / fj
    for i, (fi, _, _) in enumerate(inst['hot']):
        for j, (fj, _, _) in enumerate(inst['cold']):
            s['fh'][0][i][j] = fi * q[i][j] / rows[i] if rows[i] else fi / m
            s['fc'][0][i][j] = fj * q[i][j] / cols[j] if cols[j] else fj / n
            s['a'][0][i][j] = s['h'][i][1]
            s['b'][0][i][j] = s['c'][j][0]
    for k in range(1, inst['stages']):
        fill_idle(inst, s, k)
    return finish(inst, s)


def lift_one_stage_lp(inst, q):
    """Exact lift of the one-stage LP with active-exchanger approach constraints.
    Zero heat is allowed as bypass. Unused split flow is assigned deterministically.
    """
    assert inst['stages'] == 1
    n, m = len(inst['hot']), len(inst['cold'])
    q = [[F(x) for x in row] for row in q]
    s = blank(inst)
    s['q'][0] = q
    for i, (fi, hi, _) in enumerate(inst['hot']):
        for j, (fj, cj, _) in enumerate(inst['cold']):
            if q[i][j]:
                assert hi > cj
                s['fh'][0][i][j] = s['fc'][0][i][j] = q[i][j] / (hi - cj)
        slack = fi - sum(s['fh'][0][i])
        assert slack >= 0
        s['fh'][0][i][0] += slack
        s['h'][i][1] = hi - sum(q[i]) / fi
    for j, (fj, cj, _) in enumerate(inst['cold']):
        slack = fj - sum(s['fc'][0][i][j] for i in range(n))
        assert slack >= 0
        s['fc'][0][0][j] += slack
        s['c'][j][0] = cj + sum(q[i][j] for i in range(n)) / fj
    for i, (_, hi, _) in enumerate(inst['hot']):
        for j, (_, cj, _) in enumerate(inst['cold']):
            fh, fc = s['fh'][0][i][j], s['fc'][0][i][j]
            s['a'][0][i][j] = hi - q[i][j] / fh if fh else hi
            s['b'][0][i][j] = cj + q[i][j] / fc if fc else cj
    return finish(inst, s)


def exact_vertices_2d(constraints):
    """Enumerate all vertices a*x+b*y<=c using exact 2x2 determinants."""
    vertices = set()
    for (a, b, c), (d, e, f) in combinations(constraints, 2):
        det = a * e - b * d
        if not det:
            continue
        x, y = (c * e - b * f) / det, (a * f - c * d) / det
        if all(u * x + v * y <= w for u, v, w in constraints):
            vertices.add((x, y))
    return sorted(vertices)


def serialize(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {k: serialize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [serialize(v) for v in obj]
    return obj


def main():
    saved = {}
    # Exhaustive one-hot/one-cold examples include zero duties and crossings.
    count = 0
    for fi, fj, hi, hout, cin, cout, ell in product(range(1, 3), range(1, 3), range(3), range(3), range(3), range(3), (1, 3)):
        if hi < hout or cin > cout:
            continue
        inst = instance([(fi, hi, hout)], [(fj, cin, cout)], ell)
        s0, s1 = zero_heat(inst), no_approach_optimum(inst)
        assert not violations(inst, s0, 'printed')
        assert not violations(inst, s1)
        eh, ec = fi * (hi - hout), fj * (cout - cin)
        assert s0['cost_unit_prices'] == eh + ec
        assert s1['cost_unit_prices'] == abs(eh - ec)
        count += 1
    # Several-stream cases exercise zero rows/columns and dormant stages.
    for n, m, ell in product(range(1, 5), range(1, 5), (1, 2, 5)):
        hot = [(i + 1, (3 * i + m) % 7, (i + m) % 2) for i in range(n)]
        hot = [(f, max(a, b), min(a, b)) for f, a, b in hot]
        cold = [(j + 1, j % 3, j % 3 + (n + j) % 5) for j in range(m)]
        inst = instance(hot, cold, ell)
        s0, s1 = zero_heat(inst), no_approach_optimum(inst)
        assert not violations(inst, s0, 'printed')
        assert not violations(inst, s1)
        eh = sum(f * (a - b) for f, a, b in inst['hot'])
        ec = sum(f * (b - a) for f, a, b in inst['cold'])
        assert s0['cost_unit_prices'] == eh + ec
        assert s1['cost_unit_prices'] == abs(eh - ec)
        count += 1
    check(f'{count} exact examples: printed zero-transfer witness and sign-only energy optimum', True)

    inst = instance([(1, 1, 0)], [(1, 2, 3)])
    s = no_approach_optimum(inst)
    check('sign-only repair permits heat from inlet 1 to inlet 2 at zero utility', not violations(inst, s) and s['cost_unit_prices'] == 0)
    check('that witness violates both physical endpoint inequalities', violations(inst, s, 'physical') == ['hot inlet >= cold outlet', 'hot outlet >= cold inlet'])
    check('printed monotonicity rejects the same nonzero heat transfer', set(violations(inst, s, 'printed')) == {'5j printed', '5k printed'})
    saved['uphill_heat_counterexample'] = {'instance': inst, 'solution': s}

    # Direct 3-Partition gadget obstruction: all items even, target odd.
    sizes, capacity = [4, 4, 4, 4, 4, 6], 13
    check('3-Partition promise holds for (4,4,4,4,4,6), B=13', sum(sizes) == 2 * capacity and all(F(capacity, 4) < a < F(capacity, 2) for a in sizes))
    possible = [subset for subset in product((0, 1), repeat=6) if sum(a * take for a, take in zip(sizes, subset)) == capacity]
    check('all 64 subsets checked: no partition into two bins of capacity 13', not possible)
    inst = instance([(a, 1, 0) for a in sizes], [(capacity, 0, 1)] * 2)
    s = lift_one_stage_lp(inst, [[F(a, 2), F(a, 2)] for a in sizes])
    check('fractional heat splitting gives the failed packing gadget zero utility', not violations(inst, s, 'physical') and s['cost_unit_prices'] == 0)
    saved['packing_obstruction'] = {'instance': inst, 'solution': s, 'items': sizes, 'bin_capacity': capacity, 'subsets_checked': 64}

    inst1 = instance([(1, 10, 0)], [(1, 8, 10), (1, 0, 8)])
    constraints = [tuple(map(F, row)) for row in [(-1, 0, 0), (0, -1, 0), (1, 0, 2), (0, 1, 8), (1, 1, 10)]]
    constraints.append((F(1, 2), F(1, 10), F(1)))
    vertices = exact_vertices_2d(constraints)
    value = max(x + y for x, y in vertices)
    check('one-stage LP exact vertex enumeration: maximum recovery 42/5', value == F(42, 5))
    # Dual certificate: 2*(x/2+y/10<=1)+(4/5)*(y<=8).
    check('one-stage upper bound 42/5 also has an exact nonnegative dual certificate', 2 * F(1, 2) == 1 and 2 * F(1, 10) + F(4, 5) == 1 and 2 + F(4, 5) * 8 == value)
    s1 = lift_one_stage_lp(inst1, [[F(2, 5), F(8)]])
    check('one-stage optimum lifts to the bilinear constraints; utility 16/5', not violations(inst1, s1, 'physical') and s1['cost_unit_prices'] == F(16, 5))
    inst2 = instance(inst1['hot'], inst1['cold'], stages=2)
    s2 = blank(inst2)
    s2.update({'h': [[F(10), F(8), F(0)]], 'c': [[F(10), F(8), F(8)], [F(8), F(8), F(0)]],
               'fh': [[[F(1), F(0)]], [[F(0), F(1)]]],
               'fc': [[[F(1), F(1)]], [[F(1), F(1)]]],
               'a': [[[F(8), F(10)]], [[F(8), F(0)]]],
               'b': [[[F(10), F(8)]], [[F(8), F(8)]]],
               'q': [[[F(2), F(0)]], [[F(0), F(8)]]]})
    finish(inst2, s2)
    check('two-stage rational witness satisfies all repaired balances and active approaches', not violations(inst2, s2, 'physical'))
    check('two stages attain zero utility, strictly better than one stage', s2['cost_unit_prices'] == 0 < s1['cost_unit_prices'])
    saved['stage_separation'] = {'one_stage': {'instance': inst1, 'solution': s1, 'lp_vertices': vertices},
                                'two_stage': {'instance': inst2, 'solution': s2}}
    out = ROOT / 'certificates/exact-checks.json'
    out.write_text(json.dumps(serialize({'arithmetic': 'exact rational, fractions.Fraction',
                                       'scope': 'model audit and examples; NOT a resolution of the open problem',
                                       'checks': REPORT, 'witnesses': saved}), indent=2) + '\n', encoding='utf-8')
    print(f'PASS all {len(REPORT)} claims; wrote {out.name}')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('FAIL', str(exc), file=sys.stderr)
        raise
