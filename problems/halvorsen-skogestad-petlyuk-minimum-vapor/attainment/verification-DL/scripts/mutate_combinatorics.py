#!/usr/bin/env python3
"""mutate_combinatorics.py -- V1 (independent verifier): mutation tests of ../../limit/scripts/check_combinatorics.py.

The brief asks that the exhaustive combinatorial check of limit.md be made to FAIL on deliberately broken input.
This script imports the unmodified check_combinatorics module and replaces, one at a time, the two functions that
encode the model (supports(), ends()) or the closed-end rule by broken versions.  For every mutation it runs
check(N, pattern) over all patterns for N = 3..7 and reports how many configurations are flagged.  A mutation that
is NOT detected is reported as FAIL of this script (the checker would then be blind to that error).
The unmodified functions are run first and must give zero errors.

Usage: python mutate_combinatorics.py            (standard library only; runtime < 5 s)
Verifier: Claude Opus 5.5 (effort xhigh), fresh context.
"""
import os
import sys
import itertools

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'limit', 'scripts'))
import check_combinatorics as cc  # noqa: E402

ORIG_SUPPORTS = cc.supports
ORIG_ENDS = cc.ends


def run_all(Nmax=7):
    flagged = total = 0
    example = None
    for N in range(3, Nmax + 1):
        for pat in itertools.product('BT', repeat=N - 1):
            total += 1
            e = cc.check(N, pat)
            if e:
                flagged += 1
                example = example or (N, ''.join(pat), e[0])
    return flagged, total, example


# ---- mutations -------------------------------------------------------------------------------------------------
def sup_no_impurity(N, pat):
    """products without the designed impurity (unforced product supports)"""
    sec, node = ORIG_SUPPORTS(N, pat)
    for p in range(1, N + 1):
        node[(p, p)] = {p}
    return sec, node


def sup_unforced_binary(N, pat):
    """binary columns unforced (T-np / B-np in every binary column) while the pattern bookkeeping is kept"""
    sec, node = ORIG_SUPPORTS(N, pat)
    for p in range(1, N):
        sec[('T', p, p + 1)] = ({p}, 'T-np')
        sec[('B', p, p + 1)] = ({p + 1}, 'B-np')
    return sec, node


def sup_T_full(N, pat):
    """non-binary top sections carrying the heavy key m (d_m != 0 in the limit)"""
    sec, node = ORIG_SUPPORTS(N, pat)
    for (kind, i, m), (s, t) in list(sec.items()):
        if kind == 'T' and m - i >= 2:
            sec[(kind, i, m)] = (set(range(i, m + 1)), t)
    return sec, node


def sup_B_full(N, pat):
    """non-binary bottom sections carrying the light key i (b_i != 0 in the limit)"""
    sec, node = ORIG_SUPPORTS(N, pat)
    for (kind, i, m), (s, t) in list(sec.items()):
        if kind == 'B' and m - i >= 2:
            sec[(kind, i, m)] = (set(range(i, m + 1)), t)
    return sec, node


def sup_wrong_pattern_products(N, pat):
    """product supports computed with the pattern reversed (B <-> T)"""
    flip = tuple('T' if c == 'B' else 'B' for c in pat)
    sec, _ = ORIG_SUPPORTS(N, pat)
    _, node = ORIG_SUPPORTS(N, flip)
    return sec, node


def ends_swapped(s):
    """top and bottom end of every section exchanged (equivalently: the closed-end rule reversed)"""
    b, t = ORIG_ENDS(s)
    return t, b


def ends_T_wrong(s):
    """top end of T[i,m] put at J[i+1,m] (a wrong coupling)"""
    kind, i, m = s
    if kind == 'T':
        return (i, m), (i + 1, m)
    return ORIG_ENDS(s)


def sup_unforced_consistent(N, pat):
    """unforced binary columns with consistent (pure) product supports: the singular design of limit.md Remark 7.1"""
    sec, node = ORIG_SUPPORTS(N, pat)
    for p in range(1, N):
        sec[('T', p, p + 1)] = ({p}, 'T-np')
        sec[('B', p, p + 1)] = ({p + 1}, 'B-np')
    for p in range(1, N + 1):
        node[(p, p)] = {p}
    return sec, node


def ends_B_swapped(s):
    """closed-end rule reversed for bottom sections only"""
    b, t = ORIG_ENDS(s)
    return (t, b) if s[0] == 'B' else (b, t)


def ends_T_skip(s):
    """top end of a long T[i,m] (m - i >= 3) put two levels up, at J[i,m-2] (breaks the h-monotone moves)"""
    kind, i, m = s
    if kind == 'T' and m - i >= 3:
        return (i, m), (i, m - 2)
    return ORIG_ENDS(s)


MUTATIONS = [
    ('product supports without designed impurities', sup_no_impurity, None),
    ('unforced binary columns, forced bookkeeping', sup_unforced_binary, None),
    ('non-binary T sections carry the heavy key', sup_T_full, None),
    ('non-binary B sections carry the light key', sup_B_full, None),
    ('product supports from the reversed pattern', sup_wrong_pattern_products, None),
    ('closed-end rule reversed (ends swapped)', None, ends_swapped),
    ('wrong top end of T sections', None, ends_T_wrong),
    ('unforced binary columns, pure products (singular design)', sup_unforced_consistent, None),
    ('closed-end rule reversed in B sections only', None, ends_B_swapped),
    ('long T sections skip a level (h not monotone by 1)', None, ends_T_skip),
]


def main():
    ok = True
    f, t, ex = run_all()
    print(f'unmodified model: {f}/{t} configurations flagged -> {"PASS" if f == 0 else "FAIL"}')
    ok &= f == 0
    for name, sup, en in MUTATIONS:
        cc.supports = sup or ORIG_SUPPORTS
        cc.ends = en or ORIG_ENDS
        try:
            f, t, ex = run_all()
            msg = f'{f}/{t} flagged; e.g. {ex}'
        except Exception as exc:          # a crash also counts as detection (assertion on interval supports)
            f, t, msg = 1, 1, f'raised {type(exc).__name__}: {exc}'
        det = f > 0
        ok &= det
        print(f'mutation "{name}": {msg} -> {"PASS (detected)" if det else "FAIL (not detected)"}')
    cc.supports, cc.ends = ORIG_SUPPORTS, ORIG_ENDS
    print('ALL PASS' if ok else 'SOME FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
