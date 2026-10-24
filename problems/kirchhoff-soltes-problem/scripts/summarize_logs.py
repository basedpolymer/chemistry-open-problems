#!/usr/bin/env python3
"""summarize_logs.py -- aggregates certificates/search/*.log into the tables of result.md and checks the
counts against OEIS (A002218: 2-connected graphs; A006799: vertex-transitive graphs).
Prints markdown tables and PASS/FAIL for: counts match OEIS, no survivor except C5, no inconclusive.
"""
import glob
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(here, '..', 'certificates', 'search')
A002218 = {3: 1, 4: 3, 5: 10, 6: 56, 7: 468, 8: 7123, 9: 194066, 10: 9743542, 11: 900969091}
A006799 = {10: 22, 11: 8, 12: 74, 13: 14, 14: 56, 15: 48, 16: 286, 17: 36, 18: 380, 19: 60, 20: 1214,
           21: 240, 22: 816, 23: 188, 24: 15506, 25: 464, 26: 4236, 27: 1434, 28: 25850, 29: 1182,
           30: 46308, 31: 2192, 32: 677402, 33: 6768, 34: 132580, 35: 11150, 36: 1963202, 37: 14602,
           38: 814216, 39: 48462, 40: 13104170, 41: 52488, 42: 9462226, 43: 99880, 44: 39134640,
           45: 399420, 46: 34333800, 47: 364724}
pat = re.compile(r'graphs_read (\d+) unparsable (\d+) n<3 (\d+) disconnected (\d+) rejected (\d+) '
                 r'\(rejected_by_cut_vertex (\d+)\) survivors (\d+) inconclusive (\d+)')
ok = True


def parse(path):
    tot = [0] * 8
    for line in open(path):
        m = pat.search(line)
        if m:
            for i in range(8):
                tot[i] += int(m.group(i + 1))
    return tot


print('### 2-connected graphs (geng -C)\n\n| n | graphs | rejected | survivors | inconclusive | OEIS A002218 |\n|---|---|---|---|---|---|')
total = 0
for n in range(3, 12):
    p = os.path.join(S, f'geng_C_n{n}.log')
    if not os.path.exists(p):
        print(f'| {n} | (missing) | | | | {A002218[n]} |')
        continue
    t = parse(p)
    total += t[0]
    good = t[0] == A002218[n] and t[7] == 0 and (t[6] == (1 if n == 5 else 0))
    ok &= good
    print(f'| {n} | {t[0]:,} | {t[4]:,} | {t[6]} | {t[7]} | {A002218[n]:,} {"ok" if t[0] == A002218[n] else "MISMATCH"} |')
print(f'\ntotal 2-connected graphs tested: {total:,}\n')

print('### Holt-Royle census (vertex-transitive, all valencies)\n\n| n | graphs | disconnected | rejected | survivors | inconclusive | OEIS A006799 |\n|---|---|---|---|---|---|---|')
vt_total = vt_conn = 0
for n in range(10, 48):
    p = os.path.join(S, f'census_n{n}.log')
    if not os.path.exists(p):
        print(f'| {n} | (missing) | | | | | {A006799[n]} |')
        continue
    t = parse(p)
    vt_total += t[0]
    vt_conn += t[0] - t[3] - t[2]
    good = t[0] == A006799[n] and t[6] == 0 and t[7] == 0
    ok &= good
    print(f'| {n} | {t[0]:,} | {t[3]:,} | {t[4]:,} | {t[6]} | {t[7]} | {A006799[n]:,} {"ok" if t[0] == A006799[n] else "MISMATCH"} |')
print(f'\ncensus graphs read: {vt_total:,}; connected (tested): {vt_conn:,}\n')

print('### Kernel enumeration (Theorem 4)\n\n| c | graphs | rejected | survivors |\n|---|---|---|---|')
for c in range(2, 7):
    p = os.path.join(S, f'kernel_c{c}.log')
    if os.path.exists(p):
        t = parse(p)
        ok &= t[6] == 0 and t[7] == 0
        print(f'| {c} | {t[0]:,} | {t[4]:,} | {t[6]} |')

print('\n### Other exact searches\n\n| search | graphs | disconnected | rejected | survivors | inconclusive |\n|---|---|---|---|---|---|')
for p in sorted(glob.glob(os.path.join(S, 'regular_*.log')) + glob.glob(os.path.join(S, 'families_small_*.log'))):
    if p.endswith('.gen.log'):
        continue
    t = parse(p)
    ok &= t[6] == 0 and t[7] == 0
    print(f'| {os.path.basename(p)[:-4]} | {t[0]:,} | {t[3]:,} | {t[4]:,} | {t[6]} | {t[7]} |')
print('\nPASS search logs: counts match OEIS where applicable, only survivor C5' if ok else '\nFAIL search logs')
sys.exit(0 if ok else 1)
