#!/usr/bin/env python3
"""Collect the results of the generation and check logs into a table (Markdown), for the README and the paper.

Environment: Python 3 (standard library).  Runtime: < 1 s.
Usage: python3 collect_results.py <log directory>
Reads <prefix>.gen.log (generator, opt or pi mode), <prefix>.pt.log (generator in Ponitz-Tittmann mode, optional)
and <prefix>.check.log (latsaw_check) for prefixes dia_k<K> and fcc_k<K>.
"""
import glob
import os
import re
import sys


def grab(path, pattern, group=1):
    if not os.path.exists(path):
        return None
    m = None
    for line in open(path):
        mm = re.search(pattern, line)
        if mm:
            m = mm
    return m.group(group) if m else None


def main():
    d = sys.argv[1]
    rows = []
    for lat, pre in (('diamond', 'dia'), ('fcc', 'fcc')):
        ks = sorted(int(re.search(r'_k(\d+)\.check\.log$', f).group(1)) for f in glob.glob(os.path.join(d, pre + '_k*.check.log')))
        for K in ks:
            b = os.path.join(d, f'{pre}_k{K}')
            states = grab(b + '.gen.log', r'states=(\d+)')
            flt = grab(b + '.gen.log', r'lambda_float=\[([0-9.]+),')
            pt = grab(b + '.pt.log', r'lambda_float=\[([0-9.]+),')
            ratio = grab(b + '.check.log', r'exact max ratio = (\d+/\d+) = ([0-9.]+)', 1)
            verdict = grab(b + '.check.log', r'^(PASS|FAIL)  =>', 1)
            bound = grab(b + '.check.log', r'<= (\d+/\d+) =')
            n0 = grab(b + '.check.log', r'a_n = c_n for n <= (\d+)')
            secs = grab(b + '.check.log', r'time ([0-9.]+) s')
            rows.append((lat, K, states, pt, flt, bound, ratio, n0, verdict, secs))
    print('| lattice | K | states | P-T child (float) | optimised (float) | certified bound P/Q | exact max ratio | a_n = c_n up to n | checker | check time (s) |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for r in rows:
        print('| ' + ' | '.join('' if x is None else str(x) for x in r) + ' |')


if __name__ == '__main__':
    main()
