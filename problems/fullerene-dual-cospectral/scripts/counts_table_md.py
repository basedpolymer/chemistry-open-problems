#!/usr/bin/env python3
"""counts_table_md.py -- print certificates/counts_per_n.csv as a compact markdown table (for result.md)."""
import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(os.path.join(ROOT, 'certificates', 'counts_per_n.csv'))))
print('| n | isomers (= A007894) | A(T) classes | isomers in them | A(T6) classes | A+D classes | 2A+D classes | T5 non-unique |')
print('|---|---|---|---|---|---|---|---|')
for r in rows:
    print('| %s | %s | %s | %s | %s | %s | %s | %s |' % (
        r['n'], r['isomers_buckygen'], r['T_cospectral_classes'], r['T_isomers_in_classes'], r['T6_cospectral_classes'],
        r['signless_A_plus_D_classes'] if r['signless_A_plus_D_classes'] != '' else 'n.s.',
        r['twoA_plus_D_classes'] if r['twoA_plus_D_classes'] != '' else 'n.s.', r['T5_nonunique_spectra_hash']))
