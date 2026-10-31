#!/usr/bin/env python3
"""t6_summary.py LOGDIR -- collects the T^6 scan (dualspec mask 2, K = 16, mod 2^16) and its exact settlement into
../certificates/t6_scan.json.  PASS/FAIL:
  S1  isomer counts agree with ../../certificates/counts_per_n.csv for n <= 138 (and C150 = 14059174, ../../result.md);
  S2  for n <= NMAX every hash collision is settled; the exact T^6-cospectral classes are at n = 32, 36, 40, 52 only.
"""
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fdc import ROOT, HERE

logdir = os.path.expanduser(sys.argv[1])
known = {}
with open(os.path.join(ROOT, 'certificates', 'counts_per_n.csv')) as f:
    for r in csv.DictReader(f):
        known[int(r['n'])] = int(r['isomers_buckygen'])
known[150] = 14059174
rows = []
ok1 = ok2 = True
exact = []
for n in [20] + list(range(24, 400, 2)):
    p = os.path.join(logdir, 'n%d.txt' % n)
    if not os.path.exists(os.path.join(logdir, 'n%d.done' % n)):
        continue
    text = open(p).read()
    m = re.search(r'graphs=(\d+) distinct_keys=(\d+) keys_with_multiple_graphs=(\d+)', text)
    graphs, keys, multi = map(int, m.groups())
    if n in known:
        ok1 &= graphs == known[n]
    settled = os.path.join(logdir, 'n%d.settled.json' % n)
    rec = dict(n=n, isomers=graphs, distinct_keys=keys, keys_with_collisions=multi)
    if multi:
        if not os.path.exists(settled):
            ok2 = False
            rec['settled'] = False
        else:
            s = json.load(open(settled))
            rec['settled'] = True
            rec['exact_classes'] = len(s['exact_cospectral_classes'])
            rec['first_differing_coefficient'] = s.get('first_diff')
            if s['exact_cospectral_classes']:
                exact.append(n)
    rows.append(rec)
ok2 &= exact == [32, 36, 40, 52]
out = os.path.join(os.path.dirname(HERE), 'certificates', 't6_scan.json')
json.dump(rows, open(out, 'w'), indent=0)
print('PASS' if ok1 else 'FAIL', ' S1 isomer counts')
print('PASS' if ok2 else 'FAIL', ' S2 all collisions settled; exact T^6 classes only at', exact, '; n from', rows[0]['n'], 'to', rows[-1]['n'])
sys.exit(0 if ok1 and ok2 else 1)
