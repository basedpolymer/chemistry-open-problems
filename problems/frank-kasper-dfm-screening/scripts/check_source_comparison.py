#!/usr/bin/env python3
"""Compare the manual SI transcription against the inherited data and claims.

Python 3.12+, standard library only, < 1 second. This verifies consistency
of transcribed values; reading the source images remains a manual step.
If the local source PDF is present, its SHA-256 is checked as well.
"""
from pathlib import Path
from fractions import Fraction
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def main():
    si = json.loads((ROOT/'data'/'pnas_si_crosscheck.json').read_text())
    pdf = ROOT/'certificates'/'sources'/'reddy2018'/'supplements'/si['filename']
    if pdf.exists():
        if hashlib.sha256(pdf.read_bytes()).hexdigest() != si['sha256']:
            raise ValueError('source PDF hash mismatch')
        print('PASS local SI PDF SHA-256')
    else:
        print('INFO source PDF absent; checking the retained manual transcription')
    rows = list(csv.reader(line for line in (ROOT/'data'/'published_reddy2018.csv').read_text().splitlines()
                           if line and not line.startswith('#')))
    values = {row[0]: row[-1] for row in rows if row[4] == 'Unconstrained'}
    if len(si['rows']) != 12 or len(values) != 12:
        raise ValueError('expected twelve source structures')
    for row in si['rows']:
        if Fraction(row['F']) != Fraction(values[row['name']]):
            raise ValueError('source discrepancy for '+row['name'])
        print('PASS SI', row['figure'], row['name'], 'F =', row['F'])
    for claim in json.loads((ROOT/'certificates'/'claims.json').read_text()):
        if claim['name'] == 'sigma':
            continue
        if Fraction(claim['threshold']) != Fraction(values[claim['name']]):
            raise ValueError('wrong source threshold for '+claim['name'])
        print('PASS witness comparison threshold', claim['name'])


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('FAIL', str(exc))
        raise SystemExit(1)
