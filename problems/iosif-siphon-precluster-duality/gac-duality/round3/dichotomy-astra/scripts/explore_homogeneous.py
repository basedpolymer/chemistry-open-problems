"""A2 exploration, Python 3.14 / standard library, exact Fraction arithmetic.

Expected runtime: under 5 minutes for the default degrees 2..5 and 3 rows.
Imports earlier code unchanged. Prints PASS/FAIL for structural checks; the
certificate library's classifications are exploratory and not independently
audited here. No floating-point dynamical evidence is used.
"""
from pathlib import Path
import sys
import time
from itertools import combinations
from collections import Counter

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
sys.dont_write_bytecode = True
sys.path.insert(0, str(BASE / 'round3' / 'dichotomy' / 'scripts'))
sys.path.insert(0, str(BASE.parent / 'scripts'))
import crn
import d2cert


def from_rows(rows):
    n = len(rows)
    cols = [tuple(row[j] for row in rows) for j in range(3)]
    return crn.Network(n, [(cols[j], cols[(j+1) % 3]) for j in range(3)])


def main():
    start = time.perf_counter()
    print('Environment: Python standard library; certificates use exact rationals.')
    print('Runtime: measured at end; expected under 5 minutes.')
    degree_max = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    stats = Counter()
    examples = {}
    for d in range(2, degree_max+1):
        vectors = [(p, q, d-p-q) for p in range(d+1) for q in range(d-p+1)
                   if (p, q, d-p-q) != (d//3,)*3]
        for rows in combinations(vectors, 3):
            N = from_rows(rows)
            if crn.rank(N.Gamma) != 2:
                continue
            h = crn.nullspace(list(zip(*N.Gamma)), 3)[0]
            if not (all(v > 0 for v in h) or all(v < 0 for v in h)):
                continue
            result = d2cert.analyse(N.n, N.reactions)
            key = (d, result['verdict'])
            stats[key] += 1
            if key not in examples:
                examples[key] = rows
                print('EXAMPLE', key, rows, result['none'], flush=True)
        print('COUNTS', d, dict(stats), flush=True)
    ok = not any(key[1] == 'CONFLICT' for key in stats)
    print(('PASS' if ok else 'FAIL') + ': H0 imported certificates have no conflict.')
    print('Claim homogeneous implies ALL:', 'REFUTED' if any(k[1] == 'NONE' for k in stats) else 'UNRESOLVED')
    print('Runtime seconds:', round(time.perf_counter()-start, 3))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
