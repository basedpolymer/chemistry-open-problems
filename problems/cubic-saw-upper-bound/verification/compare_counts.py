#!/usr/bin/env python3
"""Compare a_n (paths of the automaton actually checked by sawcheck, printed with the 5th argument) with the exact
numbers c_n of self-avoiding walks (Clisby-Liang-Slade tables).  Lemma 1 predicts a_n >= c_n for every n.
Usage: compare_counts.py <sawcheck output> <saw_dX.txt>.  Prints the largest n with a_n = c_n and PASS/FAIL."""
import sys
lines = open(sys.argv[1]).read().split('a_n (paths')[1].splitlines()[1:]
a = {int(x.split()[0]): int(x.split()[1]) for x in lines if x.strip() and x.split()[0].isdigit()}
c = {}
for x in open(sys.argv[2]):
    p = x.split()
    if len(p) >= 2 and p[0].isdigit():
        c[int(p[0])] = int(p[1])
common = sorted(set(a) & set(c))
bad = [n for n in common if a[n] < c[n]]
eq = max([n for n in common if all(a[m] == c[m] for m in common if m <= n)], default=-1)
print('n <= %d compared; a_n = c_n for n <= %d; first excess a_n/c_n - 1 at n=%d: %.3e; at n=%d: %.4f; %s'
      % (common[-1], eq, eq + 1, a[eq + 1] / c[eq + 1] - 1 if eq + 1 in a else 0, common[-1],
         a[common[-1]] / c[common[-1]], 'PASS (a_n >= c_n)' if not bad else 'FAIL at n=%s' % bad))
sys.exit(1 if bad else 0)
