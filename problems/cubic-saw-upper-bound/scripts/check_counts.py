#!/usr/bin/env python3
"""Compare the automaton path counts a_n printed in a generator log ("n=.. a_n=..") with the exact SAW counts c_n
(Clisby-Liang-Slade tables, file with lines "n c_n").  Checks a_n >= c_n for every n present in both (necessary
consequence of Lemma 1) and reports the largest n0 such that a_n = c_n for all n <= n0.
Usage: check_counts.py <gen.log> <saw_dX.txt>.  Python 3 standard library.  Runtime < 1 s.
Prints PASS/FAIL; exit code 0 iff a_n >= c_n for all compared n."""
import re
import sys

log, counts = sys.argv[1], sys.argv[2]
a = {int(m.group(1)): int(m.group(2)) for m in re.finditer(r"n=\s*(\d+) a_n=(\d+)", open(log).read())}
c = {}
for line in open(counts):
    parts = line.split()
    if len(parts) >= 2 and parts[0].isdigit():
        c[int(parts[0])] = int(parts[1])
common = sorted(set(a) & set(c))
bad = [n for n in common if a[n] < c[n]]
n0 = 0
for n in common:
    if a[n] == c[n]:
        n0 = n
    else:
        break
print("%s a_n >= c_n for n = %d..%d (%d values); a_n = c_n exactly for n <= %d; a_%d/c_%d = %.6f"
      % ("PASS" if not bad else "FAIL", common[0], common[-1], len(common), n0, common[-1], common[-1],
         a[common[-1]] / c[common[-1]]))
sys.exit(0 if not bad else 1)
