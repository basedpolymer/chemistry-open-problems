#!/usr/bin/env python3
import sys, time
sys.path.insert(0, '.')
from x6_instance import *
b = int(sys.argv[1]); a0 = int(sys.argv[2]); a1 = int(sys.argv[3]); mfac = float(sys.argv[4])
res = []
for a in range(a0, a1+1):
    s = solve(a, b, int(mfac*a))
    if s:
        n = a*(b+1)+s[0]
        res.append((n, a, s)); print(n, a, s, flush=True)
