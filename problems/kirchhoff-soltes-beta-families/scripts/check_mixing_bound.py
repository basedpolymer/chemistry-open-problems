"""check_mixing_bound.py -- Proposition 5.2: the constant lambda(B_n) of Theorem 5.1 for the family B_n of all rooted
connected graphs of order <= n, computed in exact rational arithmetic from `gadgetbound -d`.

For a value rho^ > 0, a class O (of a rooted graph (B, x), equal (rho, delta), rho = rho^) has m = |B| - 1, g = |O|,
a = Rt_B(x) - m rho^, R = a + m rho^.  Theorem 5.1 bounds the proportion of good vertices of large single-hub splices
whose good vertices lie in branches from B_n by lambda(B_n) = max over rho^ of
   max( single: g/m if a <= 0, g rho^/R if a > 0 ;  pairs with a1 > 0 > a2: (|a2| g1 + a1 g2)/(|a2| m1 + a1 m2) ).
Since a pair value is a weighted mediant of g1/m1 and the single efficiency g2/m2 of the second class, only pairs with
g1/m1 > 2/3 can exceed 2/3; these are enumerated exactly.
Claims:  (M1) every single efficiency is <= 2/3 (re-derived from the exact class data);
         (M2) every balanced pair with equal rho has value <= 2/3, hence lambda(B_n) = 2/3 (attained by K4 - e).
Usage: python3 check_mixing_bound.py [NMAX] [binary]   (default NMAX = 8, binary = gadgetbound).
Runtime: NMAX = 8 about 1 min; NMAX = 9 about 20 min (Fedora 44 / WSL).
"""
import os, sys, subprocess
from fractions import Fraction
from collections import defaultdict
from kfexact import Checker

HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
GB = sys.argv[2] if len(sys.argv) > 2 else 'gadgetbound'
C = Checker()
two3 = Fraction(2, 3)
byrho = defaultdict(lambda: ([], []))   # rho -> (classes with a > 0 and g/m > 2/3, classes with a < 0)
best_single = Fraction(0)
nclasses = 0
for n in range(3, NMAX + 1):
    gen = subprocess.Popen(['geng', '-c', '-q', str(n)], stdout=subprocess.PIPE)
    gb = subprocess.Popen([os.path.join(HERE, GB), '-d'], stdin=gen.stdout, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, text=True)
    m = n - 1
    for line in gb.stdout:
        if not line.startswith('D '):
            continue
        p = line.split()
        g, rn, an, den = int(p[3]), int(p[6]), int(p[7]), int(p[8])
        if rn <= 0:
            continue
        nclasses += 1
        rho, a = Fraction(rn, den), Fraction(an, den)
        R = a + m * rho
        single = Fraction(g, m) if a <= 0 else g * rho / R
        best_single = max(best_single, single)
        if a > 0 and 3 * g > 2 * m:
            byrho[rho][0].append((g, m, a))
        elif a < 0:
            byrho[rho][1].append((g, m, a))
    gb.wait()
C.check(best_single == two3, '(M1) max single efficiency over %d classes with rho > 0, order <= %d: %s' % (nclasses, NMAX, best_single))
best_pair = Fraction(0)
npairs = 0
for rho, (P, Q) in byrho.items():
    for (g1, m1, a1) in P:
        for (g2, m2, a2) in Q:
            npairs += 1
            v = (-a2 * g1 + a1 * g2) / (-a2 * m1 + a1 * m2)
            best_pair = max(best_pair, v)
C.check(best_pair <= two3, '(M2) %d balanced pairs (a1 > 0 > a2, equal rho, 3 g1 > 2 m1): max value %s <= 2/3 (%.4f)'
        % (npairs, best_pair, float(best_pair)))
C.finish()
