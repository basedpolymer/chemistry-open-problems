"""R2: recompute, in exact arithmetic, the share of ALL verdicts in Q1's eight rank-one families.

Data transcribed from gac-duality/round2/persistence/notes.md (table "Counts (checked, logs/enum_s1_*.log ...)"),
whose numbers V3 recounted (verification/V3/verification.md, "rank-one table (8 families) CONFIRMED").
The paper (Sec. 6.1, "Counts and credits") states "persistent for every rate vector for 17--35 %".
Environment: any Python 3 (standard library). Runtime: < 1 s.
"""
import sys
from fractions import Fraction
rows = {  # (n, m, R): (networks, ALL, NONE)
    (2, 2, 4): (12, 4, 8), (2, 3, 6): (460, 115, 345), (2, 4, 6): (3544, 381, 3163), (2, 5, 5): (2742, 662, 2080),
    (3, 2, 6): (12, 4, 8), (3, 3, 6): (426, 101, 325), (3, 4, 5): (829, 245, 584), (4, 3, 4): (1, 1, 0)}
fails = 0
for fam, (n, a, z) in rows.items():
    ok = a + z == n
    fails += not ok
    print(f"{'PASS' if ok else 'FAIL'} {fam}: ALL + NONE = networks ({a} + {z} = {n}); ALL share {a}/{n} = {float(Fraction(a, n)):.4f}")
shares = [Fraction(a, n) for (n, a, z) in rows.values() if n > 1]
lo, hi = min(shares), max(shares)
print(f"families with more than one network: ALL share from {lo} = {float(lo):.4f} to {hi} = {float(hi):.4f}")
claim_lo = Fraction(17, 100)
inside = lo >= Fraction(165, 1000)
print(("PASS" if inside else "FAIL") + f" paper's lower bound 17 %: smallest share {float(lo)*100:.2f} % "
      + ("is" if inside else "is NOT") + " at least 16.5 %")
fails += not inside
print(f"FAILS: {fails}")
sys.exit(1 if fails else 0)
