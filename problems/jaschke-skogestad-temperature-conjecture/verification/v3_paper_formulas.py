#!/usr/bin/env python3
"""Independent verifier check V3: the printed formulas of Jaschke & Skogestad (2014), transcribed by the verifier
from 400-dpi renderings of journal pages 46-48 (not copied from the solver's files).

Notation: t1, t2, t3 = theta_{1,1}, theta_{2,1}, theta_{3,1}; h1, h2, h3 = theta_{h1,1}, theta_{h2,1}, theta_{h3,1}.
In the controlled-variable formulas the cold temperatures are independent measurements (for given hot inlets, any
cold profile is reached by choosing the UA_i), so identities are tested in Q(t1, t2, t3, h1, h2, h3).

Claims (PASS/FAIL):
  P1  (A.8) raw coefficient of p11 == t1^2/h1                      (a1, same for both recursions)
  P2  (A.10) raw coefficient of p21 == a2                          (same for both recursions)
  P3  (A.12) raw coefficient of p31 == corrected a3 = c3 (t3 + t2 - a1 - a2)
  P4  (A.12) raw coefficient of p31 != published a3 = c3 (t3 + t2 - a2) (= the last line of (A.12)), and the
      difference is exactly c3*a1
  P5  (A.12) middle line (the 'simplified' form) == published a3 (so the term is lost between line 1 and line 2)
Environment: WSL Fedora, Python 3, sympy 1.14. Runtime: < 5 s.
"""
import sys
import sympy as sp

FAIL = 0
def report(name, ok):
    global FAIL
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok:
        FAIL += 1

t1, t2, t3, h1, h2, h3 = sp.symbols('t1 t2 t3 h1 h2 h3')
den = h1*(-t2*h2 + h2*h3 + t1*t2 - t1*h3)          # printed denominator of (A.8), (A.10), (A.12)

A8_p11 = (-t1**2*t2*h2 + t1**3*t2 + t1**2*h2*h3 - t1**3*h3) / den
A10_p21 = (t1**2*t2*h1 + t1**2*t2**2 - t1**3*t2 - t1**2*t2*h3 + t1**3*h3 - t1**2*h1*h3 - t2**3*h1 + t2**2*h1*h3) / den
A12_p31 = (-t1**2*t2*h1 + t1**2*t2*h2 + t1**2*t3*h1 - t1**2*t3*h2 - t1**2*t2**2 + t3**2*h1*h2 + t1*t2**2*h1
           - t2**2*h1*h2 - t1*t3**2*h1 + t1**2*t2*t3 - t2**2*t3*h1 + t2**3*h1) / den
A12_mid = (t3 - t2)*(t3 + t2 - ((t2 - t1)*(t2 + t1 - t1**2/h1))/(h2 - t1)) / (h3 - t2)

c = [t1/h1, (t2 - t1)/(h2 - t1), (t3 - t2)/(h3 - t2)]
th = [0, t1, t2, t3]
a_corr, a_pap = [], []
for i in range(3):
    a_corr.append(c[i]*(th[i+1] + th[i] - sum(a_corr)))
    a_pap.append(c[i]*(th[i+1] + th[i] - (a_pap[-1] if a_pap else 0)))

z = lambda e: sp.cancel(sp.together(e)) == 0
report("P1 (A.8) raw p11 coefficient == a1 = t1^2/h1", z(A8_p11 - a_corr[0]))
report("P2 (A.10) raw p21 coefficient == a2", z(A10_p21 - a_corr[1]))
report("P3 (A.12) raw p31 coefficient == corrected a3", z(A12_p31 - a_corr[2]))
report("P4 (A.12) raw p31 coefficient != published a3, difference == c3*a1",
       (not z(A12_p31 - a_pap[2])) and z(a_pap[2] - A12_p31 - c[2]*a_corr[0]))
report("P5 (A.12) middle line == published a3 (term lost between line 1 and line 2)", z(A12_mid - a_pap[2]))
# sensitivity of P3 to the transcription: every single-term sign flip of the 12-term numerator breaks P3
num = sp.expand(A12_p31*den)
terms = sp.Add.make_args(num)
robust = all(not z((num - 2*tm)/den - a_corr[2]) for tm in terms)
report(f"P6 P3 is not an accident of transcription: each of the {len(terms)} single sign flips breaks it", robust and len(terms) == 12)
print("ALL PASS" if FAIL == 0 else f"{FAIL} FAIL")
sys.exit(1 if FAIL else 0)
