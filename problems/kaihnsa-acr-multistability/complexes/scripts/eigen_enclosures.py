"""Ball-arithmetic enclosures of the Jacobian eigenvalues of N1 and N2 at their positive steady states
(runtime: under 1 second; needs python-flint: run in the WSL sci environment with  python eigen_enclosures.py).

Complement to certify_examples.py (which already proves stability exactly by Routh-Hurwitz): the exact rational
Jacobians are taken from certify_examples.py, their characteristic polynomials are formed exactly and their roots
are enclosed with arb (python-flint, 200 bits). PASS/FAIL per point; exit code 1 on failure.
"""
import sys
import time
from fractions import Fraction as Fr

import flint
import certify_examples as ce

T0 = time.time()
FAIL = 0
flint.ctx.prec = 200
for name, N in ce.NETWORKS.items():
    f = ce.ode(N["reactions"], N["rates"])
    for p, want_stable in zip(N["points"], N["stable"]):
        p = tuple(Fr(v) for v in p)
        J = ce.jacobian(f, p)
        a1, a2, a3 = ce.rh(J)
        poly = flint.fmpq_poly([flint.fmpq(v.numerator, v.denominator) for v in (a3, a2, a1, Fr(1))])
        roots = poly.complex_roots()
        lab = "%s at (%s)" % (name, ", ".join(str(v) for v in p))
        print(lab, " characteristic polynomial: lambda^3 + (%s) lambda^2 + (%s) lambda + (%s)" % (a1, a2, a3))
        signs = []
        for r, mult in roots:
            print("    eigenvalue in", r, " multiplicity", mult)
            signs.append(-1 if r.real < 0 else (1 if r.real > 0 else 0))
        if want_stable:
            ok = sorted(signs) == [-1, -1, -1]
            print(("PASS " if ok else "FAIL ") + lab + ": all eigenvalues have certified negative real part")
        else:
            ok = sorted(signs) == [-1, -1, 1]
            print(("PASS " if ok else "FAIL ") + lab + ": one eigenvalue with certified positive real part")
        FAIL += not ok
print("runtime %.3f s" % (time.time() - T0))
if FAIL:
    print("%d check(s) FAILED" % FAIL)
    sys.exit(1)
print("ALL PASS")
