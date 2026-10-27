"""V1 symbolic certificates for Q2's witnesses (sympy; run under WSL).

Verifier V1 (Claude Opus 5.5, effort xhigh), 27 October 2026. For each witness N, V = E(lambda) with symbolic
lambda > 0 parametrises C° (relative interior of the flux cone). CF(N) is one block iff for every lambda > 0,
ker A_V = R*1, i.e. rank A_V = m - 1. Certificate: an (m-1)-minor of A_V(lambda) that is a nonzero polynomial with
all coefficients of one sign. Also prints A_V(lambda) for XY4 and Z4.
"""
import sys, time, platform
from itertools import combinations
sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
import sympy as sp
from v1lib import flux_rays, source_classes, split_choice, duplicate, Net
from examples import EX

T0 = time.time()
print("environment: Python", platform.python_version(), "sympy", sp.__version__)
FAILS = []


def claim(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" -- " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def AV_symbolic(N):
    R = flux_rays(N)
    lam = sp.symbols("l1:%d" % (len(R) + 1), positive=True)
    V = [sum(lam[k] * sp.Rational(R[k][1][t].numerator, R[k][1][t].denominator) for k in range(len(R)))
         for t in range(N.r)]
    keys, cls = source_classes(N)
    A = sp.Matrix(N.n, len(cls), lambda i, c: sum(V[t] * N.G[i][t] for t in cls[c]))
    return A, lam, V


def one_block_certificate(N):
    A, lam, V = AV_symbolic(N)
    m = A.shape[1]
    if m == 1:
        return True, "single source class"
    for rows in combinations(range(A.shape[0]), m - 1):
        for cols in combinations(range(m), m - 1):
            d = sp.expand(A.extract(list(rows), list(cols)).det())
            if d == 0:
                continue
            coeffs = sp.Poly(d, *lam).coeffs()
            if all(c > 0 for c in coeffs) or all(c < 0 for c in coeffs):
                return True, "minor rows %s cols %s = %s" % (rows, cols, sp.factor(d))
    return False, "no sign-definite (m-1)-minor"


nets = {k: EX[k] for k in ("B2dup", "K3", "PiCF2", "B2", "B3")}
nets["B2dup-split"] = split_choice(EX["B2dup"], [0, 1])
nets["B3-split"] = split_choice(EX["B3"], [1, 2])
for k, N in nets.items():
    ok, why = one_block_certificate(N)
    claim("CF(%s) is one block for every V in C° (certificate)" % k, ok, why[:300])

for k in ("XY4", "Z4"):
    A, lam, V = AV_symbolic(EX[k])
    print("   %s: V(lambda) = %s ; A_V = %s" % (k, V, A.tolist()))
# XY4: columns proportional to (1,-1) with coefficients V1-V2 and V3-V4: generic rank 1, rank 0 at V=(1,1,2,2)
A, lam, V = AV_symbolic(EX["XY4"])
sol = sp.solve([V[0] - V[1]], lam, dict=True)
claim("XY4: the set {V in C°: A_V = 0} is nonempty (V1 = V2 forces V3 = V4)", bool(sol) and
      all(sp.simplify(x.subs(sol[0])) == 0 for x in A))
print("runtime %.1f s" % (time.time() - T0))
sys.exit(1 if FAILS else 0)
