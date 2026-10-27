"""Verifier V2: a molecularity-3 counterexample to the dynamical part of Conjecture CR (found by cr_search2.py).

NC3 = {D -> B+C, C -> A, A+B -> B+C, 2A+B -> C+D}  (species A, B, C, D; n = 4, r = 4, molecularity 3).
Checks (exact): connected, (H1)-(H3) on NC3 and NC3*; s = 2, Petri clusters {D}, {C}, {A,B}: q = 3 = s + 1; TEFC, not
PC, not EC; no critical siphon on either side; a rational positive steady state with det(J|_S) < 0 (index -1 while
(-1)^s = +1), hence (degree argument) at least two more positive steady states in its class; then msolve (exact real
root isolation) counts all real solutions of the steady-state equations in that class.

Environment: WSL Fedora 44, Python 3.12, msolve.  PASS/FAIL, exit 1 on FAIL.
"""
import sys, os, time, random, subprocess, re, tempfile
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

T0 = time.time()
FAILS = []


def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + name + ((" -- " + info) if info else ""))
    if not cond:
        FAILS.append(name)


names = ["A", "B", "C", "D"]
NC3 = parse("D->B+C, C->A, A+B->B+C, 2A+B->C+D", names)
n, r = 4, 4
s = rank(Gamma(NC3))
cl = clusters(NC3)
D3 = dual(NC3)
check("connected, (H1)-(H3) on NC3 and NC3*", connected(NC3) and H123(NC3) and connected(D3) and H123(D3))
check("s = 2, q = 3 (clusters {D},{C},{A,B}): cluster rank equation", s == 2 and len(cl) == 3,
      str([(sorted(S), sorted(T)) for S, T in cl]))
tefc = all(set(i for i, x in enumerate(a) if x) == set(i for i, x in enumerate(b) if x)
           for a, _ in NC3 for b, _ in NC3 if any(x and y for x, y in zip(a, b)))
check("TEFC (A+B and 2A+B have the same support), not PC, not EC", tefc and not is_PC(NC3) and not is_EC(NC3))
check("no critical siphon in NC3 nor NC3*", not critical_siphons(NC3) and not critical_siphons(D3))
print("   conservation rays:", [[str(t) for t in z] for z in cons_rays(NC3)], " flux rays:",
      [[str(t) for t in v] for v in flux_rays(NC3)], " NC3* =", D3)

rng = random.Random(3)
rays = flux_rays(NC3)
best = None
for trial in range(40000):
    x = [F(rng.randint(1, 8), rng.randint(1, 4)) for _ in range(n)]
    wts = [F(rng.randint(1, 8), rng.randint(1, 4)) for _ in rays]
    v = [sum(w * ray[t] for w, ray in zip(wts, rays)) for t in range(r)]
    k = [v[t] / monomial(x, NC3[t][0]) for t in range(r)]
    q = strip_zero_roots(charpoly(jacobian(NC3, k, x)), n - s)
    if q[-1] < 0:  # det(J|_S) = (-1)^2 q(0)
        size = max(max(t.numerator, t.denominator) for t in x + k)
        if best is None or size < best[0]:
            best = (size, x, k, q[-1])
check("a rational positive steady state with det(J|_S) < 0 exists", best is not None)
_, x, k, dJ = best
check("f_k(x) = 0 exactly", all(e == 0 for e in field(NC3, k, x)))
Z = cons_rays(NC3)
tot = [sum(z[i] * x[i] for i in range(n)) for z in Z]
print("   x =", [str(t) for t in x], " k =", [str(t) for t in k], " det(J|_S) =", str(dJ),
      " class:", [f"{'+'.join(f'{z[i]}{names[i]}' for i in range(n) if z[i])} = {T}" for z, T in zip(Z, tot)])

# msolve: all real solutions of f = 0 together with the two conservation equations
lines = []
for i in range(n):
    terms = []
    for t, (a, b) in enumerate(NC3):
        g = b[i] - a[i]
        if g:
            mono = "*".join(f"{names[j]}^{a[j]}" for j in range(n) if a[j])
            terms.append(f"({g}*{k[t].numerator}/{k[t].denominator})*{mono}")
    if terms:
        lines.append("+".join(terms))
for z, T in zip(Z, tot):
    lines.append("+".join(f"({z[i].numerator}/{z[i].denominator})*{names[i]}" for i in range(n) if z[i])
                 + f"-({T.numerator}/{T.denominator})")
# msolve wants integer coefficients: clear denominators per polynomial with sympy
import sympy as sp
syms = sp.symbols(names)
polys = []
for L in lines:
    e = sp.together(sp.sympify(L.replace("^", "**"), locals=dict(zip(names, syms))))
    num, den = sp.fraction(e)
    polys.append(sp.expand(num))
with tempfile.TemporaryDirectory() as td:
    fin = os.path.join(td, "in.ms")
    fout = os.path.join(td, "out.ms")
    with open(fin, "w") as fh:
        fh.write(",".join(names) + "\n0\n")
        fh.write(",\n".join(str(p).replace("**", "^") for p in polys) + "\n")
    subprocess.run(["msolve", "-f", fin, "-o", fout], check=True, capture_output=True)
    out = open(fout).read()
# parse isolating intervals: numbers like  a / 2^b  or integers
txt = out.replace("\n", "")
txt = re.sub(r"(-?\d+)\s*/\s*2\^(\d+)", lambda m: f"Fraction({m.group(1)}, 2**{m.group(2)})", txt)
txt = txt.rstrip(":").strip()
data = eval(txt, {"Fraction": Fraction})
sols = data[1][1]
pos = []
for box in sols:
    lo = [iv[0] for iv in box]
    hi = [iv[1] for iv in box]
    if all(l > 0 for l in lo):
        pos.append([(float(a) + float(b)) / 2 for a, b in zip(lo, hi)])
    elif any(h < 0 for h in hi):
        pass
    else:
        print("   ambiguous box (touches 0):", [(float(a), float(b)) for a, b in zip(lo, hi)])
print(f"   msolve: {len(sols)} real solutions in the class, {len(pos)} with all coordinates > 0 (isolating boxes)")
check("exactly three positive steady states in the class (msolve)", len(pos) == 3)
import numpy as np
for p_ in pos:
    J = np.array([[float(e) for e in row] for row in jacobian(NC3, k, [Fraction(c).limit_denominator(10 ** 12) for c in p_])])
    ev = sorted(np.linalg.eigvals(J), key=lambda z: z.real)
    print("   steady state", [f"{c:.8g}" for c in p_], " eigenvalues:", [f"{e.real:.4g}" for e in ev])
print(f"runtime {time.time() - T0:.1f} s; {len(FAILS)} FAIL")
sys.exit(1 if FAILS else 0)
