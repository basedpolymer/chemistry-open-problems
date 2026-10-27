"""Verifier V2: probe of Conjecture CR (round2/ec/notes.md §7.3) on the non-PC cluster-rank network found by
random_wfec_v2.py (seed 3): is J|_S Hurwitz at every positive steady state?  If a steady state with
sign det(J|_S) != (-1)^s exists, the degree argument gives a second steady state in its class.

Environment: WSL Fedora 44, Python 3.12, fractions (exact), sympy for the final exact certificate.
"""
import sys, os, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

T0 = time.time()
N = [((0, 0, 0, 1, 0), (0, 0, 1, 0, 0)), ((0, 0, 0, 2, 0), (0, 0, 2, 0, 0)), ((0, 2, 0, 0, 0), (0, 0, 0, 0, 1)),
     ((1, 0, 2, 0, 1), (1, 2, 2, 0, 0)), ((2, 0, 2, 0, 2), (2, 2, 0, 2, 1))]
n, r = nsp(N), len(N)
s = rank(Gamma(N))
print("Gamma rows:", Gamma(N))
print("s =", s, "clusters:", [(sorted(S), sorted(T)) for S, T in clusters(N)], "PC:", is_PC(N), "EC:", is_EC(N))
print("connected", connected(N), "H1", conservative(N), "H2", consistent(N), "H3", rows_distinct(N))
print("flux rays", [[str(v) for v in ray] for ray in flux_rays(N)])
print("conservation rays", [[str(v) for v in z] for z in cons_rays(N)])
print("critical siphons N, N*:", critical_siphons(N), critical_siphons(dual(N)))
rng = random.Random(7)
rays = flux_rays(N)
found = None
cnt = 0
for trial in range(3000):
    x = [F(rng.randint(1, 40), rng.randint(1, 40)) for _ in range(n)]
    wts = [F(rng.randint(1, 30), rng.randint(1, 30)) for _ in rays]
    v = [sum(w * ray[t] for w, ray in zip(wts, rays)) for t in range(r)]
    k = [v[t] / monomial(x, N[t][0]) for t in range(r)]
    cp = charpoly(jacobian(N, k, x))
    q = strip_zero_roots(cp, n - s)
    if q is None or not is_hurwitz(q):
        cnt += 1
        if found is None:
            found = (x, k, q)
print(f"steady states sampled: 3000; non-Hurwitz on S: {cnt}")
if found:
    x, k, q = found
    print("first non-Hurwitz steady state x =", [str(t) for t in x], "k =", [str(t) for t in k])
    print("q(t) coefficients:", [str(c) for c in q], " det sign of J|_S = sign((-1)^s q(0)):",
          "(-1)^s q(0) =", str((-1) ** s * q[-1]))
print(f"runtime {time.time() - T0:.1f} s")
