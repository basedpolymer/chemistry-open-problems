"""Verifier V2: random search for small counterexamples to the dynamical part of Conjecture CR (round2/ec/notes.md
§7.3): connected, (H1), (H2), s = #clusters - 1, and a positive steady state at which J|_S has
sign det(J|_S) = (-1)^(s+1) (then, by the degree argument, the class holds at least two positive steady states,
since no critical siphon exists by Theorem C^gen and classes are compact).

Exact arithmetic for the certificate (rational x, rational v in ker Gamma, k = v / x^y).  Usage:
python cr_search.py SEED MAXN MAXMOL MAXR SECONDS
"""
import sys, os, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

T0 = time.time()
seed, MAXN, MAXMOL, MAXR, SECS = (int(a) for a in sys.argv[1:6])
rng = random.Random(seed)


def rand_cplx(n):
    while True:
        c = tuple(rng.randint(0, MAXMOL) for _ in range(n))
        if 0 < sum(c) <= MAXMOL:
            return c


best = None
tested = 0
found = 0
while time.time() - T0 < SECS:
    n = rng.randint(2, MAXN)
    r = rng.randint(3, MAXR)
    N = sorted(set((rand_cplx(n), rand_cplx(n)) for _ in range(r)))
    N = [(a, b) for a, b in N if a != b]
    if len(N) < 3 or any(all(a[i] == 0 and b[i] == 0 for a, b in N) for i in range(n)):
        continue
    if not connected(N) or not consistent(N) or not conservative(N):
        continue
    s = rank(Gamma(N))
    q = len(clusters(N))
    if s != q - 1 or s == 0 or is_PC(N):
        continue
    tested += 1
    rays = flux_rays(N)
    for trial in range(60):
        x = [F(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(n)]
        wts = [F(rng.randint(1, 30), rng.randint(1, 30)) for _ in rays]
        v = [sum(w * ray[t] for w, ray in zip(wts, rays)) for t in range(len(N))]
        k = [v[t] / monomial(x, N[t][0]) for t in range(len(N))]
        cp = charpoly(jacobian(N, k, x))
        qq = strip_zero_roots(cp, n - s)
        if qq is None:
            continue
        sgn = (-1) ** s * qq[-1]  # det(J|_S)
        if sgn != 0 and (sgn > 0) != (s % 2 == 0):
            size = (n, len(N), max(sum(a) for a, b in N + [(b, a) for a, b in N]))
            found += 1
            if best is None or size < best[0]:
                best = (size, N, x, k, sgn)
                print("found", size, "N =", N, "x =", [str(t) for t in x], "k =", [str(t) for t in k],
                      "det(J|_S) =", str(sgn), "s =", s, "H3:", rows_distinct(N), flush=True)
            break
print(f"seed {seed}: tested {tested} non-PC cluster-rank networks, {found} with a wrong-sign steady state; "
      f"runtime {time.time() - T0:.0f} s")
