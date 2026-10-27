"""Verifier V2: targeted search for small counterexamples to the dynamical part of Conjecture CR.

Generator: clusters on species blocks; in each block 1-3 sources (sharing species, not necessarily proportional);
products y' with Z y' = Z y for a random Z with n-q+1 rows (so s <= q-1); molecularity of every complex <= MAXMOL.
Kept: connected, (H1), (H2), s = q-1, not PC.  For each, 60 rational steady states x with k = v / x^Ye are sampled
(exact) and det(J|_S) is compared with (-1)^s.  Usage: python cr_search2.py SEED MAXN MAXMOL SECONDS
"""
import sys, os, time, random
from itertools import product
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

T0 = time.time()
seed, MAXN, MAXMOL, SECS = (int(a) for a in sys.argv[1:5])
rng = random.Random(seed)


def gen(n, q):
    while True:
        lab = [rng.randrange(q) for _ in range(n)]
        if len(set(lab)) == q:
            break
    B = [[p for p in range(n) if lab[p] == c] for c in range(q)]
    srcs = []
    for blk in B:
        ys = set()
        full = [0] * n
        for p in blk:
            full[p] = 1
        if sum(full) > MAXMOL:
            return None
        ys.add(tuple(full))
        for _ in range(rng.randint(0, 2)):
            y = [0] * n
            for p in blk:
                y[p] = rng.choice([0, 1, 1, 2])
            if not any(y) or sum(y) > MAXMOL:
                continue
            ys.add(tuple(y))
        srcs.append(sorted(ys))
    rows = n - q + 1
    while True:
        Z = [[rng.choice([0, 0, 1, 1, 2]) for _ in range(n)] for _ in range(rows)]
        if all(any(Z[i][p] for i in range(rows)) for p in range(n)) and rank(Z) == rows:
            break
    cand = [c for c in product(range(MAXMOL + 1), repeat=n) if 0 < sum(c) <= MAXMOL]

    def zv(y):
        return tuple(sum(Z[i][p] * y[p] for p in range(n)) for i in range(rows))
    N = []
    for ys in srcs:
        for y in ys:
            prods = [c for c in cand if c != y and zv(c) == zv(y)]
            if not prods:
                return None
            for _ in range(rng.randint(1, 2)):
                N.append((y, rng.choice(prods)))
    N = sorted(set(N))
    if not connected(N) or not conservative(N) or not consistent(N):
        return None
    if rank(Gamma(N)) != len(clusters(N)) - 1 or is_PC(N):
        return None
    return N


tested = found = 0
best = None
seen = set()
while time.time() - T0 < SECS:
    n = rng.randint(3, MAXN)
    q = rng.randint(2, min(3, n))
    N = gen(n, q)
    if N is None or tuple(N) in seen:
        continue
    seen.add(tuple(N))
    tested += 1
    s = rank(Gamma(N))
    rays = flux_rays(N)
    for trial in range(60):
        x = [F(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(n)]
        wts = [F(rng.randint(1, 30), rng.randint(1, 30)) for _ in rays]
        v = [sum(w * ray[t] for w, ray in zip(wts, rays)) for t in range(len(N))]
        k = [v[t] / monomial(x, N[t][0]) for t in range(len(N))]
        qq = strip_zero_roots(charpoly(jacobian(N, k, x)), n - s)
        if qq is None:
            continue
        detJS = (-1) ** s * qq[-1]
        if detJS != 0 and (detJS > 0) != (s % 2 == 0):
            found += 1
            size = (n, len(N), max(max(sum(a), sum(b)) for a, b in N))
            if best is None or size < best[0]:
                best = (size, N)
                print("found", size, "N =", N, "x =", [str(t) for t in x], "k =", [str(t) for t in k],
                      "det(J|_S) =", str(detJS), "s =", s, "H3:", rows_distinct(N), flush=True)
            break
print(f"seed {seed}: tested {tested} non-PC cluster-rank networks (n <= {MAXN}, molecularity <= {MAXMOL}), "
      f"{found} with a wrong-sign steady state; runtime {time.time() - T0:.0f} s")
