"""Random search for candidate K8 (agent P): two networks with (H1), (H2) and rows and columns of Gamma pairwise linearly independent (strong (H3)), the same twin classes and the same minimal
siphons but different preclusters of the dual.  Environment: Python 3.12+, standard library, exact rationals.
Usage: python k8_search.py N MAXMOL RMIN RMAX SECONDS SEED.  Prints the first collision found (or none)."""
import sys, time, random
from itertools import product
from petri import Network, twin_classes

n, m, rmin, rmax, secs, seed = (int(a) for a in sys.argv[1:7])
rng = random.Random(seed)
C = [c for c in product(range(m + 1), repeat=n) if 0 < sum(c) <= m]
seen = {}
t0 = time.time(); tried = good = 0
key = lambda part: tuple(sorted(tuple(sorted(b)) for b in part))
while time.time() - t0 < secs:
    r = rng.randint(rmin, rmax)
    rx = set()
    while len(rx) < r:
        a, b = rng.choice(C), rng.choice(C)
        if a != b:
            rx.add((a, b))
    rx = tuple(sorted(rx)); tried += 1
    N = Network(n, rx)
    if not all(any(a[i] or b[i] for a, b in rx) for i in range(n)):
        continue
    tw = [b for b in twin_classes(N) if len(b) >= 2]
    if not tw or not N.rows_pairwise_independent() or not N.columns_pairwise_independent() or not N.is_consistent() or not N.is_conservative():
        continue
    good += 1
    k = (key(tw), key(N.minimal_siphons()))
    P = key(N.dual().preclusters())
    if k in seen and seen[k][0] != P:
        print("COLLISION", k, "\n  N1 =", seen[k][1], "P(N1*) =", seen[k][0], "\n  N2 =", rx, "P(N2*) =", P)
        break
    seen.setdefault(k, (P, rx))
print(f"# tried {tried}, with (H1)-(H3) and a twin pair {good}, distinct keys {len(seen)}, {time.time()-t0:.0f} s")
