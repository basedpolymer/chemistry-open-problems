"""V1: smallest networks with pi != CF in the two-species families (exact CF, s = 1). Verifier V1 (Claude Opus 5.5)."""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from v1lib import *
import check_Q2_counts as C
for mol, R in ((2, 5), (3, 4)):
    nets = C.enumerate_n2(mol, R)
    hits = [N for N in nets if flat_closure_partition(N) != C.cf_rank_one(N)]
    print("family (2,%d,%d): pi != CF on %d networks, smallest r = %d" % (mol, R, len(hits), min(N.r for N in hits)))
    for N in sorted(hits, key=lambda N: N.r)[:3]:
        print("   ", N, "pi", flat_closure_partition(N), "CF", C.cf_rank_one(N))
