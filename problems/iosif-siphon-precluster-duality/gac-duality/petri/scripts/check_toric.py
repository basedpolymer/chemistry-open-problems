"""Cluster-throughput kernel of mass-action EC networks (agent P; conjecture K10 of petri/candidates.md).

Environment: Python 3.12, standard library only, exact rationals.  Usage: python check_toric.py N MAXMOL RMAX ORD
[NPROC] [NK].  Runtime: printed at the end.

For an EC network with distinct source complexes y_c (pairwise disjoint supports) and rate constants k, the
mass-action vector field is f(x) = A_k phi(x), phi_c(x) = x^{y_c}, A_k[:, c] = sum_{t in T_c} k_t (y'_t - y_c).
Positive steady states <=> phi(x) in ker A_k cap R_{>0}.  Claim tested (K10): if N is well-formed EC (connected,
(H1), (H2), s = #source complexes - 1), then for EVERY k > 0, ker A_k is one-dimensional and spanned by a
positive vector u(k); hence the positive steady states are exactly {x > 0 : x^{y_c} = lambda u_c(k), lambda > 0},
a nonempty log-affine (toric) set of dimension n - s.  For each network NK rate vectors are drawn with a fixed seed
(entries in {1/7, 1/3, 1/2, 1, 2, 3, 7}) plus k = (1,...,1).  PASS/FAIL on the WF class; on non-WF EC networks the
script only reports how often the property holds.
"""
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import permutations
from multiprocessing import Pool

from petri import Network, is_connected, ec_sets, rank
from crn import nullspace
from enumerate_ec import gen_networks, canonical, set_partitions

VALS = [Fraction(1, 7), Fraction(1, 3), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(3), Fraction(7)]


def throughput_kernel(N, k):
    classes = ec_sets(N)
    cols = []
    for c in classes:
        col = [Fraction(0)] * N.n
        for t in c:
            a, b = N.reactions[t]
            for i in range(N.n):
                col[i] += k[t] * (b[i] - a[i])
        cols.append(col)
    A = [[cols[c][i] for c in range(len(classes))] for i in range(N.n)]
    return A, nullspace(A, len(classes))


def check(N, rng, nk):
    s = N.rank()
    ks = [[Fraction(1)] * N.r] + [[rng.choice(VALS) for _ in range(N.r)] for _ in range(nk)]
    for k in ks:
        A, ker = throughput_kernel(N, k)
        if len(ker) != 1:
            return False, ("kernel dim", len(ker), k)
        u = ker[0]
        if not (all(x > 0 for x in u) or all(x < 0 for x in u)):
            return False, ("sign", u, k)
        if rank(A) != s:
            return False, ("rank", rank(A), s, k)
    return True, None


def worker(args):
    n, m, rmax, ordinary, pidx, nk = args
    rng = random.Random(1000 * pidx + 7)
    perms = list(permutations(range(n)))
    seen = set()
    out = []
    for rx in gen_networks(n, m, rmax, ordinary, pidx):
        rx = tuple(sorted(rx))
        can = canonical(rx, perms)
        if can != rx or can in seen:
            continue
        seen.add(can)
        N = Network(n, rx)
        if not is_connected(N) or not N.is_consistent() or not N.is_conservative():
            continue
        wf = N.rank() == len(ec_sets(N)) - 1
        ok, info = check(N, rng, nk)
        out.append((rx, wf, ok, info))
    return out


def main():
    n, m, rmax, ordinary = (int(a) for a in sys.argv[1:5])
    nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    nk = int(sys.argv[6]) if len(sys.argv) > 6 else 3
    t0 = time.time()
    print(f"# check_toric.py n={n} m={m} r<={rmax} ordinary={ordinary} NK={nk}; Python {platform.python_version()}"
          f" on {platform.system()}")
    nparts = sum(1 for _ in set_partitions(list(range(n))))
    cnt = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
    first = {}
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, [(n, m, rmax, ordinary, p, nk) for p in range(nparts)]):
            for rx, wf, ok, info in chunk:
                cnt[(wf, ok)] += 1
                if not ok and wf not in first:
                    first[wf] = (rx, info)
                elif not ok and wf in first and len(rx) < len(first[wf][0]):
                    first[wf] = (rx, info)
    bad = cnt[(True, False)]
    print(f"{'PASS' if bad == 0 else 'FAIL'}  K10 WF-EC: ker A_k one-dimensional and positive for all sampled k: "
          f"{cnt[(True, True)]}/{cnt[(True, True)] + bad}" + (f"  first failure {first[True]}" if bad else ""))
    print(f"REPORT  non-WF EC (H1,H2,connected): property holds for {cnt[(False, True)]}/"
          f"{cnt[(False, True)] + cnt[(False, False)]}" + (f"  e.g. failure {first[False]}" if False in first else ""))
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
