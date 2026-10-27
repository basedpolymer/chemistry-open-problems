"""Verifier V2: independent exhaustive enumeration of small EC networks with P's conventions (petri/scripts/
enumerate_ec.py docstring: up to species permutation; distinct sources with pairwise disjoint supports covering all
species, coefficients >= 1, molecularity <= MAXMOL; each source has a nonempty set of distinct nonzero products of
molecularity <= MAXMOL different from it; at most RMAX reactions; connected, (H1), (H2)), count EC and WF-EC, and
run the exact checks of Theorems T, C, C', U(b), U(c), U*(a) on every WF-EC network.

Environment: WSL Fedora 44, Python 3.12, standard library (fractions), multiprocessing.
Usage: python enum_wfec_v2.py N MAXMOL RMAX [NPROC]
"""
import sys, os, time, random
from itertools import product, permutations, combinations
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

n, MAXMOL, RMAX = (int(a) for a in sys.argv[1:4])
NPROC = int(sys.argv[4]) if len(sys.argv) > 4 else 8
T0 = time.time()


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        for k in range(len(p)):
            yield p[:k] + [[first] + p[k]] + p[k + 1:]
        yield [[first]] + p


cplx = [c for c in product(range(MAXMOL + 1), repeat=n) if 0 < sum(c) <= MAXMOL]
perms = list(permutations(range(n)))


def canon(N):
    best = None
    for p in perms:
        M = tuple(sorted((tuple(a[p[i]] for i in range(n)), tuple(b[p[i]] for i in range(n))) for a, b in N))
        if best is None or M < best:
            best = M
    return best


def gen_all():
    seen = set()
    for part in set_partitions(list(range(n))):
        # sources: one per block
        choices = []
        for blk in part:
            opts = []
            for coeffs in product(range(1, MAXMOL + 1), repeat=len(blk)):
                if sum(coeffs) <= MAXMOL:
                    v = [0] * n
                    for i, c in zip(blk, coeffs):
                        v[i] = c
                    opts.append(tuple(v))
            choices.append(opts)
        for srcs in product(*choices):
            m = len(srcs)
            prods = [[c for c in cplx if c != y] for y in srcs]

            def rec(i, acc, budget):
                if i == m:
                    yield list(acc)
                    return
                for size in range(1, budget - (m - i - 1) + 1):
                    for P in combinations(prods[i], size):
                        yield from rec(i + 1, acc + [(srcs[i], q) for q in P], budget - size)
            for N in rec(0, [], RMAX):
                c = canon(N)
                if c in seen:
                    continue
                seen.add(c)
                yield [tuple(map(tuple, r)) for r in c]


def work(N):
    N = [tuple(r) for r in N]
    if not connected(N) or not conservative(N) or not consistent(N):
        return None
    res = {"EC": 1}
    if not is_WFEC(N):
        res["WF"] = 0
        return res
    res["WF"] = 1
    rng = random.Random(hash(tuple(N)) & 0xffffffff)
    r = len(N)
    s = rank(Gamma(N))
    srcs, ct = ec_structure(N)
    m = len(srcs)
    res["impure"] = int(any(a[i] and b[i] for a, b in N for i in range(n)))
    res["unbalanced"] = int(not is_balanced(N))
    res["T"] = int(all(u_of_k(N, rand_pos(rng, r)) is not None for _ in range(2)))
    res["C"] = int(not critical_siphons(N))
    res["C'"] = int(not critical_siphons(dual(N)))
    ok = True
    for _ in range(2):
        k = rand_pos(rng, r)
        x = rand_pos(rng, n)
        ok = ok and hurwitz_on_S(N, k, x)
    res["U(b)"] = int(ok)
    if s > 0:
        BS = row_basis(transpose(Gamma(N)))
        BD = row_basis([[F(srcs[c][p] - srcs[0][p]) for p in range(n)] for c in range(1, m)])
        res["U(c)"] = int(maxminor_sign_condition(BS, BD))
        res["S!=D"] = int(rank(BS + BD) > s)
    else:
        res["U(c)"] = 1
        res["S!=D"] = 0
    D = dual(N)
    ok = True
    for _ in range(2):
        cp = charpoly(jacobian(D, rand_pos(rng, n), rand_pos(rng, r)))
        q = strip_zero_roots(cp, r - s)
        ok = ok and q is not None and is_hurwitz(q)
    res["U*(a)"] = int(ok)
    return res


if __name__ == "__main__":
    nets = list(gen_all())
    with Pool(NPROC) as pool:
        out = [o for o in pool.map(work, nets, chunksize=50) if o is not None]
    ec = len(out)
    wf = [o for o in out if o["WF"]]
    print(f"# enum_wfec_v2.py n={n} molecularity<={MAXMOL} r<={RMAX}: generated {len(nets)} canonical EC networks")
    print(f"connected+(H1)+(H2): EC = {ec}, WF-EC = {len(wf)}  (impure {sum(o['impure'] for o in wf)}, "
          f"unbalanced {sum(o['unbalanced'] for o in wf)}, S != D {sum(o['S!=D'] for o in wf)})")
    bad = False
    for key in ("T", "C", "C'", "U(b)", "U(c)", "U*(a)"):
        okc = sum(o[key] for o in wf)
        tag = "PASS" if okc == len(wf) else "FAIL"
        bad |= tag == "FAIL"
        print(f"{tag} {key}: {okc}/{len(wf)}")
    print(f"runtime {time.time() - T0:.1f} s")
    sys.exit(1 if bad else 0)
