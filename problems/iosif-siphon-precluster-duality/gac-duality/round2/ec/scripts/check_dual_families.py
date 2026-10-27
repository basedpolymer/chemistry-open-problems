"""Exact checks of Theorem U* (round2/ec/notes.md): for a WF-EC network N and ANY positive point k of N* (not only
steady states), with any kappa > 0, the Jacobian J* of N* has det(tI - J*) = t^(r-s) q(t) with q Hurwitz; the
m x m matrix G* = E_T^T diag(1/k) Gamma^T diag(omega a) E_S is Metzler and u(1/k) (positive kernel vector of
A_{1/k}) is a left null vector of G*.  Agent Q3, Claude Opus 5.5 (effort xhigh).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rationals.
Usage (from round2/ec/):
  python scripts/check_dual_families.py P N MAXMOL RMAX ORDINARY [NPROC]   (agent P's EC generator)
  python scripts/check_dual_families.py U N RMAX MAXMOL [NPROC]            (unbalanced generator of dual_search.py)
Claims (PASS/FAIL, exit 1 on FAIL) on every WF-EC network: Ustar (three points k, kappa per network: all ones and
two seeded random rationals); UN (Theorem U(b) and the GS identity at two random NON-steady points x of N, random k:
G = Dphi(x) A_k Metzler, G u(k) = 0, J|_S Hurwitz).
"""
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import permutations
from multiprocessing import Pool

import eclib
from eclib import Network, ECStructure


def dual_certificate(N, E, k, kappa):
    n, r = N.n, N.r
    G = N.Gamma
    Ye = N.Ye
    omega = [kappa[p] * eclib.monomial(k, Ye[p]) for p in range(n)]
    a = [E.sources[E.cl_of_species[p]][p] for p in range(n)]
    m = E.m
    Gs = [[Fraction(0)] * m for _ in range(m)]
    for t in range(r):
        c = E.cl_of_reaction[t]
        for p in range(n):
            if G[p][t]:
                Gs[c][E.cl_of_species[p]] += G[p][t] * a[p] * omega[p] / k[t]
    if any(Gs[c][d] < 0 for c in range(m) for d in range(m) if c != d):
        return False, "G* not Metzler"
    ker = eclib.nullspace(E.A([1 / x for x in k]), m)
    if len(ker) != 1:
        return False, "ker A_{1/k} not one-dimensional"
    xi = ker[0]
    if not (all(x > 0 for x in xi) or all(x < 0 for x in xi)):
        return False, "u(1/k) not positive"
    if any(sum(xi[c] * Gs[c][d] for c in range(m)) != 0 for d in range(m)):
        return False, "u(1/k) not a left null vector of G*"
    J = [[Fraction(0)] * r for _ in range(r)]          # J*[t][t'] = sum_p Gamma[p][t] omega_p Ye[p][t'] / k_t'
    for t in range(r):
        for t2 in range(r):
            J[t][t2] = sum(G[p][t] * omega[p] * Ye[p][t2] for p in range(n)) / k[t2]
    cp = eclib.charpoly(J)
    s = N.rank()
    low = r - s
    if any(cp[i] != 0 for i in range(low)):
        return False, "J* has extra zero eigenvalues"
    q = cp[low:]
    if q[0] == 0 or not eclib.is_hurwitz(q):
        return False, f"q not Hurwitz: {q}"
    return True, ""


def n_side_certificate(N, E, rng):
    """Theorem U(b) at a NON-steady point and the identity behind Theorem GS: at a random rational x > 0 and random
    rational k > 0, G = Dphi(x) A_k is Metzler, G u(k) = 0 with u(k) > 0 the kernel vector of A_k, and
    det(tI - J) = t^(n-s) q(t) with q Hurwitz."""
    x = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(N.n)]
    k = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(N.r)]
    A = E.A(k)
    ker = eclib.nullspace(A, E.m)
    if len(ker) != 1 or not (all(v > 0 for v in ker[0]) or all(v < 0 for v in ker[0])):
        return False
    u = [abs(v) for v in ker[0]]
    phi = [eclib.monomial(x, y) for y in E.sources]
    Dphi = [[phi[c] * E.sources[c][p] / x[p] for p in range(N.n)] for c in range(E.m)]
    G = eclib.matmul(Dphi, A)
    if any(G[c][d] < 0 for c in range(E.m) for d in range(E.m) if c != d):
        return False
    if any(sum(G[c][d] * u[d] for d in range(E.m)) != 0 for c in range(E.m)):
        return False
    J = eclib.jacobian(N, k, x)
    cp = eclib.charpoly(J)
    low = N.n - N.rank()
    if any(cp[i] != 0 for i in range(low)):
        return False
    q = cp[low:]
    return q[0] != 0 and eclib.is_hurwitz(q)


def evaluate(N, rng):
    E = ECStructure(N)
    pts = [([Fraction(1)] * N.r, [Fraction(1)] * N.n)]
    for _ in range(2):
        pts.append(([Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(N.r)],
                    [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(N.n)]))
    ok = all(dual_certificate(N, E, k, kap)[0] for k, kap in pts)
    okN = all(n_side_certificate(N, E, rng) for _ in range(2))
    return {"Ustar": ok, "UN": okN}, {}, eclib.is_balanced(N)


def worker_P(args):
    from enumerate_ec import gen_networks, canonical
    n, m, rmax, ordinary, pidx = args
    perms = list(permutations(range(n)))
    rng = random.Random(4242 + pidx)
    seen = set()
    out = []
    for rx in gen_networks(n, m, rmax, ordinary, pidx):
        rx = tuple(sorted(rx))
        can = canonical(rx, perms)
        if can != rx or can in seen:
            continue
        seen.add(can)
        N = Network(n, rx)
        if not eclib.is_wf_ec(N):
            continue
        out.append((rx,) + evaluate(N, rng))
    return out


def worker_U(args):
    rx, n, seed = args
    rng = random.Random(seed)
    N = Network(n, rx)
    return [(rx,) + evaluate(N, rng)]


def main():
    mode = sys.argv[1]
    t0 = time.time()
    print(f"# check_dual_families.py {' '.join(sys.argv[1:])}; Python {platform.python_version()} on "
          f"{platform.system()}")
    if mode == "P":
        from enumerate_ec import set_partitions
        n, m, rmax, ordinary = (int(a) for a in sys.argv[2:6])
        nproc = int(sys.argv[6]) if len(sys.argv) > 6 else 12
        nparts = sum(1 for _ in set_partitions(list(range(n))))
        tasks = [(n, m, rmax, ordinary, p) for p in range(nparts)]
        worker = worker_P
    else:
        import dual_search
        n, rmax, maxmol = (int(a) for a in sys.argv[2:5])
        nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 12
        nets, seen = [], set()
        for rx in dual_search.gen(n, rmax, maxmol):
            key = tuple(sorted(rx))
            if key in seen:
                continue
            seen.add(key)
            N = Network(n, key)
            if eclib.is_wf_ec(N) and not eclib.is_balanced(N):
                nets.append(key)
        tasks = [(rx, n, 99 + i) for i, rx in enumerate(nets)]
        worker = worker_U
    tally, rep, fail = {}, {}, {}
    nb = 0
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, tasks):
            for rx, res, reps, balanced in chunk:
                nb += not balanced
                for k, v in res.items():
                    h, f = tally.get(k, (0, 0))
                    tally[k] = (h + v, f + (not v))
                    if not v and k not in fail:
                        fail[k] = rx
                for k, v in reps.items():
                    h, f = rep.get(k, (0, 0))
                    rep[k] = (h + v, f + (not v))
                    if not v and k not in fail:
                        fail[k] = rx
    bad = 0
    print(f"# WF-EC networks tested: {sum(tally.get('Ustar', (0, 0)))} (unbalanced: {nb})")
    for k, (h, f) in sorted(tally.items()):
        bad += f > 0
        extra = f"  failure: {Network(n, fail[k]).describe()}" if f else ""
        print(f"{'PASS' if f == 0 else 'FAIL'}  {k}: {h}/{h + f}{extra}")
    for k, (h, f) in sorted(rep.items()):
        extra = f"  first failure: {Network(n, fail[k]).describe()}" if f else ""
        print(f"REPORT  {k}: {h}/{h + f}{extra}")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
