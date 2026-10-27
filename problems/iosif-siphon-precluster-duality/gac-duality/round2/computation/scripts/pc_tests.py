"""Final re-test of Q3-PC (round2/ec/candidates.md, 19:49): well-formed proportional-conflict networks.

Environment: Python 3.12 (WSL Fedora 44), standard library + msolve 0.10.1; exact arithmetic.
Usage (from round2/computation/): python scripts/pc_tests.py [COUNT] [SEED] [NPROC]
PC: two reactions sharing a reactant species have positively proportional sources.  Tested networks: Q3's example
NP = {A+2B -> C, 2A+4B -> C+D, C -> A+2B, D -> A+2B}, and seeded random PC networks (n in {3, 4}; species blocks;
per block a primitive vector p with entries in {1, 2}; sources a*p, a in {1, 2, 3}, molecularity <= 8; products
sums of one or two block vectors c*p', c in {1, 2}, molecularity <= 8) that are connected, (H1), (H2), s = #clusters - 1 and neither EC nor
block-homogeneous (some block has a non-constant p).  Claims (proved by Q3), on N and on N*:
  no critical siphon; absolutely consistent when the exact verdict is decided; the class Jacobian is Hurwitz at 3
  random positive points for 2 random rate vectors (det(tI - J) = t^(dim - s) q(t), q Hurwitz); msolve finds exactly
  one positive steady state in 2 random (rates, class) samples.
"""
import sys
import time
import random
import platform
from fractions import Fraction
from multiprocessing import Pool

import e2lib
from e2lib import CNet, ecore
import q3_tests
import k12_abscons
import ms_dual


def is_pc(rx):
    for (a1, _), (a2, _) in ((p, q) for i, p in enumerate(rx) for q in rx[i + 1:]):
        if any(x and y for x, y in zip(a1, a2)):
            # positively proportional: a1 * c = a2 * d for positive c, d
            i0 = next(i for i in range(len(a1)) if a1[i])
            if any(a1[i] * a2[i0] != a2[i] * a1[i0] for i in range(len(a1))):
                return False
    return True


def is_bh(rx):
    return all(len({x for x in a if x}) <= 1 for a, b in rx)


def side_ok(p, rx, rng):
    X = CNet(p, rx)
    G = X.Gamma
    r = len(rx)
    s = X.rank()
    out = {}
    mins = [ecore.to_set(S) for S in ecore.minimal_masks(ecore.siphon_masks(X))]
    sf = ecore.semiflow_supports(X)
    out["crit"] = any(not any(set(F) <= set(S) for F in sf) for S in mins)
    out["v"] = k12_abscons.verdict(p, rx)[0]
    hw = True
    for _ in range(2):
        k = [Fraction(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(r)]
        for _ in range(3):
            x = [Fraction(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(p)]
            J = [[sum(Fraction(G[i][t]) * k[t] * ms_dual.mono(x, rx[t][0]) * rx[t][0][j] / x[j] for t in range(r))
                  for j in range(p)] for i in range(p)]
            cp = q3_tests.charpoly(J)
            if any(c != 0 for c in cp[s + 1:]) or not q3_tests.hurwitz(cp[:s + 1]):
                hw = False
    out["hurwitz"] = hw
    cnt = []
    drx = e2lib.dual_rx(p, rx)
    for _ in range(2):
        k = [Fraction(rng.randint(1, 20), rng.randint(1, 5)) for _ in range(r)]
        x0 = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(p)]
        cnt.append(ms_dual.msolve_count(r, drx, k, x0)[0])
    out["counts"] = cnt
    return out


def check(args):
    n, rx, seed = args
    rng = random.Random(seed)
    return dict(n=n, rx=rx, N=side_ok(n, rx, rng), D=side_ok(len(rx), e2lib.dual_rx(n, rx), rng))


def gen(rng):
    n = rng.choice((3, 4))
    sp = list(range(n))
    rng.shuffle(sp)
    blocks = []
    while sp:
        k = rng.randint(1, len(sp))
        blocks.append(tuple(sorted(sp[:k])))
        sp = sp[k:]
    vecs = {}
    for B in blocks:
        v = [0] * n
        for i in B:
            v[i] = rng.choice((1, 2))
        if all(v[i] == 2 for i in B):
            for i in B:
                v[i] = 1
        vecs[B] = v
    rx = set()
    for _ in range(rng.randint(2, 5)):
        B = rng.choice(blocks)
        a = rng.choice((1, 2, 3))
        src = tuple(a * x for x in vecs[B])
        if sum(src) > 8:
            src = tuple(vecs[B])
        prod = [0] * n
        for _ in range(rng.randint(1, 2)):
            B2 = rng.choice(blocks)
            c = rng.choice((1, 1, 2))
            for i in range(n):
                prod[i] += c * vecs[B2][i]
        if sum(prod) > 8:
            continue
        prod = tuple(prod)
        if prod != src:
            rx.add((src, prod))
    rx = tuple(sorted(rx))
    if len(rx) < 2 or not is_pc(rx) or e2lib.is_EC(rx) or is_bh(rx):
        return None
    if not all(any(a[i] or b[i] for a, b in rx) for i in range(n)) or not e2lib.petri_connected(n, rx):
        return None
    N = CNet(n, rx)
    if not (N.is_conservative() and N.is_consistent()):
        return None
    ncl = len(ecore.petri_clusters(N)[0])
    if N.rank() != ncl - 1:
        return None
    return n, rx


def main():
    count = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    t0 = time.time()
    print(f"# pc_tests.py count={count} seed={seed}; Python {platform.python_version()} on {platform.system()}; "
          f"exact; msolve 0.10.1")
    NP = (((1, 2, 0, 0), (0, 0, 1, 0)), ((2, 4, 0, 0), (0, 0, 1, 1)), ((0, 0, 1, 0), (1, 2, 0, 0)),
          ((0, 0, 0, 1), (1, 2, 0, 0)))
    N = CNet(4, NP)
    okNP = (is_pc(NP) and not e2lib.is_EC(NP) and not is_bh(NP) and e2lib.petri_connected(4, NP)
            and N.is_conservative() and N.is_consistent() and N.rank() == len(ecore.petri_clusters(N)[0]) - 1)
    print(f"{'PASS' if okNP else 'FAIL'} NP hypotheses: PC, not EC, not BH, connected, (H1), (H2), s = #clusters - 1")
    rng = random.Random(seed)
    seen = set()
    jobs = [(4, NP, 1)]
    tries = 0
    while len(jobs) < count + 1 and tries < 3000 * count:
        tries += 1
        g = gen(rng)
        if g is None:
            continue
        key = (g[0], ecore.canonical(g[0], g[1]))
        if key in seen:
            continue
        seen.add(key)
        jobs.append((g[0], g[1], rng.randint(0, 10 ** 9)))
    print(f"# NP + {len(jobs) - 1} random PC networks (not EC, not BH) after {tries} draws")
    with Pool(nproc) as pool:
        res = pool.map(check, jobs, chunksize=2)
    ok = okNP
    for name, f in (("no critical siphon", lambda o: not o["crit"]),
                    ("absolutely consistent (decided)", lambda o: o["v"] == "OPEN" or o["v"].startswith("AC")),
                    ("class Jacobian Hurwitz at 6 random points", lambda o: o["hurwitz"]),
                    ("exactly one positive steady state (msolve, 2 samples)",
                     lambda o: all(c == 1 for c in o["counts"] if c is not None))):
        for sd in ("N", "D"):
            bad = [q for q in res if not f(q[sd])]
            ok &= not bad
            print(f"{'PASS' if not bad else 'FAIL'} {name}, side {sd}: {len(res) - len(bad)}/{len(res)}")
            for q in bad[:2]:
                print(f"     counterexample: N = {e2lib.describe(q['n'], q['rx'])}; {q[sd]}")
    print(f"# NP: {res[0]}")
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
