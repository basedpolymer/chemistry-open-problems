"""Q3's Corollary BH (round2/ec/candidates.md, 19:16): block-homogeneous networks generated directly.

Environment: Python 3.12 (WSL Fedora 44), standard library + msolve 0.10.1 (uniqueness counts); exact arithmetic.
Usage (from round2/computation/): python scripts/bh_tests.py [COUNT] [SEED] [NPROC]
Generator (seeded, not uniform): n in {3, 4}; a random partition of the species into blocks B_1..B_m; 2..5
reactions, each with source a*1_B (B a random block, a in {1, 2}, total molecularity <= 4) and a random nonzero
product of molecularity <= 4 different from the source.  Kept: every block is used by a source, connected, (H1),
(H2), s = m - 1 (m = number of blocks), distinct canonical forms.
Claims (proved by Q3): K* is WF-EC; no critical siphon in K nor K*; K absolutely consistent (exact verdict when
decided); at 3 random positive points (not steady states) and 2 random rate vectors the Jacobian restricted to S is
Hurwitz (det(tI - J) = t^(n-s) q(t), q Hurwitz, exact); msolve: exactly one positive steady state in 2 random
(k, class) samples.
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


def gen(rng):
    n = rng.choice((3, 4))
    sp = list(range(n))
    rng.shuffle(sp)
    blocks = []
    while sp:
        k = rng.randint(1, len(sp))
        blocks.append(tuple(sorted(sp[:k])))
        sp = sp[k:]
    rx = set()
    prev = None
    for _ in range(rng.randint(2, 5)):
        if prev is not None and rng.random() < 0.5:
            B = prev[0]
            a = 3 - prev[1]  # same block, the other multiplicity: a varies within the block (not EC)
        else:
            B = rng.choice(blocks)
            a = rng.choice((1, 2))
        if a * len(B) > 4:
            a = 1
        prev = (B, a)
        src = tuple(a if i in B else 0 for i in range(n))
        while True:
            prod = [0] * n
            for _ in range(rng.randint(1, 4)):
                prod[rng.randrange(n)] += 1
            prod = tuple(prod)
            if prod != src:
                break
        rx.add((src, prod))
    rx = tuple(sorted(rx))
    used = {tuple(i for i in range(n) if a[i]) for a, b in rx}
    if any(B not in used for B in blocks):
        return None
    if not e2lib.petri_connected(n, rx):
        return None
    N = CNet(n, rx)
    if not (N.is_conservative() and N.is_consistent()) or N.rank() != len(blocks) - 1:
        return None
    return n, rx, len(blocks)


def check(args):
    n, rx, m, seed = args
    rng = random.Random(seed)
    N = CNet(n, rx)
    D = N.dual()
    G = N.Gamma
    r = len(rx)
    s = N.rank()
    out = dict(n=n, rx=rx)
    drx = e2lib.dual_rx(n, rx)
    out["dual_wfec"] = (e2lib.is_EC(drx) and e2lib.petri_connected(r, drx) and D.is_conservative()
                        and D.is_consistent() and s == e2lib.n_sources(drx) - 1)
    out["crit"] = bool(ecore.minimal_masks(ecore.siphon_masks(N))) and any(
        not any(set(F) <= set(ecore.to_set(S)) for F in ecore.semiflow_supports(N))
        for S in ecore.minimal_masks(ecore.siphon_masks(N)))
    out["critD"] = any(not any(set(F) <= set(ecore.to_set(S)) for F in ecore.semiflow_supports(D))
                       for S in ecore.minimal_masks(ecore.siphon_masks(D)))
    out["v"] = k12_abscons.verdict(n, rx)[0]
    ok = True
    for _ in range(2):
        k = [Fraction(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(r)]
        for _ in range(3):
            x = [Fraction(rng.randint(1, 20), rng.randint(1, 20)) for _ in range(n)]
            J = [[sum(Fraction(G[i][t]) * k[t] * ms_dual.mono(x, rx[t][0]) * rx[t][0][j] / x[j] for t in range(r))
                  for j in range(n)] for i in range(n)]
            cp = q3_tests.charpoly(J)
            if any(c != 0 for c in cp[s + 1:]) or not q3_tests.hurwitz(cp[:s + 1]):
                ok = False
    out["hurwitz"] = ok
    cnts = []
    for _ in range(2):
        k = [Fraction(rng.randint(1, 20), rng.randint(1, 5)) for _ in range(r)]
        x0 = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(n)]
        cnts.append(ms_dual.msolve_count(r, drx, k, x0))
    out["counts"] = cnts
    return out


def main():
    count = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    t0 = time.time()
    print(f"# bh_tests.py count={count} seed={seed}; Python {platform.python_version()} on {platform.system()}; "
          f"exact; msolve 0.10.1")
    rng = random.Random(seed)
    seen = set()
    jobs = []
    tries = 0
    while len(jobs) < count and tries < 2000 * count:
        tries += 1
        g = gen(rng)
        if g is None:
            continue
        n, rx, m = g
        key = (n, ecore.canonical(n, rx))
        if key in seen:
            continue
        seen.add(key)
        jobs.append((n, rx, m, rng.randint(0, 10 ** 9)))
    print(f"# {len(jobs)} block-homogeneous networks (connected, (H1), (H2), s = #blocks - 1) after {tries} draws")
    with Pool(nproc) as pool:
        res = pool.map(check, jobs, chunksize=5)
    ok = True
    tests = [("dual WF-EC", lambda q: q["dual_wfec"]), ("no critical siphon in K", lambda q: not q["crit"]),
             ("no critical siphon in K*", lambda q: not q["critD"]),
             ("absolutely consistent (decided cases)", lambda q: q["v"] == "OPEN" or q["v"].startswith("AC")),
             ("Jacobian on S Hurwitz at 6 random points", lambda q: q["hurwitz"]),
             ("msolve: exactly one positive steady state (2 samples)",
              lambda q: all(c[0] == 1 for c in q["counts"] if c[0] is not None))]
    for name, f in tests:
        bad = [q for q in res if not f(q)]
        ok &= not bad
        print(f"{'PASS' if not bad else 'FAIL'} {name}: {len(res) - len(bad)}/{len(res)}")
        for q in sorted(bad, key=lambda q: (len(q["rx"]), q["rx"]))[:2]:
            print(f"     counterexample: {e2lib.describe(q['n'], q['rx'])}  {q.get('counts')}")
    vc = {}
    for q in res:
        vc[q["v"]] = vc.get(q["v"], 0) + 1
    print(f"# absolute-consistency verdicts: {vc}; EC among them: {sum(e2lib.is_EC(q['rx']) for q in res)}")
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
