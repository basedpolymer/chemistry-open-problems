"""Tests of agent Q1's candidates (round2/persistence/candidates.md), independent of Q1's code where possible.

Environment: Python 3.12 (WSL Fedora 44); exact rational arithmetic for every claim (standard library); numpy/scipy
only for the numerical evidence of part G (floating point, reported, never a PASS/FAIL).  Imports Q1's
s1lib.classify and drainsets.drain_sets (unchanged) only to compare with them; my own s = 1 decision is
s1_persistence.analyse.
Usage (from round2/computation/): python scripts/q1_tests.py [--random K] [--seed S]

  A  Q1-1 / Q1-9: on every first-round network with s = 1 and (H1), (H2), my decision of N* (ALL / SOME / HYPERPLANE /
     NONE) and Q1's classify(N*) agree (ALL <-> ALL, NONE <-> NONE, SOME <-> SOME, HYPERPLANE <-> TWO-SIDED); counts.
  B  Q1-2 beyond the records: K seeded random weakly reversible networks with s = 1 and (H1)-(H3) (n = 2..4, up to 3
     linkage classes, molecularity <= 5, 2..5 complexes per class): N* is never SOME or HYPERPLANE (my decision) and
     Q1's classify never returns SOME / TWO-SIDED.
  C  Q1-3: n = 2, s = 1, (H1)-(H3): N* is ALL iff every source has strictly more of its consumed species than of
     its produced species.
  D  Q1-4: s = 1, (H1), (H2): if no reaction of N has a net product in its source, N* is ALL.
  E  Q1-6 (conjecture): N weakly reversible, (H1)-(H3), pure: N* has no drain species and no drain set (Q1's
     drain_sets, first found), and for s = 1 my decision is ALL.  REFUTED if a counterexample is found.
  F  Q1-5 consequence: a drain set Z of N* contains no support of a semi-positive conservation law of N*.
  G  Q1-5 numerics (evidence): for the first 30 drain-set duals, random kappa, x_Z(0) = 1e-3, x(0) = 1 elsewhere;
     report how many have max_Z x(T) < 1e-3 at T = 1e3 (scipy LSODA).
"""
import os
import sys
import time
import random
import platform
from fractions import Fraction

import e2lib
from e2lib import CNet, endo
import s1_persistence as S1

QDIR = os.path.normpath(os.path.join(e2lib.HERE, "..", "..", "persistence", "scripts"))
sys.path.insert(0, QDIR)
import s1lib  # noqa: E402
import drainsets  # noqa: E402

MAP = {"ALL": "ALL", "NONE": "NONE", "SOME": "SOME", "HYPERPLANE": "TWO-SIDED"}


def random_s1_wr(rng):
    n = rng.randint(2, 4)
    while True:
        a = [rng.randint(-3, 3) for _ in range(n)]
        if any(x > 0 for x in a) and any(x < 0 for x in a) and len(set(a)) == n:
            break
    from math import gcd
    g = 0
    for x in a:
        g = gcd(g, abs(x))
    a = [x // g for x in a]
    rx = []
    used = set()
    for _ in range(rng.randint(1, 3)):
        y0 = [rng.randint(0, 3) for _ in range(n)]
        levels = [t for t in range(-6, 7) if all(0 <= y0[i] + t * a[i] for i in range(n))
                  and 0 < sum(y0[i] + t * a[i] for i in range(n)) <= 5]
        if len(levels) < 2:
            continue
        k = rng.randint(2, min(5, len(levels)))
        lv = sorted(rng.sample(levels, k))
        cs = [tuple(y0[i] + t * a[i] for i in range(n)) for t in lv]
        if any(c in used for c in cs):
            continue
        used.update(cs)
        order = cs[:]
        rng.shuffle(order)
        arcs = {(order[i], order[(i + 1) % k]) for i in range(k)}  # a Hamiltonian cycle: strongly connected
        for _ in range(rng.randint(0, k)):
            u, v = rng.sample(cs, 2)
            arcs.add((u, v))
        rx.extend(arcs)
    if not rx:
        return None
    rx = tuple(sorted(rx))
    if not all(any(p[i] or q[i] for p, q in rx) for i in range(n)):
        return None
    N = CNet(n, rx)
    if N.rank() != 1 or not (N.is_conservative() and N.is_consistent() and N.rows_distinct()):
        return None
    return n, rx


def main():
    t0 = time.time()
    argv = sys.argv[1:]
    K = int(argv[argv.index("--random") + 1]) if "--random" in argv else 3000
    seed = int(argv[argv.index("--seed") + 1]) if "--seed" in argv else 5
    print(f"# q1_tests.py  Python {platform.python_version()} on {platform.system()}; exact; random {K}, seed {seed}")
    ok = True
    recs = e2lib.load(e2lib.E_FILES)
    s1 = [rc for rc in recs if rc["s"] == 1 and rc["N"]["cons"] and rc["N"]["consist"]]
    # A
    agree = {}
    bad = []
    mine = {}
    for rc in s1:
        n, rx = rc["n"], rc["rx"]
        m = S1.analyse(n, rx)["cat"]
        mine[(n, rx)] = m
        q, _ = s1lib.classify(len(rx), e2lib.dual_rx(n, rx))
        agree[(m, q)] = agree.get((m, q), 0) + 1
        if MAP.get(m) != q:
            bad.append((m, q, rx))
    ok &= not bad
    print(f"{'PASS' if not bad else 'FAIL'} A (Q1-1/Q1-9): my s = 1 decision vs Q1's classify on {len(s1)} networks: "
          f"{dict(sorted(agree.items()))}")
    for b in bad[:3]:
        print(f"     disagreement {b}")
    wr = [rc for rc in s1 if e2lib.H(rc) and rc["N"]["WR"]]
    cnt = {}
    for rc in wr:
        c = mine[(rc["n"], rc["rx"])]
        cnt[c] = cnt.get(c, 0) + 1
    print(f"#   WR with (H1)-(H3), s = 1: {len(wr)} networks, my decision {cnt} (Q1 reports 150 ALL / 411 NONE)")
    # B
    rng = random.Random(seed)
    seen = set()
    cats = {}
    qcats = {}
    badB = []
    tries = 0
    while len(seen) < K and tries < 200 * K:
        tries += 1
        g = random_s1_wr(rng)
        if g is None:
            continue
        n, rx = g
        if (n, rx) in seen:
            continue
        seen.add((n, rx))
        c = S1.analyse(n, rx)["cat"]
        q, _ = s1lib.classify(len(rx), e2lib.dual_rx(n, rx))
        cats[c] = cats.get(c, 0) + 1
        qcats[q] = qcats.get(q, 0) + 1
        if c in ("SOME", "HYPERPLANE") or q in ("SOME", "TWO-SIDED") or MAP.get(c) != q:
            badB.append((c, q, n, rx))
    ok &= not badB
    print(f"{'PASS' if not badB else 'FAIL'} B (Q1-2): {len(seen)} random WR networks, s = 1, (H1)-(H3): my decision "
          f"{cats}, Q1's {qcats}")
    for b in badB[:3]:
        print(f"     counterexample {b[0]}/{b[1]}: N = {e2lib.describe(b[2], b[3])}")
    # C
    badC = []
    napp = 0
    for rc in s1:
        if rc["n"] != 2 or not e2lib.H(rc):
            continue
        napp += 1
        n, rx = rc["n"], rc["rx"]
        cond = True
        for a, b in rx:
            cons = [i for i in range(2) if b[i] < a[i]]
            prod = [i for i in range(2) if b[i] > a[i]]
            if not (cons and prod and a[cons[0]] > a[prod[0]]):
                cond = False
        if cond != (mine[(n, rx)] == "ALL"):
            badC.append(rx)
    ok &= not badC
    print(f"{'PASS' if not badC else 'FAIL'} C (Q1-3): n = 2 criterion vs exact decision ({napp} networks, "
          f"{len(badC)} failures)")
    for b in badC[:3]:
        print(f"     failure: {e2lib.describe(2, b)}")
    # D
    badD = []
    napp = 0
    for rc in s1:
        n, rx = rc["n"], rc["rx"]
        if any(a[i] > 0 and b[i] > a[i] for a, b in rx for i in range(n)):
            continue
        napp += 1
        if mine[(n, rx)] != "ALL":
            badD.append(rx)
    ok &= not badD
    print(f"{'PASS' if not badD else 'FAIL'} D (Q1-4): no net product in a source => ALL ({napp} networks, "
          f"{len(badD)} failures)")
    # E, F
    allrecs = recs + e2lib.load([os.path.join(e2lib.DATA_E2, f) for f in ("sample_5_3_5_8_s3.pkl.gz",
                                                                           "sample_6_3_5_8_s4.pkl.gz")])
    wrh = [rc for rc in allrecs if e2lib.H(rc) and rc["N"]["WR"]]
    pure = [rc for rc in wrh if not any(a[i] and b[i] for a, b in rc["rx"] for i in range(rc["n"]))]
    refE = []
    adls = 0
    for rc in pure:
        n, rx = rc["n"], rc["rx"]
        d = e2lib.dual_rx(n, rx)
        if not rc["D"]["CS"]:
            adls += 1
        bad = bool(endo.drain_species(len(rx), d)) or bool(drainsets.drain_sets(len(rx), d, first_only=True))
        if rc["s"] == 1 and mine.get((n, rx)) not in (None, "ALL"):
            bad = True
        if bad:
            refE.append(rc)
    by_s = {}
    for rc in pure:
        by_s[rc["s"]] = by_s.get(rc["s"], 0) + 1
    print(f"{'HOLDS' if not refE else 'REFUTED'} E (Q1-6, conjecture): {len(pure)} pure WR networks with (H1)-(H3) "
          f"(by s: {dict(sorted(by_s.items()))}); N* without critical siphon (ADLS => ALL): {adls}; "
          f"counterexamples (drain species / drain set / s = 1 not ALL): {len(refE)}")
    refE.sort(key=lambda rc: (rc["n"] + rc["r"], rc["r"], rc["rx"]))
    for rc in refE[:2]:
        print(f"     counterexample: N = {e2lib.describe(rc['n'], rc['rx'])}; N* = "
              f"{e2lib.describe(rc['r'], e2lib.dual_rx(rc['n'], rc['rx']), [f'k{j+1}' for j in range(rc['r'])])}")
    nds = 0
    badF = []
    witnesses = []
    for rc in wrh:
        n, rx = rc["n"], rc["rx"]
        d = e2lib.dual_rx(n, rx)
        ds = drainsets.drain_sets(len(rx), d, first_only=True)
        if not ds:
            continue
        nds += 1
        Z = set(ds[0][0])
        if any(set(F) <= Z for F in rc["D"]["SF"]):
            badF.append(rx)
        if len(witnesses) < 30:
            witnesses.append((len(rx), d, sorted(Z)))
    ok &= not badF
    print(f"{'PASS' if not badF else 'FAIL'} F (Q1-5 consequence): {nds} WR duals with a drain set; a drain set never "
          f"contains a semiflow support of N* ({len(badF)} failures)")
    # G numerics
    try:
        import numpy as np
        from scipy.integrate import solve_ivp
        dec = 0
        tot = 0
        for p, d, Z in witnesses:
            r_ = len(d)
            kap = np.exp(np.random.default_rng(len(d) + p).uniform(np.log(1e-2), np.log(1e2), r_))
            S_ = np.array([d[j][0] for j in range(r_)], dtype=float)
            D_ = np.array([np.array(d[j][1]) - np.array(d[j][0]) for j in range(r_)], dtype=float)

            def rhs(t, x):
                x = np.maximum(x, 0.0)
                rates = kap * np.prod(x ** S_, axis=1)
                return rates @ D_
            x0 = np.ones(p)
            x0[Z] = 1e-3
            sol = solve_ivp(rhs, (0, 1e3), x0, method="LSODA", rtol=1e-10, atol=1e-14)
            tot += 1
            if sol.success and max(sol.y[Z, -1]) < 1e-3:
                dec += 1
        print(f"# G (Q1-5 numerics, evidence only): {dec}/{tot} drain-set duals have max_Z x(1e3) < x_Z(0) = 1e-3")
    except Exception as e:  # noqa: BLE001
        print(f"# G skipped: {e}")
    print(f"{'PASS' if ok else 'FAIL'}: all claims (E is a conjecture test, not a claim)")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
