"""Re-derivation of the coordinator's facts 1, 2, 5 and 6 (BRIEF-exploration.md, "Facts established here").

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (crn.py via ecore.py).
Usage: python facts.py [--nproc K] [WR_RECORD_FILES...]

Part A (all networks, no hypothesis): every reaction set of r distinct reactions (source != product, nonzero
complexes of molecularity <= m), every species occurring, one representative per species-permutation orbit, for
the families (n, m, rmax) in FAMILIES_ALL.  For each network N, computed independently on N and on N*:
  F1   Gamma(N*) = Gamma(N)^T; N* has r species, n reactions, the same rank s; N** = N; dim ker Gamma(N*) = n - s.
  F5a  S is a siphon of N  <=>  S, as a set of reactions of N*, is a "T-trap": every species of N* that is a
       product of a reaction of S is a reactant of some reaction of S.  (All nonempty S.)
  F5b  S is a trap of N  <=>  S, as a set of reactions of N*, is a "T-siphon": every species of N* that is a
       reactant of a reaction of S is a product of some reaction of S.
  F5c  supports of extreme semi-positive conservation laws of N = internal cycles of N* (supports of the extreme
       rays of the flux cone of N*), so the Angeli-De Leenheer-Sontag condition "every siphon of N contains a
       semiflow support" is "every T-trap of N* contains an internal cycle of N*".
Part B (weakly reversible networks; record files written by enum_wr.py): counts for facts 2 and 6, by family:
  F2   conservative WR networks, and how many have a WR dual (all deficiencies; deficiency 0);
  F6   conservative WR networks with a critical minimal siphon (containing no semiflow support), and those whose
       minimal siphons are exactly the semiflow supports.
  Also the coordinator's families (n=2,m<=2,r<=4; n=3,m=1,r<=6; n=3,m<=2,r<=3) are regenerated here to compare
  with the coordinator's counts 32 / 9 / 24 / 8 / 9 / 23.
Part C: F0, the WR generator of ecore.py reproduces the numbers of WR networks with (H1)-(H3) found by the
  independent subset enumeration of result 17 (enumerate_small.py logs) on its five families.
"""
import sys
import time
import gzip
import pickle
import platform
from itertools import combinations, permutations
from multiprocessing import Pool

import ecore
from ecore import CNet, canonical, complexes, Network, siphon_masks, trap_masks, rank, record
from ecore import lc_instances, wr_networks_from

FAMILIES_ALL = [(2, 2, 5), (2, 3, 4), (3, 1, 6), (3, 2, 4), (4, 1, 6), (4, 2, 3)]


def _masks_sets(net):
    src = [frozenset(i for i in range(net.n) if a[i] > 0) for a, b in net.reactions]
    prd = [frozenset(i for i in range(net.n) if b[i] > 0) for a, b in net.reactions]
    return src, prd


def check_all(args):
    n, m, rmax, first = args
    C = complexes(n, m)
    RX = sorted((a, b) for a in C for b in C if a != b)
    rest = RX[first + 1:]
    cnt = 0
    fails = []
    for k in range(rmax):
        for tail in combinations(rest, k):
            rx = (RX[first],) + tail
            if not all(any(a[i] or b[i] for a, b in rx) for i in range(n)):
                continue
            if canonical(n, rx) != tuple(sorted(rx)):
                continue
            cnt += 1
            N = CNet(n, rx)
            D = N.dual()
            G = N.Gamma
            s = rank(G)
            # F1
            ok1 = (D.n == N.r and D.r == N.n and D.Gamma == [list(c) for c in zip(*G)] and rank(D.Gamma) == s
                   and Network.dual(D).reactions == N.reactions
                   and len(ecore.crn.nullspace(D.Gamma, D.r)) == n - s)
            if not ok1:
                fails.append(("F1", rx))
            # F5a/F5b, computed on N* directly (species of N* = reactions of N)
            sN = set(siphon_masks(N))
            tN = set(trap_masks(N))
            dsrc, dprd = _masks_sets(D)  # D reaction i (= species i of N): reactant/product species of N*
            for S in range(1, 1 << n):
                Sset = [i for i in range(n) if S >> i & 1]
                produced = set().union(*(dprd[i] for i in Sset))
                consumed = set().union(*(dsrc[i] for i in Sset))
                ttrap = produced <= consumed
                tsiph = consumed <= produced
                if (S in sN) != ttrap:
                    fails.append(("F5a", rx))
                    break
                if (S in tN) != tsiph:
                    fails.append(("F5b", rx))
                    break
            # F5c
            sf = set(frozenset(i for i in range(n) if z[i] > 0) for z in Network.conservation_rays(N))
            ic = set(frozenset(j for j in range(D.r) if v[j] > 0) for v in Network.flux_rays(D))
            if sf != ic:
                fails.append(("F5c", rx))
    return cnt, fails


def coordinator_family_records(nproc):
    """WR networks of the coordinator's three families, as records."""
    out = {}
    for (n, m, R) in [(2, 2, 4), (3, 1, 6), (3, 2, 3)]:
        C = complexes(n, m)
        inst = lc_instances(C, R)
        reps = set()
        for f in range(len(inst)):
            for rx in wr_networks_from(f, inst, R, n):
                if canonical(n, rx) == rx:
                    reps.add(rx)
        out[(n, m, R)] = [record(rx, n) for rx in sorted(reps)]
    return out


def counts_f2_f6(recs):
    cons = [r for r in recs if r["N"]["cons"]]
    consH = [r for r in cons if r["N"]["rows_distinct"]]
    def row(rs):
        return dict(
            nets=len(rs),
            WRdual=sum(r["D"]["WR"] for r in rs),
            def0=sum(r["N"]["deficiency"] == 0 for r in rs),
            def0_WRdual=sum(r["N"]["deficiency"] == 0 and r["D"]["WR"] for r in rs),
            crit=sum(len(r["N"]["CS"]) > 0 for r in rs),
            MS_eq_SF=sum(r["N"]["siphons_eq_SF"] for r in rs),
            no_crit=sum(len(r["N"]["CS"]) == 0 for r in rs),
        )
    return row(cons), row(consH)


def smallest(rs, pred):
    best = None
    for r in rs:
        if pred(r):
            key = (r["n"] + r["r"], r["r"], r["rx"])
            if best is None or key < best[0]:
                best = (key, r)
    return best[1] if best else None


def main():
    argv = sys.argv[1:]
    nproc = 8
    if "--nproc" in argv:
        i = argv.index("--nproc")
        nproc = int(argv[i + 1])
        del argv[i:i + 2]
    files = argv
    t0 = time.time()
    print(f"# facts.py  Python {platform.python_version()} on {platform.system()} {platform.release()}, nproc={nproc}")
    allok = True
    # ---------------- Part A
    print("# Part A: all networks (no hypothesis), facts 1 and 5")
    totA = 0
    failsA = []
    with Pool(nproc) as pool:
        for (n, m, rmax) in FAMILIES_ALL:
            C = complexes(n, m)
            nR = len(C) * (len(C) - 1)
            cnt = 0
            for c, f in pool.imap_unordered(check_all, [(n, m, rmax, i) for i in range(nR)]):
                cnt += c
                failsA.extend(f)
            totA += cnt
            print(f"#   family n={n} m<={m} r<={rmax}: {cnt} networks (orbits), cumulative failures {len(failsA)}")
            sys.stdout.flush()
    for tag, doc in [("F1", "Gamma(N*)=Gamma^T, n<->r, same s, N**=N, dim ker Gamma(N*) = n-s"),
                     ("F5a", "siphons of N = T-traps of N* (produced <= consumed)"),
                     ("F5b", "traps of N = T-siphons of N* (consumed <= produced)"),
                     ("F5c", "semiflow supports of N = internal cycles of N*")]:
        bad = [x for x in failsA if x[0] == tag]
        ok = not bad
        allok &= ok
        print(f"{'PASS' if ok else 'FAIL'} {tag}: {doc}, on all {totA} networks of Part A"
              + (f"; first failure {bad[0][1]}" if bad else ""))
    # ---------------- Part C: the WR generator of ecore against result 17
    print("# Part C: WR generator (ecore.lc_instances / wr_networks_from) against the logs of result 17")
    expectedWR = {(2, 2, 5): 16, (3, 2, 4): 112, (3, 3, 3): 123, (4, 1, 12): 84, (4, 2, 3): 4}
    okC = True
    for (n, m, R), exp in expectedWR.items():
        C = complexes(n, m)
        inst = lc_instances(C, R)
        cnt = 0
        for f in range(len(inst)):
            for rx in wr_networks_from(f, inst, R, n):
                if canonical(n, rx) == rx:
                    N = CNet(n, rx)
                    if N.is_conservative() and N.is_consistent() and N.rows_distinct():
                        cnt += 1
        okC &= cnt == exp
        print(f"#   (n,m,R)=({n},{m},{R}): WR networks with (H1)-(H3): {cnt}; result 17: {exp}")
    allok &= okC
    print(f"{'PASS' if okC else 'FAIL'} F0: the WR generator reproduces the WR counts of result 17 on its five families")
    # ---------------- Part B
    print("# Part B: weakly reversible networks, facts 2 and 6 (conservative = (H1); 'H1+H3' adds distinct rows)")
    coord = coordinator_family_records(nproc)
    union = {}
    for fam, rs in coord.items():
        for r in rs:
            union[(r["n"], r["rx"])] = r
    ru = list(union.values())
    a, b = counts_f2_f6(ru)
    sumfam = [counts_f2_f6(rs) for rs in coord.values()]
    print(f"#   coordinator's families, union of orbits ({len(ru)} WR networks):")
    print(f"#     (H1):    {a}")
    print(f"#     (H1,H3): {b}")
    tot_sum = {k: sum(x[0][k] for x in sumfam) for k in a}
    print(f"#     (H1), families summed with multiplicity: {tot_sum}")
    target = dict(nets=32, WRdual=9, def0=24, def0_WRdual=8, crit=9, MS_eq_SF=23)
    matches = [name for name, d in [("union H1", a), ("union H1H3", b), ("summed H1", tot_sum)]
               if all(d[k] == v for k, v in target.items())]
    print(f"# coordinator's counts 32/9/24/8 (fact 2) and 9/23 (fact 6) reproduced by: {matches or 'none'}")
    # smallest WR network with non-WR dual
    sm = smallest(ru, lambda r: r["N"]["cons"] and not r["D"]["WR"])
    if sm:
        print(f"#   smallest conservative WR network with non-WR dual: {ecore.describe(sm['rx'], sm['n'])}"
              f"  dual: {ecore.describe(CNet(sm['n'], sm['rx']).dual().reactions, sm['r'], [f'k{j+1}' for j in range(sm['r'])])}")
    ok = sm is not None
    allok &= ok
    print(f"{'PASS' if ok else 'FAIL'} F2: weak reversibility is not preserved by N -> N* (witness above)")
    # larger families from record files
    for fn in files:
        with gzip.open(fn, "rb") as fh:
            recs = pickle.load(fh)
        fam = recs[0]["family"] if recs else None
        a, b = counts_f2_f6(recs)
        print(f"#   {fn}: family {fam}, {len(recs)} WR networks")
        print(f"#     (H1):    {a}")
        print(f"#     (H1,H3): {b}")
        bys = {}
        for r in recs:
            if r["N"]["cons"] and r["N"]["rows_distinct"]:
                d = bys.setdefault(r["s"], [0, 0, 0, 0, 0, 0])
                d[0] += 1
                d[1] += r["D"]["WR"]
                d[2] += r["N"]["deficiency"] == 0
                d[3] += r["N"]["deficiency"] == 0 and r["D"]["WR"]
                d[4] += len(r["N"]["CS"]) > 0
                d[5] += r["N"]["siphons_eq_SF"]
        for s in sorted(bys):
            print(f"#     s={s}: (H1,H3) nets {bys[s][0]}, WR dual {bys[s][1]}, def0 {bys[s][2]}, def0&WR dual "
                  f"{bys[s][3]}, critical siphon {bys[s][4]}, minimal siphons = semiflow supports {bys[s][5]}")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
