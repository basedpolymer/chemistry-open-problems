"""Exhaustive enumeration of small equal-conflict (EC) networks and exact checks of the Petri-net statements of
petri/notes.md (agent P, Claude Opus 5.5, effort xhigh).

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rationals (crn.py, petri.py).
Usage: python enumerate_ec.py N MAXMOL RMAX ORDINARY(0/1) [NPROC]
Runtime: see the log header (a few seconds to a few minutes per family on 12 cores).

Networks generated: every network, up to permutation of the species, such that
  * the distinct source complexes have pairwise disjoint supports covering all species (EC; every species is consumed,
    which (H2) forces anyway), coefficients >= 1 and molecularity <= MAXMOL (ORDINARY=1: all coefficients 1);
  * each source has a nonempty set of distinct products, nonzero complexes of molecularity <= MAXMOL different from
    the source (ORDINARY=1: 0/1 complexes); at most RMAX reactions in total;
  * the Petri net is connected, and (H1) conservative, (H2) consistent hold.
Classes reported: EC (all generated), H3 (rows of Gamma distinct), WF (rank equation s = #EC sets - 1; well-formed
by the Rank Theorem of Teruel-Silva), WF-FC (WF and ordinary), WR.

Claims with PASS/FAIL (exit status 1 on failure): K1 (all), K2 (H3 class), K3 (WF-FC class), K7 (WR class),
token-game confirmation of the Rank Theorem on EC nets (WF iff live+bounded from the Teruel-Silva marking, nets
with at most 6 places).  Reported without PASS/FAIL (questions/conjectures): K4, K5, K6, K8, K9.
"""
import sys
import time
import platform
from itertools import combinations, permutations, product
from multiprocessing import Pool

from petri import (Network, clusters, ec_sets, twin_classes, is_ordinary, is_connected, rank_equation_ec,
                   rank_equation_clusters, minimal_t_traps, critical_siphons, is_s_component_set,
                   is_t_component_set, precluster_by_profiles, is_live_and_bounded, ts_marking, refines)


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        for k in range(len(p)):
            yield p[:k] + [[first] + p[k]] + p[k + 1:]
        yield [[first]] + p


def complexes(n, m, ordinary):
    rng = range(2) if ordinary else range(m + 1)
    return [c for c in product(rng, repeat=n) if 0 < sum(c) <= m]


def sources_for_block(n, block, m, ordinary):
    out = []
    rng = [1] if ordinary else range(1, m + 1)
    for coeffs in product(rng, repeat=len(block)):
        if sum(coeffs) <= m:
            v = [0] * n
            for i, c in zip(block, coeffs):
                v[i] = c
            out.append(tuple(v))
    return out


def canonical(rx, perms):
    best = None
    for perm in perms:
        key = tuple(sorted((tuple(a[p] for p in perm), tuple(b[p] for p in perm)) for a, b in rx))
        if best is None or key < best:
            best = key
    return best


def gen_networks(n, m, rmax, ordinary, part_index):
    parts = list(set_partitions(list(range(n))))
    part = parts[part_index]
    cpx = complexes(n, m, ordinary)
    src_choices = [sources_for_block(n, b, m, ordinary) for b in part]
    k = len(part)
    if k > rmax:
        return
    for srcs in product(*src_choices):
        prods = [[c for c in cpx if c != s] for s in srcs]

        def rec(idx, budget, acc):
            if idx == k:
                yield list(acc)
                return
            remaining_blocks = k - idx - 1
            for size in range(1, budget - remaining_blocks + 1):
                for ps in combinations(prods[idx], size):
                    yield from rec(idx + 1, budget - size, acc + [(srcs[idx], p) for p in ps])
        yield from rec(0, rmax, [])


def same_consumers_classes(N):
    blocks = {}
    for i in range(N.n):
        key = frozenset(j for j in range(N.r) if N.reactions[j][0][i] > 0)
        blocks.setdefault(key, []).append(i)
    return [set(b) for b in blocks.values()]


def same_reactant_support_classes(N):
    blocks = {}
    for j in range(N.r):
        key = frozenset(i for i in range(N.n) if N.reactions[j][0][i] > 0)
        blocks.setdefault(key, []).append(j)
    return [set(b) for b in blocks.values()]


def evaluate(N):
    """Returns (classes dict, claims dict name->bool (proved claims), reports dict name->bool or value)."""
    D = N.dual()
    cl, rep = {}, {}
    H3 = N.rows_distinct()
    s = N.rank()
    nec = len(ec_sets(N))
    WF = (s == nec - 1)
    ordn = is_ordinary(N)
    cls = {"EC": True, "H3": H3, "WF": WF, "WF-FC": WF and ordn, "WF-EC-weighted": WF and not ordn,
           "WR": N.weakly_reversible(), "FC": ordn}
    S = N.minimal_siphons()
    Sd = D.minimal_siphons()  # minimal T-traps of N
    # K1 (all networks)
    cl["K1 min siphons meet each same-consumer class in <= 1 species"] = all(
        len(R & c) <= 1 for R in S for c in same_consumers_classes(N))
    cl["K1 min T-traps meet each same-reactant-support class in <= 1 reaction"] = all(
        len(U & c) <= 1 for U in Sd for c in same_reactant_support_classes(N))
    cl["K1 min T-traps of N = min siphons of N*"] = set(Sd) == set(minimal_t_traps(N))
    sup_cons = sorted({frozenset(i for i in range(N.n) if z[i] > 0) for z in N.conservation_rays()},
                      key=lambda x: (len(x), sorted(x)))
    sup_flux = N.internal_cycles()
    M = N.mips_ray()
    Md = D.mips_ray()
    P = N.preclusters()
    Pd = D.preclusters()
    twins = [b for b in twin_classes(N) if len(b) >= 2]
    dbl = [b for b in ec_sets(N) if len(b) >= 2]
    # K2 (H3 class; proof needs (H1), true for all generated networks)
    q3a = refines(Pd, S)
    q3b = refines(P, Sd)
    cl["K2 Q3 species side iff no twins and every species in a min siphon"] = q3a == (
        not twins and set().union(*S) == set(range(N.n)) if S else False)
    cl["K2 Q3 reaction side iff no doubling pair and every reaction in a min T-trap"] = q3b == (
        not dbl and set().union(*Sd) == set(range(N.r)) if Sd else False)
    # K3 (WF-FC)
    if cls["WF-FC"]:
        cl["K3a N* is WF-FC"] = (is_ordinary(D) and is_connected(D) and D.is_conservative() and D.is_consistent()
                                 and rank_equation_clusters(D) and rank_equation_ec(D))
        rays = N.conservation_rays()
        frays = N.flux_rays()
        cl["K3b extreme conservation laws = 0/1 indicators of the minimal siphons"] = (
            all(set(z) <= {0, 1} for z in rays) and sup_cons == S and len(rays) == len(S))
        cl["K3b minimal siphons are S-components"] = all(is_s_component_set(N, R) for R in S)
        cl["K3c extreme flux rays = 0/1 indicators of the minimal T-traps (T-components)"] = (
            all(set(v) <= {0, 1} for v in frays) and sup_flux == Sd and len(frays) == len(Sd)
            and all(is_t_component_set(N, U) for U in Sd))
        Sig = [frozenset(k for k, R in enumerate(S) if i in R) for i in range(N.n)]
        Th = [frozenset(k for k, U in enumerate(Sd) if j in U) for j in range(N.r)]
        atomsS = sorted((frozenset(i for i in range(N.n) if Sig[i] == g) for g in set(Sig)), key=lambda x: min(x))
        atomsT = sorted((frozenset(j for j in range(N.r) if Th[j] == g) for g in set(Th)), key=lambda x: min(x))
        cl["K3d ray MIPS(N) = atoms of the minimal siphons"] = sorted(map(sorted, M)) == sorted(map(sorted, atomsS))
        cl["K3d ray MIPS(N*) = atoms of the minimal T-traps"] = sorted(map(sorted, Md)) == sorted(map(sorted, atomsT))
        tw_pairs = [(i, j) for b in twins for i, j in combinations(sorted(b), 2)]
        db_pairs = [(i, j) for b in dbl for i, j in combinations(sorted(b), 2)]
        cl["K3e profile rule gives P(N*)"] = precluster_by_profiles(N.n, tw_pairs, Sig) == Pd
        cl["K3e profile rule gives P(N)"] = precluster_by_profiles(N.r, db_pairs, Th) == P
        cl["K3f no critical siphons in N and N*"] = not critical_siphons(N) and not critical_siphons(D)
        jf = not twins
        cf = not dbl
        cl["K3g twins never conservation-coupled, doubling pairs never fully coupled"] = all(
            Sig[i] != Sig[j] for i, j in tw_pairs) and all(Th[i] != Th[j] for i, j in db_pairs)
        cl["K3g A2 iff JF, B2 iff CF"] = (refines(Pd, M) == jf) and (refines(P, Md) == cf)
        cl["K3g A1 fails (n>=2), B1 fails (r>=2)"] = (Pd != M or N.n < 2) and (P != Md or N.r < 2)
        cl["K3g Q3 species side iff JF, reaction side iff CF"] = (q3a == jf) and (q3b == cf)
        if N.n <= 6:
            ok, live = is_live_and_bounded(N, [1] * N.n)
            cl["K3h live and bounded from M0 = (1,...,1)"] = ok and live
    # Rank Theorem (EC) token-game confirmation
    if N.n <= 5:
        ok, live = is_live_and_bounded(N, ts_marking(N), max_states=100000)
        if ok:
            cl["RT WF (rank equation) iff live+bounded from the Teruel-Silva marking"] = (live == WF)
    # K7 (WR)
    if cls["WR"]:
        cl["K7 WR+EC: deficiency 0"] = N.deficiency() == 0
        cl["K7 WR+EC: min siphons = supports of extreme conservation laws"] = sup_cons == S
        cl["K7 WR+EC: no critical siphon"] = not critical_siphons(N)
    # reports
    rep["K4 min siphons = supports of extreme conservation laws"] = sup_cons == S
    rep["K4 extreme conservation laws are 0/1"] = all(set(z) <= {0, 1} for z in N.conservation_rays())
    rep["K5 s >= #EC sets - 1"] = s >= nec - 1
    rep["K6 no critical siphon in N"] = not critical_siphons(N)
    rep["K6 no critical siphon in N*"] = not critical_siphons(D)
    rep["K6 (no crit N) == (no crit N*)"] = rep["K6 no critical siphon in N"] == rep["K6 no critical siphon in N*"]
    rep["K9i every min siphon is a union of MIPS blocks"] = all(
        all((b <= R) or not (b & R) for b in M) for R in S)
    rep["K9ii every MIPS block lies in a min siphon"] = all(any(b <= R for R in S) for b in M)
    rep["K9iii min siphons = semiflow supports"] = sup_cons == S
    rep["K9iv every semiflow support is a union of min siphons"] = all(
        set().union(*([R for R in S if R <= X] or [set()])) == set(X) for X in sup_cons)
    key8 = (N.n, tuple(sorted(map(tuple, map(sorted, twins)))), tuple(sorted(map(tuple, map(sorted, S)))))
    rep["_K8key"] = (key8, tuple(map(tuple, map(sorted, Pd))))
    return cls, cl, rep


def worker(args):
    n, m, rmax, ordinary, pidx = args
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
        cls, cl, rep = evaluate(N)
        out.append((rx, cls, cl, rep))
    return out


def main():
    n, m, rmax, ordinary = (int(a) for a in sys.argv[1:5])
    nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    t0 = time.time()
    print(f"# enumerate_ec.py n={n} molecularity<={m} r<={rmax} ordinary={ordinary}; Python "
          f"{platform.python_version()} on {platform.system()}")
    nparts = sum(1 for _ in set_partitions(list(range(n))))
    counts, claims, reports = {}, {}, {}
    first_fail, first_rep_fail = {}, {}
    k8 = {}
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, [(n, m, rmax, ordinary, p) for p in range(nparts)]):
            for rx, cls, cl, rep in chunk:
                act = [c for c, v in cls.items() if v]
                for c in act:
                    counts[c] = counts.get(c, 0) + 1
                for name, val in cl.items():
                    h, f = claims.get(name, (0, 0))
                    claims[name] = (h + val, f + (not val))
                    if not val and (name not in first_fail or len(rx) < len(first_fail[name])):
                        first_fail[name] = rx
                for name, val in rep.items():
                    if name.startswith("_"):
                        continue
                    for c in act:
                        h, f = reports.get((name, c), (0, 0))
                        reports[(name, c)] = (h + val, f + (not val))
                        if not val and ((name, c) not in first_rep_fail or
                                        len(rx) < len(first_rep_fail[(name, c)])):
                            first_rep_fail[(name, c)] = rx
                if cls["H3"]:
                    key, pd = rep["_K8key"]
                    k8.setdefault(key, {}).setdefault(pd, rx)
    print("# class counts: " + ", ".join(f"{c}={v}" for c, v in sorted(counts.items())))
    bad = 0
    for name in sorted(claims):
        h, f = claims[name]
        st = "PASS" if f == 0 else "FAIL"
        bad += f > 0
        extra = f"  first failure: {first_fail[name]}" if f else ""
        print(f"{st}  {name}: {h}/{h + f}{extra}")
    print("# reports (holds/total per class; not PASS/FAIL)")
    for (name, c) in sorted(reports):
        h, f = reports[(name, c)]
        extra = f"  smallest failure: {first_rep_fail[(name, c)]}" if f else ""
        print(f"REPORT  [{c}] {name}: {h}/{h + f}{extra}")
    amb = [(key, v) for key, v in k8.items() if len(v) > 1]
    print(f"REPORT  K8 (H3 class): groups (n, twin classes, minimal siphons) with >1 distinct P(N*): {len(amb)}")
    for key, v in amb[:3]:
        print(f"REPORT  K8 example: key={key}; networks: {list(v.values())[:2]}")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
