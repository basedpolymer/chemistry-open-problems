"""Exhaustive test of precise formulations of Iosif's Conjecture 22 and Question 24 on small networks.

Environment: Python 3.12 (WSL Fedora 44), standard library only; exact rational arithmetic (crn.py).
Runtimes on 12 cores (see the logs in certificates/): n=2,m=2,r<=5: 5 s; n=3,m=2,r<=4: 1 min;
n=4,m=1,r<=12: 1 min; n=4,m=2,r<=3: 2 min; n=3,m=3,r<=3: 5 min.

Usage: python enumerate_small.py N_SPECIES MAX_MOLECULARITY MAX_REACTIONS [NPROC] [JSON_OUT]

Networks: sets of r distinct reactions y -> y' (y != y', both nonzero complexes of molecularity at most m; a
conservative network has no zero complex), every species occurring, one representative per orbit of the
species permutations.  Hypotheses of Conjecture 22 (H): conservative, consistent (= a positive steady state for
some rate constants), rows of Gamma pairwise distinct (no two species with exactly the same rate).

For every network satisfying H the script evaluates, on both sides of the duality, for the three precluster
variants (onestep = Definition 19 of the paper, thesis = Algorithm 4.4 of the thesis, closure) and the MIPS
readings (ray = Definition 17/Remark 18, semiflow = supports of extreme conservation laws, set = 'set of types'
reading, only when the conservation rank is <= 2):
    side A (ground set = species of N):   P = preclusters of N*, M = MIPS of N,  S = minimal siphons of N,
                                          C = internal cycles of N* (= supports of extreme conservation laws)
    side B (ground set = reactions of N): P = preclusters of N,  M = MIPS of N*, S = minimal siphons of N*,
                                          C = internal cycles of N
the statements P=M, P<=M, M<=P, |P|=|M| ('<=' = refines), C<=P (Remark 23), P=S, S<=P, P<=S,
'every block of P is a siphon', 'every block of P is a union of minimal siphons', |P|=|S|, and also
'every MIPS is a siphon' and S<=M.  It checks the exact criteria of Theorem 3 of the manuscript against the
direct computation (PASS/FAIL) and exits with status 1 on any disagreement.
"""
import sys
import json
import time
from itertools import combinations, permutations, product
from multiprocessing import Pool

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from crn import Network, refines, mips_set_q2  # noqa: E402

VARIANTS = ["onestep", "thesis", "closure"]
CLASSES = ["all", "strong", "dualH", "wr", "wr_def0"]
CLASS_DOC = {
    "all": "H: conservative, consistent, rows of Gamma distinct",
    "strong": "H + rows and columns of Gamma pairwise linearly independent",
    "dualH": "H + columns of Gamma distinct (so N* also satisfies H)",
    "wr": "H + weakly reversible",
    "wr_def0": "H + weakly reversible + deficiency zero",
}


def complexes(n, m):
    return [c for c in product(range(m + 1), repeat=n) if 0 < sum(c) <= m]


def canonical(n, rxns, perms):
    best = None
    for perm in perms:
        key = tuple(sorted((tuple(a[p] for p in perm), tuple(b[p] for p in perm)) for a, b in rxns))
        if best is None or key < best:
            best = key
    return best


def union_of(blocks_family, X):
    """X is a union of members of blocks_family."""
    U = set()
    for S in blocks_family:
        if S <= X:
            U |= S
    return U == set(X)


def evaluate(N):
    res = {}
    crit_fail = []
    D = N.dual()
    sides = {
        "A": dict(net=D, base=N, M={"ray": N.mips_ray(), "semiflow": N.mips_semiflow(), "set": mips_set_q2(N)},
                  S=N.minimal_siphons(), C=D.internal_cycles()),
        "B": dict(net=N, base=D, M={"ray": D.mips_ray(), "semiflow": D.mips_semiflow(), "set": mips_set_q2(D)},
                  S=D.minimal_siphons(), C=N.internal_cycles()),
    }
    extra = {}
    for side, d in sides.items():
        net, base, S, C = d["net"], d["base"], d["S"], d["C"]
        parts = {}
        for v in VARIANTS:
            P, edges = net.preclusters(v, return_graph=True)
            parts[v] = P
            for reading, M in d["M"].items():
                if M is None:
                    continue
                tag = f"{side}|{v}|{reading}"
                res[tag + "|P=M"] = set(P) == set(M)
                res[tag + "|P<=M"] = refines(P, M)
                res[tag + "|M<=P"] = refines(M, P)
                res[tag + "|#P=#M"] = len(P) == len(M)
            tag = f"{side}|{v}"
            res[tag + "|C<=P"] = refines(C, P)
            res[tag + "|P=S"] = set(P) == set(S)
            res[tag + "|S<=P"] = refines(S, P)
            res[tag + "|P<=S"] = refines(P, S)
            res[tag + "|P_blocks_siphons"] = all(base.is_siphon(b) for b in P)
            res[tag + "|P_blocks_unions_of_S"] = all(union_of(S, b) for b in P)
            res[tag + "|#P=#S"] = len(P) == len(S)
            # ---- criteria of Theorem 3 (ray reading of MIPS) ----
            M = d["M"]["ray"]
            blk = {i: b for b in M for i in b}
            dp = net.doubling_pairs()
            all_dp_coupled = all(blk[i] == blk[j] for i, j in dp)
            in_dp = {i for e in dp for i in e}
            nontriv_have_dp = all(any(i in in_dp for i in b) for b in M if len(b) > 1)
            nonisolated = {i for e in edges for i in e}
            touched = all(i in nonisolated for b in M if len(b) > 1 for i in b)
            if refines(P, M) != all_dp_coupled:
                crit_fail.append((side, v, "P<=M"))
            if (set(P) == set(M)) != (all_dp_coupled and nontriv_have_dp):
                crit_fail.append((side, v, "P=M"))
            if refines(M, P) != touched:
                crit_fail.append((side, v, "M<=P"))
        extra[side + "_variants_differ"] = len({tuple(sorted(map(tuple, map(sorted, parts[v])))) for v in VARIANTS}) > 1
        ms = d["M"]["set"]
        extra[side + "_set_reading_computed"] = ms is not None
        extra[side + "_set_differs_from_ray"] = ms is not None and set(ms) != set(d["M"]["ray"])
        res[f"{side}|MIPS_are_siphons"] = all(base.is_siphon(b) for b in d["M"]["ray"])
        res[f"{side}|S<=MIPS"] = refines(S, d["M"]["ray"])
    res["same-network|#P(N)=#M(N)"] = len(N.preclusters("thesis")) == len(N.mips_ray())
    return res, crit_fail, extra


def worker(args):
    n, m, rmax, first = args
    C = complexes(n, m)
    R = sorted((a, b) for a in C for b in C if a != b)
    perms = list(permutations(range(n)))
    out = []
    rest = R[first + 1:]
    for k in range(0, rmax):
        for tail in combinations(rest, k):
            rx = (R[first],) + tail
            if not all(any(a[i] or b[i] for a, b in rx) for i in range(n)):
                continue
            ok = True
            for i in range(n):
                cons = any(b[i] < a[i] for a, b in rx)
                prod_ = any(b[i] > a[i] for a, b in rx)
                if cons != prod_:
                    ok = False
                    break
            if not ok:
                continue
            if canonical(n, rx, perms) != rx:
                continue
            N = Network(n, rx)
            if not N.is_consistent() or not N.is_conservative():
                continue
            if not N.rows_distinct():
                out.append(("excluded_rows", rx, None))
                continue
            wr = N.weakly_reversible()
            flags = {
                "all": True,
                "strong": N.rows_pairwise_independent() and N.columns_pairwise_independent(),
                "dualH": N.columns_distinct(),
                "wr": wr,
                "wr_def0": wr and N.deficiency() == 0,
            }
            res, crit_fail, extra = evaluate(N)
            out.append(("H", rx, (flags, res, crit_fail, extra)))
    return out


def main():
    n, m, rmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    json_out = sys.argv[5] if len(sys.argv) > 5 else None
    t0 = time.time()
    C = complexes(n, m)
    nR = len(C) * (len(C) - 1)
    print(f"# enumerate_small.py n={n} species, molecularity <= {m}, r <= {rmax} reactions; "
          f"{len(C)} nonzero complexes, {nR} possible reactions")
    stats, first_fail = {}, {}
    count = {c: 0 for c in CLASSES}
    excluded_rows = []
    crit_failures = []
    extra_counts = {}
    examples_variants_differ = []
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, [(n, m, rmax, f) for f in range(nR)]):
            for kind, rx, data in chunk:
                if kind == "excluded_rows":
                    excluded_rows.append(rx)
                    continue
                flags, res, crit_fail, extra = data
                cls = [c for c in CLASSES if flags[c]]
                for c in cls:
                    count[c] += 1
                for k, v in extra.items():
                    extra_counts[k] = extra_counts.get(k, 0) + (1 if v else 0)
                if extra["A_variants_differ"] or extra["B_variants_differ"]:
                    examples_variants_differ.append(rx)
                if crit_fail:
                    crit_failures.append((rx, crit_fail))
                for st, val in res.items():
                    for c in cls:
                        h, f = stats.get((st, c), (0, 0))
                        stats[(st, c)] = (h + val, f + (not val))
                        if not val:
                            cand = (len(rx), rx)
                            if (st, c) not in first_fail or cand < first_fail[(st, c)]:
                                first_fail[(st, c)] = cand
    dt = time.time() - t0
    print("# classes:")
    for c in CLASSES:
        print(f"#   {c:8s} {count[c]:7d}  ({CLASS_DOC[c]})")
    print(f"# conservative + consistent networks excluded only by equal rows of Gamma: {len(excluded_rows)}")
    for k in sorted(extra_counts):
        print(f"# {k}: {extra_counts[k]}")
    print("#")
    print("# statement (side|precluster variant|MIPS reading|relation): holds/total per class")
    sts = sorted({st for st, _ in stats})
    # rows of the 'onestep' and 'closure' variants are printed only if they differ from the 'thesis' rows
    same_as_thesis = True
    for st in sts:
        for v in ("onestep", "closure"):
            if f"|{v}|" in st:
                twin = st.replace(f"|{v}|", "|thesis|")
                if any(stats.get((st, c)) != stats.get((twin, c)) for c in CLASSES):
                    same_as_thesis = False
    print(f"# the three precluster variants give identical statistics: {same_as_thesis}")
    if same_as_thesis:
        sts = [st for st in sts if "|onestep|" not in st and "|closure|" not in st]
    header = f"{'statement':40s}" + "".join(f"{c:>16s}" for c in CLASSES)
    print(header)
    for st in sts:
        row = f"{st:40s}"
        for c in CLASSES:
            h, f = stats.get((st, c), (0, 0))
            row += f"{f'{h}/{h + f}':>16s}"
        print(row)
    print("#")
    print("# smallest failing network per statement (class 'all'; thesis variant and ray reading shown)")
    for st in sts:
        if ("thesis" in st and "ray" in st) or ("thesis" in st and st.count("|") == 2) or st.count("|") <= 1:
            if (st, "all") in first_fail:
                print(f"{st:40s} {Network(n, first_fail[(st, 'all')][1]).describe()}")
    ok = not crit_failures
    print(f"{'PASS' if ok else 'FAIL'}: criteria of Theorem 3 agree with the direct computation on all "
          f"{count['all']} networks (both sides, three precluster variants)")
    for rx, cf in crit_failures[:20]:
        print("  criterion failure:", Network(n, rx).describe(), cf)
    if examples_variants_differ:
        print("# example where the precluster variants differ:", Network(n, examples_variants_differ[0]).describe())
    print(f"# runtime {dt:.1f} s")
    if json_out:
        with open(json_out, "w") as fh:
            json.dump({
                "n": n, "m": m, "rmax": rmax, "counts": count,
                "stats": {f"{st}||{c}": stats[(st, c)] for st, c in stats},
                "smallest_failing": {f"{st}||{c}": Network(n, first_fail[(st, c)][1]).describe()
                                     for st, c in first_fail},
                "variants_differ_examples": [Network(n, rx).describe() for rx in examples_variants_differ[:50]],
            }, fh, indent=0)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
