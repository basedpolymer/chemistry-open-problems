"""Block-partition lower bound for the minimum number of matches, by exhaustive exact enumeration.

Lemma (proved in result.md).  Call a nonempty set K of streams a *block* if it is heat balanced
(sum of hot heat in K = sum of cold heat in K) and cascade feasible (for every interval t, the hot heat of K in
intervals <= t is >= the cold demand of K in intervals <= t).  In every feasible solution, the vertex set of every
connected component of the match graph is a block; hence #matches >= (n+m) - c_max, where c_max is the largest
number of blocks in a partition of all streams into blocks.

This script computes, by meet-in-the-middle enumeration of ALL zero-sum subsets (exact integers: heat x lcm of
the decimal denominators), the set F of all blocks, and decides whether a partition into 2 and into 3 blocks
exists.  It prints the resulting lower bound and PASS/FAIL against the claimed value.

usage: python3 block_bound.py <instance.dat> <claimed_LB> [--list2]
  --list2 : also print every 2-block partition (needed for the UB/LB analysis of instances with LB = n+m-2)
Environment: Python 3.13 + numpy.  Runtime: < 1 min for N <= 27 streams, a few minutes for N = 33-38.
Memory: < 1.5 GB for N <= 38.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hen_instance import parse_dat


def half_tables(vals, profs):
    """all subsets of a half: sums (int64 / python int) and cumulative profiles (int64 matrix)."""
    h = len(vals)
    size = 1 << h
    sums = np.zeros(size, dtype=object)
    cum = np.zeros((size, profs.shape[1]), dtype=np.int64)
    for idx in range(h):
        step = 1 << idx
        sums[step:2 * step] = sums[:step] + vals[idx]
        cum[step:2 * step] = cum[:step] + profs[idx]
    return sums, cum


def blocks(inst):
    n, m, k = inst.n, inst.m, inst.k
    L = inst.denominators_lcm()
    vals = [int(inst.hot_total(i) * L) for i in range(n)] + [-int(inst.cold_total(j) * L) for j in range(m)]
    rows = [[int(x * L) for x in inst.QH[i]] for i in range(n)] + [[-int(x * L) for x in inst.QC[j]] for j in range(m)]
    for v, r in zip(vals, rows):
        assert sum(r) == v and v != 0
    # exactness of the int64 profile arithmetic: all partial sums bounded by sum |values| < 2^62
    assert sum(abs(v) for v in vals) < 2 ** 62
    profs = np.cumsum(np.array(rows, dtype=np.int64), axis=1)  # cumulative profile of each stream
    N = n + m
    h = N // 2
    sA, cA = half_tables(vals[:h], profs[:h])
    sB, cB = half_tables(vals[h:], profs[h:])
    # hash join on sums (exact python ints)
    idxB = {}
    for b, s in enumerate(sB):
        idxB.setdefault(int(s), []).append(b)
    pa, pb = [], []
    for a, s in enumerate(sA):
        lst = idxB.get(-int(s))
        if lst:
            pa.extend([a] * len(lst))
            pb.extend(lst)
    pa = np.array(pa, dtype=np.int64)
    pb = np.array(pb, dtype=np.int64)
    nzero = len(pa) - 1  # minus the empty set
    feas = np.zeros(len(pa), dtype=bool)
    CH = 200000
    for s in range(0, len(pa), CH):
        tot = cA[pa[s:s + CH]] + cB[pb[s:s + CH]]
        feas[s:s + CH] = (tot >= 0).all(axis=1)
    masks = set()
    for a, b, f in zip(pa.tolist(), pb.tolist(), feas.tolist()):
        if f:
            msk = a | (b << h)
            if msk:
                masks.add(msk)
    return N, nzero, masks


def analyse(inst, list2=False):
    t0 = time.time()
    N, nzero, F = blocks(inst)
    allm = (1 << N) - 1
    assert allm in F, "the whole instance must be a block (balanced, cascade feasible)"
    # 2-partitions {A, all\A} with A containing stream 0
    two = sorted(A for A in F if A != allm and (A & 1) and (allm & ~A) in F)
    # 3-partitions: A contains stream 0; X = all\A in F; B contains lowest(X), B subset X, X\B in F
    bylow = {}
    for S in F:
        low = (S & -S).bit_length() - 1
        bylow.setdefault(low, []).append(S)
    three = None
    for A in F:
        if A == allm or not (A & 1):
            continue
        X = allm & ~A
        if X not in F:
            continue
        low = (X & -X).bit_length() - 1
        for B in bylow.get(low, []):
            if B & A == 0 and B != X and (X & ~B) in F:
                three = (A, B, X & ~B)
                break
        if three:
            break
    cmax_ub = 2 if (two and not three) else (1 if not two else None)
    return dict(N=N, nzero=nzero, nblocks=len(F), two=two, three=three, cmax_ub=cmax_ub, time=time.time() - t0)


def name_of(inst, v):
    return f"H{v}" if v < inst.n else f"C{v - inst.n}"


if __name__ == "__main__":
    inst = parse_dat(sys.argv[1])
    claimed = int(sys.argv[2])
    res = analyse(inst, "--list2" in sys.argv)
    N = res["N"]
    print(f"instance {inst.name}: N = n+m = {N}; zero-sum nonempty subsets: {res['nzero']}; blocks: {res['nblocks']}; "
          f"2-block partitions: {len(res['two'])}; 3-block partition exists: {res['three'] is not None} ({res['time']:.1f}s)")
    if res["three"] is not None:
        print("  a 3-block partition exists -> this script only certifies LB = N-3 or less; example:",
              [[name_of(inst, v) for v in range(N) if S >> v & 1] for S in res["three"]])
        lb = None
    else:
        lb = N - (2 if res["two"] else 1)
    if "--list2" in sys.argv:
        for A in res["two"]:
            B = ((1 << N) - 1) & ~A
            print("  2-partition:", [name_of(inst, v) for v in range(N) if A >> v & 1], "|",
                  [name_of(inst, v) for v in range(N) if B >> v & 1])
    ok = lb is not None and lb >= claimed
    print(f"block-partition lower bound: {lb}  (claimed {claimed}) -> {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)
