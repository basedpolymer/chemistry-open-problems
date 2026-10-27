"""V1 checks of agent E's generalised readings of "refine" (computation/notes.md section 3).

Verifier V1 (Claude Opus 5.5, effort xhigh), 27 October 2026. Exact, standard library only, independent of ecore.py.
Checks: Prop E.1, E.2 (on random WR networks and the small exhaustive family), Prop E.3 (Galois connection) by brute
force on all partitions and union-closed families of small ground sets, and E's smallest counterexamples for the
refuted readings E-a3..a8, E-b1..b4, E-c1..c4, E-d1..d4 (re-checked exactly).
"""
import sys, time, random, platform
from itertools import combinations, product
sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from v1lib import *
from check_P import family
from check_Q2 import random_wr

T0 = time.time()
print("environment: Python", platform.python_version(), platform.system())
FAILS, RESULTS = [], []


def claim(name, ok, detail=""):
    RESULTS.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name + ((" -- " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


# ------------------------------------------------------------------ Prop E.3 (Galois connection), brute force

def set_partitions(xs):
    if not xs:
        yield []
        return
    first, rest = xs[0], xs[1:]
    for p in set_partitions(rest):
        for k in range(len(p)):
            yield p[:k] + [[first] + p[k]] + p[k + 1:]
        yield [[first]] + p


def galois():
    ok_all, ok_iso, cnt = True, True, 0
    for nX in (1, 2, 3):
        X = list(range(nX))
        allsets = [frozenset(c) for k in range(nX + 1) for c in combinations(X, k)]
        nonempty = [s for s in allsets if s]
        # all union-closed families containing the empty set
        fams = []
        for mask in range(1 << len(nonempty)):
            F = {frozenset()} | {nonempty[i] for i in range(len(nonempty)) if mask >> i & 1}
            if all(a | b in F for a in F for b in F):
                fams.append(F)
        parts = list(set_partitions(X))
        for P in parts:
            LP = {frozenset().union(*[frozenset(P[i]) for i in range(len(P)) if m >> i & 1]) for m in range(1 << len(P))}
            cP = lambda U: frozenset().union(*[frozenset(b) for b in P if set(b) & U]) if U else frozenset()
            for F in fams:
                iF = lambda V: frozenset().union(*[U for U in F if U <= V])
                for U in F:
                    for V in LP:
                        cnt += 1
                        a, b, c = cP(U) <= V, U <= V, U <= iF(V)
                        if not (a == b == c):
                            ok_all = False
            for Q in parts:
                LQ = {frozenset().union(*[frozenset(Q[i]) for i in range(len(Q)) if m >> i & 1]) for m in range(1 << len(Q))}
                iQ = lambda V: frozenset().union(*[U for U in LQ if U <= V])
                inverse = all(iQ(cP(U)) == U for U in LQ) and all(cP(iQ(V)) == V for V in LP)
                if inverse != (canon(P) == canon(Q)):
                    ok_iso = False
    claim("Prop E.3: c_P(U) <= V <=> U <= V <=> U <= i_F(V) for all partitions P and union-closed F (|X| <= 3; %d triples)"
          % cnt, ok_all)
    claim("Prop E.3: for F = L(Q), c_P and i_Q are inverse isomorphisms iff P = Q (|X| <= 3)", ok_iso)


# ------------------------------------------------------------------ Props E.1, E.2

def props_E12():
    rng = random.Random(99)
    bad1, bad2, n1, n2 = [], [], 0, 0
    for _ in range(3000):
        N = random_wr(rng, rng.randint(2, 4), 2)
        if N is None or not is_WR(N) or any(N.src(j) == N.prod(j) for j in range(N.r)):
            continue
        n1 += 1
        twoout = len({N.src(j) for j in range(N.r)}) < N.r
        P, M = preclusters_def19(N), full_coupling(N)
        if refines(P, M) == twoout:
            bad1.append(N)
        if not twoout and P != [[j] for j in range(N.r)]:
            bad1.append(N)
        if is_H1(N) and is_H3(N):
            n2 += 1
            D = N.dual()
            PD, MN = preclusters_def19(D), ray_mips(N)
            disc = PD == [[i] for i in range(N.n)]
            two = False
            for i, j in combinations(range(N.n), 2):
                rest = [N.G[k] for k in range(N.n) if k not in (i, j)]
                if srank(N) - (rank(rest) if rest else 0) >= 1:
                    two = True
            if not disc or not refines(PD, MN) or refines(MN, PD) != (MN == [[i] for i in range(N.n)]) or \
                    (MN == [[i] for i in range(N.n)]) == two:
                bad2.append(N)
    claim("Prop E.1 on %d random WR networks (P(N) <= M(N*) iff no complex has two outgoing reactions)" % n1,
          not bad1, repr(bad1[:1]))
    claim("Prop E.2 on %d random WR networks with (H1),(H3)" % n2, not bad2, repr(bad2[:1]))


# ------------------------------------------------------------------ smallest counterexamples of E's readings

def readings(N):
    """Evaluate E's readings on the side network N (which should satisfy (H1),(H2)); P(dual), M = ray MIPS."""
    D = N.dual()
    PD, M = preclusters_def19(D), ray_mips(N)
    S = [set(x) for x in min_siphons(N)]
    def union_of(block, fam):
        sub = [f for f in fam if f <= set(block)]
        return bool(sub) and set().union(*sub) == set(block)
    r = {}
    r["b1 every MIPS is a union of minimal siphons"] = all(union_of(m, S) for m in M)
    r["b2 every minimal siphon is a union of MIPS"] = all(union_of(s, [set(m) for m in M]) for s in S)
    r["b3 every precluster of the dual is a union of minimal siphons"] = all(union_of(b, S) for b in PD)
    r["b4 every minimal siphon is a union of preclusters of the dual"] = all(union_of(s, [set(b) for b in PD]) for s in S)
    ns = lambda P: sorted(sorted(b) for b in P if len(b) >= 2)
    r["d1 non-singleton P(dual) = non-singleton M"] = ns(PD) == ns(M)
    r["d2 non-singleton P(dual) blocks refine M"] = refines([b for b in PD if len(b) >= 2], M)
    r["d3 non-singleton M blocks lie in blocks of P(dual)"] = refines([m for m in M if len(m) >= 2], PD)
    r["d4 non-singleton min siphons = non-singleton P(dual) blocks"] = sorted(sorted(s) for s in S if len(s) >= 2) == ns(PD)
    r["c4 join(P(dual),M) = join(P(dual),M,min siphons)"] = join(PD, M, N.n) == components(
        N.n, [(a, b) for P_ in (PD, M) for blk in P_ for a, b in zip(blk, blk[1:])] +
        [(a, b) for s in S for a, b in zip(sorted(s), sorted(s)[1:])])
    return r


def smallest():
    cases = [
        ("E-b1/b2 (side D)", net_from_str("X2 -> X1, X1 + X2 -> 2X2", species_order=["X1", "X2"]).dual(),
         ["b1 every MIPS is a union of minimal siphons", "b2 every minimal siphon is a union of MIPS"]),
        ("E-b3, E-d1, E-d3, E-d4 (side D of R)", net_from_str("X2 -> X1, X1 -> X2", species_order=["X1", "X2"]).dual(),
         ["b3 every precluster of the dual is a union of minimal siphons", "d1 non-singleton P(dual) = non-singleton M",
          "d3 non-singleton M blocks lie in blocks of P(dual)", "d4 non-singleton min siphons = non-singleton P(dual) blocks"]),
        ("E-b4 (side D)", net_from_str("X1 + X2 -> 2X2, X1 + X2 -> 2X1", species_order=["X1", "X2"]).dual(),
         ["b4 every minimal siphon is a union of preclusters of the dual"]),
        ("E-d2 (side N)", net_from_str("X3 -> X2, X1 + X2 -> X1 + X3", species_order=["X1", "X2", "X3"]),
         ["d2 non-singleton P(dual) blocks refine M"]),
        ("E-c4 (side N)", net_from_str("X3 -> X1 + X2, X1 + 2X2 -> X2 + X3", species_order=["X1", "X2", "X3"]),
         ["c4 join(P(dual),M) = join(P(dual),M,min siphons)"]),
    ]
    for label, N, keys in cases:
        r = readings(N)
        hyp = is_H1(N) and is_H2(N)
        claim("%s: %r satisfies (H1),(H2) [(H3): %s] and refutes %s" % (label, N, is_H3(N),
                                                         "; ".join(k.split()[0] for k in keys)),
              hyp and all(not r[k] for k in keys), str({k.split()[0]: r[k] for k in keys}))
    # E-a3..a8 on R (and E-a5 / E-a7 examples)
    R = net_from_str("X2 -> X1, X1 -> X2", species_order=["X1", "X2"])
    P, Ms, PD, M = preclusters_def19(R), full_coupling(R), preclusters_def19(R.dual()), ray_mips(R)
    claim("E-a3, a4, a6, a8 refuted by R: P(R) != M(R*), P(R*) != M(R), M(R) not<= P(R*), M(R*) not<= P(R)",
          P != Ms and PD != M and not refines(M, PD) and not refines(Ms, P))
    N = net_from_str("X3 -> X2, X1 + X2 -> X1 + X3", species_order=["X1", "X2", "X3"])
    claim("E-a5 refuted: X3 -> X2, X1+X2 -> X1+X3 has (H1)-(H3), P(N*) not<= M(N) (dual has the loop k2 -> k2)",
          is_H1(N) and is_H2(N) and is_H3(N) and not refines(preclusters_def19(N.dual()), ray_mips(N)) and
          any(N.dual().src(j) == N.dual().prod(j) for j in range(N.n)))
    N = net_from_str("X2 -> X1, X1 + X2 -> 2X2, X1 + X2 -> 2X1", species_order=["X1", "X2"])
    claim("E-a7 refuted: X2 -> X1, X1+X2 -> 2X2, X1+X2 -> 2X1 has (H1)-(H3), P(N) not<= M(N*)",
          is_H1(N) and is_H2(N) and is_H3(N) and not refines(preclusters_def19(N), full_coupling(N)))
    # E-c1, E-c1', E-c2, E-c3
    N = net_from_str("X2 -> X1 + X2, X1 + X2 -> 2X1 + X2, 2X1 + X2 -> X2", species_order=["X1", "X2"])
    lc = linkage_classes(N)[0]
    claim("E-c1 refuted: WR %r, P(N) v M(N*) != linkage classes" % N,
          is_WR(N) and join(preclusters_def19(N), full_coupling(N), N.r) != lc)
    N = net_from_str("0 -> X2, 0 -> X1, X2 -> 0, X2 -> X1, X1 -> 0", species_order=["X1", "X2"])
    claim("E-c1' refuted: WR deficiency 0 %r, P(N) v M(N*) != linkage classes" % N,
          is_WR(N) and deficiency(N) == 0 and join(preclusters_def19(N), full_coupling(N), N.r) != linkage_classes(N)[0])
    N = net_from_str("X2 -> X1, X1 + X2 -> 2X2", species_order=["X1", "X2"])
    claim("E-c2 refuted: %r (H), P(N) v M(N*) does not refine the linkage classes" % N,
          is_H1(N) and is_H2(N) and is_H3(N) and not refines(join(preclusters_def19(N), full_coupling(N), N.r),
                                                           linkage_classes(N)[0]))
    N = net_from_str("X3 -> X2, X1 + 2X2 -> X1 + 2X3", species_order=["X1", "X2", "X3"])
    claim("E-c3 refuted: %r connected with (H), P(N*) v M(N) not one block" % N,
          is_connected(N) and is_H1(N) and is_H2(N) and is_H3(N) and
          len(join(preclusters_def19(N.dual()), ray_mips(N), N.n)) > 1)
    # Y5 illustration
    Y5 = net_from_str("A -> B, B -> C, B -> D, B -> E, C -> A, D -> A, E -> A")
    claim("E1 illustration Y5: P(Y5) = {A->B},{six others}; M(Y5*) = {A->B},{B->C,C->A},{B->D,D->A},{B->E,E->A}; "
          "join has two blocks; WR, deficiency 0, (H1)-(H3)",
          preclusters_def19(Y5) == [[0], [1, 2, 3, 4, 5, 6]] and full_coupling(Y5) == [[0], [1, 4], [2, 5], [3, 6]]
          and len(join(preclusters_def19(Y5), full_coupling(Y5), Y5.r)) == 2 and is_WR(Y5) and deficiency(Y5) == 0
          and is_H1(Y5) and is_H2(Y5) and is_H3(Y5))


if __name__ == "__main__":
    galois()
    props_E12()
    smallest()
    print("runtime %.1f s; %d claims, %d FAIL" % (time.time() - T0, len(RESULTS), len(FAILS)))
    if FAILS:
        print("FAILED:", FAILS)
        sys.exit(1)
