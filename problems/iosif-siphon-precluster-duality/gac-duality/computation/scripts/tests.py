"""Tests of candidate relations (agents P and G, and E's generalised readings of "refine") on the enumerated records.

Environment: Python 3.12 (WSL Fedora 44), standard library only, exact rational arithmetic (crn.py via ecore.py).
Usage: python tests.py RECORD_FILES... [--only PREFIX[,PREFIX...]] [--json OUT] [--nproc K]
Candidates whose id starts with 'X-' (endotacticity, drain species; slower) run only when selected with --only.

Records are produced by enum_wr.py (weakly reversible families), enum_all.py (all networks with (H1),(H2)) and
sample_wr.py (seeded samples).  Duplicates across files (same species count and canonical reaction set) are
counted once.  Every candidate is a function rec -> None (not applicable) / True / False.  For each candidate the
script prints the number of applicable networks, how many satisfy it, the same restricted to the classes
  H     : (H1)-(H3) for N
  strong: H and rows and columns of Gamma pairwise linearly independent
  WR    : H and N weakly reversible;  WR0: WR and deficiency 0
and the smallest counterexample (by n + r, then r).  Status of a candidate: 'proved' (claimed proved by its author:
a counterexample is a FAIL of the script), 'test' (conjecture, question or reading: PASS = holds on every tested
network, REFUTED = counterexample found; neither outcome is a failure of the script).
Side conventions: a 'side' is (network, its dual).  Side N: ground sets species of N (P(N*), M(N), S(N)); side D:
ground sets reactions of N (P(N), M(N*), S(N*)).  Symmetric candidates are tested on both sides.
"""
import sys
import time
import gzip
import json
import pickle
import platform
from fractions import Fraction
from itertools import combinations

import ecore
from ecore import CNet, refines, rank

CANDS = []


def cand(cid, doc, status="test", sides=False):
    def deco(f):
        CANDS.append(dict(id=cid, doc=doc, status=status, f=f, sides=sides))
        return f
    return deco


# ----------------------------------------------------------------------------------------------
# side abstraction
# ----------------------------------------------------------------------------------------------
class Side:
    """One orientation: the network X (species ground set) together with its dual Xd."""
    def __init__(self, rec, which):
        n, rx = rec["n"], rec["rx"]
        if which == "N":
            self.n, self.rx, self.X, self.Xd = n, rx, rec["N"], rec["D"]
        else:
            Ye = [tuple(a[i] for a, b in rx) for i in range(n)]
            Yp = [tuple(b[i] for a, b in rx) for i in range(n)]
            self.n, self.rx, self.X, self.Xd = len(rx), tuple(zip(Ye, Yp)), rec["D"], rec["N"]
        self.r = len(self.rx)
        self.which = which
        self.rec = rec
        self._c = {}

    def net(self):
        if "net" not in self._c:
            self._c["net"] = CNet(self.n, self.rx)
        return self._c["net"]

    def ye_rows(self):
        return [tuple(a[i] for a, b in self.rx) for i in range(self.n)]

    def H(self):
        X = self.X
        return X["cons"] and X["consist"] and X["rows_distinct"]

    def twins(self):
        rows = self.ye_rows()
        return [(i, j) for i, j in combinations(range(self.n), 2) if rows[i] == rows[j] and any(rows[i])]

    def doubling(self):
        return [(i, j) for i, j in combinations(range(self.r), 2)
                if self.rx[i][0] == self.rx[j][0] and any(self.rx[i][0])]

    def src_sets(self):
        return [frozenset(i for i in range(self.n) if a[i] > 0) for a, b in self.rx]

    def ordinary(self):
        return all(x in (0, 1) for a, b in self.rx for x in a + b)

    def EC(self):
        ss = self.src_sets()
        return all(self.rx[j][0] == self.rx[k][0] for j, k in combinations(range(self.r), 2) if ss[j] & ss[k])

    def FC(self):  # agent P's definition: ordinary + equal conflict
        return self.ordinary() and self.EC()

    def balanced_joins(self):
        return all(len({x for x in a if x > 0}) <= 1 for a, b in self.rx)

    def JF(self):
        return all(sum(1 for x in a if x > 0) <= 1 for a, b in self.rx)

    def CF(self):
        ss = self.src_sets()
        return all(sum(1 for s in ss if i in s) <= 1 for i in range(self.n))

    def connected(self):
        n, r = self.n, self.r
        parent = list(range(n + r))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for j, (a, b) in enumerate(self.rx):
            for i in range(n):
                if a[i] > 0 or b[i] > 0:
                    parent[find(i)] = find(n + j)
        return len({find(x) for x in range(n + r)}) == 1

    def nclusters(self):
        CLr, CLs = ecore.petri_clusters(self.net())
        # clusters = components of the species-reaction graph of input arcs; count all n + r nodes
        n, r = self.n, self.r
        parent = list(range(n + r))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for j, (a, b) in enumerate(self.rx):
            for i in range(n):
                if a[i] > 0:
                    parent[find(i)] = find(n + j)
        return len({find(x) for x in range(n + r)})

    def ndoubling_classes(self):
        return len({a for a, b in self.rx if any(a)})

    def s(self):
        return self.rec["s"]

    def rank_eq(self):
        return self.s() == self.nclusters() - 1

    def WFFC(self):
        X = self.X
        return self.connected() and self.FC() and X["cons"] and X["consist"] and self.rank_eq()


def sides(rec):
    return [Side(rec, "N"), Side(rec, "D")]


def is_union_of(family, X):
    U = set()
    for S in family:
        if S <= X:
            U |= S
    return U == set(X)


def join(P, Q, m):
    parent = list(range(m))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for B in list(P) + list(Q):
        B = sorted(B)
        for x in B[1:]:
            parent[find(x)] = find(B[0])
    blocks = {}
    for x in range(m):
        blocks.setdefault(find(x), set()).add(x)
    return set(frozenset(b) for b in blocks.values())


def lc_partition(rx):
    """Partition of the reactions by linkage class."""
    cs = []
    for a, b in rx:
        for c in (a, b):
            if c not in cs:
                cs.append(c)
    idx = {c: i for i, c in enumerate(cs)}
    parent = list(range(len(cs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in rx:
        parent[find(idx[a])] = find(idx[b])
    blocks = {}
    for j, (a, b) in enumerate(rx):
        blocks.setdefault(find(idx[a]), set()).add(j)
    return set(frozenset(b) for b in blocks.values())


def ADLS(X):
    """Angeli-De Leenheer-Sontag criterion: every siphon contains a semiflow support (no critical minimal siphon)."""
    return len(X["CS"]) == 0


# ==============================================================================================
# E: generalised readings of "refine" (task 3 of BRIEF-E)
# ==============================================================================================
@cand("E-a1", "WR: P(N) refines M(N*) iff no complex is the source of two reactions (every linkage class is a "
      "directed cycle)", status="proved")
def _(rec):
    N = rec["N"]
    if not (N["WR"] and N["consist"] and rec["D"]["cons"]):
        return None
    branching = len(rec["N"]["doubling"]) > 0
    return refines(N["P"], rec["D"]["M"]) == (not branching)


@cand("E-a2", "WR+H: M(N) refines P(N*) iff M(N) is discrete iff im Gamma has no vector supported on exactly two "
      "species", status="proved")
def _(rec):
    N = rec["N"]
    if not (N["WR"] and Side(rec, "N").H()):
        return None
    M = N["M"]
    disc = all(len(b) == 1 for b in M)
    # im Gamma contains a vector supported on {i,j} iff rank of Gamma with rows i,j deleted drops by 1
    G = CNet(rec["n"], rec["rx"]).Gamma
    s = rec["s"]
    two = False
    for i, j in combinations(range(rec["n"]), 2):
        rows = [G[k] for k in range(rec["n"]) if k not in (i, j)]
        # vectors of im Gamma vanishing outside {i,j}: dim = s - rank(rows of the other species)
        if s - (rank(rows) if rows else 0) >= 1:
            two = True
            break
    return refines(M, rec["D"]["P"]) == disc and disc == (not two)


@cand("E-a3", "H: P(N) = M(N*) (reaction side, equality)")
def _(rec):
    if not Side(rec, "N").H():
        return None
    return set(rec["N"]["P"]) == set(rec["D"]["M"])


@cand("E-a4", "H: P(N*) = M(N) (species side, equality)")
def _(rec):
    if not Side(rec, "N").H():
        return None
    return set(rec["D"]["P"]) == set(rec["N"]["M"])


@cand("E-a5", "H: P(N*) refines M(N) (species side; exceptional class = some twin pair not conservation-coupled)")
def _(rec):
    if not Side(rec, "N").H():
        return None
    return refines(rec["D"]["P"], rec["N"]["M"])


@cand("E-a6", "H: M(N) refines P(N*) (species side)")
def _(rec):
    if not Side(rec, "N").H():
        return None
    return refines(rec["N"]["M"], rec["D"]["P"])


@cand("E-a7", "H: P(N) refines M(N*) (reaction side; exceptional class = some doubling pair not fully coupled)")
def _(rec):
    if not Side(rec, "N").H():
        return None
    return refines(rec["N"]["P"], rec["D"]["M"])


@cand("E-a8", "H: M(N*) refines P(N) (reaction side)")
def _(rec):
    if not Side(rec, "N").H():
        return None
    return refines(rec["D"]["M"], rec["N"]["P"])


@cand("E-b1", "(H1),(H2) for the side's network (both sides): every MIPS block is a union of minimal siphons", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    return all(is_union_of(sd.X["MS"], B) for B in sd.X["M"])


@cand("E-b2", "(H1),(H2) for the side's network (both sides): every minimal siphon is a union of MIPS blocks", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    return all(is_union_of(sd.X["M"], S) for S in sd.X["MS"])


@cand("E-b3", "(H1),(H2) for the side's network (both sides): every precluster of the dual is a union of minimal siphons", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    return all(is_union_of(sd.X["MS"], B) for B in sd.Xd["P"])


@cand("E-b4", "(H1),(H2) for the side's network (both sides): every minimal siphon is a union of preclusters of the dual", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    return all(is_union_of(sd.Xd["P"], S) for S in sd.X["MS"])


@cand("E-c1", "WR: P(N) v M(N*) (join of partitions = image lattice of the Galois connection) = linkage classes")
def _(rec):
    N = rec["N"]
    if not (N["WR"] and N["consist"] and rec["D"]["cons"]):
        return None
    return join(N["P"], rec["D"]["M"], rec["r"]) == lc_partition(rec["rx"])


@cand("E-c1'", "WR, deficiency 0: P(N) v M(N*) = linkage classes")
def _(rec):
    N = rec["N"]
    if not (N["WR"] and N["deficiency"] == 0 and N["consist"] and rec["D"]["cons"]):
        return None
    return join(N["P"], rec["D"]["M"], rec["r"]) == lc_partition(rec["rx"])


@cand("E-c2", "H: P(N) v M(N*) refines the linkage-class partition")
def _(rec):
    if not Side(rec, "N").H():
        return None
    return refines(join(rec["N"]["P"], rec["D"]["M"], rec["r"]), lc_partition(rec["rx"]))


@cand("E-c3", "H: P(N*) v M(N) is the one-block partition when N is connected (species side join)")
def _(rec):
    sd = Side(rec, "N")
    if not (sd.H() and sd.connected()):
        return None
    return len(join(rec["D"]["P"], rec["N"]["M"], rec["n"])) == 1


@cand("E-c4", "(H1),(H2) for the side's network (both sides): join(P(dual), M) = join(P(dual), M, minimal siphons) (siphons add nothing to the "
      "join)", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    J = join(sd.Xd["P"], sd.X["M"], sd.n)
    return J == join(J, [S for S in sd.X["MS"]], sd.n)


@cand("E-d1", "(H1),(H2) for the side's network (both sides): non-singleton blocks of P(dual) = non-singleton blocks of M", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    return {b for b in sd.Xd["P"] if len(b) > 1} == {b for b in sd.X["M"] if len(b) > 1}


@cand("E-d2", "(H1),(H2) for the side's network (both sides): non-singleton blocks of P(dual) refine M (singletons removed)", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    return refines([b for b in sd.Xd["P"] if len(b) > 1], sd.X["M"])


@cand("E-d3", "(H1),(H2) for the side's network (both sides): non-singleton blocks of M lie in blocks of P(dual)", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    return refines([b for b in sd.X["M"] if len(b) > 1], sd.Xd["P"])


@cand("E-d4", "(H1),(H2) for the side's network (both sides): non-singleton minimal siphons = non-singleton blocks of P(dual)", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    return {b for b in sd.X["MS"] if len(b) > 1} == {b for b in sd.Xd["P"] if len(b) > 1}


@cand("E-e1", "WR+H: N* has (H1)-(H3) as well (dual inherits the hypotheses)")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H()):
        return None
    return Side(rec, "D").H()


@cand("E-e2", "WR+H: ADLS criterion (no critical siphon) for N implies it for N*")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H()):
        return None
    return (not ADLS(rec["N"])) or ADLS(rec["D"])


@cand("E-e3", "WR+H: ADLS criterion for N* implies it for N")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H()):
        return None
    return (not ADLS(rec["D"])) or ADLS(rec["N"])


@cand("E-e4", "WR+H and N* WR: ADLS(N) iff ADLS(N*)")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H() and rec["D"]["WR"]):
        return None
    return ADLS(rec["N"]) == ADLS(rec["D"])


@cand("E-e4b", "N and N* both WR with (H1)-(H3): ADLS(N) iff ADLS(N*)")
def _(rec):
    if not (rec["N"]["WR"] and rec["D"]["WR"] and Side(rec, "N").H() and Side(rec, "D").H()):
        return None
    return ADLS(rec["N"]) == ADLS(rec["D"])


@cand("E-e5", "WR+H and N* WR: N* satisfies (H1)-(H3)")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H() and rec["D"]["WR"]):
        return None
    return Side(rec, "D").H()


@cand("E-e6", "WR+H, deficiency 0: N* has deficiency 0")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H() and rec["N"]["deficiency"] == 0):
        return None
    return rec["D"]["deficiency"] == 0


@cand("E-e7", "WR+H: the dual satisfies (H1)-(H3) and has no critical siphon whenever N has none (ADLS transfer "
      "under H*)")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H() and Side(rec, "D").H() and ADLS(rec["N"])):
        return None
    return ADLS(rec["D"])


# ==============================================================================================
# agent P (petri/candidates.md)
# ==============================================================================================
@cand("P-K1", "no hypothesis: a minimal siphon contains at most one species of each class of species with the same "
      "set of consuming reactions (both sides)", status="proved", sides=True)
def _(sd):
    ss = sd.src_sets()
    cons = [frozenset(j for j in range(sd.r) if i in ss[j]) for i in range(sd.n)]
    for S in sd.X["MS"]:
        seen = set()
        for i in S:
            if cons[i] in seen:
                return False
            seen.add(cons[i])
    return True


@cand("P-K2", "(H1) side N / (H2) side D: P(dual) refines S iff no twin pair and every element lies in a minimal "
      "siphon", status="proved", sides=True)
def _(sd):
    if not (sd.X["cons"] and sd.X["consist"]):
        return None
    lhs = refines(sd.Xd["P"], sd.X["MS"])
    cover = set().union(*sd.X["MS"]) if sd.X["MS"] else set()
    rhs = (not sd.twins()) and cover == set(range(sd.n))
    return lhs == rhs


@cand("P-K2n", "no hypothesis (P's latest version of K2): P(dual) refines S iff no twin pair and every element lies "
      "in a minimal siphon (preclusters computed by crn.py even when the dual is not consistent)", status="proved",
      sides=True)
def _(sd):
    lhs = refines(sd.Xd["P"], sd.X["MS"])
    cover = set().union(*sd.X["MS"]) if sd.X["MS"] else set()
    rows = sd.ye_rows()
    twins_any = any(rows[i] == rows[j] for i, j in combinations(range(sd.n), 2))  # crn's doubling pairs of the
    # dual include equal ZERO rows (species consumed by no reaction); P's twins exclude them
    rhs = (not twins_any) and cover == set(range(sd.n))
    return lhs == rhs


def _wffc_sides(rec):
    return [sd for sd in sides(rec) if sd.WFFC()]


@cand("P-K3a", "WF-FC: the dual is WF-FC (Hack)", status="proved", sides=True)
def _(sd):
    if not sd.WFFC():
        return None
    other = Side(sd.rec, "D" if sd.which == "N" else "N")
    return other.WFFC()


@cand("P-K3b", "WF-FC: extreme conservation rays = 0/1 indicators of minimal siphons", status="proved", sides=True)
def _(sd):
    if not sd.WFFC():
        return None
    rays = set(sd.net().conservation_rays())
    ind = set(tuple(1 if i in S else 0 for i in range(sd.n)) for S in sd.X["MS"])
    return rays == ind


@cand("P-K3c", "WF-FC: extreme flux rays = 0/1 indicators of minimal T-traps (minimal siphons of the dual)",
      status="proved", sides=True)
def _(sd):
    if not sd.WFFC():
        return None
    rays = set(sd.net().flux_rays())
    ind = set(tuple(1 if j in S else 0 for j in range(sd.r)) for S in sd.Xd["MS"])
    return rays == ind


@cand("P-K3d", "WF-FC: ray MIPS = classes of species lying in the same minimal siphons", status="proved", sides=True)
def _(sd):
    if not sd.WFFC():
        return None
    prof = {}
    for i in range(sd.n):
        prof.setdefault(frozenset(k for k, S in enumerate(sd.X["MS"]) if i in S), set()).add(i)
    return set(sd.X["M"]) == set(frozenset(b) for b in prof.values())


def profile_rule(n, twins, MS):
    Sig = [frozenset(k for k, S in enumerate(MS) if i in S) for i in range(n)]
    edges = list(twins)
    for i, j in twins:
        for s in range(n):
            if s in (i, j):
                continue
            if Sig[s] in (Sig[i], Sig[j], Sig[i] | Sig[j]):
                edges += [(i, s), (j, s)]
    return set(ecore.partition_from_components(n, edges))


@cand("P-K3e", "WF-FC: profile rule gives P(dual)", status="proved", sides=True)
def _(sd):
    if not sd.WFFC():
        return None
    return profile_rule(sd.n, sd.twins(), sd.X["MS"]) == set(sd.Xd["P"])


@cand("P-K3f", "WF-FC: N and N* have no critical siphon", status="proved", sides=True)
def _(sd):
    if not sd.WFFC():
        return None
    return ADLS(sd.X) and ADLS(sd.Xd)


@cand("P-K3g", "WF-FC: no twin pair conservation-coupled; P(dual)<=M iff JF; P(dual)=M fails if n>=2; "
      "P(dual)<=S iff JF", status="proved", sides=True)
def _(sd):
    if not sd.WFFC():
        return None
    blk = {i: b for b in sd.X["M"] for i in b}
    ok = all(blk[i] != blk[j] for i, j in sd.twins())
    ok &= refines(sd.Xd["P"], sd.X["M"]) == sd.JF()
    if sd.n >= 2:
        ok &= set(sd.Xd["P"]) != set(sd.X["M"])
    ok &= refines(sd.Xd["P"], sd.X["MS"]) == sd.JF()
    return ok


@cand("P-K3g'", "WF-FC: reaction side: no doubling pair fully coupled; P<=M(dual) iff CF; P=M(dual) fails if r>=2; "
      "P<=S(dual) iff CF", status="proved", sides=True)
def _(sd):
    if not sd.WFFC():
        return None
    blk = {i: b for b in sd.Xd["M"] for i in b}
    ok = all(blk[i] != blk[j] for i, j in sd.doubling())
    ok &= refines(sd.X["P"], sd.Xd["M"]) == sd.CF()
    if sd.r >= 2:
        ok &= set(sd.X["P"]) != set(sd.Xd["M"])
    ok &= refines(sd.X["P"], sd.Xd["MS"]) == sd.CF()
    return ok


@cand("P-K3-need", "FC+(H1)+(H2)+connected but rank equation violated: (b)-(e) all hold (expected to FAIL "
      "somewhere: the rank equation is needed)", sides=True)
def _(sd):
    if not (sd.FC() and sd.X["cons"] and sd.X["consist"] and sd.connected()) or sd.rank_eq():
        return None
    rays = set(sd.net().conservation_rays())
    ind = set(tuple(1 if i in S else 0 for i in range(sd.n)) for S in sd.X["MS"])
    b = rays == ind
    frays = set(sd.net().flux_rays())
    find = set(tuple(1 if j in S else 0 for j in range(sd.r)) for S in sd.Xd["MS"])
    c = frays == find
    prof = {}
    for i in range(sd.n):
        prof.setdefault(frozenset(k for k, S in enumerate(sd.X["MS"]) if i in S), set()).add(i)
    d = set(sd.X["M"]) == set(frozenset(x) for x in prof.values())
    e = profile_rule(sd.n, sd.twins(), sd.X["MS"]) == set(sd.Xd["P"])
    return b and c and d and e


@cand("P-K4i", "well-formed EC (EC, connected, (H1), (H2), s = #doubling classes - 1) with balanced joins: "
      "semiflow supports = minimal siphons", sides=True)
def _(sd):
    if not (sd.EC() and sd.connected() and sd.X["cons"] and sd.X["consist"]
            and sd.s() == sd.ndoubling_classes() - 1 and sd.balanced_joins()):
        return None
    return set(sd.X["SF"]) == set(sd.X["MS"])


@cand("P-K4ii", "well-formed EC without balanced joins: semiflow supports = minimal siphons", sides=True)
def _(sd):
    if not (sd.EC() and sd.connected() and sd.X["cons"] and sd.X["consist"]
            and sd.s() == sd.ndoubling_classes() - 1 and not sd.balanced_joins()):
        return None
    return set(sd.X["SF"]) == set(sd.X["MS"])


@cand("P-K5", "EC, connected, (H1), (H2): s >= #doubling classes - 1", sides=True)
def _(sd):
    if not (sd.EC() and sd.connected() and sd.X["cons"] and sd.X["consist"]):
        return None
    return sd.s() >= sd.ndoubling_classes() - 1


@cand("P-K6", "(H1),(H2): N has no critical siphon iff N* has no critical siphon")
def _(rec):
    if not (rec["N"]["cons"] and rec["N"]["consist"]):
        return None
    return ADLS(rec["N"]) == ADLS(rec["D"])


@cand("P-K7", "WR+EC: distinct complexes have disjoint supports, deficiency 0, minimal siphons = semiflow supports, "
      "no critical siphon; with (H3)+ordinary every complex is one species", status="proved")
def _K7(rec):
    sd = Side(rec, "N")
    if not (rec["N"]["WR"] and sd.EC()):
        return None
    cs = CNet(rec["n"], rec["rx"]).complexes()
    ok = all(not any(a[i] and b[i] for i in range(rec["n"])) for a, b in combinations(cs, 2))
    ok &= rec["N"]["deficiency"] == 0
    ok &= set(rec["N"]["MS"]) == set(rec["N"]["SF"])
    ok &= ADLS(rec["N"])
    if rec["N"]["rows_distinct"] and sd.ordinary():
        ok &= all(sum(c) == 1 for c in cs)
    return ok


@cand("P-K7r", "K7 restricted (P, 18:15): WR + EC + no zero complex => all clauses of K7", status="proved")
def _(rec):
    if any(sum(c) == 0 for a, b in rec["rx"] for c in (a, b)):
        return None
    return _K7(rec)


for _k, _doc, _f in [
        ("i", "every minimal siphon is a union of ray-MIPS blocks", lambda sd: all(is_union_of(sd.X["M"], S) for S in sd.X["MS"])),
        ("ii", "every ray-MIPS block lies in a minimal siphon", lambda sd: all(any(B <= S for S in sd.X["MS"]) for B in sd.X["M"])),
        ("iii", "minimal siphons = semiflow supports", lambda sd: set(sd.X["MS"]) == set(sd.X["SF"])),
        ("iv", "every semiflow support is a union of minimal siphons", lambda sd: all(is_union_of(sd.X["MS"], F) for F in sd.X["SF"]))]:
    def _mk(f):
        def g(sd):
            if not sd.H():
                return None
            return f(sd)
        return g
    cand(f"P-K9{_k}", f"(H1)-(H3) on the side's network: {_doc}", sides=True)(_mk(_f))


@cand("P-H", "(H1) for the side's network (so its dual is consistent): every block of P(dual) with >= 2 elements "
      "is a union of MIPS blocks", status="proved", sides=True)
def _(sd):
    if not sd.X["cons"]:
        return None
    return all(is_union_of(sd.X["M"], B) for B in sd.Xd["P"] if len(B) > 1)


def _wfec(sd):
    return (sd.EC() and sd.connected() and sd.X["cons"] and sd.X["consist"]
            and all(any(a) for a, b in sd.rx) and sd.s() == sd.ndoubling_classes() - 1)


@cand("P-C", "well-formed EC (EC, connected, (H1), (H2), s = m-1): no critical siphon and minimal siphons = supports "
      "of extreme conservation laws", status="proved", sides=True)
def _(sd):
    if not _wfec(sd):
        return None
    return ADLS(sd.X) and set(sd.X["MS"]) == set(sd.X["SF"])


@cand("P-C'", "well-formed EC: the dual has no critical siphon (every T-trap contains an internal cycle)",
      status="proved", sides=True)
def _(sd):
    if not _wfec(sd):
        return None
    return ADLS(sd.Xd)


@cand("P-K14a", "well-formed EC: minimal T-traps (minimal siphons of the dual) = internal cycles", status="proved",
      sides=True)
def _(sd):
    if not _wfec(sd):
        return None
    return set(sd.Xd["MS"]) == set(sd.X["IC"])


def _pure(sd):
    return all(not any(a[i] and b[i] for i in range(sd.n)) for a, b in sd.rx)


@cand("P-K14b", "well-formed EC and pure: minimal siphons = semiflow supports = internal cycles of the dual",
      status="proved", sides=True)
def _(sd):
    if not (_wfec(sd) and _pure(sd)):
        return None
    return set(sd.X["MS"]) == set(sd.X["SF"]) == set(sd.Xd["IC"])


@cand("P-K14c", "well-formed EC and NOT pure: minimal siphons = semiflow supports (is purity needed?)", sides=True)
def _(sd):
    if not (_wfec(sd) and not _pure(sd)):
        return None
    return set(sd.X["MS"]) == set(sd.X["SF"])


@cand("P-G", "(H1)-(H3): minimal siphons = semiflow supports implies every minimal siphon is a union of MIPS and "
      "every MIPS lies in a minimal siphon", status="proved", sides=True)
def _(sd):
    if not sd.H():
        return None
    if set(sd.X["MS"]) != set(sd.X["SF"]):
        return True
    return (all(is_union_of(sd.X["M"], S) for S in sd.X["MS"])
            and all(any(B <= S for S in sd.X["MS"]) for B in sd.X["M"]))


def _rand_rates(r, seed):
    import random
    rng = random.Random(seed)
    return [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(r)]


def _positive_in_kernel(rows, m):
    """Exact: does {x in R^m : rows x = 0} contain a vector with all coordinates > 0?"""
    if m == 0:
        return False
    rays = ecore.crn.extreme_rays_nonneg_kernel([list(map(Fraction, row)) for row in rows], m) if rows else         [tuple(int(i == j) for i in range(m)) for j in range(m)]
    cover = set()
    for v in rays:
        cover |= {i for i in range(m) if v[i] > 0}
    return len(cover) == m


def _A_k(sd, k):
    srcs = sorted({a for a, b in sd.rx if any(a)})
    cols = []
    for c in srcs:
        col = [Fraction(0)] * sd.n
        for t, (a, b) in enumerate(sd.rx):
            if a == c:
                for i in range(sd.n):
                    col[i] += k[t] * (b[i] - a[i])
        cols.append(col)
    return [[cols[c][i] for c in range(len(srcs))] for i in range(sd.n)], srcs


@cand("P-K10", "EC, connected, (H1), (H2), m = #source complexes: if s = m-1 then for k = 1 and 3 seeded random "
      "rational k, ker A_k is a line spanned by a positive vector; if s >= m, some tested k has no positive "
      "steady state", status="proved", sides=True)
def _(sd):
    if not (sd.EC() and sd.connected() and sd.X["cons"] and sd.X["consist"]):
        return None
    if any(not any(a) for a, b in sd.rx):
        return None
    ks = [[Fraction(1)] * sd.r] + [_rand_rates(sd.r, 1000 + j) for j in range(3)]
    m = sd.ndoubling_classes()
    s_ = sd.s()
    if s_ == m - 1:
        for k in ks:
            A, srcs = _A_k(sd, k)
            ker = ecore.crn.nullspace(A, m)
            if len(ker) != 1:
                return False
            v = ker[0]
            if not (all(x > 0 for x in v) or all(x < 0 for x in v)):
                return False
        return True
    if s_ >= m:
        return any(not _positive_in_kernel(_A_k(sd, k)[0], m) for k in ks[1:])
    return None  # s < m-1 is excluded by K5 (tested separately)


def _dual_pss_balanced_ec(sd, kappa):
    """N* (dual of the side's network, which is EC with balanced joins and every species consumed) has a positive
    steady state for rate constants kappa (one per species of the side's network) iff there is lambda > 0 (one per
    source complex) with w_i = kappa_i * lambda_{c(i)} a conservation law of the side's network: linear, exact."""
    srcs = sorted({a for a, b in sd.rx if any(a)})
    cof = {}
    for c, y in enumerate(srcs):
        for i in range(sd.n):
            if y[i] > 0:
                cof[i] = c
    if len(cof) != sd.n:
        return None
    m = len(srcs)
    rows = []
    for a, b in sd.rx:
        row = [Fraction(0)] * m
        for i in range(sd.n):
            row[cof[i]] += kappa[i] * (b[i] - a[i])
        rows.append(row)
    return _positive_in_kernel(rows, m)


@cand("P-K11", "EC with balanced joins, connected, (H1), (H2): [N has a positive steady state for all tested k] = "
      "[N* has one for all tested kappa] = [s = #clusters - 1] (k, kappa: all ones and 3 seeded random rationals)",
      status="proved", sides=True)
def _(sd):
    if not (sd.EC() and sd.balanced_joins() and sd.connected() and sd.X["cons"] and sd.X["consist"]):
        return None
    if any(not any(a) for a, b in sd.rx):
        return None
    m = sd.ndoubling_classes()
    ks = [[Fraction(1)] * sd.r] + [_rand_rates(sd.r, 2000 + j) for j in range(3)]
    kap = [[Fraction(1)] * sd.n] + [_rand_rates(sd.n, 3000 + j) for j in range(3)]
    a1 = all(_positive_in_kernel(_A_k(sd, k)[0], m) for k in ks)
    vals = [_dual_pss_balanced_ec(sd, kp) for kp in kap]
    if any(v is None for v in vals):
        return None
    a2 = all(vals)
    a3 = sd.s() == sd.nclusters() - 1
    return a1 == a2 == a3


_K9 = {"i": lambda sd: all(is_union_of(sd.X["M"], S) for S in sd.X["MS"]),
       "ii": lambda sd: all(any(B <= S for S in sd.X["MS"]) for B in sd.X["M"]),
       "iii": lambda sd: set(sd.X["MS"]) == set(sd.X["SF"]),
       "iv": lambda sd: all(is_union_of(sd.X["MS"], F) for F in sd.X["SF"])}
_SUB = {"EC": lambda sd: sd.EC(), "FC": lambda sd: sd.FC(), "WFFC": lambda sd: sd.WFFC()}
for _k, _f in _K9.items():
    for _cn, _cf in _SUB.items():
        def _mk2(f, cf):
            def g(sd):
                if not (sd.H() and cf(sd)):
                    return None
                return f(sd)
            return g
        cand(f"P-K9{_k}-{_cn}", f"(H1)-(H3) and {_cn} on the side's network: K9({_k})", sides=True)(_mk2(_f, _cf))


# ==============================================================================================
# agent G (gac/notes.md)
# ==============================================================================================
def dual_graph_counts(rec):
    D = CNet(rec["n"], rec["rx"]).dual()
    return len(D.complexes()), D.linkage_classes(), D.deficiency(), D


@cand("G-2.1", "all networks: q = beta* + delta* (q = n - s, beta* = n - m* + l*); if N* WR then l* + delta* <= q; "
      "q = 1 and N* WR => l* = 1, delta* = 0, dual graph one directed cycle", status="proved")
def _(rec):
    n, s = rec["n"], rec["s"]
    q = n - s
    m_, l_, d_, D = dual_graph_counts(rec)
    beta = n - m_ + l_
    ok = q == beta + d_
    if rec["D"]["WR"]:
        ok &= l_ + d_ <= q
        if q == 1:
            ok &= l_ == 1 and d_ == 0 and m_ == D.r
    return ok


@cand("G-2.2", "WR+(H3): source rows distinct, product rows distinct; N* WR iff {source rows} = {product rows}; then "
      "m* = n, l* = c(pi), delta* = q - c(pi), each pi-orbit sums to a conservation law; q = 1 => pi one n-cycle",
      status="proved")
def _(rec):
    N = rec["N"]
    if not (N["WR"] and N["rows_distinct"]):
        return None
    n, rx = rec["n"], rec["rx"]
    A = [tuple(a[i] for a, b in rx) for i in range(n)]
    B = [tuple(b[i] for a, b in rx) for i in range(n)]
    ok = len(set(A)) == n and len(set(B)) == n
    eq = set(A) == set(B)
    ok &= (rec["D"]["WR"] == eq)
    if eq:
        pi = [A.index(B[i]) for i in range(n)]
        seen, cyc = set(), []
        for i in range(n):
            if i not in seen:
                c, x = [], i
                while x not in seen:
                    seen.add(x)
                    c.append(x)
                    x = pi[x]
                cyc.append(c)
        m_, l_, d_, D = dual_graph_counts(rec)
        q = n - rec["s"]
        ok &= m_ == n and l_ == len(cyc) and d_ == q - len(cyc)
        G = CNet(n, rx).Gamma
        for c in cyc:
            ok &= all(sum(G[i][j] for i in c) == 0 for j in range(len(rx)))
        if q == 1:
            ok &= len(cyc) == 1
    return ok


@cand("G-7.1", "WR + (H1) + every row of Ye has sum <= 1: disjoint cycles (each complex source of one reaction), "
      "deficiency 0, no critical siphon; with (H3) every complex is one species", status="proved")
def _(rec):
    N = rec["N"]
    n, rx = rec["n"], rec["rx"]
    if not (N["WR"] and N["cons"]):
        return None
    if any(sum(a[i] for a, b in rx) > 1 for i in range(n)):
        return None
    srcs = [a for a, b in rx]
    ok = len(set(srcs)) == len(srcs)
    ok &= N["deficiency"] == 0 and ADLS(N)
    if N["rows_distinct"]:
        ok &= all(sum(c) == 1 for c in CNet(n, rx).complexes())
    return ok


@cand("G-req1", "(H1)-(H3) for N, N WR: the ADLS criterion holds for N but fails for N* (count of witnesses; "
      "'holds' = no witness)")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H()):
        return None
    return not (ADLS(rec["N"]) and not ADLS(rec["D"]))


@cand("G-req1'", "(H1)-(H3) for N and N*, N WR: ADLS(N) and not ADLS(N*) (witnesses)")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H() and Side(rec, "D").H()):
        return None
    return not (ADLS(rec["N"]) and not ADLS(rec["D"]))


# ==============================================================================================
# transfers through the duality (coordinator, 27 Oct): endotacticity, rate-independent non-persistence
# ==============================================================================================
import endo  # noqa: E402


def _drx(rec):
    n, rx = rec["n"], rec["rx"]
    return tuple((tuple(a[i] for a, b in rx), tuple(b[i] for a, b in rx)) for i in range(n))


@cand("X-endoN", "WR+H, r <= 6 => N endotactic (Craciun-Nazarov-Pantea 2013; consistency check of endo.py)", status="proved")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H()) or rec["r"] > 6:
        return None  # validation of endo.py only; skipped for r > 6 (cost 2^(r-1) r eliminations)
    return endo.endotactic(rec["n"], rec["rx"])


@cand("X-endoD", "WR+H: N* is endotactic")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H()):
        return None
    return endo.endotactic(rec["r"], _drx(rec))


@cand("X-sendoD", "WR+H: N* is strongly endotactic")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H()):
        return None
    return endo.strongly_endotactic(rec["r"], _drx(rec))


@cand("X-sendoN", "WR+H, one linkage class, r <= 6: N strongly endotactic (Gopalkrishnan-Miller-Shiu 2014; check of endo.py)",
      status="proved")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H() and rec["N"]["nlc"] == 1) or rec["r"] > 6:
        return None  # validation of endo.py only; skipped for r > 6
    return endo.strongly_endotactic(rec["n"], rec["rx"])


@cand("X-drainN", "WR+(H1): N has no drain species (else N would be non-persistent for all rates, against the GAC "
      "in its proved cases)", status="proved")
def _(rec):
    if not (rec["N"]["WR"] and rec["N"]["cons"]):
        return None
    return not endo.drain_species(rec["n"], rec["rx"])


@cand("X-drainD", "WR+H: N* has no drain species (a witness is a WR network with (H1)-(H3) whose dual is "
      "non-persistent for EVERY choice of rate constants)")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H()):
        return None
    return not endo.drain_species(rec["r"], _drx(rec))


@cand("X-drainD1", "WR+H, one linkage class, s <= 2: N* has no drain species (witnesses: N globally attracting by "
      "Anderson 2011 / Pantea 2012 for complex-balanced rates, N* non-persistent for all rates)")
def _(rec):
    if not (rec["N"]["WR"] and Side(rec, "N").H() and rec["N"]["nlc"] == 1 and rec["s"] <= 2):
        return None
    return not endo.drain_species(rec["r"], _drx(rec))


@cand("X-L5", "manuscript Lemma 5: the three precluster variants (onestep, thesis, closure) agree, for N under (H2) "
      "and for N* under (H1)", status="proved", sides=True)
def _(sd):
    if not sd.X["consist"]:
        return None
    net = sd.net()
    p1 = set(net.preclusters("onestep"))
    return p1 == set(net.preclusters("thesis")) == set(net.preclusters("closure")) == set(sd.X["P"])


def _is_prime(k):
    return k >= 2 and all(k % d for d in range(2, int(k ** 0.5) + 1))


@cand("X-G4.1", "G Prop. 4.1: N WR, (H3), no zero reaction vector, q = 1, N* WR, n prime => N strongly endotactic",
      status="proved")
def _(rec):
    N = rec["N"]
    if not (N["WR"] and N["rows_distinct"] and rec["n"] - rec["s"] == 1 and rec["D"]["WR"] and _is_prime(rec["n"])):
        return None
    return endo.strongly_endotactic(rec["n"], rec["rx"])


@cand("X-G4.1np", "same hypotheses without primality (n not prime): N strongly endotactic (expected to fail)")
def _(rec):
    N = rec["N"]
    if not (N["WR"] and N["rows_distinct"] and rec["n"] - rec["s"] == 1 and rec["D"]["WR"]
            and not _is_prime(rec["n"])):
        return None
    return endo.strongly_endotactic(rec["n"], rec["rx"])


# ----------------------------------------------------------------------------------------------
# runner
# ----------------------------------------------------------------------------------------------
def load(files):
    seen = {}
    for fn in files:
        with gzip.open(fn, "rb") as fh:
            for r in pickle.load(fh):
                key = (r["n"], r["rx"])
                if key not in seen:
                    seen[key] = r
    return list(seen.values())


def classes(rec):
    sd = Side(rec, "N")
    h = sd.H()
    out = ["all"]
    if h:
        out.append("H")
        if rec["rows_pi"] and rec["cols_pi"]:
            out.append("strong")
        if rec["N"]["WR"]:
            out.append("WR")
            if rec["N"]["deficiency"] == 0:
                out.append("WR0")
            if rec["s"] in (2, 3, 4):
                out.append(f"WR,s={rec['s']}")
    return out


CLS = ["all", "H", "strong", "WR", "WR0", "WR,s=2", "WR,s=3", "WR,s=4"]


RECS = None
CAND_SEL = None


def eval_range(bounds):
    lo, hi = bounds
    out = []
    for idx in range(lo, hi):
        r = RECS[idx]
        res = []
        for ci, c in enumerate(CAND_SEL):
            if c["sides"]:
                for sd in sides(r):
                    v = c["f"](sd)
                    if v is not None:
                        res.append((ci, bool(v), sd.which))
            else:
                v = c["f"](r)
                if v is not None:
                    res.append((ci, bool(v), "N"))
        out.append((idx, classes(r), res))
    return out


def main():
    global RECS, CAND_SEL
    argv = sys.argv[1:]
    only = None
    jout = None
    nproc = 1
    for flag in ("--only", "--json", "--nproc"):
        if flag in argv:
            i = argv.index(flag)
            val = argv[i + 1]
            del argv[i:i + 2]
            if flag == "--only":
                only = val
            elif flag == "--json":
                jout = val
            else:
                nproc = int(val)
    t0 = time.time()
    print(f"# tests.py  Python {platform.python_version()} on {platform.system()} {platform.release()}, nproc={nproc}")
    RECS = load(argv)
    recs = RECS
    print(f"# {len(recs)} distinct networks from {len(argv)} record files: {', '.join(a.rsplit('/', 1)[-1] for a in argv)}")
    fams = {}
    for r in recs:
        fams[str(r.get("family"))] = fams.get(str(r.get("family")), 0) + 1
    print(f"# families: {fams}")
    if only is None:
        cands = [c for c in CANDS if not c["id"].startswith("X-")]
    else:
        pre = only.split(",")
        cands = [c for c in CANDS if any(c["id"].startswith(p) for p in pre)]
    CAND_SEL = cands
    print(f"# candidates: {' '.join(c['id'] for c in cands)}")
    sys.stdout.flush()
    stats = {c["id"]: {k: [0, 0] for k in CLS} for c in cands}
    first = {}
    step = 200
    ranges = [(i, min(i + step, len(recs))) for i in range(0, len(recs), step)]
    if nproc > 1:
        from multiprocessing import Pool
        pool = Pool(nproc)
        it = pool.imap_unordered(eval_range, ranges)
    else:
        it = map(eval_range, ranges)
    for chunk in it:
        for idx, cl, res in chunk:
            r = recs[idx]
            for ci, v, which in res:
                c = cands[ci]
                for k in cl:
                    stats[c["id"]][k][0] += 1
                    stats[c["id"]][k][1] += v
                if not v:
                    key = (r["n"] + r["r"], r["r"], r["rx"], which)
                    if c["id"] not in first or key < first[c["id"]][0]:
                        first[c["id"]] = (key, r, which)
    if nproc > 1:
        pool.close()
    ok_all = True
    out = {}
    print()
    print(f"{'candidate':14s} {'status':7s} " + " ".join(f"{k:>15s}" for k in CLS) + "  verdict")
    for c in cands:
        st = stats[c["id"]]
        cells = " ".join(f"{st[k][1]:>7d}/{st[k][0]:<7d}" for k in CLS)
        fail = c["id"] in first
        if c["status"] == "proved":
            verdict = "FAIL" if fail else ("PASS" if st["all"][0] else "PASS (vacuous)")
            ok_all &= not fail
        else:
            verdict = "REFUTED" if fail else ("holds on all tested" if st["all"][0] else "no applicable network")
        print(f"{c['id']:14s} {c['status']:7s} {cells}  {verdict}")
        out[c["id"]] = dict(doc=c["doc"], status=c["status"], stats=st, verdict=verdict)
    print()
    print("# statements and smallest counterexamples (side N: the network itself; side D: its dual)")
    for c in cands:
        print(f"{c['id']}: {c['doc']}")
        if c["id"] in first:
            key, r, which = first[c["id"]]
            N = CNet(r["n"], r["rx"])
            D = N.dual()
            print(f"   counterexample (side {which}; n={r['n']}, r={r['r']}, s={r['s']}, H={Side(r, 'N').H()}, "
                  f"H*={Side(r, 'D').H()}, WR={r['N']['WR']}, WR*={r['D']['WR']}, def={r['N']['deficiency']}, "
                  f"family={r.get('family')}):")
            print(f"     N : {N.describe()}")
            print(f"     N*: {D.describe()}")
            out[c["id"]]["counterexample"] = dict(side=which, N=N.describe(), Nstar=D.describe())
    if jout:
        with open(jout, "w") as fh:
            json.dump(out, fh, indent=1)
    print()
    print(f"{'PASS' if ok_all else 'FAIL'}: every candidate marked 'proved' holds on all applicable networks")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
