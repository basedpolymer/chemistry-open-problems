"""Core helpers of agent E (computation and red team) for the gac-duality round of result 17.

Environment: Python 3.12 (WSL Fedora 44), standard library only.  Exact rational arithmetic throughout (the linear
algebra is that of ../../../scripts/crn.py, imported, not edited; fractions.Fraction, no floating point).

Contents
  * CNet: crn.Network with cached cones; the dual shares the caches (flux rays of N* = conservation rays of N).
  * canonical(n, rx): canonical form of a reaction set under permutations of the species (complete invariant).
  * WR generator: weakly reversible networks as disjoint unions of strongly connected linkage classes.
  * all-network generator (every reaction set), as in enumerate_small.py of result 17.
  * record(rx, n): every object used by the tests, for N and for its dual N*.
Conventions: species 0..n-1, reactions 0..r-1; a reaction is a pair (source, product) of n-tuples of nonnegative
integers with source != product.  Subsets are frozensets of indices; partitions are tuples of frozensets.
"""
import os
import sys
import math
from itertools import combinations, permutations, product

HERE = os.path.dirname(os.path.abspath(__file__))
CRN_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "scripts"))
sys.path.insert(0, CRN_DIR)
import crn  # noqa: E402
from crn import Network, refines, rank, normalize_partition, partition_from_components  # noqa: E402,F401


# ----------------------------------------------------------------------------------------------
# cached network
# ----------------------------------------------------------------------------------------------
class CNet(Network):
    def __init__(self, n, reactions, names=None):
        super().__init__(n, reactions, names)
        self._flux = None
        self._cons = None
        self._dual = None
        self._G = None

    @property
    def Gamma(self):
        if self._G is None:
            self._G = Network.Gamma.fget(self)
        return self._G

    def flux_rays(self):
        if self._flux is None:
            if self._dual is not None and self._dual._cons is not None:
                self._flux = self._dual._cons
            else:
                self._flux = Network.flux_rays(self)
        return self._flux

    def conservation_rays(self):
        if self._cons is None:
            if self._dual is not None and self._dual._flux is not None:
                self._cons = self._dual._flux
            else:
                self._cons = Network.conservation_rays(self)
        return self._cons

    def dual(self):
        if self._dual is None:
            Ye, Yp = self.Ye, self.Yp
            D = CNet(self.r, [(tuple(Ye[i]), tuple(Yp[i])) for i in range(self.n)],
                     names=[f"k{j+1}" for j in range(self.r)])
            D._dual = self
            self._dual = D
        return self._dual


# ----------------------------------------------------------------------------------------------
# canonical form under species permutations
# ----------------------------------------------------------------------------------------------
def _signature(n, rx):
    out = []
    for i in range(n):
        out.append(tuple(sorted((a[i], b[i], tuple(sorted(a)), tuple(sorted(b))) for a, b in rx)))
    return out


def canonical(n, rx):
    """Minimum, over the orderings of the species that are nondecreasing in a permutation-invariant signature, of
    the sorted reaction list.  The set of such orderings is equivariant, so this is a complete invariant of the
    orbit of rx under the symmetric group on species."""
    sig = _signature(n, rx)
    order = sorted(range(n), key=lambda i: sig[i])
    groups = []
    for i in order:
        if groups and sig[groups[-1][0]] == sig[i]:
            groups[-1].append(i)
        else:
            groups.append([i])
    best = None
    for choice in product(*[permutations(g) for g in groups]):
        perm = [i for g in choice for i in g]
        key = tuple(sorted((tuple(a[p] for p in perm), tuple(b[p] for p in perm)) for a, b in rx))
        if best is None or key < best:
            best = key
    return best


# ----------------------------------------------------------------------------------------------
# generators
# ----------------------------------------------------------------------------------------------
def complexes(n, m, zero=False):
    cs = [c for c in product(range(m + 1), repeat=n) if sum(c) <= m and (zero or sum(c) > 0)]
    return sorted(cs)


def _strongly_connected(k, arcs):
    adj = [[] for _ in range(k)]
    radj = [[] for _ in range(k)]
    for i, j in arcs:
        adj[i].append(j)
        radj[j].append(i)
    for g in (adj, radj):
        seen = {0}
        st = [0]
        while st:
            x = st.pop()
            for y in g[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        if len(seen) != k:
            return False
    return True


def sc_patterns(k, emax, limit=4_000_000):
    """Labelled strongly connected spanning digraphs on {0..k-1} (no loops) with k <= arcs <= emax."""
    pats = []
    if k < 2:
        return pats
    arcs = [(i, j) for i in range(k) for j in range(k) if i != j]
    for e in range(k, min(emax, len(arcs)) + 1):
        if e == k:  # exactly the Hamiltonian cycles
            seen = set()
            for p in permutations(range(1, k)):
                cyc = (0,) + p
                s = tuple(sorted((cyc[t], cyc[(t + 1) % k]) for t in range(k)))
                if s not in seen:
                    seen.add(s)
                    pats.append(s)
        else:
            if math.comb(len(arcs), e) > limit:
                raise ValueError(f"sc_patterns: k={k}, e={e} too large for brute force")
            for S in combinations(arcs, e):
                if _strongly_connected(k, S):
                    pats.append(S)
    return pats


def lc_instances(C, R):
    """All linkage classes (complex subset, strongly connected arc set) on the complexes C with <= R reactions.
    Returns a list of (mask, reactions) with reactions a tuple of (source, product) pairs."""
    out = []
    for k in range(2, min(R, len(C)) + 1):
        pats = sc_patterns(k, R)
        for sub in combinations(range(len(C)), k):
            mask = 0
            for t in sub:
                mask |= 1 << t
            for pat in pats:
                out.append((mask, tuple((C[sub[i]], C[sub[j]]) for i, j in pat)))
    out.sort(key=lambda t: len(t[1]))  # by number of reactions (used to stop the search early)
    return out


def wr_networks_from(first, inst, R, n):
    """All WR reaction sets whose linkage classes are inst[first] and instances of larger index (disjoint
    complexes, <= R reactions in total).  Yields sorted reaction tuples in which every species occurs."""
    m0, r0 = inst[first]

    def rec(start, mask, rx):
        if all(any(a[i] or b[i] for a, b in rx) for i in range(n)):
            yield tuple(sorted(rx))
        for t in range(start, len(inst)):
            mt, rt = inst[t]
            if len(rx) + len(rt) > R:
                break  # inst is sorted by number of reactions
            if mt & mask:
                continue
            yield from rec(t + 1, mask | mt, rx + list(rt))
    yield from rec(first + 1, m0, list(r0))


# ----------------------------------------------------------------------------------------------
# bitmask combinatorics: siphons, traps, clusters, free choice
# ----------------------------------------------------------------------------------------------
def _masks(net):
    src = []
    prd = []
    for a, b in net.reactions:
        sm = 0
        pm = 0
        for i in range(net.n):
            if a[i] > 0:
                sm |= 1 << i
            if b[i] > 0:
                pm |= 1 << i
        src.append(sm)
        prd.append(pm)
    return src, prd


def to_set(mask):
    out = []
    i = 0
    while mask:
        if mask & 1:
            out.append(i)
        mask >>= 1
        i += 1
    return frozenset(out)


def siphon_masks(net):
    src, prd = _masks(net)
    out = []
    for S in range(1, 1 << net.n):
        if all(not (p & S) or (s & S) for s, p in zip(src, prd)):
            out.append(S)
    return out


def trap_masks(net):
    src, prd = _masks(net)
    out = []
    for S in range(1, 1 << net.n):
        if all(not (s & S) or (p & S) for s, p in zip(src, prd)):
            out.append(S)
    return out


def minimal_masks(fam):
    fam = sorted(fam, key=lambda x: bin(x).count("1"))
    mins = []
    for S in fam:
        if not any((T & S) == T for T in mins):
            mins.append(S)
    return mins


def is_siphon(net, S):
    S = set(S)
    if not S:
        return False
    for a, b in net.reactions:
        if any(b[i] > 0 for i in S) and not any(a[i] > 0 for i in S):
            return False
    return True


def is_trap(net, S):
    S = set(S)
    if not S:
        return False
    for a, b in net.reactions:
        if any(a[i] > 0 for i in S) and not any(b[i] > 0 for i in S):
            return False
    return True


def petri_clusters(net):
    """Petri-net clusters: connected components of the bipartite graph species i -- reaction j iff (Ye)_ij > 0.
    Returns (partition of the reactions, partition of the species) induced by the components."""
    n, r = net.n, net.r
    parent = list(range(n + r))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for j, (a, b) in enumerate(net.reactions):
        for i in range(n):
            if a[i] > 0:
                parent[find(i)] = find(n + j)
    comp = [find(x) for x in range(n + r)]
    rb, sb = {}, {}
    for j in range(r):
        rb.setdefault(comp[n + j], []).append(j)
    for i in range(n):
        sb.setdefault(comp[i], []).append(i)
    return normalize_partition(rb.values()), normalize_partition(sb.values())


def free_choice_flags(net):
    """(FC, EFC, EC): free choice (every arc s->t: s* = {t} or *t = {s}), extended free choice (reactions sharing a
    reactant species have the same reactant set), equal conflict (... the same source complex)."""
    src, _ = _masks(net)
    n, r = net.n, net.r
    post = [[j for j in range(r) if src[j] >> i & 1] for i in range(n)]
    fc = True
    for i in range(n):
        for j in post[i]:
            if not (len(post[i]) == 1 or src[j] == (1 << i)):
                fc = False
    efc = True
    ec = True
    for j, k in combinations(range(r), 2):
        if src[j] & src[k]:
            if src[j] != src[k]:
                efc = False
            if net.reactions[j][0] != net.reactions[k][0]:
                ec = False
    return fc, efc, ec


def semiflow_supports(net):
    return sorted({frozenset(i for i in range(net.n) if z[i] > 0) for z in net.conservation_rays()},
                  key=lambda s: (len(s), sorted(s)))


# ----------------------------------------------------------------------------------------------
# the record of a network
# ----------------------------------------------------------------------------------------------
def part(P):
    return tuple(sorted((frozenset(b) for b in P), key=lambda s: (min(s), len(s))))


def side_record(N):
    """Objects of one network (no dual-dependent items)."""
    d = {}
    d["n"], d["r"] = N.n, N.r
    d["cons"] = N.is_conservative()
    d["consist"] = N.is_consistent()
    d["rows_distinct"] = N.rows_distinct()
    d["WR"] = N.weakly_reversible()
    d["deficiency"] = N.deficiency()
    d["nlc"] = N.linkage_classes()
    d["ncomplexes"] = len(N.complexes())
    d["has_loop"] = any(a == b for a, b in N.reactions)
    d["has_zero_complex"] = any(sum(c) == 0 for c in N.complexes())
    d["P"] = part(N.preclusters("onestep"))  # meaningful under consistency (H2); see Theorem 7 of result 17
    d["M"] = part(N.mips_ray()) if d["cons"] else None
    d["SF"] = tuple(semiflow_supports(N))
    d["IC"] = tuple(N.internal_cycles())
    sy = siphon_masks(N)
    tr = trap_masks(N)
    d["nsiphons"] = len(sy)
    d["MS"] = tuple(to_set(S) for S in minimal_masks(sy))
    d["MT"] = tuple(to_set(S) for S in minimal_masks(tr))
    sfm = [sum(1 << i for i in F) for F in d["SF"]]
    d["CS"] = tuple(to_set(S) for S in minimal_masks(sy) if not any((F & S) == F for F in sfm))
    d["siphons_eq_SF"] = set(d["MS"]) == set(d["SF"])
    d["CLr"], d["CLs"] = petri_clusters(N)
    d["FC"], d["EFC"], d["EC"] = free_choice_flags(N)
    d["doubling"] = tuple(N.doubling_pairs())
    return d


def record(rx, n):
    """Record of N (key 'N') and of its dual (key 'D'), plus global items."""
    N = CNet(n, rx)
    D = N.dual()
    rec = {"rx": tuple(rx), "n": n, "r": len(rx), "s": N.rank()}
    rec["N"] = side_record(N)
    rec["D"] = side_record(D)
    rec["rows_pi"] = N.rows_pairwise_independent()
    rec["cols_pi"] = N.columns_pairwise_independent()
    return rec


def describe(rx, n, names=None):
    return CNet(n, rx, names).describe()


def H(rec):
    """Hypotheses (H1)-(H3) of Conjecture 22 for N."""
    N = rec["N"]
    return N["cons"] and N["consist"] and N["rows_distinct"]


def Hd(rec):
    """(H1)-(H3) for N*."""
    D = rec["D"]
    return D["cons"] and D["consist"] and D["rows_distinct"]
