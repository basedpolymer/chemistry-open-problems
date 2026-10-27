"""Petri-net view of reaction networks and of Iosif's network duality (agent P, Claude Opus 5.5, effort xhigh).

Environment: Python 3.12 (WSL Fedora 44 or Windows), standard library only.  Exact rational arithmetic through
crn.py (fractions.Fraction); no floating point anywhere.

A network N (crn.Network) is read as a weighted Petri net: places = species, transitions = reactions,
Pre(X_i, R_j) = (Ye)_ij, Post(X_i, R_j) = (Yp)_ij.  Iosif's dual N* is the weighted reverse-dual
rd(N) = (T, S, F^-1) (checked by `is_reverse_dual`).

Objects implemented here (crn.py is imported, never edited):
  * clusters (Desel-Esparza Def. 2.?): connected components of the undirected bipartite graph with an edge
    {X_i, R_j} iff (Ye)_ij > 0; returned as pairs (species set, reaction set);
  * equal-conflict sets (Teruel-Silva Def. 5.3): classes of reactions with the same NONZERO source complex
    (reactions with the zero source complex are singletons);  twin classes = the same for species (equal nonzero
    rows of Ye) = equal-conflict sets of N*;
  * ordinary, free choice (FC, ordinary nets: reactions sharing a reactant species have the same source),
    equal conflict (EC, weighted: same definition), balanced joins (every source complex has equal nonzero
    coefficients);
  * connectedness of the Petri net, the rank equation rank(Gamma) = #clusters - 1 (Desel-Esparza Thm 6.14) and
    rank(Gamma) = #equal-conflict sets - 1 (Teruel-Silva Thm 28);
  * T-traps (sets U of reactions with U-dot subset dot-U: every species produced by U is consumed by U) = siphons
    of N*;
  * a behavioural (token game) check of liveness for bounded nets, used to test the Rank Theorem independently.
"""
import os
import sys
from itertools import combinations
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CRN_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "scripts"))
sys.path.insert(0, CRN_DIR)
from crn import (Network, rank, refines, normalize_partition, partition_from_components,  # noqa: E402,F401
                 components, direction_key)


# ------------------------------------------------------------------------------------------------
# basic structure
# ------------------------------------------------------------------------------------------------
def is_reverse_dual(N, D):
    """D is the weighted reverse-dual of N: places of D = transitions of N, transitions of D = places of N,
    Pre_D(R_j, X_i) = Pre_N(X_i, R_j) and Post_D(R_j, X_i) = Post_N(X_i, R_j)."""
    if D.n != N.r or D.r != N.n:
        return False
    YeN, YpN, YeD, YpD = N.Ye, N.Yp, D.Ye, D.Yp
    return all(YeD[j][i] == YeN[i][j] and YpD[j][i] == YpN[i][j] for i in range(N.n) for j in range(N.r))


def src_support(N, j):
    return frozenset(i for i in range(N.n) if N.reactions[j][0][i] > 0)


def prod_support(N, j):
    return frozenset(i for i in range(N.n) if N.reactions[j][1][i] > 0)


def consumers(N, i):
    return frozenset(j for j in range(N.r) if N.reactions[j][0][i] > 0)


def producers(N, i):
    return frozenset(j for j in range(N.r) if N.reactions[j][1][i] > 0)


def is_ordinary(N):
    return all(c in (0, 1) for a, b in N.reactions for c in a + b)


def clusters(N):
    """Clusters: components of the bipartite graph species--reactions, edge iff the species is a reactant.
    Nodes 0..n-1 are species, n..n+r-1 reactions.  Returns a sorted list of (species frozenset, reaction frozenset).
    The definition (x in C, x place => x-dot in C; x transition => dot-x in C) is the same in N and rd(N)."""
    n, r = N.n, N.r
    edges = [(i, n + j) for j in range(r) for i in range(n) if N.reactions[j][0][i] > 0]
    comp = components(n + r, edges)
    blocks = {}
    for v in range(n + r):
        blocks.setdefault(comp[v], set()).add(v)
    out = []
    for b in blocks.values():
        out.append((frozenset(v for v in b if v < n), frozenset(v - n for v in b if v >= n)))
    return sorted(out, key=lambda c: (sorted(c[0]), sorted(c[1])))


def ec_sets(N):
    """Equal-conflict sets: reactions with the same nonzero source complex; zero-source reactions are singletons."""
    blocks = {}
    for j, (a, b) in enumerate(N.reactions):
        key = a if any(a) else ("zero", j)
        blocks.setdefault(key, []).append(j)
    return normalize_partition(blocks.values())


def twin_classes(N):
    """Species with equal nonzero rows of Ye (equal-conflict sets of N*); species never consumed are singletons."""
    Ye = N.Ye
    blocks = {}
    for i in range(N.n):
        key = tuple(Ye[i]) if any(Ye[i]) else ("zero", i)
        blocks.setdefault(key, []).append(i)
    return normalize_partition(blocks.values())


def coupled_conflict_sets(N):
    """Transitive closure of 'share a reactant species' (transition parts of clusters with transitions)."""
    return normalize_partition([c[1] for c in clusters(N) if c[1]])


def is_EC(N):
    """Equal conflict: two reactions sharing a reactant species have the same source complex."""
    for j, k in combinations(range(N.r), 2):
        if src_support(N, j) & src_support(N, k) and N.reactions[j][0] != N.reactions[k][0]:
            return False
    return True


def is_free_choice(N):
    """(Extended) free choice in the sense of Desel-Esparza: ordinary and EC."""
    return is_ordinary(N) and is_EC(N)


def is_balanced_joins(N):
    """Every source complex has all its nonzero coefficients equal (all input arcs of a transition equal)."""
    for a, _ in N.reactions:
        nz = {c for c in a if c > 0}
        if len(nz) > 1:
            return False
    return True


def is_connected(N):
    """The Petri net (bipartite graph with an edge for each nonzero entry of Ye or Yp) is connected."""
    n, r = N.n, N.r
    edges = [(i, n + j) for j in range(r) for i in range(n)
             if N.reactions[j][0][i] > 0 or N.reactions[j][1][i] > 0]
    comp = components(n + r, edges)
    return len(set(comp)) == 1


def rank_equation_clusters(N):
    return N.rank() == len(clusters(N)) - 1


def rank_equation_ec(N):
    return N.rank() == len(ec_sets(N)) - 1


def well_formed_by_rank_theorem(N):
    """Rank Theorem (Desel-Esparza Thm 6.14 for FC; Teruel-Silva Thm 28 for EC): for an EC net, well-formed iff
    connected (with a place and a transition), conservative, consistent and rank = #EC sets - 1.  For FC nets the
    clusters with transitions are exactly the EC sets, and every cluster has a transition when the net is
    consistent, conservative and every species occurs; we return None if N is not EC."""
    if not is_EC(N):
        return None
    return (N.n > 0 and N.r > 0 and is_connected(N) and N.is_conservative() and N.is_consistent()
            and rank_equation_ec(N))


# ------------------------------------------------------------------------------------------------
# T-traps, T-siphons, components
# ------------------------------------------------------------------------------------------------
def is_t_trap(N, U):
    """U-dot subset dot-U: every species produced by a reaction of U is consumed by a reaction of U."""
    U = set(U)
    if not U:
        return False
    prod = set().union(*(prod_support(N, j) for j in U))
    cons = set().union(*(src_support(N, j) for j in U))
    return prod <= cons


def is_t_siphon(N, U):
    """dot-U subset U-dot: every species consumed by U is produced by U (reverse-dual of a trap)."""
    U = set(U)
    if not U:
        return False
    prod = set().union(*(prod_support(N, j) for j in U))
    cons = set().union(*(src_support(N, j) for j in U))
    return cons <= prod


def is_trap(N, S):
    S = set(S)
    if not S:
        return False
    for a, b in N.reactions:
        if any(a[i] > 0 for i in S) and not any(b[i] > 0 for i in S):
            return False
    return True


def minimal_sets(family):
    fam = list(family)
    return sorted([S for S in fam if not any(T < S for T in fam)], key=lambda s: (len(s), sorted(s)))


def t_traps(N):
    return [frozenset(U) for k in range(1, N.r + 1) for U in combinations(range(N.r), k) if is_t_trap(N, U)]


def minimal_t_traps(N):
    return minimal_sets(t_traps(N))


def critical_siphons(N):
    """Siphons containing no support of a semi-positive conservation law (Angeli-De Leenheer-Sontag)."""
    sup = [frozenset(i for i in range(N.n) if z[i] > 0) for z in N.conservation_rays()]
    return [S for S in N.siphons() if not any(T <= S for T in sup)]


def is_s_component_set(N, S):
    """Ordinary S-component (Desel-Esparza Def. 5.1) with place set S: the subnet generated by S and
    dot-S cup S-dot is a strongly connected S-net, i.e. every reaction touching S has exactly one reactant and one
    product species in S, with coefficient 1, and the subnet is strongly connected."""
    S = set(S)
    T = [j for j in range(N.r) if src_support(N, j) & S or prod_support(N, j) & S]
    for j in T:
        a, b = N.reactions[j]
        ins = [i for i in S if a[i] > 0]
        outs = [i for i in S if b[i] > 0]
        if len(ins) != 1 or len(outs) != 1 or a[ins[0]] != 1 or b[outs[0]] != 1:
            return False
    # strong connectivity of the S-net (species graph i -> k if a reaction of T has i in source and k in product)
    adj = {i: set() for i in S}
    for j in T:
        a, b = N.reactions[j]
        i = next(i for i in S if a[i] > 0)
        k = next(k for k in S if b[k] > 0)
        adj[i].add(k)
    return _strongly_connected(S, adj)


def is_t_component_set(N, U):
    """Ordinary T-component with transition set U: every species touched by U is consumed by exactly one and
    produced by exactly one reaction of U, with coefficient 1, and the T-net is strongly connected."""
    U = set(U)
    if not U:
        return False
    P = set().union(*(src_support(N, j) | prod_support(N, j) for j in U))
    for i in P:
        c = [j for j in U if N.reactions[j][0][i] > 0]
        p = [j for j in U if N.reactions[j][1][i] > 0]
        if len(c) != 1 or len(p) != 1 or N.reactions[c[0]][0][i] != 1 or N.reactions[p[0]][1][i] != 1:
            return False
    adj = {j: set() for j in U}
    for i in P:
        c = next(j for j in U if N.reactions[j][0][i] > 0)
        p = next(j for j in U if N.reactions[j][1][i] > 0)
        adj[p].add(c)
    return _strongly_connected(U, adj)


def _strongly_connected(V, adj):
    V = list(V)
    if not V:
        return False

    def reach(start, nbr):
        seen = {start}
        stack = [start]
        while stack:
            x = stack.pop()
            for y in nbr[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        return seen
    radj = {v: set() for v in V}
    for v in V:
        for w in adj[v]:
            radj[w].add(v)
    return reach(V[0], adj) == set(V) and reach(V[0], radj) == set(V)


# ------------------------------------------------------------------------------------------------
# token game (behavioural check, bounded nets only)
# ------------------------------------------------------------------------------------------------
def is_live_and_bounded(N, M0, max_states=200000):
    """Explores the reachability graph from M0.  Returns (bounded_and_explored, live).  Conservative nets are
    bounded; exploration stops (returns (False, None)) beyond max_states markings."""
    Ye = [N.reactions[j][0] for j in range(N.r)]
    G = [tuple(N.reactions[j][1][i] - N.reactions[j][0][i] for i in range(N.n)) for j in range(N.r)]
    start = tuple(M0)
    index = {start: 0}
    states = [start]
    succ = [[]]
    k = 0
    while k < len(states):
        M = states[k]
        for j in range(N.r):
            if all(M[i] >= Ye[j][i] for i in range(N.n)):
                M2 = tuple(M[i] + G[j][i] for i in range(N.n))
                if M2 not in index:
                    if len(states) >= max_states:
                        return (False, None)
                    index[M2] = len(states)
                    states.append(M2)
                    succ.append([])
                succ[k].append((j, index[M2]))
        k += 1
    # Tarjan SCC (iterative)
    nS = len(states)
    idx = [None] * nS
    low = [0] * nS
    onst = [False] * nS
    st = []
    comp = [None] * nS
    counter = 0
    ncomp = 0
    for v0 in range(nS):
        if idx[v0] is not None:
            continue
        work = [(v0, 0)]
        idx[v0] = low[v0] = counter
        counter += 1
        st.append(v0)
        onst[v0] = True
        while work:
            v, pi = work[-1]
            if pi < len(succ[v]):
                work[-1] = (v, pi + 1)
                w = succ[v][pi][1]
                if idx[w] is None:
                    idx[w] = low[w] = counter
                    counter += 1
                    st.append(w)
                    onst[w] = True
                    work.append((w, 0))
                elif onst[w]:
                    low[v] = min(low[v], idx[w])
            else:
                work.pop()
                if work:
                    u = work[-1][0]
                    low[u] = min(low[u], low[v])
                if low[v] == idx[v]:
                    while True:
                        w = st.pop()
                        onst[w] = False
                        comp[w] = ncomp
                        if w == v:
                            break
                    ncomp += 1
    bottom = [True] * ncomp
    labels = [set() for _ in range(ncomp)]
    for v in range(nS):
        for j, w in succ[v]:
            if comp[w] != comp[v]:
                bottom[comp[v]] = False
            else:
                labels[comp[v]].add(j)
    live = all(labels[c] == set(range(N.r)) for c in range(ncomp) if bottom[c])
    return (True, live)


def ts_marking(N):
    """Marking of Teruel-Silva's proof of Thm 28 (3)=>(1): M0[p] = W(p, t) for some t in p-dot (here the maximum
    over consumers, which is the same number in an EC net)."""
    Ye = N.Ye
    return [max(Ye[i]) if any(Ye[i]) else 0 for i in range(N.n)]


# ------------------------------------------------------------------------------------------------
# precluster helpers
# ------------------------------------------------------------------------------------------------
def precluster_by_profiles(r, pairs, profile):
    """Precluster graph G+ of Definition 19 computed from 0/1 'profiles' (frozensets): for a doubling pair {i,j}
    with disjoint nonempty profiles, s is adjoined iff profile(s) is profile(i), profile(j) or their union.
    Valid when the vectors n_s are the indicator vectors of the profiles (well-formed free-choice nets)."""
    edges = set(pairs)
    for (i, j) in pairs:
        Pi, Pj = profile[i], profile[j]
        for s in range(r):
            if s in (i, j):
                continue
            if profile[s] in (Pi, Pj, Pi | Pj):
                edges.add(tuple(sorted((i, s))))
                edges.add(tuple(sorted((j, s))))
    return partition_from_components(r, edges)


def span_closure_preclusters(N):
    """Variant (not Iosif's): the span rule applied to whole classes.  Start from the doubling-pair components
    (equal-conflict sets); a reaction s joins a class B when n_s lies in span{n_t : t in B}; iterate on the merged
    classes.  This is the closure suggested by the ratio argument of Conradi-Flockerzi (Fact 2(a))."""
    nv = N.n_vectors()
    r = N.r
    part = [set(b) for b in ec_sets(N)]
    changed = True
    while changed:
        changed = False
        for B in part:
            if len(B) < 2:
                continue
            base = [list(nv[t]) for t in B]
            rb = rank(base)
            for C in part:
                if C is B:
                    continue
                if any(rank(base + [list(nv[s])]) == rb for s in C):
                    B |= C
                    part.remove(C)
                    changed = True
                    break
            if changed:
                break
    return normalize_partition(part)
