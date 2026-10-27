"""Splitting-invariant preclusters: flat closure, Conradi-Flockerzi clusters (agent Q2, Claude Opus 5.5, xhigh).

Environment: Python 3.12, standard library only for everything except `one_block_certificate` (sympy, WSL).
Exact rational arithmetic (fractions.Fraction, via crn.py); no floating point.

Imports crn.py (problem scripts/) and petri.py (gac-duality/petri/scripts/) without editing them.

Notation (see ../notes.md).  N has species 0..n-1 and reactions 0..r-1; Gamma = Yp - Ye.
  * source classes (doubling classes): reactions with the same source complex;  Q = r x m indicator matrix;
  * C = ker Gamma cap R^r_{>=0} (flux cone), C0 = C cap R^r_{>0} (nonempty under (H2));
  * for V in C0:  A_V = Gamma diag(V) Q  (n x m), column c = net production of the source complex c under V;
  * kappa_V = partition of the reactions: t ~ t' iff e_{c(t)} - e_{c(t')} lies in the row space of A_V
              (= the clusters of Iosif's thesis, Def. 4.6, for the rate constants k = V and steady state x = 1);
  * CF = meet of kappa_V over all V in C0 (= Conradi-Flockerzi clusters, Fact 2(b) / Algorithm 1, rate-free);
  * pi = flat closure: the finest partition coarser than the source classes whose blocks are flats of the dual
         matroid M*(Gamma) (B is a flat iff every v in ker Gamma vanishing on B vanishes on n_s in span(B) ...);
  * P  = Iosif's preclusters (Definition 19); Pcs = agent P's class-span closure; M* = full-coupling classes.
Proved in notes.md: P <= Pcs <= pi <= CF <= kappa_V (refinement order), pi = Pcs v M*.
"""
import os
import sys
import random
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
CRN_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "scripts"))
PETRI_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "petri", "scripts"))
sys.path.insert(0, CRN_DIR)
sys.path.insert(0, PETRI_DIR)
from crn import Network, rank, rref, nullspace, normalize_partition, refines, direction_key  # noqa: E402,F401
import petri  # noqa: E402


# ------------------------------------------------------------------------------------------------
# partitions
# ------------------------------------------------------------------------------------------------
def meet(P, Q):
    out = []
    for p in P:
        for q in Q:
            x = p & q
            if x:
                out.append(x)
    return normalize_partition(out)


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
    for i in range(m):
        blocks.setdefault(find(i), []).append(i)
    return normalize_partition(blocks.values())


def part_str(P, names):
    return "{" + ", ".join("{" + ",".join(names[i] for i in sorted(b)) + "}" for b in P) + "}"


# ------------------------------------------------------------------------------------------------
# basic objects
# ------------------------------------------------------------------------------------------------
def source_classes(N):
    """Reactions grouped by source complex (the zero source is treated like any other complex)."""
    blocks = {}
    for j, (a, b) in enumerate(N.reactions):
        blocks.setdefault(a, []).append(j)
    return normalize_partition(blocks.values())


def nvec(N):
    """n_t = rows of the extreme-ray matrix of C (coordinate functionals of ker Gamma under (H2))."""
    return [list(map(Fraction, v)) for v in N.n_vectors()]


def full_coupling(N):
    """Classes of fully coupled reactions (n_s positively proportional) = MIPS(N*) (manuscript Prop. 9)."""
    nv = N.n_vectors()
    key = {}
    for t, v in enumerate(nv):
        key.setdefault(direction_key(v), []).append(t)
    return normalize_partition(key.values())


def iosif_P(N):
    return N.preclusters("onestep")


def P_cs(N):
    return petri.span_closure_preclusters(N)


def flat_closure(N, start=None):
    """pi(N): finest partition coarser than `start` (default: source classes) whose blocks are flats of M*(Gamma),
    i.e. closed under: n_s in span{n_t : t in B}  =>  s in B."""
    nv = nvec(N)
    part = [set(b) for b in (start if start is not None else source_classes(N))]
    changed = True
    while changed:
        changed = False
        for B in part:
            base = [nv[t] for t in B]
            rb = rank(base) if base and base[0] else 0
            for C in part:
                if C is B:
                    continue
                if any((rank(base + [nv[s]]) if nv[s] else 0) == rb for s in C):
                    B |= C
                    part.remove(C)
                    changed = True
                    break
            if changed:
                break
    return normalize_partition(part)


# ------------------------------------------------------------------------------------------------
# A_V and kappa_V
# ------------------------------------------------------------------------------------------------
def A_matrix(N, V, classes=None):
    classes = classes if classes is not None else source_classes(N)
    G = N.Gamma
    return [[sum(Fraction(V[t]) * G[p][t] for t in c) for c in classes] for p in range(N.n)]


def kappa(N, V, blocks=None):
    """kappa_V restricted to unions of `blocks` (default: source classes): blocks c, c' are merged iff
    e_c - e_c' lies in the row space of X(V) = [Gamma_B V_B]_B.  Returns a partition of the reactions."""
    blocks = blocks if blocks is not None else source_classes(N)
    X = A_matrix(N, V, blocks)
    q = len(blocks)
    R, piv = rref(X) if X else ([], [])
    rk = len(piv)
    parent = list(range(q))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in combinations(range(q), 2):
        if find(a) == find(b):
            continue
        e = [Fraction(0)] * q
        e[a], e[b] = Fraction(1), Fraction(-1)
        if rank(R + [e]) == rk:
            parent[find(a)] = find(b)
    groups = {}
    for i in range(q):
        groups.setdefault(find(i), set()).update(blocks[i])
    return normalize_partition(groups.values())


def kernel_X(N, V, blocks):
    X = A_matrix(N, V, blocks)
    return nullspace(X, len(blocks))


def random_flux(N, rng, lo=1, hi=10**6):
    rays = N.flux_rays()
    V = [Fraction(0)] * N.r
    for v in rays:
        w = rng.randint(lo, hi)
        for t in range(N.r):
            V[t] += w * v[t]
    return V


def sub_consistent_flux(N, U, rng):
    """A positive flux of the subnetwork with reactions U (as a vector on all r reactions, zero off U), or None."""
    U = sorted(U)
    if not U:
        return [Fraction(0)] * N.r
    sub = Network(N.n, [N.reactions[t] for t in U])
    rays = sub.flux_rays()
    if not all(any(v[i] > 0 for v in rays) for i in range(len(U))):
        return None
    V = [Fraction(0)] * N.r
    for v in rays:
        w = rng.randint(1, 30)
        for i, t in enumerate(U):
            V[t] += w * v[i]
    return V


# ------------------------------------------------------------------------------------------------
# lower bound for CF: flat closure + two-block species rule (sound: proved in notes.md, Lemma 3.4)
# ------------------------------------------------------------------------------------------------
def _sign_on_cone(N, coeffs):
    """Sign of the linear functional V -> sum_t coeffs[t] V_t on C0: +1 / -1 if definite, 0 if identically zero on
    ker Gamma (rays span it under (H2)), None if indefinite (may vanish on C0)."""
    vals = [sum(Fraction(coeffs[t]) * v[t] for t in range(N.r)) for v in N.flux_rays()]
    if all(x == 0 for x in vals):
        return 0
    if all(x >= 0 for x in vals):
        return 1
    if all(x <= 0 for x in vals):
        return -1
    return None


def _ec_component_merge(N, part):
    """Max-principle rule (Lemma 3.5 of notes.md, generalising agent P's Lemma CF).  Blocks are unions of source
    classes.  Call a species exclusive if all its consumers lie in one block.  Graph on blocks: B -> beta(p) when B
    produces an exclusive species p consumed by beta(p).  A strongly connected set S of >= 2 blocks such that every
    species consumed by a block of S is exclusive and every producer of such a species lies in S is merged.
    Returns a coarser partition or None."""
    Ye, Yp = N.Ye, N.Yp
    blocks = list(part)
    q = len(blocks)
    bidx = {}
    for i, B in enumerate(blocks):
        for t in B:
            bidx[t] = i
    cons = [set() for _ in range(N.n)]   # blocks consuming p
    prod = [set() for _ in range(N.n)]   # blocks producing p
    for p in range(N.n):
        for t in range(N.r):
            if Ye[p][t] > 0:
                cons[p].add(bidx[t])
            if Yp[p][t] > 0:
                prod[p].add(bidx[t])
    consumed_by = [set() for _ in range(q)]
    for p in range(N.n):
        for b in cons[p]:
            consumed_by[b].add(p)
    good = [all(len(cons[p]) == 1 for p in consumed_by[b]) and consumed_by[b] for b in range(q)]
    adj = [set() for _ in range(q)]
    for p in range(N.n):
        if len(cons[p]) == 1:
            (b,) = tuple(cons[p])
            for a in prod[p]:
                adj[a].add(b)
    # strongly connected components among good blocks (Kosaraju)
    nodes = [b for b in range(q) if good[b]]
    order, seen = [], set()
    for s0 in nodes:
        if s0 in seen:
            continue
        stack = [(s0, iter(sorted(x for x in adj[s0] if good[x])))]
        seen.add(s0)
        while stack:
            v, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                order.append(v)
                stack.pop()
            elif nxt not in seen:
                seen.add(nxt)
                stack.append((nxt, iter(sorted(x for x in adj[nxt] if good[x]))))
    radj = [set() for _ in range(q)]
    for a in nodes:
        for b in adj[a]:
            if good[b]:
                radj[b].add(a)
    comp = {}
    for s0 in reversed(order):
        if s0 in comp:
            continue
        stack = [s0]
        comp[s0] = s0
        while stack:
            v = stack.pop()
            for w in radj[v]:
                if w not in comp:
                    comp[w] = s0
                    stack.append(w)
    sccs = {}
    for v_, c in comp.items():
        sccs.setdefault(c, set()).add(v_)
    for S in sccs.values():
        if len(S) < 2:
            continue
        ok = all(prod[p] <= S for b in S for p in consumed_by[b])
        if ok:
            newpart = [blocks[b] for b in range(q) if b not in S] + [set().union(*[blocks[b] for b in S])]
            return normalize_partition(newpart)
    return None


def pi_plus(N, start=None):
    part = flat_closure(N, start)
    G = N.Gamma
    while True:
        merged = False
        np_ = _ec_component_merge(N, part)
        if np_ is not None:
            part = flat_closure(N, np_)
            continue
        for p in range(N.n):
            active = []
            ok = True
            for B in part:
                coeffs = [G[p][t] if t in B else 0 for t in range(N.r)]
                if not any(coeffs):
                    continue
                sg = _sign_on_cone(N, coeffs)
                if sg == 0:
                    continue
                if sg is None:
                    ok = False
                active.append(B)
            if ok and len(active) == 2:
                newpart = [B for B in part if B not in active] + [active[0] | active[1]]
                part = flat_closure(N, normalize_partition(newpart))
                merged = True
                break
        if not merged:
            return part


# ------------------------------------------------------------------------------------------------
# upper bound for CF: meet of kappa_V over generic V and over 'split' special V
# ------------------------------------------------------------------------------------------------
def cf_upper(N, lower, rng, ntrials=3, max_split_blocks=10):
    up = None
    witnesses = []
    for _ in range(ntrials):
        V = random_flux(N, rng)
        K = kappa(N, V)
        up = K if up is None else meet(up, K)
        witnesses.append(("generic", V))
    q = len(lower)
    if 2 <= q <= max_split_blocks:
        blocks = list(lower)
        for k in range(1, q // 2 + 1):
            for S in combinations(range(q), k):
                U = set().union(*[blocks[i] for i in S])
                Uc = set(range(N.r)) - U
                V1 = sub_consistent_flux(N, U, rng)
                if V1 is None:
                    continue
                V2 = sub_consistent_flux(N, Uc, rng)
                if V2 is None:
                    continue
                V = [a + b for a, b in zip(V1, V2)]
                K = kappa(N, V)
                if not refines(up, K):
                    up = meet(up, K)
                    witnesses.append(("split", V))
                else:
                    up = meet(up, K)
    return up, witnesses


def cf_bounds(N, rng=None, ntrials=3):
    """(lower, upper, exact): lower <= CF <= upper in the refinement order; exact iff they coincide."""
    rng = rng or random.Random(12345)
    lo = pi_plus(N)
    up, _ = cf_upper(N, lo, rng, ntrials)
    return lo, up, lo == up


# ------------------------------------------------------------------------------------------------
# harmless operations
# ------------------------------------------------------------------------------------------------
def duplicate(N, j):
    """Add a copy of reaction j (new index r)."""
    return Network(N.n, list(N.reactions) + [N.reactions[j]], names=list(N.names))


def split_choice(N, A):
    """Split off the reactions A (all with the same source y) through a new species W (index n):
    y -> W (new reaction, index r) and W -> y'_t for t in A.  |A| = 1 subdivides a reaction."""
    A = set(A)
    y = N.reactions[min(A)][0]
    assert all(N.reactions[t][0] == y for t in A)
    n1 = N.n + 1
    rx = []
    for t, (a, b) in enumerate(N.reactions):
        a1, b1 = tuple(a) + (0,), tuple(b) + (0,)
        if t in A:
            a1 = (0,) * N.n + (1,)
        rx.append((a1, b1))
    rx.append((tuple(y) + (0,), (0,) * N.n + (1,)))
    return Network(n1, rx, names=list(N.names) + ["W"])


def merge_identical(N):
    """Merge reactions with identical source and product; returns (N', map old index -> new index)."""
    seen = {}
    rx = []
    mp = []
    for a, b in N.reactions:
        if (a, b) not in seen:
            seen[(a, b)] = len(rx)
            rx.append((a, b))
        mp.append(seen[(a, b)])
    return Network(N.n, rx, names=list(N.names)), mp


def pullback(P, mp, r_old):
    """Partition of the old reactions induced by a partition P of the new ones through mp (old -> new)."""
    inv = {}
    for t in range(r_old):
        inv.setdefault(mp[t], []).append(t)
    out = []
    for B in P:
        blk = set()
        for x in B:
            blk.update(inv.get(x, []))
        if blk:
            out.append(blk)
    return normalize_partition(out)


def restrict(P, keep):
    keep = set(keep)
    return normalize_partition([B & keep for B in P if B & keep])


# ------------------------------------------------------------------------------------------------
# exact certificate that CF is one block (sympy; WSL)
# ------------------------------------------------------------------------------------------------
def one_block_certificate(N, blocks):
    """Sufficient test that kappa_V merges all `blocks` for EVERY V in C0: some (q-1)-minor of X(V) without its
    last column, written as a polynomial in the ray weights lambda (V = sum lambda_i e_i, lambda > 0), is nonzero
    with all coefficients of one sign.  Returns (True, minor rows) or (False, None)."""
    import sympy as sp
    rays = N.flux_rays()
    lam = sp.symbols(f"l0:{len(rays)}", positive=True)
    V = [sum(lam[i] * rays[i][t] for i in range(len(rays))) for t in range(N.r)]
    G = N.Gamma
    q = len(blocks)
    if q == 1:
        return True, ()
    X = sp.Matrix([[sum(V[t] * G[p][t] for t in B) for B in blocks[:-1]] for p in range(N.n)])
    for rows in combinations(range(N.n), q - 1):
        d = sp.expand(X.extract(list(rows), list(range(q - 1))).det())
        if d == 0:
            continue
        coeffs = sp.Poly(d, *lam).coeffs()
        if all(c > 0 for c in coeffs) or all(c < 0 for c in coeffs):
            return True, rows
    return False, None
