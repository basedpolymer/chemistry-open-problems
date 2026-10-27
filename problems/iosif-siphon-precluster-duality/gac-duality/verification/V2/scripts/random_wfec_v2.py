"""Verifier V2: independent exact checks of Theorems T, C, C', U(b), U(c), U*(a) and of the G-matrix identities used in
Theorem GS, on randomly generated well-formed equal-conflict (WF-EC) networks; negative controls on EC networks with
s >= m; and Theorems C^gen / K5gen / PC(c) on random cluster-rank networks (PC and non-PC).

Generator (independent of the authors' enumerations): species split into blocks; sources on the blocks with random
coefficients 1..3 (unbalanced allowed); a random nonnegative integer matrix Z with n-q+1 rows and every column
positive; products are random nonnegative vectors y' (entries <= 3) with Z y' = Z y (so every reaction vector lies in
ker Z, rank <= q-1); networks kept if connected, (H1), (H2) (then s = q-1 by K5gen).  Self-loops (impure reactions)
are allowed.

Environment: WSL Fedora 44, Python 3.12, standard library only (fractions).  Usage: python random_wfec_v2.py SEED COUNT
PASS/FAIL per claim, exit 1 on FAIL.
"""
import sys, os, time, random
from itertools import product
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from v2lib import *  # noqa

T0 = time.time()
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
COUNT = int(sys.argv[2]) if len(sys.argv) > 2 else 200
rng = random.Random(seed)
stats = {}
fails = {}


def tally(name, ok):
    s = stats.setdefault(name, [0, 0])
    s[0] += 1
    if ok:
        s[1] += 1
    else:
        fails.setdefault(name, 0)
        fails[name] += 1


def blocks(n, q):
    while True:
        lab = [rng.randrange(q) for _ in range(n)]
        if len(set(lab)) == q:
            return [[p for p in range(n) if lab[p] == c] for c in range(q)]


def gen(n, q, kind):
    """kind: 'EC' (one source per block), 'PC' (sources a*b on a block), 'GEN' (sharing sources, not proportional)."""
    B = blocks(n, q)
    srcs_by_cluster = []
    for blk in B:
        if kind == "EC":
            y = [0] * n
            for p in blk:
                y[p] = rng.choice([1, 1, 2, 3])
            srcs_by_cluster.append([tuple(y)])
        elif kind == "PC":
            b = [0] * n
            for p in blk:
                b[p] = rng.choice([1, 1, 2])
            mult = sorted(set(rng.choice([1, 2]) for _ in range(rng.randint(1, 2))))
            srcs_by_cluster.append([tuple(a * x for x in b) for a in mult])
        else:
            # a few sources on the block whose sharing graph is connected and whose union is the block
            ys = []
            for _ in range(rng.randint(1, 3)):
                y = [0] * n
                for p in blk:
                    y[p] = rng.choice([0, 1, 1, 2])
                if not any(y):
                    y[rng.choice(blk)] = 1
                ys.append(tuple(y))
            full = [0] * n
            for p in blk:
                full[p] = rng.choice([1, 2])
            ys.append(tuple(full))  # shares a species with every other source of the block
            srcs_by_cluster.append(sorted(set(ys)))
    rows = n - q + 1
    while True:
        Z = [[rng.choice([0, 0, 1, 1, 2]) for _ in range(n)] for _ in range(rows)]
        if all(any(Z[i][p] for i in range(rows)) for p in range(n)) and rank(Z) == rows:
            break
    cand = list(product(range(4), repeat=n))

    def zval(y):
        return tuple(sum(Z[i][p] * y[p] for p in range(n)) for i in range(rows))
    N = []
    for ys in srcs_by_cluster:
        for y in ys:
            zy = zval(y)
            prods = [c for c in cand if c != y and zval(c) == zy]
            if not prods:
                return None
            for _ in range(rng.randint(1, 2)):
                N.append((y, rng.choice(prods)))
    N = sorted(set(N))
    if not connected(N) or not conservative(N) or not consistent(N):
        return None
    return N


def check_wfec(N):
    n, r = nsp(N), len(N)
    s = rank(Gamma(N))
    srcs, ct = ec_structure(N)
    m = len(srcs)
    tally("generator: WF-EC (EC, connected, H1, H2, s = m-1)", is_WFEC(N))
    # Theorem T: ker A_k = R u(k), u(k) > 0
    for _ in range(2):
        k = rand_pos(rng, r)
        tally("T: ker A_k one-dimensional and positive (random k)", u_of_k(N, k) is not None)
    # Theorem C and C'
    tally("C: N has no critical siphon", not critical_siphons(N))
    tally("C': N* has no critical siphon", not critical_siphons(dual(N)))
    # U(b) at random non-steady points, and G identities
    for _ in range(2):
        k = rand_pos(rng, r)
        x = rand_pos(rng, n)
        u = u_of_k(N, k)
        G = G_matrix(N, k, x)
        metz = all(G[c][d] >= 0 for c in range(m) for d in range(m) if c != d)
        # irreducible: G[c][d] > 0 iff d -> c ; check strong connectivity of that digraph
        adj = {c: [d for d in range(m) if d != c and G[c][d] > 0] for c in range(m)}
        reach = {0}
        stack = [0]
        while stack:
            c = stack.pop()
            for d in adj[c]:
                if d not in reach:
                    reach.add(d)
                    stack.append(d)
        radj = {c: [d for d in range(m) if d != c and G[d][c] > 0] for c in range(m)}
        reach2 = {0}
        stack = [0]
        while stack:
            c = stack.pop()
            for d in radj[c]:
                if d not in reach2:
                    reach2.add(d)
                    stack.append(d)
        irr = len(reach) == m and len(reach2) == m
        Gu = [sum(G[c][d] * u[d] for d in range(m)) for c in range(m)]
        tally("GS/U(b): G(x) Metzler, irreducible, G(x)u(k) = 0 at random (x,k)", metz and irr and all(v == 0 for v in Gu))
        nonsteady = any(v != 0 for v in field(N, k, x))
        tally("U(b): det(tI-J) = t^(n-s) q(t), q Hurwitz, at random non-steady x", hurwitz_on_S(N, k, x)
              if nonsteady else True)
    # U(c): sign condition via maximal minors, S = im Gamma vs D = span{y_c - y_d}
    if s > 0:
        BS = row_basis(transpose(Gamma(N)))
        Dv = [[F(srcs[c][p] - srcs[0][p]) for p in range(n)] for c in range(1, m)]
        BD = row_basis(Dv)
        ok = len(BD) == s and maxminor_sign_condition(BS, BD)
        tally("U(c): sigma(S) cap sigma(D^perp) = {0} (maximal minors)", ok)
        tally("info: S != D (Birch does not apply)", rank(BS + BD) > s)
    tally("info: impure", any(a[i] and b[i] for a, b in N for i in range(n)))
    tally("info: unbalanced joins", not is_balanced(N))
    # U*(a): dual Jacobian Hurwitz on S* at random positive k with random kappa
    D = dual(N)
    for _ in range(2):
        kap = rand_pos(rng, n)
        kk = rand_pos(rng, r)
        cp = charpoly(jacobian(D, kap, kk))
        q = strip_zero_roots(cp, r - s)
        tally("U*(a): det(tI-J*) = t^(r-s) q(t), q Hurwitz, at random k", q is not None and is_hurwitz(q))


def check_negative(N):
    """EC, connected, (H1), (H2), s >= m: at a positive steady state (rates k_t = v_t / x^{y_c(t)}, v > 0 in ker
    Gamma) the Hurwitz test must fail (0 has multiplicity >= n-m+1 > n-s).  Sanity check of the test itself."""
    n, r = nsp(N), len(N)
    x = rand_pos(rng, n)
    rays = flux_rays(N)
    wts = [F(rng.randint(1, 5)) for _ in rays]
    v = [sum(w * ray[t] for w, ray in zip(wts, rays)) for t in range(r)]
    k = [v[t] / monomial(x, N[t][0]) for t in range(r)]
    assert all(e == 0 for e in field(N, k, x))
    tally("control: EC with s >= m -> Hurwitz test fails at a positive steady state", not hurwitz_on_S(N, k, x))


def check_cluster(N, kind):
    n, r = nsp(N), len(N)
    s = rank(Gamma(N))
    q = len(clusters(N))
    tally(f"K5gen [{kind}]: s >= q - 1", s >= q - 1)
    if s != q - 1:
        return
    tally(f"Cgen [{kind}]: N has no critical siphon", not critical_siphons(N))
    tally(f"Cgen [{kind}]: N* has no critical siphon", not critical_siphons(dual(N)))
    if kind == "PC":
        tally("PC: N is PC and N* is PC with the same cluster count and rank", is_PC(N) and is_PC(dual(N))
              and len(clusters(dual(N))) == q and WF_cluster(dual(N)))
        for _ in range(2):
            k = rand_pos(rng, r)
            x = rand_pos(rng, n)
            tally("PC(c): J|_S Hurwitz at random point (N)", hurwitz_on_S(N, k, x))
            kap = rand_pos(rng, n)
            kk = rand_pos(rng, r)
            cp = charpoly(jacobian(dual(N), kap, kk))
            qq = strip_zero_roots(cp, r - s)
            tally("PC(c): J*|_S* Hurwitz at random point (N*)", qq is not None and is_hurwitz(qq))
        bh = all(len(set(x for x in a if x)) == 1 for a, _ in N)
        tally("info: PC network neither EC nor BH", not is_EC(N) and not bh)
    else:
        tally("info[GEN]: rank-equation network that is not PC", not is_PC(N))
        if not is_PC(N):
            for _ in range(2):
                k = rand_pos(rng, r)
                x = rand_pos(rng, n)
                hw = hurwitz_on_S(N, k, x)
                tally("CR probe [GEN, not PC]: J|_S Hurwitz at random point", hw)
                if not hw:
                    print("CR-probe non-Hurwitz point: N =", N, "k =", [str(t) for t in k], "x =", [str(t) for t in x])


done = {"EC": 0, "PC": 0, "GEN": 0, "NEG": 0}
tries = 0
seen = set()
while min(done["EC"], done["PC"], done["GEN"]) < COUNT and tries < 400 * COUNT:
    tries += 1
    kind = rng.choice(["EC", "PC", "GEN"])
    if done[kind] >= COUNT:
        continue
    n = rng.randint(3, 6)
    q = rng.randint(2, min(4, n))
    N = gen(n, q, kind)
    if N is None or tuple(N) in seen:
        continue
    seen.add(tuple(N))
    if kind == "EC":
        if is_WFEC(N):
            check_wfec(N)
            done["EC"] += 1
    else:
        check_cluster(N, kind)
        done[kind] += 1
# negative controls: EC networks with s >= m (random products, no Z constraint)
neg = 0
while neg < COUNT // 4 and tries < 800 * COUNT:
    tries += 1
    n = rng.randint(2, 4)
    q = rng.randint(1, min(2, n))
    B = blocks(n, q)
    N = []
    for blk in B:
        y = [0] * n
        for p in blk:
            y[p] = rng.choice([1, 2])
        for _ in range(rng.randint(1, 2)):
            N.append((tuple(y), tuple(rng.randint(0, 2) for _ in range(n))))
    N = sorted(set(N))
    if is_EC(N) and connected(N) and conservative(N) and consistent(N) and rank(Gamma(N)) >= len({a for a, _ in N}):
        check_negative(N)
        neg += 1
# unconstrained random networks: K5gen (s >= q-1 under connected, H1, H2) and C^gen when s = q-1
gen_ok = 0
gen_tries = 0
while gen_ok < 2 * COUNT and gen_tries < 2000 * COUNT:
    gen_tries += 1
    n = rng.randint(2, 4)
    r = rng.randint(2, 5)
    N = []
    for _ in range(r):
        a = tuple(rng.choice([0, 0, 1, 1, 2]) for _ in range(n))
        b = tuple(rng.choice([0, 0, 1, 1, 2]) for _ in range(n))
        if any(a) and a != b and sum(a) <= 3 and sum(b) <= 3:
            N.append((a, b))
    N = sorted(set(N))
    if len(N) < 2 or not connected(N) or not conservative(N) or not consistent(N):
        continue
    gen_ok += 1
    s_ = rank(Gamma(N))
    q_ = len(clusters(N))
    tally("K5gen [unconstrained random]: s >= q - 1", s_ >= q_ - 1)
    if s_ == q_ - 1:
        tally("Cgen [unconstrained random]: no critical siphon in N and N*",
              not critical_siphons(N) and not critical_siphons(dual(N)))
        tally("info[unconstrained]: cluster-rank network not EC", not is_EC(N))
print(f"seed {seed}: generated EC {done['EC']}, PC {done['PC']}, GEN {done['GEN']}, negative controls {neg}, tries {tries}")
bad = False
for name, (tot, ok) in sorted(stats.items()):
    tag = "INFO" if name.startswith("info") or name.startswith("CR probe") else ("PASS" if ok == tot else "FAIL")
    if tag == "FAIL":
        bad = True
    print(f"{tag} {name}: {ok}/{tot}")
print(f"runtime {time.time() - T0:.1f} s")
sys.exit(1 if bad else 0)
