"""V1 checks of agent Q2's corrected preclusters (round2/preclusters/notes.md) and of E2's Prop. E2.6.

Verifier V1 (Claude Opus 5.5, effort xhigh), 27 October 2026. Exact arithmetic (fractions), standard library only;
independent of crn.py / pclust.py. Symbolic one-block certificates are in check_Q2_sym.py (sympy).
Prints PASS/FAIL per claim; exits 1 on any FAIL of a claim marked proved.
"""
import sys, time, random, platform
from itertools import combinations, product
sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from v1lib import *
from examples import EX
from check_P import family

T0 = time.time()
print("environment: Python", platform.python_version(), platform.system())
QUICK = "--quick" in sys.argv
FAILS, RESULTS = [], []


def claim(name, ok, detail=""):
    RESULTS.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name + ((" -- " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


RNG = random.Random(20261027)


def rand_V(N, k=None):
    R = flux_rays(N)
    return interior_flux(N, [RNG.randint(1, 10 ** 6) for _ in R])


def kmeet(N, Vs):
    P = kappa_V(N, Vs[0])
    for V in Vs[1:]:
        P = meet(P, kappa_V(N, V))
    return P


def dimker(N, V):
    keys, cls = source_classes(N)
    return len(nullspace(A_V(N, V), ncols=len(cls)))


# ------------------------------------------------------------------ hierarchy table and named examples

def hierarchy():
    rows = []
    for k in ("B3", "R", "B2dup", "K3", "PiCF2", "XY4", "Z4", "B2", "N412", "T", "W"):
        N = EX[k]
        P, Pcs, pi = preclusters_def19(N), class_span_closure(N), flat_closure_partition(N)
        Vs = [rand_V(N) for _ in range(3)]
        kg = kmeet(N, Vs)
        rows.append((k, P, Pcs, pi, kg))
        print("   %-6s P=%s Pcs=%s pi=%s kappa(generic, meet of 3)=%s M*=%s" % (k, P, Pcs, pi, kg, full_coupling(N)))
        claim("Prop 2.1 chain P <= Pcs <= pi <= kappa_V (generic) on %s" % k,
              refines(P, Pcs) and refines(Pcs, pi) and refines(pi, kg))
    d = {r[0]: r for r in rows}
    claim("P < Pcs strictly on B3: P two blocks, Pcs one block", len(d["B3"][1]) == 2 and len(d["B3"][2]) == 1)
    claim("Pcs < pi strictly on R: Pcs singletons, pi = M* one block", d["R"][2] == [[0], [1]] and d["R"][3] == [[0, 1]])
    claim("B2dup: pi = {R1-R4},{R5,R6}; kappa_V one block at 3 generic V", d["B2dup"][3] == [[0, 1, 2, 3], [4, 5]]
          and d["B2dup"][4] == [[0, 1, 2, 3, 4, 5]])
    claim("K3: pi = the 3 source classes; kappa_V one block at 3 generic V", len(d["K3"][3]) == 3 and len(d["K3"][4]) == 1)
    claim("PiCF2 (n=2 smallest pi<CF): pi = 2 source classes; kappa_V one block", d["PiCF2"][3] == [[0, 1], [2, 3]] and
          d["PiCF2"][4] == [[0, 1, 2, 3]])
    # XY4: generic one block, special V=(1,1,2,2) two blocks; pi two blocks => CF = pi = two blocks
    N = EX["XY4"]
    Vs = [Fr(1), Fr(1), Fr(2), Fr(2)]
    claim("XY4: V=(1,1,2,2) in C°, A_V = 0, kappa_V = {R1,R2},{R3,R4}; generic kappa one block; pi two blocks "
          "=> CF = pi = two blocks (pi <= CF <= kappa_V(special))",
          all(sum(N.G[i][t] * Vs[t] for t in range(4)) == 0 for i in range(2)) and
          all(x == 0 for row in A_V(N, Vs) for x in row) and kappa_V(N, Vs) == [[0, 1], [2, 3]] and
          d["XY4"][4] == [[0, 1, 2, 3]] and d["XY4"][3] == [[0, 1], [2, 3]])
    N = EX["Z4"]
    claim("Z4: generic kappa = CF = pi = {R1},{R2},{R3,R4}; kappa at V=(1,1,2,2) = {R1,R2},{R3,R4} (coarser)",
          d["Z4"][4] == [[0], [1], [2, 3]] and d["Z4"][3] == [[0], [1], [2, 3]] and
          kappa_V(N, [Fr(1), Fr(1), Fr(2), Fr(2)]) == [[0, 1], [2, 3]])
    # (H1)-(H3) of the witnesses
    for k in ("XY4", "Z4", "PiCF2", "B2dup", "K3"):
        N = EX[k]
        claim("%s satisfies (H1),(H2),(H3)" % k, is_H1(N) and is_H2(N) and is_H3(N))
    claim("K3 is WF-FC, WR, deficiency 0", is_WFFC(EX["K3"]) and is_WR(EX["K3"]) and deficiency(EX["K3"]) == 0)


def prop32():
    N = EX["N412"]; N2 = duplicate(N, 0)
    claim("Prop 3.2(i): P(N4.12) = {R1},{R2}; with R1 duplicated P is one block",
          preclusters_def19(N) == [[0], [1]] and preclusters_def19(N2) == [[0, 1, 2]])
    B3 = EX["B3"]
    S = split_choice(B3, [1, 2])  # X->Y2, X->Y3 behind W
    claim("Prop 3.2(i): P(B3) two blocks, P(B3 with X->Y2,X->Y3 split behind W) one block",
          len(preclusters_def19(B3)) == 2 and len(preclusters_def19(S)) == 1, "split net: %r" % S)
    B2, B2d = EX["B2"], duplicate(EX["B2"], 4)
    claim("Prop 3.2(ii): Pcs and pi of B2 one block, of B2dup two blocks",
          len(class_span_closure(B2)) == 1 and len(flat_closure_partition(B2)) == 1 and
          len(class_span_closure(B2d)) == 2 and len(flat_closure_partition(B2d)) == 2)
    Sd = split_choice(B2d, [0, 1])
    claim("Prop 3.2(ii): splitting {X->Y1,X->Y2} of B2dup behind W makes Pcs and pi one block",
          len(class_span_closure(Sd)) == 1 and len(flat_closure_partition(Sd)) == 1, "%r" % Sd)
    # CF of all these is one block (EC + connected + H1 + H2, Theorem 4.1); kappa at random V
    for X in (B2, B2d, Sd, B3, S):
        claim("kappa_V one block at 2 random V for %r (EC=%s)" % (X, is_EC(X)),
              len(kmeet(X, [rand_V(X), rand_V(X)])) == 1)


# ------------------------------------------------------------------ Theorem 3.1 invariance (pointwise in V)

def lift_dup(N, t, Vp):
    """image of V' in C°(N): v_t = v'_t + v'_t'"""
    V = list(Vp[:N.r])
    V[t] = V[t] + Vp[N.r]
    return V


def check_invariance(nets, label):
    bad = []
    cnt = 0
    for N in nets:
        if not is_H2(N):
            continue
        # duplication of each reaction
        for t in range(N.r):
            Np = duplicate(N, t)
            for _ in range(2):
                Vp = rand_V(Np)
                V = lift_dup(N, t, Vp)
                K, Kp = kappa_V(N, V), kappa_V(Np, Vp)
                # phi^{-1}(K): t' (index N.r) goes with t
                lift = canon([b + ([N.r] if t in b else []) for b in K])
                cnt += 1
                if lift != Kp:
                    bad.append(("dup", N, t))
        # splitting: every source class c and every nonempty A subset of c
        keys, cls = source_classes(N)
        for c in cls:
            for a in range(1, len(c) + 1):
                for A in combinations(c, a):
                    Np = split_choice(N, A)
                    if not is_H2(Np):
                        bad.append(("split H2", N, A)); continue
                    for _ in range(2):
                        Vp = rand_V(Np)
                        V = Vp[:N.r]  # old reactions keep their fluxes
                        K, Kp = kappa_V(N, V), kappa_V(Np, Vp)
                        tA = N.r
                        lift = canon([b + ([tA] if A[0] in b else []) for b in K])
                        cnt += 1
                        if lift != Kp:
                            bad.append(("split", N, A))
    claim("Thm 3.1 (O2),(O3) pointwise: kappa'_{V'} = phi^{-1}(kappa_V) on %s (%d checks)" % (label, cnt), not bad,
          repr(bad[:2]))


def check_prop33(nets, label):
    """Prop 3.3 mechanism: if s, t fully coupled in N (H2), then after duplicating t, s and t lie in one block of P."""
    bad, cnt = [], 0
    for N in nets:
        if not is_H2(N):
            continue
        FC = full_coupling(N)
        for b in FC:
            for s, t in combinations(b, 2):
                Np = duplicate(N, t)
                P = preclusters_def19(Np)
                cnt += 1
                if not any(s in x and t in x for x in P):
                    bad.append((N, s, t))
    claim("Prop 3.3 mechanism on %s (%d fully coupled pairs)" % (label, cnt), not bad, repr(bad[:1]))


# ------------------------------------------------------------------ Theorems 4.1, 4.2, Prop 4.3, Thm 6.1, Cor 6.2-6.3

def family_tests(nets, label):
    t0 = time.time()
    c = {"H2": 0, "H1": 0, "eqiff": 0}
    bad = {}
    for N in nets:
        if is_H2(N):
            c["H2"] += 1
            Ms, pi = full_coupling(N), flat_closure_partition(N)
            Vs = [rand_V(N) for _ in range(2)]
            ks = [kappa_V(N, V) for V in Vs]
            if not (refines(Ms, pi) and all(refines(pi, k) for k in ks)):
                bad.setdefault("6.1 chain M* <= pi <= kappa_V", N)
            dp = doubling_pairs(N)
            nv = nvec(N)
            allc = all(positively_proportional(nv[i], nv[j]) for i, j in dp)
            for P_ in [pi] + ks:
                if (P_ == Ms) != allc or (len(P_) == len(Ms)) != allc or refines(P_, Ms) != allc:
                    bad.setdefault("6.1 equivalences", N)
            if not refines(preclusters_def19(N), class_span_closure(N)) or not refines(class_span_closure(N), pi):
                bad.setdefault("2.1 P <= Pcs <= pi", N)
            # Prop 2.2: pi = Pcs v M*
            if pi != join(class_span_closure(N), Ms, N.r):
                bad.setdefault("2.2 pi = Pcs v M*", N)
            # Prop 5.1 at k = V
            for V, k in zip(Vs, ks):
                if dimker(N, V) > len([b for b in k]):
                    bad.setdefault("5.1 dim ker <= #kappa", N)
            # Cor 5.3 with pi: internal cycles inside pi-blocks => dim ker A_V = #pi = #kappa_V
            cyc = [set(S) for S, _ in flux_rays(N)]
            if all(any(cc <= set(b) for b in pi) for cc in cyc):
                for V, k in zip(Vs, ks):
                    if not (dimker(N, V) == len(pi) == len(k)):
                        bad.setdefault("5.3 isolation => dim ker = #blocks", N)
            # Theorem 4.1
            if is_EC(N) and is_connected(N) and is_H1(N):
                c["EC"] = c.get("EC", 0) + 1
                for V in Vs:
                    K = nullspace(A_V(N, V), ncols=n_sources(N))
                    if len(K) != 1 or len(set(K[0])) != 1:
                        bad.setdefault("4.1 EC => ker A_V = R1", N)
            # E2.6: each linkage class consistent on its own => kappa_V at V = sum of class fluxes refines linkage
            lc, m, l = linkage_classes(N)
            okL = True
            Vsum = [Fr(0)] * N.r
            for L in lc:
                sub = Net([[N.Ye[i][t] for t in L] for i in range(N.n)], [[N.Yp[i][t] for t in L] for i in range(N.n)])
                if not is_H2(sub):
                    okL = False; break
                VL = interior_flux(sub, [RNG.randint(1, 1000) for _ in flux_rays(sub)])
                for t, v in zip(L, VL):
                    Vsum[t] = v
            if okL:
                c["E26"] = c.get("E26", 0) + 1
                if not refines(kappa_V(N, Vsum), lc):
                    bad.setdefault("E2.6", N)
        if is_H1(N):
            c["H1"] += 1
            D = N.dual()
            M = ray_mips(N)
            piD = flat_closure_partition(D)
            kD = kappa_V(D, rand_V(D))
            tw = twin_pairs(N)
            w = wvec(N)
            allc = all(positively_proportional(w[i], w[j]) for i, j in tw)
            for P_ in (piD, kD):
                if not refines(M, P_) or (P_ == M) != allc:
                    bad.setdefault("6.2 species side", N)
            # Prop 6.4(b) for pi(N*) and kappa_V(N*)
            ms = min_siphons(N)
            rhs = (not tw) and all(any(set(m_) <= S for S in ms) for m_ in M)
            for P_ in (piD, kD):
                lhs = all(any(set(b) <= S for S in ms) for b in P_)
                if lhs != rhs:
                    bad.setdefault("6.4(b)", N)
            # 6.2 last clause
            prods_are_sources = all(N.prod(j) in {N.src(i) for i in range(N.r)} for j in range(N.r))
            if prods_are_sources and is_H3(N) and (tw or piD != M):
                bad.setdefault("6.2 product complexes sources + H3 => no twins, P^(N*) = M(N)", N)
    for key in ("6.1 chain M* <= pi <= kappa_V", "6.1 equivalences", "2.1 P <= Pcs <= pi", "2.2 pi = Pcs v M*",
                "5.1 dim ker <= #kappa", "5.3 isolation => dim ker = #blocks", "4.1 EC => ker A_V = R1", "E2.6",
                "6.2 species side", "6.4(b)", "6.2 product complexes sources + H3 => no twins, P^(N*) = M(N)"):
        claim("%s on %s" % (key, label), key not in bad, ("counterexample %r" % bad[key]) if key in bad else "")
    print("   %s: (H2) %d, (H1) %d, EC-connected-(H1)(H2) %d, E2.6-applicable %d; %.1f s" %
          (label, c["H2"], c["H1"], c.get("EC", 0), c.get("E26", 0), time.time() - t0))


# ------------------------------------------------------------------ random WR networks (Thm 4.2, Cor 6.3, E2.6)

def random_wr(rng, n, mol):
    cpx_all = [y for y in product(range(mol + 1), repeat=n) if 0 < sum(y) <= mol]
    l = rng.randint(1, 2)
    cpx = rng.sample(cpx_all, min(len(cpx_all), rng.randint(2, 5)))
    # split complexes into l classes
    rng.shuffle(cpx)
    if l == 2 and len(cpx) >= 4:
        k = rng.randint(2, len(cpx) - 2)
        classes = [cpx[:k], cpx[k:]]
    else:
        classes = [cpx]
    edges = set()
    for C in classes:
        order = C[:]
        rng.shuffle(order)
        for a, b in zip(order, order[1:] + order[:1]):
            edges.add((a, b))
        for _ in range(rng.randint(0, 3)):
            a, b = rng.sample(C, 2)
            edges.add((a, b))
    edges = sorted(edges)
    Ye = [[e[0][i] for e in edges] for i in range(n)]
    Yp = [[e[1][i] for e in edges] for i in range(n)]
    if any(not any(Ye[i]) for i in range(n)):
        return None
    return Net(Ye, Yp)


def wr_tests():
    rng = random.Random(4242)
    cnt = {"WR": 0, "d0": 0}
    bad = {}
    tries = 600 if QUICK else 4000
    for _ in range(tries):
        N = random_wr(rng, rng.randint(2, 4), 2)
        if N is None or not is_WR(N):
            continue
        cnt["WR"] += 1
        V = rand_V(N)
        lc, m, l = linkage_classes(N)
        k = kappa_V(N, V)
        if not refines(k, lc):
            bad.setdefault("E2.6 WR: kappa_V refines linkage classes", N)
        if deficiency(N) == 0:
            cnt["d0"] += 1
            K = nullspace(A_V(N, V), ncols=n_sources(N))
            if len(K) != l or k != lc:
                bad.setdefault("4.2 WR d0: ker A_V = span 1_L, kappa_V = linkage classes", N)
        # Cor 6.3
        twoout = len({N.src(j) for j in range(N.r)}) < N.r
        if (k == full_coupling(N)) == twoout:
            bad.setdefault("6.3 WR: kappa_V = M* iff no complex with two outgoing reactions", N)
    for key in ("E2.6 WR: kappa_V refines linkage classes", "4.2 WR d0: ker A_V = span 1_L, kappa_V = linkage classes",
                "6.3 WR: kappa_V = M* iff no complex with two outgoing reactions"):
        claim("%s (%d random WR networks, %d with deficiency 0)" % (key, cnt["WR"], cnt["d0"]), key not in bad,
              ("counterexample %r" % bad[key]) if key in bad else "")


if __name__ == "__main__":
    hierarchy()
    prop32()
    if QUICK:
        nets = list(family(2, 2, 3))
        check_invariance([N for N in nets if is_H2(N)][:60], "first 60 (H2) networks n=2,mol<=2,r<=3")
        check_prop33(nets, "n=2,mol<=2,r<=3")
        family_tests(nets, "n=2,mol<=2,r<=3 (zero complex allowed)")
    else:
        nets = list(family(2, 2, 4))
        check_invariance([N for N in nets if is_H2(N)], "all (H2) networks n=2,mol<=2,r<=4")
        check_prop33(nets, "n=2,mol<=2,r<=4")
        family_tests(nets, "n=2,mol<=2,r<=4 (zero complex allowed)")
        nets3 = [N for N in family(3, 2, 3, zero=False)]
        check_invariance([N for N in nets3 if is_H2(N)][:400], "first 400 (H2) networks n=3,mol<=2,r<=3")
        family_tests(nets3, "n=3,mol<=2,r<=3")
        family_tests(list(family(2, 3, 3, zero=False)), "n=2,mol<=3,r<=3")
    wr_tests()
    print("runtime %.1f s; %d claims, %d FAIL" % (time.time() - T0, len(RESULTS), len(FAILS)))
    if FAILS:
        print("FAILED:", FAILS)
        sys.exit(1)
