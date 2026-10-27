"""V1 checks of agent P's Petri-net notes (petri/notes.md, all sections except 7).

Verifier V1 (Claude Opus 5.5, effort xhigh), 27 October 2026. Exact arithmetic (fractions), standard library only;
independent of crn.py / petri.py. Prints PASS/FAIL per claim; exits 1 on any FAIL of a claim marked proved.
Usage: python check_P.py [--quick]
"""
import sys, time, random, platform
from itertools import combinations, product
sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from v1lib import *
from examples import EX, K8

T0 = time.time()
print("environment: Python", platform.python_version(), platform.system())
QUICK = "--quick" in sys.argv
FAILS = []
RESULTS = []


def claim(name, ok, detail=""):
    RESULTS.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name + ((" -- " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


# ------------------------------------------------------------------ helpers specific to P's claims

def consumers(N, i):
    return frozenset(j for j in range(N.r) if N.Ye[i][j] > 0)


def reactant_set(N, j):
    return frozenset(i for i in range(N.n) if N.Ye[i][j] > 0)


def thmA(N):
    """Theorem A: no minimal siphon contains i != j with equal consumer sets; no minimal T-trap contains two
    reactions with equal reactant sets."""
    for R in min_siphons(N):
        for i, j in combinations(sorted(R), 2):
            if consumers(N, i) == consumers(N, j):
                return False
    mt = minimal(t_traps(N))
    for U in mt:
        for a, b in combinations(sorted(U), 2):
            if reactant_set(N, a) == reactant_set(N, b):
                return False
    return True


def thmB_species(N):
    """Theorem B, species side: [every precluster of N* lies in a minimal siphon of N] iff [no twin pair (equal
    rows of Ye, zero rows included) and every species lies in a minimal siphon]."""
    D = N.dual()
    P = preclusters_def19(D)
    ms = min_siphons(N)
    lhs = all(any(set(b) <= S for S in ms) for b in P)
    rhs = (len(twin_pairs(N)) == 0) and all(any(i in S for S in ms) for i in range(N.n))
    return lhs == rhs


def propH(N):
    """Prop H: (H1) => every precluster of N* with >= 2 elements is a union of ray-MIPS of N."""
    D = N.dual()
    P = preclusters_def19(D)
    M = ray_mips(N)
    for b in P:
        if len(b) < 2:
            continue
        for m in M:
            if set(m) & set(b) and not set(m) <= set(b):
                return False
    return True


def propG(N):
    """Prop G: under (H1), (iii) min siphons = supports of extreme laws implies (i) every min siphon is a union of
    MIPS and (ii) every MIPS lies in a min siphon. Returns (iii, i, ii, iv)."""
    ms = set(min_siphons(N))
    sups = {frozenset(S) for S, _ in cons_rays(N)}
    M = ray_mips(N)
    iii = ms == sups
    i_ = all(all(set(m) <= S or not (set(m) & S) for m in M) for S in ms)
    ii = all(any(set(m) <= S for S in ms) for m in M)
    iv = all(Z == frozenset().union(*[S for S in ms if S <= Z]) if any(S <= Z for S in ms) else False for Z in sups)
    return iii, i_, ii, iv


def dictionary(N):
    """Prop 1.1/1.2 entries that can be checked mechanically."""
    D = N.dual()
    out = {}
    # 1.1 weighted reverse dual: build Petri arcs and reverse them
    arcs = {}
    for i in range(N.n):
        for j in range(N.r):
            if N.Ye[i][j]:
                arcs[("p", i, "t", j)] = N.Ye[i][j]
            if N.Yp[i][j]:
                arcs[("t", j, "p", i)] = N.Yp[i][j]
    rd = {}
    for (a, x, b, y), w in arcs.items():  # reverse arc and swap kinds: place<->transition
        swap = {"p": "t", "t": "p"}
        rd[(swap[b], y, swap[a], x)] = w
    arcsD = {}
    for i in range(D.n):
        for j in range(D.r):
            if D.Ye[i][j]:
                arcsD[("p", i, "t", j)] = D.Ye[i][j]
            if D.Yp[i][j]:
                arcsD[("t", j, "p", i)] = D.Yp[i][j]
    out["1.1 reverse-dual"] = rd == arcsD
    out["siphons=T-traps(N*)"] = set(siphons(N)) == set(t_traps(D))
    out["traps=T-siphons(N*)"] = set(traps(N)) == set(t_siphons(D))
    cN = {(S, T) for S, T in clusters(N)}
    cD = {(T, S) for S, T in clusters(D)}
    out["clusters swap"] = cN == cD
    out["FC self-dual"] = (is_ordinary(N) and is_EC(N)) == (is_ordinary(D) and is_EC(D))
    out["EC+balanced self-dual"] = (is_EC(N) and balanced_joins(N)) == (is_EC(D) and balanced_joins(D))
    out["connected"] = is_connected(N) == is_connected(D)
    out["H1<->H2"] = is_H1(N) == is_H2(D) and is_H2(N) == is_H1(D)
    out["H3<->cols"] = is_H3(D) == (len({tuple(c) for c in transpose(N.G)}) == N.r)
    out["rank"] = srank(N) == srank(D)
    out["cons supports = internal cycles of N*"] = ({frozenset(S) for S, _ in cons_rays(N)} ==
                                                   {frozenset(S) for S, _ in flux_rays(D)})
    if is_H1(N):
        out["MIPS(N)=fullcoupling(N*)"] = ray_mips(N) == full_coupling(D)
    # twin classes of N = doubling classes of N*
    out["twins=doubling(N*)"] = set(twin_pairs(N)) == set(doubling_pairs(D))
    return out


def prop22(N):
    """Prop 2.2: if every reaction has nonzero source: doubling classes refine T_c, equality iff EC.
    If every species is consumed: twin classes refine S_c, equality iff N* EC."""
    res = []
    cl = clusters(N)
    if all(any(N.src(j)) for j in range(N.r)):
        dc = {}
        for j in range(N.r):
            dc.setdefault(N.src(j), set()).add(j)
        dcl = [frozenset(v) for v in dc.values()]
        Tc = [T for S, T in cl if T]
        ref = all(any(d <= T for T in Tc) for d in dcl)
        eq = set(dcl) == set(Tc)
        res.append(ref and (eq == is_EC(N)))
    if all(any(N.Ye[i]) for i in range(N.n)):
        tc = {}
        for i in range(N.n):
            tc.setdefault(tuple(N.Ye[i]), set()).add(i)
        tcl = [frozenset(v) for v in tc.values()]
        Sc = [S for S, T in cl if S]
        ref = all(any(d <= S for S in Sc) for d in tcl)
        eq = set(tcl) == set(Sc)
        res.append(ref and (eq == is_EC(N.dual())))
    return all(res)


# ------------------------------------------------------------------ enumeration

def complexes(n, mol, zero=True):
    out = []
    for y in product(range(mol + 1), repeat=n):
        if sum(y) <= mol and (zero or sum(y) > 0):
            out.append(y)
    return out


def family(n, mol, rmax, zero=True, require_all_species=True):
    C = complexes(n, mol, zero)
    reacs = [(a, b) for a in C for b in C if a != b]
    for r in range(1, rmax + 1):
        for S in combinations(reacs, r):
            Ye = [[S[j][0][i] for j in range(r)] for i in range(n)]
            Yp = [[S[j][1][i] for j in range(r)] for i in range(n)]
            if require_all_species and any(not any(Ye[i]) and not any(Yp[i]) for i in range(n)):
                continue
            yield Net(Ye, Yp)


def run_family(label, gen, limit=None):
    t = time.time()
    cnt = {"nets": 0, "A": 0, "B": 0, "dict": 0, "22": 0, "H": 0, "Happl": 0, "G": 0, "Gappl": 0, "iv_fail": 0,
           "Hred": 0}
    bad = {}
    for k, N in enumerate(gen):
        if limit and k >= limit:
            break
        cnt["nets"] += 1
        D = N.dual()
        for side in (N, D):
            if not thmA(side):
                bad.setdefault("A", side)
            if not thmB_species(side):
                bad.setdefault("B", side)
            if not prop22(side):
                bad.setdefault("2.2", side)
        dct = dictionary(N)
        if not all(dct.values()):
            bad.setdefault("dict " + ",".join(k for k, v in dct.items() if not v), N)
        if is_H1(N):
            cnt["Happl"] += 1
            if not propH(N):
                bad.setdefault("H", N)
            iii, i_, ii, iv = propG(N)
            cnt["Gappl"] += 1
            if iii and not (i_ and ii):
                bad.setdefault("G", N)
            if not iv:
                cnt["iv_fail"] += 1
        if is_H2(N):
            # reaction side of H: every precluster of N with >= 2 elements is a union of fully coupled classes
            P = preclusters_def19(N)
            M = full_coupling(N)
            okH = all(not (len(b) >= 2 and any(set(m) & set(b) and not set(m) <= set(b) for m in M)) for b in P)
            cnt["Hred"] += 1
            if not okH:
                bad.setdefault("H reaction side", N)
    for k in ("A", "B", "2.2", "H", "G", "H reaction side"):
        claim("%s family %s (%d networks, both sides)" % (k, label, cnt["nets"]), k not in bad,
              ("counterexample: %r" % bad[k]) if k in bad else "")
    dkeys = [k for k in bad if k.startswith("dict")]
    claim("Prop 1.1/1.2 dictionary, family %s" % label, not dkeys, "; ".join("%s: %r" % (k, bad[k]) for k in dkeys))
    print("   family %s: %d networks, (H1) %d, (H2) %d, Prop G(iv) fails on %d (H1)-networks; %.1f s" %
          (label, cnt["nets"], cnt["Happl"], cnt["Hred"], cnt["iv_fail"], time.time() - t))


# ------------------------------------------------------------------ examples

def examples():
    # Remark 2.5
    for k, want in ((2, 1), (3, 2), (4, 2)):
        N = EX["B%d" % k]
        P19, Pth = preclusters_def19(N), preclusters_thesis(N)
        Pcs = class_span_closure(N)
        claim("Rem 2.5: P(B%d) has %d block(s) (Def 19), thesis Alg 4.4 agrees; class-span closure one block" % (k, want),
              len(P19) == want and P19 == Pth and len(Pcs) == 1, "P=%s Pcs=%s" % (P19, Pcs))
        claim("Rem 2.5: B%d is WF-FC with (H1)-(H3)" % k, is_WFFC(N) and is_H3(N))
    # B3 explicit blocks: {all but Z->X}, {Z->X}
    claim("Rem 2.5: P(B3) = {X->Ya, Ya->Z}, {Z->X}", preclusters_def19(EX["B3"]) == [[0, 1, 2, 3, 4, 5], [6]])

    # Theorem E (i): R, T, W, W*, F are WF-FC with (H1)-(H3) and WF-FC duals
    for k in ("R", "T", "W", "Wstar", "F"):
        N = EX[k]
        claim("Thm E(i): %s WF-FC, (H1)-(H3), dual WF-FC with (H3)" % k,
              is_WFFC(N) and is_H3(N) and is_WFFC(N.dual()) and is_H3(N.dual()))
    # Table 1 formulations on the five networks
    table = {}
    for k in ("R", "T", "W", "Wstar", "F"):
        N = EX[k]; D = N.dual()
        PD, M = preclusters_def19(D), ray_mips(N)
        PN, MD = preclusters_def19(N), full_coupling(N)
        S, SD = [sorted(x) for x in min_siphons(N)], [sorted(x) for x in min_siphons(D)]
        f = {}
        f["A1"] = PD == M; f["A2"] = refines(PD, M); f["A3"] = refines(M, PD); f["A4"] = len(PD) == len(M)
        f["B1"] = PN == MD; f["B2"] = refines(PN, MD); f["B3"] = refines(MD, PN); f["B4"] = len(PN) == len(MD)
        f["C"] = len(PN) == len(M)
        for side, (P, SS) in (("s", (PD, S)), ("r", (PN, SD))):
            f["Q1" + side] = canon(SS) == P
            f["Q2" + side] = all(any(set(x) <= set(b) for b in P) for x in SS)
            f["Q3" + side] = all(any(set(b) <= set(x) for x in SS) for b in P)
            f["Q4" + side] = all(frozenset(b) in set(siphons(N if side == "s" else D)) for b in P)
            f["Q5" + side] = all(set(b) == set().union(*[set(x) for x in SS if set(x) <= set(b)]) for b in P)
            f["Q6" + side] = len(P) == len(SS)
        cyc = [set(Sx) for Sx, _ in flux_rays(N)]
        f["R23"] = all(any(c <= set(b) for b in PN) for c in cyc)
        table[k] = f
    fails = {form: [k for k in table if not table[k][form]] for form in table["R"]}
    print("   Table 1 on WF-FC networks (formulation: networks where it FAILS):")
    for form, ks in fails.items():
        print("     %-4s %s" % (form, ks))
    claim("Thm E(i): every formulation of Table 1 fails on one of R, T, W, W*, F (at least one side for Q's)",
          all(fails[f] or fails.get(f[:-1] + ("r" if f.endswith("s") else "s")) for f in fails))
    stated = {"A1": "R T W Wstar F", "A4": "R T W Wstar F", "B1": "R T W Wstar F", "B4": "R T W Wstar F",
              "A2": "Wstar F", "A3": "R T W", "B2": "W F", "B3": "R T Wstar", "C": "R T Wstar F", "R23": "R T Wstar"}
    for form, s in stated.items():
        claim("Thm E(i) list: %s fails exactly on %s" % (form, s), sorted(fails[form]) == sorted(s.split()),
              "computed %s" % fails[form])

    # Theorem E (ii) on WF-FC examples
    for k in ("R", "T", "W", "Wstar", "F", "B2", "B3", "B4", "K3", "Y5"):
        N = EX[k]; D = N.dual()
        if not is_WFFC(N):
            claim("Thm E(ii) %s is WF-FC" % k, False); continue
        JF = all(sum(1 for x in N.src(j) if x) <= 1 for j in range(N.r))
        CF = all(sum(1 for j in range(N.r) if N.Ye[i][j] > 0) <= 1 for i in range(N.n))
        PD, M, PN, MD = preclusters_def19(D), ray_mips(N), preclusters_def19(N), full_coupling(N)
        Q3s = all(any(set(b) <= x for x in min_siphons(N)) for b in PD)
        Q3r = all(any(set(b) <= x for x in min_siphons(D)) for b in PN)
        ok = (PD != M and PN != MD and refines(PD, M) == JF and refines(PN, MD) == CF and Q3s == JF and Q3r == CF)
        claim("Thm E(ii) on %s (JF=%s, CF=%s)" % (k, JF, CF), ok)

    # Theorem D on WF-FC examples
    for k in ("R", "T", "W", "Wstar", "F", "B2", "B3", "B4", "K3", "Y5"):
        claim("Thm D(a)-(h) on %s" % k, thmD(EX[k]))

    # Prop G and (iv)
    D_ = EX["D"]
    iii, i_, ii, iv = propG(D_)
    claim("Prop G: D refutes (i), (ii), (iii); M(D) one block; S(D) = {X1,X2},{X1,X3}",
          (not iii) and (not i_) and (not ii) and ray_mips(D_) == [[0, 1, 2]] and
          sorted(sorted(s) for s in min_siphons(D_)) == [[0, 1], [0, 2]])
    N = EX["N412"]
    iii, i_, ii, iv = propG(N)
    claim("Prop G: N4.12 refutes (iv): law X1+X2, only minimal siphon {X1}", (not iv) and
          [sorted(s) for s in min_siphons(N)] == [[0]] and [S for S, _ in cons_rays(N)] == [(0, 1)])

    # Theorem C' / D counterexample remark: D (not WF) minimal siphons {X1,X2},{X1,X3}; only internal cycle of D*
    Dd = EX["D"].dual()
    claim("Rem after C': internal cycles of D* = {{X1,X2,X3}} (as reactions of D*)",
          [S for S, _ in flux_rays(Dd)] == [(0, 1, 2)])

    # P's table of the manuscript's counterexamples (EC/balanced/ordinary/WF/crit columns; abs. cons. = V2)
    def crit_all(N):
        sups = [frozenset(S) for S, _ in cons_rays(N)]
        return len([R for R in siphons(N) if not any(Z <= R for Z in sups)])
    rows = {"R": (1, 1, 1, 1), "T": (1, 1, 1, 1), "W": (1, 1, 1, 1), "Wstar": (1, 1, 1, 1), "F": (1, 1, 1, 1),
            "D": (1, 1, 0, 0), "U": (0, 1, 0, 0), "N3": (0, 1, 0, 0)}
    for k, (ec, bal, ordi, wf) in rows.items():
        N = EX[k]
        got = (int(is_EC(N)), int(balanced_joins(N)), int(is_ordinary(N)), int(is_WFEC(N)))
        print("   table row %s: EC/bal/ord/WF = %s ; critical siphons (all / minimal): N %d / %d, N* %d / %d" %
              (k, got, crit_all(N), len(critical_siphons(N)), crit_all(N.dual()), len(critical_siphons(N.dual()))))
        claim("table row %s (N side EC/bal/ord/WF)" % k, got == (ec, bal, ordi, wf))
    Ns = EX["N3"].dual()
    print("   N3*: EC %s bal %s ord %s WF %s" % (is_EC(Ns), balanced_joins(Ns), is_ordinary(Ns), is_WFEC(Ns)))
    Us = EX["U"].dual()
    print("   U*:  EC %s bal %s ord %s WF %s" % (is_EC(Us), balanced_joins(Us), is_ordinary(Us), is_WFEC(Us)))

    # EC not self-dual example
    N = net_from_str("2X1+X2 -> X3, X3 -> 2X1+X2")
    claim("EC alone not self-dual: 2X1+X2 <=> X3 is EC, dual is not EC", is_EC(N) and not is_EC(N.dual()))

    # K8 pairs
    k8()

    # Prop F on examples (WR EC)
    for txt in ("X1 -> X2, X2 -> X3, X3 -> X1", "X1 + X2 -> X3, X3 -> X1 + X2", "X1 -> X2, X2 -> X1, X3 -> X4, X4 -> X3",
                "2X1 -> X2, X2 -> X3 + X4, X3 + X4 -> 2X1"):
        N = net_from_str(txt)
        claim("Prop F on %r" % txt, propF(N))
    N = net_from_str("0 -> X2, X2 -> X1, X1 -> 0", species_order=["X1", "X2"])
    claim("Prop F needs 'no zero complex': 0 -> X2 -> X1 -> 0 has no siphon at all (clause fails as stated)",
          len(siphons(N)) == 0 and is_EC(N) and is_WR(N))


def propF(N):
    """Prop F: WR + EC, no zero complex => complexes pairwise support-disjoint, deficiency 0, minimal siphons per
    linkage class = one species from each complex, = supports of extreme laws, no critical siphon; if ordinary and
    (H3): every complex a single species."""
    if not (is_WR(N) and is_EC(N)):
        return None
    cpx = sorted({N.src(j) for j in range(N.r)} | {N.prod(j) for j in range(N.r)})
    disjoint = all(not any(a[i] and b[i] for i in range(N.n)) for a, b in combinations(cpx, 2))
    ok = disjoint and deficiency(N) == 0
    # predicted minimal siphons
    lc, m, l = linkage_classes(N)
    pred = set()
    for block in lc:
        cs = sorted({N.src(j) for j in block} | {N.prod(j) for j in block})
        for choice in product(*[[i for i in range(N.n) if c[i]] for c in cs]):
            pred.add(frozenset(choice))
    ms = set(min_siphons(N))
    sups = {frozenset(S) for S, _ in cons_rays(N)}
    ok = ok and ms == pred and ms == sups and not critical_siphons(N)
    if is_ordinary(N) and is_H3(N):
        ok = ok and all(sum(c) == 1 for c in cpx)
    return ok


def thmD(N):
    D = N.dual()
    ok = is_WFFC(D)  # (a)
    ms = min_siphons(N)
    cr = cons_rays(N)
    ind = lambda S, n: [Fr(int(i in S)) for i in range(n)]
    ok = ok and sorted(tuple(v) for _, v in cr) == sorted(tuple(ind(S, N.n)) for S in ms)  # (b)
    mt = minimal(t_traps(N))
    ok = ok and sorted(tuple(v) for _, v in flux_rays(N)) == sorted(tuple(ind(U, N.r)) for U in mt)  # (c)
    # (d) MIPS = atoms of minimal siphons
    prof = {}
    for i in range(N.n):
        prof.setdefault(frozenset(k for k, S in enumerate(ms) if i in S), []).append(i)
    ok = ok and canon(prof.values()) == ray_mips(N)
    # (e) profile rule for P(N*)
    Sig = {i: frozenset(k for k, S in enumerate(ms) if i in S) for i in range(N.n)}
    edges = []
    for i, j in twin_pairs(N):
        if not (Sig[i] and Sig[j]) or (Sig[i] & Sig[j]):
            return False
        edges.append((i, j))
        for s in range(N.n):
            if s not in (i, j) and Sig[s] in (Sig[i], Sig[j], Sig[i] | Sig[j]):
                edges += [(i, s), (j, s)]
    ok = ok and components(N.n, edges) == preclusters_def19(D)
    # dual profile rule for P(N)
    Th = {t: frozenset(k for k, U in enumerate(mt) if t in U) for t in range(N.r)}
    edges = []
    for i, j in doubling_pairs(N):
        edges.append((i, j))
        for s in range(N.r):
            if s not in (i, j) and Th[s] in (Th[i], Th[j], Th[i] | Th[j]):
                edges += [(i, s), (j, s)]
    ok = ok and components(N.r, edges) == preclusters_def19(N)
    ok = ok and not critical_siphons(N) and not critical_siphons(D)  # (f)
    w, nv = wvec(N), nvec(N)
    ok = ok and all(not positively_proportional(w[i], w[j]) for i, j in twin_pairs(N))  # (g)
    ok = ok and all(not positively_proportional(nv[i], nv[j]) for i, j in doubling_pairs(N))
    ok = ok and live_bounded_all_ones(N)  # (h)
    return ok


def live_bounded_all_ones(N, cap=20000):
    """Token game from M0 = (1,...,1): explore the reachability graph (bounded under (H1)); live iff every bottom
    SCC contains an occurrence of every transition."""
    M0 = tuple([1] * N.n)
    idx = {M0: 0}
    marks = [M0]
    succ = [[]]
    k = 0
    while k < len(marks):
        M = marks[k]
        for t in range(N.r):
            if all(M[i] >= N.Ye[i][t] for i in range(N.n)):
                M2 = tuple(M[i] + N.G[i][t] for i in range(N.n))
                if M2 not in idx:
                    idx[M2] = len(marks); marks.append(M2); succ.append([])
                    if len(marks) > cap:
                        return False
                succ[k].append((t, idx[M2]))
        k += 1
    # Tarjan SCC (iterative)
    nV = len(marks)
    index = [None] * nV; low = [0] * nV; onst = [False] * nV; st = []; comp = [None] * nV; c = 0; counter = 0
    for v0 in range(nV):
        if index[v0] is not None:
            continue
        stack = [(v0, 0)]
        index[v0] = low[v0] = counter; counter += 1; st.append(v0); onst[v0] = True
        while stack:
            v, pi = stack[-1]
            if pi < len(succ[v]):
                stack[-1] = (v, pi + 1)
                w = succ[v][pi][1]
                if index[w] is None:
                    index[w] = low[w] = counter; counter += 1; st.append(w); onst[w] = True
                    stack.append((w, 0))
                elif onst[w]:
                    low[v] = min(low[v], index[w])
            else:
                stack.pop()
                if stack:
                    u = stack[-1][0]; low[u] = min(low[u], low[v])
                if low[v] == index[v]:
                    while True:
                        w = st.pop(); onst[w] = False; comp[w] = c
                        if w == v:
                            break
                    c += 1
    bottom = [True] * c
    fired = [set() for _ in range(c)]
    for v in range(nV):
        for t, w in succ[v]:
            if comp[w] != comp[v]:
                bottom[comp[v]] = False
            else:
                fired[comp[v]].add(t)
    return all(len(fired[k]) == N.r for k in range(c) if bottom[k])


def k8():
    for a, b in (("N1", "N2"), ("E1", "E2")):
        A, B = K8[a], K8[b]
        h = all(is_H1(X) and is_H2(X) and is_H3(X) for X in (A, B))
        def twinclasses(N):
            d = {}
            for i in range(N.n):
                if any(N.Ye[i]):
                    d.setdefault(tuple(N.Ye[i]), []).append(i)
            return canon(d.values())
        same_tw = twinclasses(A) == twinclasses(B)
        same_cl = sorted(sorted(S) for S, T in clusters(A) if S) == sorted(sorted(S) for S, T in clusters(B) if S)
        same_ms = set(min_siphons(A)) == set(min_siphons(B))
        PA, PB = preclusters_def19(A.dual()), preclusters_def19(B.dual())
        print("   K8 %s/%s: twins %s / %s, min siphons %s, P(N*) %s vs %s, clusters same: %s, zero rows of Gamma: %s/%s"
              % (a, b, twinclasses(A), twinclasses(B), sorted(sorted(s) for s in min_siphons(A)), PA, PB, same_cl,
                 [i for i in range(A.n) if not any(A.G[i])], [i for i in range(B.n) if not any(B.G[i])]))
        claim("K8 species-side pair %s/%s: (H1)-(H3), same twins and minimal siphons, different P(N*)" % (a, b),
              h and same_tw and same_ms and PA != PB)
        if a == "N1":
            claim("K8 N1/N2 as stated: P(N1*) = {X1,X3,X4},{X2}; P(N2*) = one block; same species parts of clusters",
                  PA == [[0, 2, 3], [1]] and PB == [[0, 1, 2, 3]] and same_cl)
    A, B = K8["M1"], K8["M2"]
    h = all(is_H1(X) and is_H2(X) and is_H3(X) for X in (A, B))
    dA = canon(components(A.r, doubling_pairs(A))); dB = canon(components(B.r, doubling_pairs(B)))
    tA, tB = set(minimal(t_traps(A))), set(minimal(t_traps(B)))
    PA, PB = preclusters_def19(A), preclusters_def19(B)
    print("   K8 M1/M2: doubling classes %s / %s, min T-traps %s / %s, P %s / %s" %
          (dA, dB, sorted(sorted(u) for u in tA), sorted(sorted(u) for u in tB), PA, PB))
    claim("K8 reaction-side pair M1/M2: (H1)-(H3), same doubling classes and minimal T-traps; P(M1) one block, "
          "P(M2) = {R1},{R2},{R3,R4}", h and dA == dB and tA == tB and PA == [[0, 1, 2, 3]] and PB == [[0], [1], [2, 3]])
    claim("K8 M1/M2 have no zero row of Gamma (P: 'a pair without zero rows')",
          all(all(any(r) for r in X.G) for X in (A, B)))
    print("   K8 M1/M2: dual (H3) (columns of Gamma distinct): M1 %s, M2 %s" % (is_H3(A.dual()), is_H3(B.dual())))


# ------------------------------------------------------------------ random WF-EC networks (Theorems C, C', D)

def random_ec(rng, n, ordinary):
    species = list(range(n))
    rng.shuffle(species)
    m = rng.randint(1, min(n, 4))
    cuts = sorted(rng.sample(range(1, n), m - 1)) if m > 1 else []
    parts = [species[a:b] for a, b in zip([0] + cuts, cuts + [n])]
    # every species occurs in a source (so it is consumed)
    sources = []
    for p in parts:
        k = rng.randint(1, len(p))
        sup = rng.sample(p, k)
        y = [0] * n
        for i in sup:
            y[i] = 1 if ordinary else rng.randint(1, 2)
        sources.append(y)
    # species not in a chosen source: put them into sources anyway to keep every species consumed
    for i in range(n):
        if not any(s[i] for s in sources):
            # add to the source of its part
            for p, s in zip(parts, sources):
                if i in p:
                    s[i] = 1 if ordinary else rng.randint(1, 2)
    Ye, Yp = [[] for _ in range(n)], [[] for _ in range(n)]
    for s in sources:
        for _ in range(rng.randint(1, 3)):
            while True:
                k = rng.randint(1, 2)
                p = [0] * n
                for i in rng.sample(range(n), k):
                    p[i] = 1 if ordinary else rng.randint(1, 2)
                if p != s:
                    break
            for i in range(n):
                Ye[i].append(s[i]); Yp[i].append(p[i])
    return Net(Ye, Yp)


def random_wfec(seed, count, ordinary, nmin=3, nmax=5, tries=400000):
    rng = random.Random(seed)
    out, seen = [], set()
    for _ in range(tries):
        N = random_ec(rng, rng.randint(nmin, nmax), ordinary)
        key = (tuple(map(tuple, N.Ye)), tuple(map(tuple, N.Yp)))
        if key in seen:
            continue
        seen.add(key)
        if N.r > 9:
            continue
        if is_WFEC(N):
            out.append(N)
            if len(out) >= count:
                break
    return out


def thmC_Cp(N):
    D = N.dual()
    noC = not critical_siphons(N)
    noCp = not critical_siphons(D)
    ms = set(min_siphons(N)) == {frozenset(S) for S, _ in cons_rays(N)}
    mt = set(minimal(t_traps(N))) == {frozenset(S) for S, _ in flux_rays(N)}
    return noC, noCp, ms, mt


def wfec_tests():
    t = time.time()
    cnt = 25 if QUICK else 150
    for ordinary in (True, False):
        nets = random_wfec(7 if ordinary else 11, cnt, ordinary)
        pure = [N for N in nets if is_pure(N)]
        impure = [N for N in nets if not is_pure(N)]
        badC = [N for N in pure if not thmC_Cp(N)[0] or not thmC_Cp(N)[2]]
        badCp = [N for N in nets if not thmC_Cp(N)[1] or not thmC_Cp(N)[3]]
        imp = [N for N in impure if not thmC_Cp(N)[0]]
        lab = "ordinary (WF-FC)" if ordinary else "weighted"
        claim("Thm C on %d random pure WF-EC %s networks (no critical siphon; min siphons = extreme-law supports)"
              % (len(pure), lab), not badC, repr(badC[:1]))
        claim("Thm C' on %d random WF-EC %s networks (N* no critical siphon; min T-traps = internal cycles)"
              % (len(nets), lab), not badCp, repr(badCp[:1]))
        print("   impure WF-EC %s networks: %d, of which with a critical siphon: %d (Theorem C without purity)" %
              (lab, len(impure), len(imp)))
        bal = [N for N in nets if balanced_joins(N)]
        claim("balanced-joins WF-EC => N* WF-EC (%d %s networks)" % (len(bal), lab),
              all(is_WFEC(N.dual()) and is_EC(N.dual()) for N in bal))
        if ordinary:
            badD = [N for N in nets if not thmD(N)]
            claim("Thm D(a)-(h) on %d random WF-FC networks" % len(nets), not badD, repr(badD[:1]))
    print("   random WF-EC tests: %.1f s" % (time.time() - t))


if __name__ == "__main__":
    examples()
    wfec_tests()
    if QUICK:
        run_family("n=2,mol<=2,r<=3 (zero complex allowed)", family(2, 2, 3))
    else:
        run_family("n=2,mol<=2,r<=4 (zero complex allowed)", family(2, 2, 4))
        run_family("n=3,mol<=2,r<=3 (no zero complex)", family(3, 2, 3, zero=False))
        run_family("n=2,mol<=3,r<=3 (no zero complex)", family(2, 3, 3, zero=False))
    print("runtime %.1f s; %d claims, %d FAIL" % (time.time() - T0, len(RESULTS), len(FAILS)))
    if FAILS:
        print("FAILED:", FAILS)
        sys.exit(1)
