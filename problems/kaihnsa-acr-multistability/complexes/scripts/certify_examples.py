"""Exact certificates for three 6-complex examples N1, N2, N3 (runtime: well under 1 second; Python standard
library only: fractions).  Usage:  python certify_examples.py      (exit code 1 if any check fails)

All three networks are at most bimolecular, have species X, Y, Z, exactly 6 complexes (all of them reactant
complexes) and dimension 3 (full-dimensional, so the only compatibility class is R^3_{>0}).

N1 (main example; one linkage class, deficiency 2, 6 reactions):
    Z -> Y+Z (5),  2Z -> Z (1),  Y -> Y+Z (1),  Y+Z -> 2X (4),  X+Z -> 2Z (10),  2X -> Y (7)
    weight w = (1,2,0):  f_X + 2 f_Y = 2 k1 z - k5 x z = 10 z (1 - x).
    positive steady states (exactly): (1,19/4,1/2) stable, (1,3,1) unstable, (1,7/4,7/2) stable.

N2 (two linkage classes, weakly reversible, deficiency 1, linkage-class deficiencies 0 and 0, 7 reactions):
    2Z -> X+Z (1),  X+Z -> 2Z (5),
    Y+Z -> X (2),  Y+Z -> X+Y (5),  X -> 2X (3),  X+Y -> X (1),  2X -> Y+Z (3)
    weight w = (1,0,1):  f_X + f_Z = k4 x - k7 x^2 = 3 x (1 - x).
    positive steady states (exactly): (1,3/2,1/2) stable, (1,1,1) unstable, (1,3/7,3) stable.

N3 (one linkage class, deficiency 2, 8 reactions):
    2Z -> Y+Z (1),  Y -> X (6),  Y+Z -> 2Z (4),  Y+Z -> 2X (1),
    X -> X+Z (1),  X+Z -> Y+Z (10),  X+Z -> 2X (3),  2X -> X (1)
    weight w = (1,1,1):  f_X + f_Y + f_Z = k5 x - k8 x^2 = x (1 - x).
    positive steady states (exactly): (1,1,1) stable, (1,23/18,3/2) unstable, (1,3/2,2) stable.

Checks for each network (PASS/FAIL):
  S1 counts: species, reactions, complexes, reactant complexes, molecularity <= 2
  S2 dim S = 3 (full-dimensional: the only stoichiometric compatibility class is R^3_{>0})
  S3 linkage classes, terminal strong linkage classes, deficiency, linkage-class deficiencies, weak reversibility
  A1 symbolic (all rate constants): every reaction except the two reactions out of the pair of complexes
     {y, y + X} keeps w.(complex) unchanged, and w.f = x^y (a - b x) with a, b positive multiples of the two
     corresponding rate constants; hence for EVERY rate vector every positive steady state has x = a/b (ACR in X)
  A2 with the chosen rates w.f = c0 x^y (1 - x), so x = 1 at every positive steady state
  P1 on x = 1 one equation reads A(z) - D(z) y with A, D having positive coefficients (so y = A/D > 0 is forced)
  P2 D(z) * (other equation)(1, A/D, z) equals the claimed cubic, whose exact factorisation gives three positive
     roots; hence the positive steady states are exactly the three listed points
  P3 f vanishes exactly at the three points
  J  exact Jacobians; Routh-Hurwitz (a1 > 0, a3 > 0, a1 a2 > a3) at the two stable points; det J > 0 at the middle
     point (so it has a positive real eigenvalue); det J != 0 at all three (nondegenerate)
  E  for every rate vector there is a positive steady state (cubic with sign change between 0 and infinity,
     checked symbolically), so the network has unconditional ACR in X
"""
import sys
import time
from fractions import Fraction as Fr

SPECIES = ["X", "Y", "Z"]


def parse(s):
    v = [0, 0, 0]
    if s == "0":
        return tuple(v)
    for part in s.split("+"):
        part = part.strip()
        coef = 1
        if part[0].isdigit():
            coef = int(part[0])
            part = part[1:]
        v[SPECIES.index(part)] += coef
    return tuple(v)


def cname(c):
    parts = []
    for k in range(3):
        if c[k] == 1:
            parts.append(SPECIES[k])
        elif c[k] > 1:
            parts.append(str(c[k]) + SPECIES[k])
    return "+".join(parts) if parts else "0"


NETWORKS = {
    "N1": dict(
        reactions=[("Z", "Y+Z", "k1"), ("2Z", "Z", "k2"), ("Y", "Y+Z", "k3"), ("Y+Z", "2X", "k4"),
                   ("X+Z", "2Z", "k5"), ("2X", "Y", "k6")],
        rates=dict(k1=5, k2=1, k3=1, k4=4, k5=10, k6=7),
        w=(1, 2, 0), pair=("k1", "k5"),
        wf_sym={("k1", (0, 0, 1)): 2, ("k5", (1, 0, 1)): -1},
        wf_num={(0, 0, 1): 10, (1, 0, 1): -10},
        expect=dict(species=3, reactions=6, complexes=6, reactants=6, ell=1, t=1, delta=2, lc_deficiencies=[2],
                    weakly_reversible=False),
        ysolve=1, zeq=2,
        cubic={3: -4, 2: 20, 1: -23, 0: 7},  # = -(2z - 1)(z - 1)(2z - 7)
        factors=[(-1, 2), (-1, 1), (-7, 2)],
        factor_sign=-1,
        points=[(1, Fr(19, 4), Fr(1, 2)), (1, 3, 1), (1, Fr(7, 4), Fr(7, 2))],
        stable=[True, False, True],
        symY={("k1", (0, 0, 1)): 1, ("k4", (0, 1, 1)): -1, ("k6", (2, 0, 0)): 1},
        symZ={("k2", (0, 0, 2)): -1, ("k3", (0, 1, 0)): 1, ("k4", (0, 1, 1)): -1, ("k5", (1, 0, 1)): 1},
        existence_note="on x = c = 2 k1/k5: y = (k1 z + k6 c^2)/(k4 z); the cubic "
                       "(-k2 z^2 + k5 c z)(k4 z) + (k3 - k4 z)(k1 z + k6 c^2) = -k2 k4 z^3 + k1 k4 z^2 "
                       "+ (k1 k3 - k4 k6 c^2) z + k3 k6 c^2 has constant term > 0 and leading coefficient < 0",
    ),
    "N2": dict(
        reactions=[("2Z", "X+Z", "k1"), ("Y+Z", "X", "k2"), ("Y+Z", "X+Y", "k3"), ("X", "2X", "k4"),
                   ("X+Z", "2Z", "k5"), ("X+Y", "X", "k6"), ("2X", "Y+Z", "k7")],
        rates=dict(k1=1, k2=2, k3=5, k4=3, k5=5, k6=1, k7=3),
        w=(1, 0, 1), pair=("k4", "k7"),
        wf_sym={("k4", (1, 0, 0)): 1, ("k7", (2, 0, 0)): -1},
        wf_num={(1, 0, 0): 3, (2, 0, 0): -3},
        expect=dict(species=3, reactions=7, complexes=6, reactants=6, ell=2, t=2, delta=1, lc_deficiencies=[0, 0],
                    weakly_reversible=True),
        ysolve=1, zeq=2,
        cubic={3: -2, 2: 9, 1: -10, 0: 3},  # = -(z - 1)(2z - 1)(z - 3)
        factors=[(-1, 1), (-1, 2), (-3, 1)],
        factor_sign=-1,
        points=[(1, Fr(3, 2), Fr(1, 2)), (1, 1, 1), (1, Fr(3, 7), 3)],
        stable=[True, False, True],
        symY={("k2", (0, 1, 1)): -1, ("k6", (1, 1, 0)): -1, ("k7", (2, 0, 0)): 1},
        symZ={("k1", (0, 0, 2)): -1, ("k2", (0, 1, 1)): -1, ("k3", (0, 1, 1)): -1, ("k5", (1, 0, 1)): 1,
              ("k7", (2, 0, 0)): 1},
        existence_note="on x = c = k4/k7: y = k7 c^2/(k2 z + k6 c); the cubic "
                       "(-k1 z^2 + k5 c z + k7 c^2)(k2 z + k6 c) - (k2+k3) k7 c^2 z has constant term "
                       "k6 k7 c^3 > 0 and leading coefficient -k1 k2 < 0",
    ),
    "N3": dict(
        reactions=[("2Z", "Y+Z", "k1"), ("Y", "X", "k2"), ("Y+Z", "2Z", "k3"), ("Y+Z", "2X", "k4"),
                   ("X", "X+Z", "k5"), ("X+Z", "Y+Z", "k6"), ("X+Z", "2X", "k7"), ("2X", "X", "k8")],
        rates=dict(k1=1, k2=6, k3=4, k4=1, k5=1, k6=10, k7=3, k8=1),
        w=(1, 1, 1), pair=("k5", "k8"),
        wf_sym={("k5", (1, 0, 0)): 1, ("k8", (2, 0, 0)): -1},
        wf_num={(1, 0, 0): 1, (2, 0, 0): -1},
        expect=dict(species=3, reactions=8, complexes=6, reactants=6, ell=1, t=1, delta=2, lc_deficiencies=[2],
                    weakly_reversible=False),
        ysolve=1, zeq=2,  # solve f_Y for y, substitute in f_Z
        cubic={3: -2, 2: 9, 1: -13, 0: 6},  # = -2 (z - 1)(z - 3/2)(z - 2)
        factors=[(-2, 2), (Fr(-3, 2), 1), (-2, 1)],  # (-2 + 2 z)... see check: product of (a + b z)
        factor_sign=-1,
        points=[(1, 1, 1), (1, Fr(23, 18), Fr(3, 2)), (1, Fr(3, 2), 2)],
        stable=[True, False, True],
        # symbolic f_Y, f_Z (rate name, monomial (x,y,z)) -> coefficient, used for the existence argument E
        symY={("k1", (0, 0, 2)): 1, ("k2", (0, 1, 0)): -1, ("k3", (0, 1, 1)): -1, ("k4", (0, 1, 1)): -1,
              ("k6", (1, 0, 1)): 1},
        symZ={("k1", (0, 0, 2)): -1, ("k3", (0, 1, 1)): 1, ("k4", (0, 1, 1)): -1, ("k5", (1, 0, 0)): 1,
              ("k7", (1, 0, 1)): -1},
        existence_note="on x = c = k5/k8: y = (k1 z^2 + k6 c z)/(k2 + (k3+k4) z); the cubic "
                       "(-k1 z^2 + k5 c - k7 c z)(k2 + (k3+k4) z) + (k3-k4) z (k1 z^2 + k6 c z) has constant term "
                       "k2 k5 c > 0 and leading coefficient -2 k1 k4 < 0",
    ),
}

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("PASS " if cond else "FAIL ") + name + (("  [" + detail + "]") if detail else ""))
    if not cond:
        FAIL += 1


# ---------------- polynomial helpers: dict exponent-tuple -> Fraction
def padd(p, q, s=1):
    r = dict(p)
    for m, c in q.items():
        r[m] = r.get(m, 0) + s * c
        if r[m] == 0:
            del r[m]
    return r


def pmul(p, q):
    r = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(a + b for a, b in zip(m1, m2))
            r[m] = r.get(m, 0) + c1 * c2
    return {m: c for m, c in r.items() if c != 0}


def peval(p, pt):
    s = Fr(0)
    for m, c in p.items():
        t = Fr(c)
        for a, v in zip(m, pt):
            t *= Fr(v) ** a
        s += t
    return s


def pdiff(p, i):
    r = {}
    for m, c in p.items():
        if m[i] > 0:
            mm = list(m)
            mm[i] -= 1
            r[tuple(mm)] = r.get(tuple(mm), 0) + c * m[i]
    return {m: c for m, c in r.items() if c != 0}


def rank(rows):
    M = [list(map(Fr, r)) for r in rows]
    rk, col = 0, 0
    ncol = len(M[0]) if M else 0
    while rk < len(M) and col < ncol:
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None:
            col += 1
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                fct = M[i][col] / M[rk][col]
                M[i] = [a - fct * b for a, b in zip(M[i], M[rk])]
        rk += 1
        col += 1
    return rk


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def rh(M):
    a1 = -(M[0][0] + M[1][1] + M[2][2])
    a2 = ((M[0][0] * M[1][1] - M[0][1] * M[1][0]) + (M[0][0] * M[2][2] - M[0][2] * M[2][0])
          + (M[1][1] * M[2][2] - M[1][2] * M[2][1]))
    a3 = -det3(M)
    return a1, a2, a3


def structure(reactions):
    complexes, reactants = [], []
    for s, t, _ in reactions:
        for c in (s, t):
            if parse(c) not in complexes:
                complexes.append(parse(c))
        if parse(s) not in reactants:
            reactants.append(parse(s))
    m = len(complexes)
    idx = {c: i for i, c in enumerate(complexes)}
    adj = [[False] * m for _ in range(m)]
    for s, t, _ in reactions:
        adj[idx[parse(s)]][idx[parse(t)]] = True
    reach = [[i == j or adj[i][j] for j in range(m)] for i in range(m)]
    for k in range(m):
        for i in range(m):
            for j in range(m):
                if reach[i][k] and reach[k][j]:
                    reach[i][j] = True
    comp = list(range(m))
    changed = True
    while changed:
        changed = False
        for i in range(m):
            for j in range(m):
                if (adj[i][j] or adj[j][i]) and comp[i] != comp[j]:
                    a, b = min(comp[i], comp[j]), max(comp[i], comp[j])
                    comp = [a if c == b else c for c in comp]
                    changed = True
    classes = sorted(set(comp))
    slcs = sorted(set(frozenset(j for j in range(m) if reach[i][j] and reach[j][i]) for i in range(m)), key=min)
    terminal = [K for K in slcs if all((not reach[i][j]) or (j in K) for i in K for j in range(m))]
    vecs = [tuple(b - a for a, b in zip(parse(s), parse(t))) for s, t, _ in reactions]
    s_dim = rank(vecs)
    lc_def = []
    for L in classes:
        members = [complexes[i] for i in range(m) if comp[i] == L]
        diffs = [tuple(a - b for a, b in zip(c, members[0])) for c in members[1:]]
        lc_def.append(len(members) - 1 - (rank(diffs) if diffs else 0))
    weakly_rev = all(reach[j][i] for i in range(m) for j in range(m) if adj[i][j])
    return dict(complexes=complexes, reactants=reactants, ell=len(classes), t=len(terminal), s=s_dim,
                delta=m - len(classes) - s_dim, lc_def=lc_def, weakly_reversible=weakly_rev,
                slcs=[[cname(complexes[i]) for i in sorted(K)] for K in slcs],
                terminal=[[cname(complexes[i]) for i in sorted(K)] for K in terminal])


def symbolic_ode(reactions):
    """f as dict over (rate-index-vector + x-exponent) monomials"""
    NR = len(reactions)
    f = [dict(), dict(), dict()]
    for e, (s, t, _) in enumerate(reactions):
        ys, yt = parse(s), parse(t)
        mono = tuple([1 if i == e else 0 for i in range(NR)] + list(ys))
        for sp in range(3):
            d = yt[sp] - ys[sp]
            if d:
                f[sp] = padd(f[sp], {mono: Fr(d)})
    return f


def specialise(p, reactions, rates):
    NR = len(reactions)
    r = {}
    for mono, c in p.items():
        coef = Fr(c)
        for e in range(NR):
            coef *= Fr(rates[reactions[e][2]]) ** mono[e]
        r[mono[NR:]] = r.get(mono[NR:], 0) + coef
    return {a: b for a, b in r.items() if b != 0}


def ode(reactions, rates):
    return [specialise(p, reactions, rates) for p in symbolic_ode(reactions)]


def jacobian(f, p):
    return [[peval(pdiff(f[i], j), p) for j in range(3)] for i in range(3)]


def poly_str(p):
    terms = []
    for m, c in sorted(p.items(), key=lambda t: (-sum(t[0]), t[0])):
        mon = "*".join((SPECIES[i].lower() + ("^%d" % m[i] if m[i] > 1 else "")) for i in range(3) if m[i])
        terms.append("%s%s" % (str(c), ("*" + mon) if mon else ""))
    return " + ".join(terms).replace("+ -", "- ")


def certify(name, N):
    print("=" * 100)
    print("Network", name)
    R = N["reactions"]
    for s, t, kn in R:
        print("    %-4s -> %-4s  %s = %s" % (s, t, kn, N["rates"][kn]))
    st = structure(R)
    ex = N["expect"]
    allc = st["complexes"]
    check(name + " S1 counts: 3 species, %d reactions, 6 complexes, %d reactant complexes, molecularity <= 2"
          % (ex["reactions"], ex["reactants"]),
          len(R) == ex["reactions"] and len(allc) == ex["complexes"] and len(st["reactants"]) == ex["reactants"]
          and all(sum(c) <= 2 for c in allc) and all(any(c[i] > 0 for c in allc) for i in range(3)),
          "complexes: " + ", ".join(cname(c) for c in allc))
    check(name + " S2 dim S = 3 (full-dimensional)", st["s"] == 3)
    check(name + " S3 linkage classes l=%d, terminal strong linkage classes t=%d, deficiency %d, "
                 "linkage-class deficiencies %s, weakly reversible = %s"
          % (ex["ell"], ex["t"], ex["delta"], ex["lc_deficiencies"], ex["weakly_reversible"]),
          st["ell"] == ex["ell"] and st["t"] == ex["t"] and st["delta"] == ex["delta"]
          and sorted(st["lc_def"]) == sorted(ex["lc_deficiencies"])
          and st["weakly_reversible"] == ex["weakly_reversible"],
          "strong linkage classes %s; terminal %s" % (st["slcs"], st["terminal"]))

    # ---- ACR, symbolic in the rate constants
    w = N["w"]
    ka, kb = N["pair"]
    neutral = all(sum(wi * (b - a) for wi, a, b in zip(w, parse(s), parse(t))) == 0
                  for s, t, kn in R if kn not in (ka, kb))
    fsym = symbolic_ode(R)
    NR = len(R)
    wf = {}
    for i in range(3):
        if w[i]:
            wf = padd(wf, {m: c * w[i] for m, c in fsym[i].items()})
    names = [kn for _, _, kn in R]
    want = {tuple([1 if names[e] == rn else 0 for e in range(NR)] + list(mono)): Fr(cf)
            for (rn, mono), cf in N["wf_sym"].items()}
    check(name + " A1 for all rate constants: w.f = %s with w = %s (all reactions other than %s keep w.y)"
          % (" + ".join("%s*%s*x^%s" % (cf, rn, mono) for (rn, mono), cf in N["wf_sym"].items()), w, (ka, kb)),
          neutral and wf == want)
    f = ode(R, N["rates"])
    for i in range(3):
        print("    d%s/dt = %s" % (SPECIES[i].lower(), poly_str(f[i])))
    wfn = {}
    for i in range(3):
        if w[i]:
            wfn = padd(wfn, {m: c * w[i] for m, c in f[i].items()})
    check(name + " A2 with the chosen rates w.f = %s: since it is a positive monomial times (1 - x), x = 1 at "
                 "every positive steady state" % poly_str({m: Fr(c) for m, c in N["wf_num"].items()}),
          wfn == {m: Fr(c) for m, c in N["wf_num"].items()}
          and sum(N["wf_num"].values()) == 0 and len(N["wf_num"]) == 2
          and [m for m in N["wf_num"]][1][0] == [m for m in N["wf_num"]][0][0] + 1
          and [m for m in N["wf_num"]][1][1:] == [m for m in N["wf_num"]][0][1:])

    # ---- plane x = 1
    def on_plane(p):
        r = {}
        for (a, b, c), v in p.items():
            r[(b, c)] = r.get((b, c), 0) + v
        return {q: v for q, v in r.items() if v != 0}

    fy, fz = on_plane(f[N["ysolve"]]), on_plane(f[N["zeq"]])
    aff = all(b <= 1 for (b, c) in fy) and all(b <= 1 for (b, c) in fz)
    A = {c: v for (b, c), v in fy.items() if b == 0}
    D = {c: -v for (b, c), v in fy.items() if b == 1}  # D(z) > 0 for z > 0 (positive coefficients)
    Az = {c: v for (b, c), v in fz.items() if b == 0}
    Bz = {c: v for (b, c), v in fz.items() if b == 1}
    check(name + " P1 on x = 1: f_%s = A(z) - D(z) y with A = %s, D = %s (positive coefficients), f_%s affine in y"
          % (SPECIES[N["ysolve"]], sorted(A.items()), sorted(D.items()), SPECIES[N["zeq"]]),
          aff and len(A) > 0 and len(D) > 0 and all(v > 0 for v in A.values()) and all(v > 0 for v in D.values()))

    def umul(p, q):
        r = {}
        for a, u in p.items():
            for b, v in q.items():
                r[a + b] = r.get(a + b, 0) + u * v
        return {e: c for e, c in r.items() if c != 0}

    def uadd(p, q):
        r = dict(p)
        for e, c in q.items():
            r[e] = r.get(e, 0) + c
        return {e: c for e, c in r.items() if c != 0}

    Ncub = uadd(umul(Az, D), umul(Bz, A))
    target = {e: Fr(c) for e, c in N["cubic"].items()}
    prod = {0: Fr(N["factor_sign"])}
    for a, b in N["factors"]:
        prod = umul(prod, {0: Fr(a), 1: Fr(b)})
    roots = sorted(Fr(-a) / Fr(b) for a, b in N["factors"])
    check(name + " P2 D(z) f_%s(1, A/D, z) = %s = %s * prod (a + b z) over %s; positive roots z = %s"
          % (SPECIES[N["zeq"]], sorted(Ncub.items(), reverse=True), N["factor_sign"], N["factors"],
             [str(r) for r in roots]),
          Ncub == target and prod == target and all(r > 0 for r in roots) and len(set(roots)) == 3)
    pts = [(Fr(1), peval({(c,): v for c, v in A.items()}, (z,)) / peval({(c,): v for c, v in D.items()}, (z,)), z)
           for z in roots]
    expected_pts = [tuple(Fr(v) for v in p) for p in N["points"]]
    check(name + " P3 the positive steady states are exactly " + ", ".join(
        "(" + ", ".join(str(v) for v in p) + ")" for p in expected_pts) + " and f vanishes there exactly",
          sorted(pts) == sorted(expected_pts) and all(peval(fi, p) == 0 for fi in f for p in expected_pts))

    # ---- stability
    for p, want_st in zip(expected_pts, N["stable"]):
        J = jacobian(f, p)
        a1, a2, a3 = rh(J)
        lab = "(" + ", ".join(str(v) for v in p) + ")"
        print("    at %s: J = %s" % (lab, [[str(v) for v in row] for row in J]))
        print("      det(lambda I - J) = lambda^3 + (%s) lambda^2 + (%s) lambda + (%s);  a1 a2 - a3 = %s"
              % (a1, a2, a3, a1 * a2 - a3))
        check(name + " J  %s nondegenerate (det J = %s)" % (lab, det3(J)), det3(J) != 0)
        if want_st:
            check(name + " J  %s exponentially stable (Routh-Hurwitz)" % lab, a1 > 0 and a3 > 0 and a1 * a2 > a3)
        else:
            check(name + " J  %s unstable (det J > 0)" % lab, det3(J) > 0)

    # ---- existence of a positive steady state for every rate vector (symbolic structure)
    def symdict(p):
        out = {}
        for mono, c in p.items():
            out[(R[mono[:NR].index(1)][2], mono[NR:])] = c
        return out

    check(name + " E  symbolic f_Y, f_Z as used in: " + N["existence_note"]
          + " => a positive steady state exists for every rate vector (unconditional ACR in X)",
          symdict(fsym[1]) == N["symY"] and symdict(fsym[2]) == N["symZ"])
    return f, expected_pts


if __name__ == "__main__":
    T0 = time.time()
    for nm, net in NETWORKS.items():
        certify(nm, net)
    print("runtime %.3f s" % (time.time() - T0))
    if FAIL:
        print("%d check(s) FAILED" % FAIL)
        sys.exit(1)
    print("ALL PASS")
