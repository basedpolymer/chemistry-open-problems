"""s = 2: certificates for persistence of the dual N*, facet analysis, and numerical evidence (task 2 of BRIEF-E2;
Q1-7, Q1-8).

Environment: Python 3.12 (WSL Fedora 44).  Certificates are exact (standard library, Fractions); the numerical part
(numpy/scipy LSODA) is evidence only and never decides a PASS/FAIL.  Imports Q1's drainsets.drain_sets (unchanged).
Usage (from round2/computation/): python scripts/s2_persistence.py [--limit K] [--nproc P] [--numeric M]

Networks: weakly reversible N with (H1)-(H3) and s = 2 from every first-round record file (at most K per file,
seeded).  For N*:
  ALL  certified if N* has no critical siphon (Angeli-De Leenheer-Sontag Thm 2: persistent for every kappa);
  NONE certified if N* has a drain species (E.4) or a drain set (Q1's Theorem D): non-persistent for every kappa;
  facets: Z is a realizable facet zero set if Z is a siphon of N*, rank Gamma*[Z,:] = 1 with spanning vector g of one
       strict sign on Z (then every point with zero set exactly Z is in the relative interior of the Z-facet of its
       class).  Near such a point the Z-coordinates move along g, and the normal speed has the s = 1 expansion of
       Prop. E2.1 with J = Z: the facet is IN (good order strictly lower), OUT (bad strictly lower), ESC (tie with a bad
       Newton vertex) or DOM (tie, bad exponents in the good hull).  Q1's Conjecture F: an OUT facet forces NONE.
Per network (first M undecided, by size): 4 sampled rational kappa (10^k * p/q, k in -2..2).  For each kappa:
  'cert'    exact: in one of the classes T*u (u random positive rational, T = 10^-6, 10^-4, ..., 10^6) a vertex of the
            class polygon is an equilibrium whose Jacobian restricted to S* has trace < 0 and determinant > 0, so it
            attracts interior points: N* is not persistent for this kappa (certificate);
  'num-ext' numerical only: LSODA from T*(1,..,1) with one coordinate 1e-3*T (T = 1e-2, 1, 1e2), some coordinate ends
            below 1e-6*T and decreasing at time 1e3 (evidence of extinction, not a proof);
  'none'    neither.  A network extinct for some kappa and not for another is a
candidate against Q1-7 (dichotomy at s = 2), to be certified separately.
"""
import os
import sys
import time
import random
import platform
from fractions import Fraction
from itertools import combinations
from multiprocessing import Pool

import e2lib
from e2lib import CNet, crn, ecore, endo

QDIR = os.path.normpath(os.path.join(e2lib.HERE, "..", "..", "persistence", "scripts"))
sys.path.insert(0, QDIR)
import drainsets  # noqa: E402


def facet_types(p, rx):
    """Facet analysis of the network (p species, reactions rx).  Returns list of (Z, type)."""
    X = CNet(p, rx)
    G = X.Gamma  # p x nreac
    nre = len(rx)
    out = []
    sy = ecore.siphon_masks(X)
    for Sm in sy:
        Z = sorted(ecore.to_set(Sm))
        rows = [G[i] for i in Z]
        if crn.rank(rows) != 1:
            continue
        # spanning vector g on Z: Gamma[Z, j] = g * c_j
        j0 = next(j for j in range(nre) if any(G[i][j] for i in Z))
        g = [Fraction(G[i][j0]) for i in Z]
        if not (all(x > 0 for x in g) or all(x < 0 for x in g)):
            continue
        if g[0] > 0:
            g = [-x for x in g]  # orient: moving along +g decreases x_Z (towards the facet)
        c = []
        for j in range(nre):
            col = [Fraction(G[i][j]) for i in Z]
            k = next((t for t in range(len(Z)) if g[t] != 0), 0)
            c.append(col[k] / g[k])
        # reaction j pushes towards the facet iff c_j > 0
        d = {j: sum(rx[j][0][i] for i in Z) for j in range(nre) if c[j] != 0}
        bad = [j for j in d if c[j] > 0]
        good = [j for j in d if c[j] < 0]
        if not bad:
            out.append((Z, "IN"))
            continue
        if not good:
            out.append((Z, "OUT"))
            continue
        db, dg = min(d[j] for j in bad), min(d[j] for j in good)
        if dg < db:
            out.append((Z, "IN"))
        elif db < dg:
            out.append((Z, "OUT"))
        else:
            Zc = [i for i in range(p) if i not in Z]
            ge = sorted({tuple(rx[j][0][i] for i in Zc) for j in good if d[j] == dg})
            esc = False
            for j in bad:
                if d[j] != db:
                    continue
                e = tuple(rx[j][0][i] for i in Zc)
                if e in ge:
                    continue
                if e2lib.in_convex_hull(list(e), [list(x) for x in ge]) is None:
                    esc = True
                    break
            out.append((Z, "ESC" if esc else "DOM"))
    return out


SCALES = [Fraction(10) ** k for k in range(-6, 7, 2)]
NSCALES = [1e-2, 1.0, 1e2]


def stable_vertex(p, rx, kappa, x0):
    """Exact: a vertex of the class of x0 that is an equilibrium with the 2x2 restricted Jacobian of trace < 0 and
    determinant > 0 (asymptotically stable, hence attracting interior points: non-persistence for kappa)."""
    import g6_props
    G = CNet(p, rx).Gamma
    B = g6_props.basis([[Fraction(G[i][j]) for i in range(p)] for j in range(len(rx))])
    if len(B) != 2:
        return None
    for v in g6_props.vertices(x0, B, p):
        f = [Fraction(0)] * p
        for k, (a, b) in enumerate(rx):
            rate = kappa[k]
            for i in range(p):
                rate *= v[i] ** a[i]
            for i in range(p):
                f[i] += rate * (b[i] - a[i])
        if any(f):
            continue
        J = [[Fraction(0)] * p for _ in range(p)]
        for k, (a, b) in enumerate(rx):
            for j in range(p):
                if a[j] == 0:
                    continue
                der = kappa[k] * a[j]
                for i in range(p):
                    e = a[i] - (1 if i == j else 0)
                    der *= v[i] ** e
                for i in range(p):
                    J[i][j] += der * (b[i] - a[i])
        JB = [[sum(J[i][l] * B[c][l] for l in range(p)) for i in range(p)] for c in range(2)]
        # coordinates of JB[c] in the basis B (least squares exact: Gram system)
        g11 = sum(x * x for x in B[0])
        g12 = sum(x * y for x, y in zip(B[0], B[1]))
        g22 = sum(x * x for x in B[1])
        det = g11 * g22 - g12 * g12
        M = []
        for c in range(2):
            r1 = sum(x * y for x, y in zip(B[0], JB[c]))
            r2 = sum(x * y for x, y in zip(B[1], JB[c]))
            M.append(((g22 * r1 - g12 * r2) / det, (-g12 * r1 + g11 * r2) / det))
        tr = M[0][0] + M[1][1]
        dt = M[0][0] * M[1][1] - M[0][1] * M[1][0]
        if tr < 0 and dt > 0:
            return v
    return None


def numeric(p, rx, facets, seed):
    """Exact stable-vertex certificates for 4 sampled kappa over scaled classes, and LSODA evidence over the same
    scaled classes for the kappa without a certificate."""
    import warnings
    warnings.filterwarnings("ignore")
    import numpy as np
    from scipy.integrate import solve_ivp
    rng = random.Random(seed)
    S_ = np.array([a for a, b in rx], dtype=float)
    D_ = np.array([[bb - aa for aa, bb in zip(a, b)] for a, b in rx], dtype=float)
    res = []
    for _ in range(4):
        kappa = [Fraction(10) ** rng.randint(-2, 2) * Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in rx]
        cert = False
        u = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(p)]
        for T in SCALES:
            if stable_vertex(p, rx, kappa, [T * x for x in u]) is not None:
                cert = True
                break
        if cert:
            res.append("cert")
            continue
        kap = np.array([float(k) for k in kappa])

        def rhs(t, x):
            x = np.maximum(x, 0.0)
            return (kap * np.prod(x ** S_, axis=1)) @ D_
        ext = False
        for T in NSCALES:
            for i0 in range(p):
                x0 = float(T) * np.ones(p)
                x0[i0] = float(T) * 1e-3
                sol = solve_ivp(rhs, (0, 1e3), x0, method="LSODA", rtol=1e-9, atol=1e-15)
                if not sol.success:
                    continue
                xe = sol.y[:, -1]
                f = rhs(0, xe)
                if any(xe[i] < 1e-6 * float(T) and f[i] <= 0 for i in range(p)):
                    ext = True
                    break
            if ext:
                break
        res.append("num-ext" if ext else "none")
    return res


def vertex_scaling(p, rx):
    """Rate-independent NONE certificate (Prop. E2.7 of notes.md), s = 2.  Z: set of species with
    rank Gamma[Z,:] = 2, every reaction source meets Z (so the vertex v with zero set Z is an equilibrium for every
    kappa), W = projection of S onto R^Z contains a strictly negative vector (Z is the zero set of a vertex of some
    class).  The linearisation at v acts on W by L = sum_k kappa_k v^{w_k} M_k over the reactions with source
    e_j + w_k (j in Z, w_k off Z), M_k = Delta_k|_Z e_j^T.  Under the scaling of the class by T, v^{w_k} scales by
    T^{|w_k|}; if, in the lowest (or highest) degree group, every trace(M_k|W) <= 0 with one < 0 and every mixed
    determinant coefficient m(M_k|W, M_l|W) >= 0 with one > 0 (det(M_k|W) = 0: rank one), then for every kappa
    the vertex is linearly asymptotically stable in the classes T*u for T small (large): non-persistent for every
    kappa.  Returns (Z, 'T->0' or 'T->inf') or None."""
    G = CNet(p, rx).Gamma
    nre = len(rx)
    cols = [[Fraction(G[i][j]) for i in range(p)] for j in range(nre)]
    import g6_props
    B = g6_props.basis(cols)
    if len(B) != 2:
        return None
    for size in range(2, p + 1):
        for Z in combinations(range(p), size):
            BZ = [[b[i] for i in Z] for b in B]
            if crn.rank(BZ) != 2:
                continue
            if not all(any(rx[k][0][i] for i in Z) for k in range(nre)):
                continue
            # realizability: some t with t1*BZ0 + t2*BZ1 <= -1 componentwise
            ineqs = [([BZ[0][a], BZ[1][a]], -1) for a in range(len(Z))]
            if e2lib.fm_point(ineqs, 2) is None:
                continue
            # coordinates in the basis BZ of a vector of W
            g11 = sum(x * x for x in BZ[0])
            g12 = sum(x * y for x, y in zip(BZ[0], BZ[1]))
            g22 = sum(x * x for x in BZ[1])
            det = g11 * g22 - g12 * g12

            def coords(vec):
                r1 = sum(x * y for x, y in zip(BZ[0], vec))
                r2 = sum(x * y for x, y in zip(BZ[1], vec))
                return ((g22 * r1 - g12 * r2) / det, (-g12 * r1 + g11 * r2) / det)
            groups = {}
            for k, (a, b) in enumerate(rx):
                if sum(a[i] for i in Z) != 1:
                    continue
                j = next(i for i in Z if a[i] == 1)
                jz = Z.index(j)
                al, be = coords([Fraction(b[i] - a[i]) for i in Z])
                # M_k|W = [al; be] [BZ0_j, BZ1_j]
                Mk = ((al * BZ[0][jz], al * BZ[1][jz]), (be * BZ[0][jz], be * BZ[1][jz]))
                deg = sum(a[i] for i in range(p) if i not in Z)
                groups.setdefault(deg, []).append(Mk)
            if not groups:
                continue
            for lab, deg in (("T->0", min(groups)), ("T->inf", max(groups))):
                Ms = groups[deg]
                tr = [M[0][0] + M[1][1] for M in Ms]
                mix = [A[0][0] * Bm[1][1] + A[1][1] * Bm[0][0] - A[0][1] * Bm[1][0] - A[1][0] * Bm[0][1]
                       for A, Bm in combinations(Ms, 2)]
                if all(t <= 0 for t in tr) and any(t < 0 for t in tr) and all(m >= 0 for m in mix) \
                        and any(m > 0 for m in mix):
                    return (Z, lab)
    return None


def classify(args):
    n, rx, cs_empty, seed, numeric_on = args
    r = len(rx)
    d = e2lib.dual_rx(n, rx)
    out = {"n": n, "rx": rx}
    if cs_empty:
        out["cat"] = "ALL"
        return out
    if endo.drain_species(r, d):
        out["cat"] = "NONE-drain-species"
        return out
    if drainsets.drain_sets(r, d, first_only=True):
        out["cat"] = "NONE-drain-set"
        return out
    vs = vertex_scaling(r, d)
    if vs is not None:
        out["cat"] = "NONE-vertex-scaling"
        out["vs"] = vs
        return out
    ft = facet_types(r, d)
    out["facets"] = ft
    types = {t for Z, t in ft}
    out["cat"] = "UNDECIDED"
    out["ftypes"] = types
    if numeric_on:
        try:
            out["num"] = numeric(r, d, ft, seed)
        except Exception as e:  # noqa: BLE001
            out["num"] = str(e)
    return out


def main():
    t0 = time.time()
    argv = sys.argv[1:]
    limit = int(argv[argv.index("--limit") + 1]) if "--limit" in argv else 100000
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 4
    M = int(argv[argv.index("--numeric") + 1]) if "--numeric" in argv else 400
    print(f"# s2_persistence.py  Python {platform.python_version()} on {platform.system()}; certificates exact; "
          f"numerics scipy LSODA (evidence); limit {limit}/file, numeric on first {M} undecided")
    rng = random.Random(22)
    jobs = []
    for fn in e2lib.E_FILES:
        recs = [rc for rc in e2lib.load([fn]) if rc["s"] == 2 and e2lib.H(rc) and rc["N"]["WR"]]
        rng.shuffle(recs)
        for rc in recs[:limit]:
            jobs.append((rc["n"], rc["rx"], not rc["D"]["CS"], rng.randint(0, 10 ** 9)))
    # dedupe
    seen = set()
    uj = []
    for j in jobs:
        if (j[0], j[1]) not in seen:
            seen.add((j[0], j[1]))
            uj.append(j)
    uj.sort(key=lambda j: (j[0] + len(j[1]), len(j[1]), j[1]))
    print(f"# {len(uj)} WR networks with (H1)-(H3) and s = 2")
    # numeric only for the first M undecided (in size order): two passes
    with Pool(nproc) as pool:
        res = pool.map(classify, [(n, rx, ce, sd, False) for n, rx, ce, sd in uj], chunksize=20)
        und = [q for q in res if q["cat"] == "UNDECIDED"][:M]
        num = pool.map(classify, [(q["n"], q["rx"], False, 7, True) for q in und], chunksize=2)
    numd = {(q["n"], q["rx"]): q.get("num") for q in num}
    cnt = {}
    for q in res:
        cnt[q["cat"]] = cnt.get(q["cat"], 0) + 1
    print(f"# certificates: {dict(sorted(cnt.items()))}")
    ft = {}
    for q in res:
        if q["cat"] == "UNDECIDED":
            key = tuple(sorted(q["ftypes"]))
            ft[key] = ft.get(key, 0) + 1
    print(f"# undecided networks by the set of facet types present: {dict(sorted(ft.items()))}")
    # stable-vertex certificates and numerics (per kappa: 'cert' = exact stable boundary vertex in some scaled class,
    # 'num-ext' = numerical extinction only, 'none' = neither)
    stats = {"cert-all": 0, "cert-some": 0, "cert-none": 0, "error": 0}
    detail = {}
    mixed = []
    out_never = []
    for q in res:
        k = (q["n"], q["rx"])
        if k not in numd:
            continue
        v = numd[k]
        if not isinstance(v, list):
            stats["error"] += 1
            continue
        key = tuple(sorted(v))
        detail[key] = detail.get(key, 0) + 1
        nc = sum(x == "cert" for x in v)
        if nc == len(v):
            stats["cert-all"] += 1
        elif nc == 0:
            stats["cert-none"] += 1
        else:
            stats["cert-some"] += 1
        if "none" in v and any(x != "none" for x in v):
            mixed.append(q)
        if all(x == "none" for x in v) and ("OUT" in q["ftypes"] or "ESC" in q["ftypes"]):
            out_never.append(q)
    print(f"# {len(numd)} undecided networks, 4 sampled kappa each, 7 scaled classes (10^-6..10^6): stable boundary vertex "
          f"(exact) for all kappa {stats['cert-all']}, for some {stats['cert-some']}, for none {stats['cert-none']}, "
          f"errors {stats['error']}")
    print(f"# per-network outcome multisets (cert / num-ext / none): {dict(sorted(detail.items()))}")
    print(f"# OUT/ESC facet present but neither certificate nor numerical extinction for any kappa (against Q1-8): "
          f"{len(out_never)}")
    for q in out_never[:3]:
        print(f"     N = {e2lib.describe(q['n'], q['rx'])}; facets {q['facets']}")
    print(f"# networks with 'none' for some kappa and extinction (cert or numerical) for another (candidates against "
          f"Q1-7, to be examined): {len(mixed)}")
    for q in mixed[:5]:
        print(f"  mixed: N = {e2lib.describe(q['n'], q['rx'])}; facets {q['facets']}; {numd[(q['n'], q['rx'])]}")
    sm = {}
    for q in res:
        if q["cat"] not in sm:
            sm[q["cat"]] = q
    for c, q in sorted(sm.items()):
        print(f"  smallest {c}: N = {e2lib.describe(q['n'], q['rx'])}; N* = "
              f"{e2lib.describe(len(q['rx']), e2lib.dual_rx(q['n'], q['rx']), [f'k{j+1}' for j in range(len(q['rx']))])}")
    print("PASS: (no claim; certificates are proved statements applied exactly)")
    print(f"# runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
