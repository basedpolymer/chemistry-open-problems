"""Network-level test of the dimension-two theorem (floating point; a search for counterexamples, not a proof).

Runtime: about 15 minutes on 11 cores (WSL, Python 3.12, numpy, scipy).
Usage:   python dim2_network_scan.py [nproc] [n5_sample] [seed]

Enumerates at-most-bimolecular networks on 3 species X1, X2, X3 whose stoichiometric subspace S has dimension
exactly 2: S = ker(w) for one of the 61 primitive normal vectors w spanned by bimolecular reaction vectors, all
networks with 3 or 4 reactions inside ker(w) (up to permutation of species), and a random sample of networks with
5 reactions.  Filters (all justified by standard results): a positive flux vector exists (otherwise there is no
positive steady state), at least 3 distinct reactant complexes (KNS Thm 4.9(1)).

For each network and 3 random rate vectors (log-uniform in [1/10, 10]):
  * the positive steady states in 8 random compatibility classes {w.x = T} are computed by a resultant
    (two quadratics in two unknowns -> quartic) and polished by Newton's method;
  * if some species has the same value at all of them (relative spread < 1e-7, at least 2 classes with steady
    states), the system is flagged as an ACR candidate and examined in 60 further classes: ACR is re-checked, and in
    every class the nondegenerate positive steady states are classified by the eigenvalues of Df restricted to S.
Claims (PASS/FAIL):
  * KNS Example 4.10 is found (sanity check of the search);
  * no ACR system has a class with two exponentially stable positive steady states (the theorem);
  * whenever an ACR system has two nondegenerate positive steady states in a class, their determinants of Df|_S have
    opposite signs (what the proof actually shows: one of the two is a saddle).
"""
import itertools
import sys
import time
from multiprocessing import Pool
import numpy as np
from scipy.optimize import linprog

T0 = time.time()
NPROC = int(sys.argv[1]) if len(sys.argv) > 1 else 12
N5 = int(sys.argv[2]) if len(sys.argv) > 2 else 60000
ONLY3 = len(sys.argv) > 4 and sys.argv[4] == 'only3'
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 7

COMPLEXES = [c for c in itertools.product(range(3), repeat=3) if sum(c) <= 2]
REACTIONS = [(a, b) for a in COMPLEXES for b in COMPLEXES if a != b]
PERMS = list(itertools.permutations(range(3)))
MONOS = COMPLEXES                      # monomials of degree <= 2 = the complexes


def vec(r):
    return tuple(r[1][i] - r[0][i] for i in range(3))


def normals():
    D = sorted({vec(r) for r in REACTIONS})
    W = set()
    for a, b in itertools.combinations(D, 2):
        w = np.cross(a, b)
        if not w.any():
            continue
        w = w // abs(np.gcd.reduce(w))
        if w[np.nonzero(w)[0][0]] < 0:
            w = -w
        W.add(tuple(int(t) for t in w))
    return sorted(W)


def perm_c(c, p):
    return tuple(c[p[i]] for i in range(3))


def canon(net):
    return min(tuple(sorted((perm_c(a, p), perm_c(b, p)) for a, b in net)) for p in PERMS)


def networks():
    seen = set()
    out3, out4, pool5 = [], [], []
    for w in normals():
        Rw = [r for r in REACTIONS if np.dot(w, vec(r)) == 0]
        for k in (3, 4):
            for net in itertools.combinations(Rw, k):
                V = np.array([vec(r) for r in net])
                if np.linalg.matrix_rank(V) != 2:
                    continue
                key = canon(net)
                if key in seen:
                    continue
                seen.add(key)
                (out3 if k == 3 else out4).append((w, key))
        pool5.append((w, Rw))
    return out3, out4, pool5


def sample5(pool5, nsamp, rng):
    seen, out = set(), []
    weights = np.array([len(R) ** 5 for w, R in pool5], dtype=float)
    weights /= weights.sum()
    tries = 0
    while len(out) < nsamp and tries < 20 * nsamp:
        tries += 1
        w, Rw = pool5[rng.choice(len(pool5), p=weights)]
        if len(Rw) < 5:
            continue
        idx = rng.choice(len(Rw), 5, replace=False)
        net = [Rw[i] for i in idx]
        V = np.array([vec(r) for r in net])
        if np.linalg.matrix_rank(V) != 2:
            continue
        key = canon(net)
        if key in seen:
            continue
        seen.add(key)
        out.append((w, key))
    return out


def positive_flux(net):
    G = np.array([vec(r) for r in net], dtype=float).T
    r = len(net)
    res = linprog(np.zeros(r), A_eq=G, b_eq=np.zeros(3), bounds=[(1, None)] * r, method='highs')
    return res.status == 0


# ---------------------------------------------------------------- polynomial machinery
def field_coeffs(net, k):
    """F[i, mono] = coefficient of monomial in f_i."""
    F = np.zeros((3, len(MONOS)))
    for (a, b), kk in zip(net, k):
        j = MONOS.index(a)
        for i in range(3):
            F[i, j] += kk * (b[i] - a[i])
    return F


def f_eval(F, x):
    m = np.array([x[0] ** e[0] * x[1] ** e[1] * x[2] ** e[2] for e in MONOS])
    return F @ m


def jac(F, x):
    J = np.zeros((3, 3))
    for j, e in enumerate(MONOS):
        for v in range(3):
            if e[v] > 0:
                ee = list(e)
                ee[v] -= 1
                J[:, v] += F[:, j] * e[v] * x[0] ** ee[0] * x[1] ** ee[1] * x[2] ** ee[2]
    return J


def affine_product(L1, L2):
    """product of affine forms (c0, cu, cv) -> quadratic coefficients dict {(i,j): coef} in u^i v^j."""
    out = {}
    for (i1, c1) in ((0, L1[0]), (1, L1[1]), (2, L1[2])):
        for (i2, c2) in ((0, L2[0]), (1, L2[1]), (2, L2[2])):
            e = [0, 0]
            for t in (i1, i2):
                if t == 1:
                    e[0] += 1
                elif t == 2:
                    e[1] += 1
            out[tuple(e)] = out.get(tuple(e), 0.0) + c1 * c2
    return out


def class_polys(F, w, T, m, ab):
    """restrict f_a, f_b to the class w.x = T, with x_m eliminated; return 2 dicts {(i,j): coef}."""
    a, b = ab
    L = [None] * 3
    L[a] = (0.0, 1.0, 0.0)
    L[b] = (0.0, 0.0, 1.0)
    L[m] = (T / w[m], -w[a] / w[m], -w[b] / w[m])
    one = (1.0, 0.0, 0.0)
    polys = []
    for i in ab:
        g = {}
        for j, e in enumerate(MONOS):
            if F[i, j] == 0:
                continue
            fs = [L[v] for v in range(3) for _ in range(e[v])]
            if len(fs) == 0:
                prod = {(0, 0): 1.0}
            elif len(fs) == 1:
                prod = affine_product(fs[0], one)
            else:
                prod = affine_product(fs[0], fs[1])
            for key, cc in prod.items():
                g[key] = g.get(key, 0.0) + F[i, j] * cc
        polys.append(g)
    return polys


def P(g, key):
    return g.get(key, 0.0)


def solve_class(F, w, T):
    m = int(np.argmax(np.abs(w)))
    ab = tuple(i for i in range(3) if i != m)
    g1, g2 = class_polys(F, w, T, m, ab)
    # as polynomials in v: A v^2 + B(u) v + C(u); numpy poly1d in u
    def parts(g):
        A = np.poly1d([P(g, (0, 2))])
        B = np.poly1d([P(g, (1, 1)), P(g, (0, 1))])
        C = np.poly1d([P(g, (2, 0)), P(g, (1, 0)), P(g, (0, 0))])
        return A, B, C
    A1, B1, C1 = parts(g1)
    A2, B2, C2 = parts(g2)
    if abs(A1.c[0]) < 1e-14 and abs(A2.c[0]) < 1e-14:
        R = B1 * C2 - B2 * C1
    else:
        R = (A1 * C2 - A2 * C1) ** 2 - (A1 * B2 - A2 * B1) * (B1 * C2 - B2 * C1)
    scale = max(1e-300, max(abs(R.c)))
    if max(abs(R.c)) < 1e-10 * (1 + max(abs(np.concatenate([A1.c, B1.c, C1.c, A2.c, B2.c, C2.c]))) ** 4):
        return None                      # resultant vanishes: continuum or common factor
    us = [z.real for z in np.roots(R.c) if abs(z.imag) < 1e-7 * (1 + abs(z)) and z.real > 0]
    sols = []
    for u in us:
        cand_v = []
        for (A, B, C) in ((A1, B1, C1), (A2, B2, C2)):
            co = [A(u), B(u), C(u)]
            while co and abs(co[0]) < 1e-14:
                co = co[1:]
            if len(co) >= 2:
                cand_v += [z.real for z in np.roots(co) if abs(z.imag) < 1e-6 * (1 + abs(z))]
        for v in cand_v:
            if v <= 0:
                continue
            x = np.zeros(3)
            x[ab[0]], x[ab[1]] = u, v
            x[m] = (T - w[ab[0]] * u - w[ab[1]] * v) / w[m]
            if x[m] <= 0:
                continue
            # Newton polish on (f_a, f_b) inside the class
            Bm = np.zeros((3, 2))
            Bm[ab[0], 0] = 1
            Bm[ab[1], 1] = 1
            Bm[m, 0] = -w[ab[0]] / w[m]
            Bm[m, 1] = -w[ab[1]] / w[m]
            ok = True
            for _ in range(30):
                fx = f_eval(F, x)[list(ab)]
                Jx = jac(F, x)[list(ab)] @ Bm
                try:
                    dx = np.linalg.solve(Jx, -fx)
                except np.linalg.LinAlgError:
                    break
                x = x + Bm @ dx
                if np.any(x <= 0):
                    ok = False
                    break
                if np.linalg.norm(dx) < 1e-14 * (1 + np.linalg.norm(x)):
                    break
            if not ok:
                continue
            if np.linalg.norm(f_eval(F, x)) > 1e-8 * (1 + np.abs(F).max() * (1 + np.linalg.norm(x)) ** 2):
                continue
            if x.min() < 1e-9 * max(1.0, x.max()):
                continue                 # numerically on the boundary, not a positive steady state
            if all(np.linalg.norm(x - y) > 1e-7 * (1 + np.linalg.norm(x)) for y in sols):
                sols.append(x)
    return sols


def S_basis(w):
    w = np.array(w, dtype=float)
    _, _, Vt = np.linalg.svd(w[None, :])
    return Vt[1:].T                      # 3 x 2 orthonormal basis of ker w


def classify(F, w, x):
    B = S_basis(w)
    JS = B.T @ jac(F, x) @ B
    ev = np.linalg.eigvals(JS)
    d = np.linalg.det(JS)
    scale = 1 + np.abs(JS).max() ** 2
    if abs(d) < 1e-9 * scale:
        return 'degenerate', d
    if np.all(ev.real < -1e-10 * (1 + np.abs(JS).max())):
        return 'stable', d
    return 'unstable', d


def random_class(w, rng):
    x0 = np.exp(rng.uniform(np.log(0.05), np.log(20), 3))
    return float(np.dot(w, x0))


def normal_of(net):
    """primitive integer normal vector of the (2-dimensional) span of the reaction vectors of net."""
    V = [np.array(vec(r)) for r in net]
    for a, b in itertools.combinations(V, 2):
        w = np.cross(a, b)
        if w.any():
            return w // abs(np.gcd.reduce(w))
    raise ValueError("rank < 2")


def test_network(args):
    _, net, seed = args
    rng = np.random.default_rng(seed)
    # the normal vector is recomputed from the (canonically permuted) network itself
    w = normal_of(net)
    assert all(np.dot(w, vec(r)) == 0 for r in net)
    res = dict(acr=0, acr_multi=0, two_stable=0, opposite_ok=0, opposite_bad=0, examples=[])
    for trial in range(3):
        k = np.exp(rng.uniform(np.log(0.1), np.log(10), len(net)))
        F = field_coeffs(net, k)
        pts = []
        ncls = 0
        for _ in range(8):
            T = random_class(w, rng)
            s = solve_class(F, w, T)
            if s is None:
                pts = None
                break
            if s:
                ncls += 1
                pts += s
        if pts is None or ncls < 2:
            continue
        X = np.array(pts)
        acr_sp = [i for i in range(3) if (X[:, i].max() - X[:, i].min()) <= 1e-7 * X[:, i].mean()]
        if not acr_sp:
            continue
        # closer look; the class data count only if ACR is confirmed
        allpts, multi, still_acr = [], False, True
        loc = dict(two_stable=0, opposite_ok=0, opposite_bad=0)
        loc_ex = []
        for _ in range(60):
            T = random_class(w, rng)
            s = solve_class(F, w, T)
            if s is None:
                still_acr = False
                break
            allpts += s
            kinds = [classify(F, w, x) for x in s]
            nd = [(kd, d) for kd, d in kinds if kd != 'degenerate']
            if len(nd) >= 2:
                multi = True
                if sum(1 for kd, d in nd if kd == 'stable') >= 2:
                    loc['two_stable'] += 1
                    loc_ex.append(('two_stable', net, list(k), T))
                if len(nd) == 2:
                    if np.sign(nd[0][1]) != np.sign(nd[1][1]):
                        loc['opposite_ok'] += 1
                    else:
                        loc['opposite_bad'] += 1
                        loc_ex.append(('same_sign', net, list(k), T))
        if not still_acr or not allpts:
            continue
        X = np.array(allpts + pts)
        acr_sp = [i for i in acr_sp if (X[:, i].max() - X[:, i].min()) <= 1e-7 * X[:, i].mean()]
        if not acr_sp:
            continue
        res['acr'] += 1
        for key in loc:
            res[key] += loc[key]
        res['examples'] += loc_ex
        if multi:
            res['acr_multi'] += 1
            res['examples'].append(('acr_multi', net, list(k), None))
    return res


def main():
    rng = np.random.default_rng(SEED)
    out3, out4, pool5 = networks()
    out5 = sample5(pool5, N5, rng)
    allnets = [(w, net) for w, net in (out3 if ONLY3 else out3 + out4 + out5)]
    print("networks up to species permutation: 3 reactions %d, 4 reactions %d, 5 reactions (sample) %d  (%.0f s)" % (
        len(out3), len(out4), len(out5), time.time() - T0), flush=True)
    kept = [(w, net) for w, net in allnets
            if len({a for a, b in net}) >= 3 and positive_flux(net)]
    print("after filters (>= 3 reactant complexes, positive flux): %d  (%.0f s)" % (len(kept), time.time() - T0),
          flush=True)
    jobs = [(w, net, SEED * 1000003 + i) for i, (w, net) in enumerate(kept)]
    tot = dict(acr=0, acr_multi=0, two_stable=0, opposite_ok=0, opposite_bad=0)
    examples = []
    acr_nets = []
    with Pool(NPROC) as pool:
        for (w, net, _), r in zip(jobs, pool.imap(test_network, jobs, chunksize=50)):
            for key in tot:
                tot[key] += r[key]
            examples += r['examples']
            if r['acr_multi']:
                acr_nets.append(net)
    print("totals over (network, rate vector) pairs:", tot, " (%.0f s)" % (time.time() - T0), flush=True)
    print("ACR + nondegenerate multistationarity found in %d networks, e.g.:" % len(acr_nets))
    for net in acr_nets[:12]:
        print("   ", " , ".join("%s->%s" % (a, b) for a, b in net))
    FAIL = []
    kns = canon([((1, 1, 0), (0, 0, 2)), ((0, 0, 1), (1, 0, 0)), ((0, 0, 2), (0, 2, 0))])
    found = kns in {canon(n) for n in acr_nets}
    print(("PASS " if found else "FAIL ") + "KNS Example 4.10 is detected as ACR + nondegenerately multistationary")
    if not found:
        FAIL.append("kns")
    ok = tot['two_stable'] == 0
    print(("PASS " if ok else "FAIL ") + "no ACR system with two exponentially stable positive steady states in a class")
    if not ok:
        FAIL.append("two stable")
        for e in examples[:5]:
            print("   ", e)
    ok = tot['opposite_bad'] == 0 and tot['opposite_ok'] > 0
    print(("PASS " if ok else "FAIL ") + "two nondegenerate positive steady states in a class of an ACR system always "
          "have determinants of opposite signs (%d classes)" % tot['opposite_ok'])
    if not ok:
        FAIL.append("opposite")
        for e in examples[:5]:
            print("   ", e)
    print("runtime %.0f s" % (time.time() - T0))
    print("ALL PASS" if not FAIL else "FAILURES")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
