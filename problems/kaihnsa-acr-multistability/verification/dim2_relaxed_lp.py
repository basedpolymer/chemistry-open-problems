"""Adversarial test of the dimension-two theorem on a RELAXED class of vector fields (floating point; a search for
counterexamples, not a proof).

Runtime: about 14 minutes (WSL, Python 3.12, numpy, scipy/HiGHS; single core).
Usage:   python dim2_relaxed_lp.py [samples_per_setting] [seed]

The proof of the dimension-two theorem uses, besides f(S-valued) and deg f <= 2, only the two sign properties of the
coefficient vectors v_y of the monomials x^y (the lemma "coefficient cones"):
   (C1) (v_y)_i >= 0 for i not in supp(y);   (C2) sum_i (v_y)_i <= 0 when |y| = 2.
So the theorem should hold for every polynomial vector field in normal form
   f = psi * alpha + (x1 - c) * P * ell,   alpha_1 = 1, ell_1 = 0,  deg psi <= 2,  deg P <= 1,
satisfying (C1)-(C2): if all positive zeros of f lie in H = {x1 = c} (global ACR), two zeros p, q = p + tau*ell
cannot both be exponentially stable for Df|_S, S = span(alpha, ell).

For fixed (alpha, ell, c, p, q), the LOCAL conditions
   psi(p) = psi(q) = 0;  P(p), P(q) of opposite signs;  det Df|_S > 0 and trace Df|_S < 0 at p and q;  (C1), (C2)
are LINEAR in the coefficients of (psi, P) (det = -P grad(psi).ell, trace = grad(psi).alpha on H).  We sample
(alpha, ell, c, p, q) at random, find feasible (psi, P) by linear programming with random objectives, and then
test the GLOBAL condition: by the theorem, psi must vanish somewhere on Pi^+ or Pi^- (Pi^pm = {P = 0, x > 0,
pm(x1 - c) > 0}), i.e. there must be a positive steady state off H, so ACR fails.  A sample where no such zero is
found would be a counterexample candidate to the proof (it is printed and re-examined with a finer search).

Note: P(p) P(q) < 0 is necessary for two stable states (det formula), so this covers every potential
counterexample with alpha_1 != 0; alpha_1 = 0 is excluded by an eigenvector argument checked by hand.

As a control, the same search is run WITHOUT (C2) (products of any molecularity): the remark of the manuscript
says the argument then breaks down; we report whether the global condition can then hold.
"""
import sys
import time
import itertools
import numpy as np
from scipy.optimize import linprog

T0 = time.time()
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 2026
rng = np.random.default_rng(SEED)


def monomials(n):
    out = []
    for d in range(3):
        for e in itertools.product(range(3), repeat=n):
            if sum(e) == d:
                out.append(e)
    return out


def build(n, alpha, ell, c, p, q, sp_sign, use_c2, stab=True, qsign=1, delta=1e-2, bound=50.0):
    M = monomials(n)
    nm = len(M)
    nv = nm + n + 1                    # psi coefficients, then P = pi0 + sum pi_j x_j
    idx = {e: i for i, e in enumerate(M)}

    def R_coeff(e):
        """coefficient vector (over the variables) of monomial e in (x1 - c) * P."""
        v = np.zeros(nv)
        deg = sum(e)
        if deg == 0:
            v[nm + 0] = -c
        elif deg == 1:
            j = e.index(1)
            v[nm + 1 + j] = -c
            if j == 0:
                v[nm + 0] += 1.0
        else:
            if e[0] >= 1:
                # e = x1 * x_jj (jj = 0 for x1^2): coefficient pi_jj
                rest = list(e)
                rest[0] -= 1
                jj = [k for k in range(n) if rest[k] == 1][0]
                v[nm + 1 + jj] = 1.0
        return v

    def psi_coeff(e):
        v = np.zeros(nv)
        v[idx[e]] = 1.0
        return v

    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    # cone conditions
    for e in M:
        supp = [i for i in range(n) if e[i] > 0]
        rows = [alpha[i] * psi_coeff(e) + ell[i] * R_coeff(e) for i in range(n)]
        for i in range(n):
            if i not in supp:
                A_ub.append(-rows[i]); b_ub.append(0.0)           # rows[i] >= 0
        if sum(e) == 2 and use_c2:
            A_ub.append(sum(rows)); b_ub.append(0.0)              # sum <= 0

    def eval_row(pt):
        return np.array([np.prod([pt[i] ** e[i] for i in range(n)]) for e in M])

    def grad_rows(pt):
        G = np.zeros((n, nm))
        for k, e in enumerate(M):
            for i in range(n):
                if e[i] > 0:
                    ee = list(e); ee[i] -= 1
                    G[i, k] = e[i] * np.prod([pt[j] ** ee[j] for j in range(n)])
        return G

    def P_row(pt):
        v = np.zeros(nv)
        v[nm] = 1.0
        v[nm + 1:] = pt
        return v

    for pt, sgn in ((p, sp_sign), (q, -sp_sign)):
        r = np.zeros(nv); r[:nm] = eval_row(pt)
        A_eq.append(r); b_eq.append(0.0)                          # psi(pt) = 0
        A_ub.append(-sgn * P_row(pt)); b_ub.append(-1.0)          # sgn * P(pt) >= 1
        G = grad_rows(pt)
        gl = np.zeros(nv); gl[:nm] = ell @ G                      # grad psi . ell
        ga = np.zeros(nv); ga[:nm] = alpha @ G                    # grad psi . alpha
        if stab:
            # det = -P * (grad psi . ell) > 0  <=>  sgn * (grad psi . ell) <= -delta
            A_ub.append(sgn * gl); b_ub.append(-delta)
            A_ub.append(ga); b_ub.append(-delta)                  # trace < 0
    if not stab:
        # h_L(t) = psi(p + t ell) must have degree 2: qsign * Q(ell') >= 1, Q = quadratic part of psi
        r = np.zeros(nv)
        for k, e in enumerate(M):
            if sum(e) == 2:
                r[k] = np.prod([ell[i] ** e[i] for i in range(n)])
        A_ub.append(-qsign * r); b_ub.append(-1.0)
    return M, nm, nv, np.array(A_ub), np.array(b_ub), np.array(A_eq), np.array(b_eq), bound


def psi_eval(M, coef, X):
    """X: (N, n) points."""
    out = np.zeros(X.shape[0])
    for k, e in enumerate(M):
        if coef[k] != 0:
            out += coef[k] * np.prod(X ** np.array(e), axis=1)
    return out


def zero_off_H(n, M, psi, P, c, nsamp=40000, fine=False):
    """Search for a zero of psi on Pi^+ or Pi^- (positive points of {P=0} with x1 != c).
    Returns ('zero', point) if psi takes both signs on Pi^+ or on Pi^- (a zero exists by convexity),
    or ('none', info)."""
    pi0, pi = P[0], P[1:]
    ks = [k for k in range(n) if abs(pi[k]) > 1e-9]
    if not ks:
        return ('none', 'P constant')
    res = {}
    for k in ks:
        lo, hi = (-4, 4) if not fine else (-7, 7)
        Xs = np.exp(rng.uniform(lo, hi, size=(nsamp, n)) * np.log(10) / 1.0 * 0.5)
        Xs[:, 0] = np.where(rng.random(nsamp) < 0.5, c * np.exp(rng.uniform(-3, 3, nsamp)),
                            c + rng.uniform(-1, 1, nsamp) * c * 0.999)
        if k == 0:
            Xs[:, 0] = -(pi0 + Xs[:, 1:] @ pi[1:]) / pi[0]
        else:
            others = [j for j in range(n) if j != k]
            Xs[:, k] = -(pi0 + Xs[:, others] @ pi[others]) / pi[k]
        keep = np.all(Xs > 0, axis=1) & (np.abs(Xs[:, 0] - c) > 1e-9)
        Xs = Xs[keep]
        if len(Xs) == 0:
            continue
        v = psi_eval(M, psi, Xs)
        for side, mask in (('+', Xs[:, 0] > c), ('-', Xs[:, 0] < c)):
            vv = v[mask]
            if len(vv) == 0:
                continue
            res.setdefault(side, [np.inf, -np.inf])
            res[side][0] = min(res[side][0], vv.min())
            res[side][1] = max(res[side][1], vv.max())
    for side, (mn, mx) in res.items():
        if mn < 0 < mx:
            return ('zero', side)
    return ('none', res)


def realizable_subspaces(n):
    """2-dimensional subspaces S of R^n spanned by two bimolecular reaction vectors with e1 not in S^perp,
    returned as (ell, alpha0): ell a unit vector spanning S cap e1^perp, alpha0 in S with alpha0_1 = 1."""
    C = [np.array(e) for e in itertools.product(range(3), repeat=n) if sum(e) <= 2]
    D = {tuple(b - a) for a in C for b in C if not np.array_equal(a, b)}
    D = [np.array(d, dtype=float) for d in sorted(D)]
    seen, out = set(), []
    order = rng.permutation(len(D) * len(D)) if n > 3 else range(len(D) * len(D))
    for t in order:
        i, j = divmod(int(t), len(D))
        if i >= j:
            continue
        B = np.vstack([D[i], D[j]])
        if np.linalg.matrix_rank(B) < 2:
            continue
        # canonical key: projector onto S, rounded
        Qm, _ = np.linalg.qr(B.T)
        key = tuple(np.round((Qm @ Qm.T).ravel(), 8))
        if key in seen:
            continue
        seen.add(key)
        # ell spans S cap e1^perp
        a, b = B[0], B[1]
        ell = b[0] * a - a[0] * b if (a[0] != 0 or b[0] != 0) else None
        if ell is None or np.linalg.norm(ell) < 1e-12:
            continue                     # S inside e1^perp: e1 in S^perp, excluded by the lemma
        ell = ell / np.linalg.norm(ell)
        alpha0 = a / a[0] if a[0] != 0 else b / b[0]
        out.append((ell, alpha0))
        if n > 3 and len(out) >= 400:
            break
    return out


SUBSPACES = {}


def run(n, use_c2, ns, stab=True):
    stats = dict(tried=0, feasible=0, acr_fails=0, candidates=0)
    cands = []
    if n not in SUBSPACES:
        SUBSPACES[n] = realizable_subspaces(n)
    subs = SUBSPACES[n]
    for it in range(ns):
        if rng.random() < 0.25:
            # generic (alpha, ell)
            alpha = np.concatenate([[1.0], rng.choice([-1, 1], n - 1) * rng.uniform(0, 3, n - 1)])
            if rng.random() < 0.3:          # some zero entries
                alpha[1 + rng.integers(n - 1)] = 0.0
            if rng.random() < 0.7 and np.any(alpha[1:] < 0):
                # bias towards sum(alpha) <= 0, without which no quadratic monomial x_j x_k (j,k >= 2) is allowed
                neg = alpha < 0
                target = -(1.0 + alpha[~neg].sum()) - rng.uniform(0, 1.5)
                alpha[neg] *= target / alpha[neg].sum()
            ell = np.concatenate([[0.0], rng.normal(size=n - 1)])
            ell /= np.linalg.norm(ell)
        else:
            # S spanned by bimolecular reaction vectors; alpha = alpha0 + lam * ell, with lam random or special
            if rng.random() < 0.4:
                # molecularity-preserving subspaces (S inside {sum = 0}), where (C2) holds with equality
                pool = [sa for sa in subs if abs(sa[0].sum()) < 1e-9 and abs(sa[1].sum()) < 1e-9] or subs
            else:
                pool = subs
            ell, alpha0 = pool[rng.integers(len(pool))]
            specials = [-alpha0[i] / ell[i] for i in range(1, n) if abs(ell[i]) > 1e-12]
            if abs(ell[1:].sum()) > 1e-12:
                specials.append(-(alpha0.sum()) / ell.sum())
            u = rng.random()
            if specials and u < 0.4:
                lam = rng.choice(specials)
            elif len(specials) >= 2 and u < 0.8:
                lam = rng.uniform(min(specials), max(specials))
            else:
                lam = rng.normal() * 2
            alpha = alpha0 + lam * ell
        for rep in range(4):
            c = float(np.exp(rng.uniform(-1, 1)))
            p = np.concatenate([[c], np.exp(rng.uniform(-1.5, 1.5, n - 1))])
            taus = []
            for _ in range(50):
                tau = rng.normal() * 2
                if np.all(p + tau * ell > 0) and abs(tau) > 0.05:
                    taus.append(tau)
            if not taus:
                continue
            q = p + taus[0] * ell
            for sp_sign in (1, -1):
                stats['tried'] += 1
                M, nm, nv, A_ub, b_ub, A_eq, b_eq, bound = build(n, alpha, ell, c, p, q, sp_sign, use_c2, stab,
                                                                 rng.choice([-1, 1]))
                for trial in range(2):
                    obj = rng.normal(size=nv)
                    r = linprog(obj, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(-bound, bound)] * nv,
                                method='highs')
                    if r.status != 0:
                        break
                    stats['feasible'] += 1
                    sol = r.x
                    psi, P = sol[:nm], sol[nm:]
                    verdict, info = zero_off_H(n, M, psi, P, c)
                    if verdict != 'zero':
                        verdict, info = zero_off_H(n, M, psi, P, c, nsamp=400000, fine=True)
                    if verdict == 'zero':
                        stats['acr_fails'] += 1
                    else:
                        stats['candidates'] += 1
                        cands.append(dict(alpha=alpha, ell=ell, c=c, p=p, q=q, psi=psi, P=P, info=info))
    return stats, cands


FAIL = []
MODES = [("theorem", True, True), ("global", True, False), ("control", False, True)]
for n in (3, 4):
    for mode, use_c2, stab in MODES:
        st, cands = run(n, use_c2, NS, stab)
        label = "n=%d, mode %s" % (n, mode)
        print("%s: %s  (%.0f s)" % (label, st, time.time() - T0), flush=True)
        if mode == "theorem":
            ok = st['candidates'] == 0
            print(("PASS " if ok else "FAIL ") + "%s: no configuration with two exponentially stable zeros on a class "
                  "line and global ACR (%d locally feasible LP solutions, all violating global ACR)" % (
                      label, st['feasible']))
        elif mode == "global":
            ok = st['candidates'] == 0 and st['feasible'] > 0
            print(("PASS " if ok else "FAIL ") + "%s: whenever P changes sign between two zeros p, q of a class line "
                  "(cones with (C2), deg h_L = 2), a positive zero off H exists, i.e. ACR fails "
                  "(%d feasible LP solutions)" % (label, st['feasible']))
        else:
            ok = True
            print("     control without (C2): %d feasible, %d with a zero off H found, %d without" % (
                st['feasible'], st['acr_fails'], st['candidates']))
        if not ok:
            FAIL.append(label)
        for cd in cands[:2]:
            print("   configuration where no zero off H was found:",
                  {k: np.round(v, 5) if isinstance(v, np.ndarray) else v for k, v in cd.items()})
print("runtime %.0f s" % (time.time() - T0))
print("ALL PASS" if not FAIL else "FAILURES")
sys.exit(1 if FAIL else 0)
