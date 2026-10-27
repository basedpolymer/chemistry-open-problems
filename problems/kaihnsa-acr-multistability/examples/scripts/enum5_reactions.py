#!/usr/bin/env python3
"""
enum5_reactions.py -- exhaustive enumeration of 3-species, at-most-bimolecular networks with exactly 5 reactions
(5 distinct reactant complexes; full-dimensional networks need at least 5 = n+2 reactant complexes, KNS Thm 5.1).

Runtime: about 2 minutes on 12 cores (WSL, Python 3.12, sympy 1.14; measured 112 s).  Usage:  python enum5_reactions.py [nproc]

All decisions are exact (integer / rational arithmetic, sympy exact real-root counting).  For a network with
stoichiometric matrix Gamma (3x5, rank 3) and reactant matrix Y (3x5), the theory of report.md, Section 2, gives:

  (a) positive steady states exist for some kappa  iff  ker Gamma meets R^5_{>0};
  (b) if rank [Y; 1^T] < 4, every positive steady state is degenerate;
  (c) otherwise positive steady states <-> roots s in (0,1) of g(s) = K, g(s) = sum_i c_i log mu_i(s), where
      mu(s) = (1-s) r1 + s r2 parametrizes the cone ker Gamma cap R^5_{>0} (r1, r2 its extreme rays) and c is the
      affine dependence of the reactant complexes; sign det J = (fixed sign) * sign g'(s);
  (d) with P = number of distinct zeros (in P^1) of the linear forms lambda_i on ker Gamma:
        P <= 3  -> at most 2 positive steady states, so no multistability;
        P  = 4  -> ACR together with >= 3 positive steady states forces the two reactions of the unique pole class
                   {a,b} to satisfy y_a - y_b = +-e_j (then ACR in X_j holds for every kappa);
        P  = 5  -> ACR with >= 3 positive steady states needs a triple point of the curve (g, h_j); such networks
                   would be written to enum5_out/p5_candidates.json (none occur).
  The decisive count is simpler: for EVERY consistent network the number of sign changes of g' on (0,1) is
  computed exactly (P <= 3 gives at most one by degree); it never reaches 2, so no full-dimensional 5-reaction
  network has three nondegenerate positive steady states, with or without ACR.
  (e) multistability needs >= 3 positive steady states with two of them having det J < 0, hence g' must change
      sign at least twice on (0,1).

Output: counts per stage, and the lists of surviving networks (P=4 structural candidates; P=5 candidates) in
enum5_out/.  Prints PASS lines for the summary claims; exits non-zero if the enumeration is inconsistent.
"""
import itertools, sys, os, time, json
from fractions import Fraction
from multiprocessing import Pool
import sympy as sp

T0 = time.time()
COMPLEXES = [c for c in itertools.product(range(3), repeat=3) if sum(c) <= 2]
PERMS = list(itertools.permutations(range(3)))
S_ = sp.Symbol('s')


def pc(c, p):
    return tuple(c[p[i]] for i in range(3))


def canonical_subsets():
    out = []
    for sub in itertools.combinations(COMPLEXES, 5):
        key = tuple(sorted(sub))
        can = min(tuple(sorted(pc(c, p) for c in sub)) for p in PERMS)
        if key == can:
            out.append(key)
    return out


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def frac_rank(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    ncol = len(M[0])
    for col in range(ncol):
        piv = None
        for r in range(rk, len(M)):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for r in range(len(M)):
            if r != rk and M[r][col] != 0:
                f = M[r][col] / M[rk][col]
                M[r] = [M[r][k] - f * M[rk][k] for k in range(ncol)]
        rk += 1
    return rk


def nullspace(rows, ncol):
    M = sp.Matrix(rows)
    ns = M.nullspace()
    out = []
    for v in ns:
        den = sp.ilcm(*[sp.fraction(x)[1] for x in v])
        out.append([int(x * den) for x in v])
    return out


def solve_affine(Y, j):
    """d with Y d = e_j and sum d = 0 (rational)."""
    A = sp.Matrix([list(Y[0]), list(Y[1]), list(Y[2]), [1] * 5])
    b = sp.Matrix([1 if i == j else 0 for i in range(3)] + [0])
    sol, params = A.gauss_jordan_solve(b)
    sol = sol.subs({p: 0 for p in params})
    return [sp.Rational(x) for x in sol]


def reactant_data(R):
    Y = [[R[i][k] for i in range(5)] for k in range(3)]
    rkY = frac_rank(Y)
    rkY1 = frac_rank(Y + [[1] * 5])
    if rkY1 < 4:
        return dict(Y=Y, ok=False, rkY=rkY, rkY1=rkY1)
    c = nullspace(Y + [[1] * 5], 5)
    assert len(c) == 1
    c = c[0]
    ds = [solve_affine(Y, j) for j in range(3)]
    return dict(Y=Y, ok=True, c=c, d=ds)


def kernel_basis(G):
    """G: list of 5 columns (3-vectors). Returns 5x2 integer kernel basis (rows K_i) or None if rank < 3."""
    for S in itertools.combinations(range(5), 3):
        D = det3(G[S[0]], G[S[1]], G[S[2]])
        if D != 0:
            break
    else:
        return None
    a, b, cc = G[S[0]], G[S[1]], G[S[2]]
    out = [[0, 0] for _ in range(5)]
    others = [i for i in range(5) if i not in S]
    for col, i in enumerate(others):
        gi = G[i]
        # Cramer: solve [a b c] z = -gi * D  -> z_k = -det(replace k by gi)
        z0 = -det3(gi, b, cc)
        z1 = -det3(a, gi, cc)
        z2 = -det3(a, b, gi)
        out[S[0]][col] = z0
        out[S[1]][col] = z1
        out[S[2]][col] = z2
        out[i][col] = D
    return out


def cone_rays(K):
    """Return (r1, r2) extreme rays of {theta : K theta >= 0} if its interior is nonempty, else None."""
    cands = []
    for k in K:
        rho = (-k[1], k[0])
        for r in (rho, (-rho[0], -rho[1])):
            if r == (0, 0):
                continue
            if all(kk[0] * r[0] + kk[1] * r[1] >= 0 for kk in K):
                # dedupe by direction
                if not any(r[0] * q[1] - r[1] * q[0] == 0 and r[0] * q[0] + r[1] * q[1] > 0 for q in cands):
                    cands.append(r)
    if len(cands) != 2:
        return None
    r1, r2 = cands
    th = (r1[0] + r2[0], r1[1] + r2[1])
    if all(kk[0] * th[0] + kk[1] * th[1] > 0 for kk in K):
        return r1, r2
    return None


def pole_classes(K):
    cls = []
    for i, k in enumerate(K):
        for cl in cls:
            q = K[cl[0]]
            if k[0] * q[1] - k[1] * q[0] == 0:
                cl.append(i)
                break
        else:
            cls.append([i])
    return cls


def sign_changes_01(K, r1, r2, w, cls):
    """Number of roots of odd multiplicity in (0,1) of the numerator of sum_i w_i mu_i'/mu_i."""
    u = [k[0] * r1[0] + k[1] * r1[1] for k in K]
    v = [k[0] * r2[0] + k[1] * r2[1] - u[i] for i, k in enumerate(K)]
    reps = [cl[0] for cl in cls]
    W = [sum(sp.Rational(w[i]) for i in cl) for cl in cls]
    mus = [u[q] + S_ * v[q] for q in reps]
    N = 0
    for a, q in enumerate(reps):
        term = W[a] * v[q]
        for b, q2 in enumerate(reps):
            if b != a:
                term = term * mus[b]
        N += term
    N = sp.expand(N)
    if N == 0:
        return 0, N
    P = sp.Poly(N, S_)
    cnt = 0
    for fac, mult in P.sqf_list()[1]:
        if mult % 2 == 1:
            n = fac.count_roots(0, 1)
            n -= (1 if fac.eval(0) == 0 else 0) + (1 if fac.eval(1) == 0 else 0)
            cnt += n
    return cnt, N


def work(R):
    rd = reactant_data(R)
    stats = dict(total=0, rank3=0, consistent=0, degenerateY=0, P=[0] * 6, P4_pair=0, P4_pair_sc2=0,
                 P5_sc2=0, sc2_any=0, sc_hist=[0, 0, 0, 0])
    p4 = []
    many = []
    p5 = []
    prods = [[c for c in COMPLEXES if c != R[i]] for i in range(5)]
    for P_ in itertools.product(*prods):
        stats['total'] += 1
        G = [tuple(P_[i][k] - R[i][k] for k in range(3)) for i in range(5)]
        # quick necessary condition for consistency: each species row has + and - entries
        okq = True
        for k in range(3):
            row = [g[k] for g in G]
            if not (max(row) > 0 and min(row) < 0):
                okq = False
                break
        if not okq:
            continue
        K = kernel_basis(G)
        if K is None:
            continue
        stats['rank3'] += 1
        rays = cone_rays(K)
        if rays is None:
            continue
        stats['consistent'] += 1
        if not rd['ok']:
            stats['degenerateY'] += 1
            continue
        cls = pole_classes(K)
        Pn = len(cls)
        stats['P'][Pn] += 1
        if Pn <= 3:
            continue
        r1, r2 = rays
        net = [[list(R[i]), list(P_[i])] for i in range(5)]
        sc, N = sign_changes_01(K, r1, r2, rd['c'], cls)
        stats['sc_hist'][min(sc, 3)] += 1
        if sc >= 2:
            stats['sc2_any'] += 1
            many.append(dict(net=net, P=Pn, sc=sc, N=str(N)))
        if Pn == 4:
            pair = [cl for cl in cls if len(cl) == 2][0]
            a, b = pair
            diff = [R[a][k] - R[b][k] for k in range(3)]
            js = [j for j in range(3) if sorted(map(abs, diff)) == [0, 0, 1] and diff[j] != 0]
            if not js:
                continue
            stats['P4_pair'] += 1
            if sc >= 2:
                stats['P4_pair_sc2'] += 1
                p4.append(dict(net=net, K=K, r1=r1, r2=r2, pair=pair, acr_species=js[0], N=str(N)))
        else:
            if sc >= 2:
                stats['P5_sc2'] += 1
                p5.append(dict(net=net, K=K, r1=r1, r2=r2, c=rd['c'], d=[[str(x) for x in dj] for dj in rd['d']],
                               N=str(N)))
    return R, stats, p4, p5, many


def main():
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    subs = canonical_subsets()
    print("canonical reactant 5-subsets:", len(subs), flush=True)
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "enum5_out")
    os.makedirs(outdir, exist_ok=True)
    tot = dict(total=0, rank3=0, consistent=0, degenerateY=0, P=[0] * 6, P4_pair=0, P4_pair_sc2=0, P5_sc2=0,
               sc2_any=0, sc_hist=[0, 0, 0, 0])
    P4, P5, MANY = [], [], []
    with Pool(nproc) as pool:
        for R, st, p4, p5, many in pool.imap_unordered(work, subs):
            MANY += many
            for k in tot:
                if k in ('P', 'sc_hist'):
                    tot[k] = [x + y for x, y in zip(tot[k], st[k])]
                else:
                    tot[k] += st[k]
            P4 += p4
            P5 += p5
    print(json.dumps(tot))
    with open(os.path.join(outdir, "p4_structural.json"), "w") as f:
        json.dump(P4, f)
    with open(os.path.join(outdir, "p5_candidates.json"), "w") as f:
        json.dump(P5, f)
    with open(os.path.join(outdir, "three_or_more.json"), "w") as f:
        json.dump(MANY, f)
    ok = tot['total'] == len(subs) * 9 ** 5
    print("PASS" if ok else "FAIL", "enumerated", tot['total'], "networks (", len(subs), "reactant sets x 9^5 )")
    print("P4 structural-ACR candidates with >=2 sign changes of g':", len(P4))
    print("P5 candidates with >=2 sign changes of g':", len(P5))
    print("networks (any P>=4) whose g' changes sign >= 2 times on the cone (needed for 3 nondegenerate positive "
          "steady states):", len(MANY), " histogram of sign changes (0,1,2,>=3):", tot['sc_hist'])
    claim = (len(P4) == 0 and len(P5) == 0)
    print("PASS" if claim else "FAIL", "no full-dimensional 5-reaction network has ACR together with >= 3 nondegenerate positive steady "
          "states for some kappa (hence none has ACR + multistability)")
    claim2 = (len(MANY) == 0)
    print("PASS" if claim2 else "INFO", "no full-dimensional 5-reaction network has >= 3 nondegenerate positive steady states at all"
          if claim2 else "some 5-reaction networks have >= 3 nondegenerate positive steady states")
    ok = ok and claim
    print("runtime %.1f s" % (time.time() - T0))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
