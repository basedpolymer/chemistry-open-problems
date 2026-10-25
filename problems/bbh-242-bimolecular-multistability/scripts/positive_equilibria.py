#!/usr/bin/env python3
"""Certificate for Proposition prop:four(2) of paper.tex (positive equilibria of four-source networks).

For every dynamical-equivalence class of dynamically nontrivial bimolecular (2,4,2) networks with FOUR
distinct sources (838 classes), exact rational computation of:

  * the extreme rays u, v of ker(Gamma) cap R^4_{>=0}  (h(t) = (1-t) u + t v, t in I = (0,1));
  * the affine relation c (left kernel of [A | 1]), when the sources are affinely spanning;
  * Q(t) = det(Gamma diag(h(t)) A), the sign of the Jacobian determinant at the positive equilibrium
    attached to t (det J = lambda^2 Q(t) / (x y));
  * N(t) = (prod_i h_i(t)) * Phi'(t) with Phi(t) = sum_i c_i ln h_i(t)  (positive equilibria are the roots
    of Phi(t) = c . ln(kappa), BBH 2024, Lemma 11);

and checks:
  C1  the four sources are affinely spanning (rank [A|1] = 3)  or else flags the class as "degenerate";
  C2  Q * N has no sign change on I (so Q and Phi' change sign at the same points of I);
  C3  Q has at most one sign change on I.
Consequence (paper.tex, Corollary cor:pos): for every kappa at most one positive equilibrium has det J > 0, hence
at most one positive equilibrium is asymptotically stable; classes failing C1 have no isolated positive
equilibrium.

Also lists the classes for which Q (equivalently Phi') changes sign on I: these are exactly the classes with
two positive equilibria for some kappa (the 30 networks of Table 1 of BBH 2024 must appear).

Runtime: about 10 s.  Environment: Python 3.12, sympy 1.14 (exact rationals), WSL Fedora 44.
Prints PASS/FAIL per check; exit status 1 on failure.
"""
import sys
from itertools import combinations
from math import gcd
import sympy as sp
from networks import *

t = sp.symbols('t')
FAILS = []


def kernel_cone_rays(G):
    """extreme rays of ker(G) cap R^m_{>=0} for a 2 x m integer matrix of rank 2 (m=4), assuming a
    strictly positive kernel vector exists. Returns two primitive integer vectors."""
    M = sp.Matrix(G)
    K = M.nullspace()
    assert len(K) == 2
    b1, b2 = [sp.Matrix(k) for k in K]
    m = M.shape[1]
    cands = []
    for i in range(m):
        # vectors s b1 + r b2 with coordinate i equal to 0
        a, b = b1[i], b2[i]
        if a == 0 and b == 0:
            continue
        w = b * b1 - a * b2
        w = w * sp.lcm([sp.fraction(x)[1] for x in w])
        wi = [int(x) for x in w]
        g = 0
        for x in wi:
            g = gcd(g, abs(x))
        wi = [x // g for x in wi]
        for s in (wi, [-x for x in wi]):
            if all(x >= 0 for x in s):
                cands.append(tuple(s))
    cands = sorted(set(cands))
    # keep the two extreme ones (all candidates lie on the boundary; in a 2-d pointed cone there are 2)
    assert len(cands) == 2, cands
    return cands


def sign_changes_on_I(p):
    """number of sign changes of the polynomial p(t) on (0,1) (p not identically zero)"""
    P = sp.Poly(p, t)
    sq = sp.sqf_list(P)[1]
    n = 0
    for fac, mult in sq:
        if mult % 2 == 1:
            n += fac.count_roots(0, 1) - (1 if fac.eval(0) == 0 else 0) - (1 if fac.eval(1) == 0 else 0)
    return n


def analyse(net):
    srcs = [r[0] for r in net]
    assert len(set(srcs)) == 4
    Gam = [[vec(r)[0] for r in net], [vec(r)[1] for r in net]]
    A = sp.Matrix([[s[0], s[1]] for s in srcs])
    A1 = A.row_join(sp.ones(4, 1))
    if A1.rank() < 3:
        return dict(degenerate=True)
    c = A1.T.nullspace()[0]
    c = c * sp.lcm([sp.fraction(x)[1] for x in c])
    u, v = kernel_cone_rays(Gam)
    h = [(1 - t) * u[i] + t * v[i] for i in range(4)]
    Q = sp.expand((sp.Matrix(Gam) * sp.diag(*h) * A).det())
    N = sp.expand(sum(c[i] * (v[i] - u[i]) * sp.prod([h[j] for j in range(4) if j != i]) for i in range(4)))
    return dict(degenerate=False, u=u, v=v, c=list(c), Q=Q, N=N)


def main():
    cl = classes()
    reps = []
    for key, mem in cl.items():
        net = mem[0]
        if len(key) == 4 and cone_canon([vec(r) for r in net]) == ('plane',):
            reps.append(net)
    reps.sort(key=net_str)
    print("4-source dynamically nontrivial classes:", len(reps))
    ok = len(reps) == 838
    print(("PASS" if ok else "FAIL"), "count equals 838 (BBH 2024, Lemma 24)")
    if not ok:
        FAILS.append("count")
    ndeg = 0; nfold = []; bad = []; nconst = []
    for net in reps:
        d = analyse(net)
        if d['degenerate']:
            ndeg += 1
            continue
        Q, N = d['Q'], d['N']
        if N == 0:
            # Phi is constant on I: for each kappa the positive equilibria are either none or a curve
            nconst.append(net_str(net))
            if Q != 0:
                bad.append((net_str(net), "N = 0 but Q != 0"))
            continue
        if Q == 0:
            bad.append((net_str(net), "Q = 0 but N != 0"))
            continue
        c2 = sign_changes_on_I(sp.expand(Q * N))
        cq = sign_changes_on_I(Q)
        cn = sign_changes_on_I(N)
        if c2 != 0 or cq != cn:
            bad.append((net_str(net), "C2 fails", c2, cq, cn))
        if cq > 1:
            bad.append((net_str(net), "C3 fails", cq))
        PQ = sp.Poly(Q, t)
        nroots = PQ.count_roots(0, 1) - (1 if PQ.eval(0) == 0 else 0) - (1 if PQ.eval(1) == 0 else 0)
        if nroots > 0:
            nfold.append(net)
            if cq != 1:
                bad.append((net_str(net), "Q has a root in I without sign change"))
    print("classes with affinely dependent sources (no isolated positive equilibrium):", ndeg)
    print("classes with Phi constant (N = Q = 0; positive equilibria never isolated):", len(nconst))
    for s_ in nconst:
        print("    ", s_)
    ok = len(nconst) == 9
    print(("PASS" if ok else "FAIL"), "there are 9 such classes (BBH 2024, Remark 26(f): 838 - 829 = 9)")
    if not ok:
        FAILS.append("nconst")
    print(("PASS" if not bad else "FAIL"), "C2 and C3 hold for every other class (%d)" % (len(reps) - ndeg - len(nconst)))
    for b in bad:
        print("   ", b)
    if bad:
        FAILS.append("C2C3")
    print("classes where Q vanishes on I (zero eigenvalue at a positive equilibrium for some kappa):", len(nfold))
    t1 = set(canon_class(parse_net(s)) for s in TABLE1)
    f1 = set(canon_class(n) for n in nfold)
    ok = (t1 == f1)
    print(("PASS" if ok else "FAIL"), "these are exactly the 30 networks of Table 1 of BBH 2024")
    if not ok:
        FAILS.append("table1")
        print("   only in ours:", [net_str(cl[k][0]) for k in f1 - t1])
        print("   only in Table 1:", [net_str(cl[k][0]) for k in t1 - f1])
    print("SUMMARY:", "ALL PASS" if not FAILS else "FAILURES: %s" % FAILS)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
