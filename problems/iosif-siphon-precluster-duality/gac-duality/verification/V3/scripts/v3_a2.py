"""Verifier V3 (Claude Opus 5.5, effort xhigh), 27 October 2026: independent exact re-checks of agent A2's claims
(round3/dichotomy-astra/notes.md, GPT-6-Astra).  Does not import A2's or D2's code (own v3lib + sympy only).

  A    Theorem A: network certificates, dual equations, boundary vector-field signs (Sec. 2.2), integrating-factor
       identity (A.1), and the margin bound of Sec. 2.3 on random exact points of U (random rates and T).
  AF   Theorem AF: u = 3..7, R = 1..7: hypotheses, dual, identity (AF.1), facet c = 0 is OUT; R >= 3 margin bound
       (sampled exactly); R = 2 and R = 1 extinction inequalities (symbolic reduction + sampled exactly).
  E    Theorem E: the source-span hypothesis is equivalent to E2's Cor. E2.5' hypothesis (im Yp in im Ye); checked
       on random WR networks; the equilibrium curve of the Theorem-A dual meets each class once (dimension r - s = 1).
  L    Proposition L: random single-linkage rank-two networks: no siphon edge of N* with opposite active vectors.
  O    Theorem O: uniqueness parametrisation (identities, monotonicity), equilibrium, Jacobian, eigenvalues, H.
  S7   Section 7 candidate: exact network, vertex p realizable, Jacobian [[-k2,0],[k2,-k4]] in (k1,k3).
Usage: python v3_a2.py   (exit 1 on FAIL).  Environment: Python 3.12, sympy 1.14; exact arithmetic.
"""
import sys
import time
import random
import platform
from fractions import Fraction as Fr
from itertools import combinations
import sympy as sp

import v3lib as L

FAILS = 0


def claim(name, ok, msg=""):
    global FAILS
    print(("PASS" if ok else "FAIL") + ": " + name + ((" -- " + msg) if msg else ""))
    FAILS += 0 if ok else 1


def net_uR(u, R):
    # N_(u,R): u X3 -> X1+X2 -> X2+2X3+R X4 -> u X3
    c1 = (0, 0, u, 0)
    c2 = (1, 1, 0, 0)
    c3 = (0, 1, 2, R)
    return [(c1, c2), (c2, c3), (c3, c1)]


def dual_field(nd, Rd, kap, x):
    f = [0] * nd
    for i, (s, p) in enumerate(Rd):
        rate = kap[i]
        for j in range(nd):
            rate *= x[j] ** s[j]
        for j in range(nd):
            f[j] += rate * (p[j] - s[j])
    return [sp.expand(e) for e in f]


def main():
    t0 = time.time()
    rnd = random.Random(20261027)
    print("# v3_a2.py  Python %s on %s; sympy %s; exact" % (platform.python_version(), platform.system(),
                                                            sp.__version__))
    al, be, ga, de, T = sp.symbols('alpha beta gamma delta T', positive=True)
    a, b, c = sp.symbols('a b c', positive=True)
    # ------------------------------------------------------------------ A
    NA = net_uR(3, 3)
    G = L.Gamma(4, NA)
    claim("A.net  Gamma = [[1,-1,0],[1,0,-1],[-3,2,1],[0,3,-3]]; WR, (m,l,s,delta) = (3,1,2,0); (H1) with "
          "(11,1,4,1) (value 12 on every complex), (H2), (H3)",
          G == [[1, -1, 0], [1, 0, -1], [-3, 2, 1], [0, 3, -3]] and L.is_wr(NA) and L.counts(4, NA) == (3, 1, 2, 0)
          and all(sum(w * y for w, y in zip((11, 1, 4, 1), cx)) == 12 for cx in L.complexes(NA)) and
          L.consistent(4, NA) and L.H3(4, NA))
    nd, Rd = L.dual(4, NA)
    claim("A.dual  N* = b -> a, b+c -> a+b, 3a+2c -> 2b+3c, 3c -> 3b",
          Rd == [((0, 1, 0), (1, 0, 0)), ((0, 1, 1), (1, 1, 0)), ((3, 0, 2), (0, 2, 3)), ((0, 0, 3), (0, 3, 0))])
    f = dual_field(3, Rd, [al, be, ga, de], [a, b, c])
    claim("A.ode  a' = alpha b + beta b c - 3 gamma a^3 c^2, b' = -alpha b + 2 gamma a^3 c^2 + 3 delta c^3, "
          "c' = -beta b c + gamma a^3 c^2 - 3 delta c^3",
          f[0] == sp.expand(al * b + be * b * c - 3 * ga * a**3 * c**2) and
          f[1] == sp.expand(-al * b + 2 * ga * a**3 * c**2 + 3 * de * c**3) and
          f[2] == sp.expand(-be * b * c + ga * a**3 * c**2 - 3 * de * c**3) and sp.expand(sum(f)) == 0)
    # boundary signs
    claim("A.bdry  a' = alpha b + beta b c > 0 on open edge a=0; b' = 2 gamma a^3 c^2 + 3 delta c^3 > 0 on open edge "
          "b=0; c=0 invariant with b' = -alpha b; field nonzero at (0,T,0) and (0,0,T), zero at (T,0,0)",
          f[0].subs(a, 0) == sp.expand(al * b + be * b * c) and f[1].subs(b, 0) == sp.expand(2 * ga * a**3 * c**2 + 3 * de * c**3)
          and f[2].subs(c, 0) == 0 and f[1].subs(c, 0) == -al * b and
          [e.subs({a: 0, b: T, c: 0}) for e in f] == [al * T, -al * T, 0] and
          [e.subs({a: 0, b: 0, c: T}) for e in f][1] == 3 * de * T**3 and
          [e.subs({a: T, b: 0, c: 0}) for e in f] == [0, 0, 0])
    fts = L.facet_types_any_rank(4, NA)
    sips = [Z for Z, t in fts if L.is_siphon(3, Rd, Z)]
    claim("A.OUT  the only siphon facet of N* is {c} (Z = {R3}) and its type is OUT",
          sips == [(2,)] and dict(fts)[(2,)] == 'OUT', str(fts))
    logV_dot = sp.expand(f[2] / c - be / al * f[1])
    claim("A.1  (log V)' = c[gamma a^3 - (3 delta + 2 beta gamma a^3/alpha) c - (3 beta delta/alpha) c^2] for "
          "V = c exp(-beta b/alpha) (exact)",
          sp.simplify(logV_dot - c * (ga * a**3 - (3 * de + 2 * be * ga * a**3 / al) * c - 3 * be * de / al * c**2)) == 0)
    # margin bound sampled exactly
    ok = True
    for _ in range(3000):
        A_, B_, G_, D_ = [Fr(rnd.randint(1, 50), rnd.randint(1, 50)) for _ in range(4)]
        TT = Fr(rnd.randint(1, 40), rnd.randint(1, 20))
        E = 3 * D_ + 2 * B_ * G_ * TT**3 / A_ + 3 * B_ * D_ / A_
        eps = min(TT / 2, Fr(1), G_ * TT**3 / (16 * E))
        bb = eps * Fr(rnd.randint(0, 999), 1000)
        cc = (eps - bb) * Fr(rnd.randint(1, 999), 1000)
        aa = TT - bb - cc
        br = G_ * aa**3 - (3 * D_ + 2 * B_ * G_ * aa**3 / A_) * cc - 3 * B_ * D_ / A_ * cc**2
        if not (aa > TT / 2 and br >= G_ * TT**3 / 16):
            ok = False
            break
    claim("A.2.3  bracket >= gamma T^3/16 on U = {b + c < epsilon} (3000 random exact points, random rates and T)", ok)
    claim("A.unit  rates = T = 1: epsilon = min(1/2, 1, 1/(16*8)) = 1/128",
          min(Fr(1, 2), Fr(1), Fr(1, 16 * 8)) == Fr(1, 128))
    # ------------------------------------------------------------------ AF
    okH = okD = okI = okO = okB = ok2 = ok1 = True
    for u in range(3, 8):
        for R in range(1, 8):
            N = net_uR(u, R)
            w = (u * R + 2, u - 2, R + 1, u - 2)
            if not (L.is_wr(N) and L.counts(4, N) == (3, 1, 2, 0) and L.H3(4, N) and L.consistent(4, N) and
                    all(sum(x * y for x, y in zip(w, cx)) == u * (R + 1) for cx in L.complexes(N))):
                okH = False
            ndd, Rdd = L.dual(4, N)
            if Rdd != [((0, 1, 0), (1, 0, 0)), ((0, 1, 1), (1, 1, 0)), ((u, 0, 2), (0, 2, u)), ((0, 0, R), (0, R, 0))]:
                okD = False
            ff = dual_field(3, Rdd, [al, be, ga, de], [a, b, c])
            if sp.simplify(ff[2] / c - be / al * ff[1] - c * ((u - 2) * ga * a**u - R * de * c**(R - 2)
                                                            - 2 * be * ga / al * a**u * c - R * be * de / al * c**(R - 1))) != 0:
                okI = False
            if dict(L.facet_types_any_rank(4, N)).get((2,)) != 'OUT':
                okO = False
            for _ in range(60):
                A_, B_, G_, D_ = [Fr(rnd.randint(1, 50), rnd.randint(1, 50)) for _ in range(4)]
                TT = Fr(rnd.randint(1, 40), rnd.randint(1, 20))
                if R >= 3:
                    Lm = (u - 2) * G_ * (TT / 2)**u
                    E = R * D_ + 2 * B_ * G_ * TT**u / A_ + R * B_ * D_ / A_
                    eps = min(TT / 2, Fr(1), Lm / (2 * E))
                    bb = eps * Fr(rnd.randint(0, 999), 1000)
                    cc = (eps - bb) * Fr(rnd.randint(1, 999), 1000)
                    aa = TT - bb - cc
                    br = ((u - 2) * G_ * aa**u - R * D_ * cc**(R - 2) - 2 * B_ * G_ / A_ * aa**u * cc
                          - R * B_ * D_ / A_ * cc**(R - 1))
                    if not br >= Lm / 2:
                        okB = False
                else:
                    if R == 2:
                        TT = min(Fr(1), D_ / ((u - 2) * G_))
                    else:
                        TT = min(Fr(1), D_ / (2 * (u - 2) * G_))
                    bb = TT * Fr(rnd.randint(0, 1000), 1001)
                    cc = (TT - bb) * Fr(rnd.randint(1, 1000), 1001)
                    aa = TT - bb - cc
                    cdot = -B_ * bb * cc + (u - 2) * G_ * aa**u * cc**2 - R * D_ * cc**R
                    if R == 2 and not cdot <= -D_ * cc**2:
                        ok2 = False
                    if R == 1 and not cdot <= -D_ / 2 * cc:
                        ok1 = False
    claim("AF.net  u=3..7, R=1..7: WR 3-cycle, (m,l,s,delta)=(3,1,2,0), (H1) with (uR+2,u-2,R+1,u-2) (value u(R+1)), "
          "(H2), (H3)", okH)
    claim("AF.dual  N* = b -> a, b+c -> a+b, ua+2c -> 2b+uc, Rc -> Rb", okD)
    claim("AF.1  (log V)' identity (exact, sympy) for all (u,R) tested", okI)
    claim("AF.OUT  facet {c} is OUT for all (u,R) tested", okO)
    claim("AF.R>=3  bracket >= L/2 on {b + c < epsilon} (60 random exact points per (u,R))", okB)
    claim("AF.R=2  with T = min(1, delta/((u-2)gamma)): c' <= -delta c^2 on the class (sampled exactly)", ok2)
    claim("AF.R=1  with T = min(1, delta/(2(u-2)gamma)): c' <= -(delta/2) c on the class (sampled exactly)", ok1)
    # symbolic reductions for R = 2, 1 (the inequality chain): c' + delta c^2 = -beta b c + c^2((u-2) gamma a^u - delta)
    uu = sp.Symbol('u', positive=True)
    claim("AF.R=2,1  c' + delta c^2 = -beta b c + c^2((u-2)gamma a^u - delta) (R=2); c' + (delta/2) c = "
          "-beta b c + c((u-2) gamma a^u c - delta/2) (R=1)",
          sp.simplify((-be * b * c + (uu - 2) * ga * a**uu * c**2 - 2 * de * c**2) + de * c**2
                      - (-be * b * c + c**2 * ((uu - 2) * ga * a**uu - de))) == 0 and
          sp.simplify((-be * b * c + (uu - 2) * ga * a**uu * c**2 - de * c) + de / 2 * c
                      - (-be * b * c + c * ((uu - 2) * ga * a**uu * c - de / 2))) == 0)
    # ------------------------------------------------------------------ E
    okE = True
    cntE = 0
    for _ in range(300):
        n = rnd.randint(2, 5)
        Rn = []
        used = set()
        for _cl in range(rnd.randint(1, 3)):
            k = rnd.randint(2, 4)
            pts = []
            while len(pts) < k:
                y = tuple(rnd.randint(0, 3) for _ in range(n))
                if sum(y) > 0 and y not in used and y not in pts:
                    pts.append(y)
            used.update(pts)
            order = list(range(k))
            rnd.shuffle(order)
            Rn += [(pts[order[i]], pts[order[(i + 1) % k]]) for i in range(k)]
        Ye = L.Ye(n, Rn)
        Gm = L.Gamma(n, Rn)
        cntE += 1
        if L.rank([list(r) for r in zip(*(Ye + Gm))]) != L.rank([list(r) for r in zip(*Ye)]) and False:
            okE = False
        # im Gamma subset im Ye  <=>  rank [Ye | Gamma] = rank Ye  (columns)
        if L.rank([Ye[i] + Gm[i] for i in range(n)]) != L.rank(Ye):
            okE = False
    claim("E.span  im(Gamma) in im(Ye) for %d random WR networks (WR-Sigma), i.e. E2's Cor. E2.5' hypothesis" % cntE,
          okE)
    # Theorem A dual: equilibria of N* for fixed rates form a curve meeting each class a+b+c = T exactly once (O)
    # ------------------------------------------------------------------ O
    cc_ = sp.symbols('c', positive=True)
    bfun = 9 * de * cc_**3 / (al - 2 * be * cc_)
    a3fun = 3 * de * cc_ * (al + be * cc_) / (ga * (al - 2 * be * cc_))
    sub = {b: bfun, c: cc_}
    eqs = [sp.simplify(fi.subs(sub).subs(a**3, a3fun)) for fi in f]
    claim("O.uniq  b' - 2c' = (2 beta c - alpha) b + 9 delta c^3; b(c), a(c)^3 as stated solve all three equilibrium "
          "equations", sp.expand(f[1] - 2 * f[2] - ((2 * be * c - al) * b + 9 * de * c**3)) == 0 and
          all(sp.simplify(e) == 0 for e in eqs), str(eqs))
    db = sp.simplify(sp.diff(bfun, cc_))
    da3 = sp.simplify(sp.diff(a3fun, cc_))
    claim("O.mono  b'(c) = %s > 0 and (a^3)'(c) = %s > 0 on 0 < c < alpha/(2 beta)" % (sp.factor(db), sp.factor(da3)),
          sp.simplify(db - 9 * de * cc_**2 * (3 * al - 4 * be * cc_) / (al - 2 * be * cc_)**2) == 0 and
          sp.simplify(da3 - 3 * de * (al**2 + 2 * al * be * cc_ - 2 * be**2 * cc_**2) / (ga * (al - 2 * be * cc_)**2)) == 0)
    # positivity of numerators on the interval: 3 alpha - 4 beta c > alpha > 0; alpha^2 + 4 alpha beta c - 2 beta^2 c^2
    # >= alpha^2 + 4 alpha beta c - alpha beta c > 0 since beta c < alpha/2
    x = sp.Symbol('x', positive=True)  # x = beta c / alpha in (0, 1/2)
    # numerators: 9 delta c^2 (3 alpha - 4 beta c) and 3 delta alpha^2 (1 + 2x - 2x^2) with x = beta c/alpha in (0,1/2);
    # 3 - 4x is linear and 1 + 2x - 2x^2 concave, so positivity at both endpoints suffices
    qx = 1 + 2 * x - 2 * x**2
    claim("O.mono2  on x = beta c/alpha in [0,1/2]: 3 - 4x > 0 (linear) and 1 + 2x - 2x^2 > 0 (concave): endpoint "
          "values positive; numerator identity alpha^2+2 alpha beta c-2 beta^2 c^2 = alpha^2 (1+2x-2x^2)",
          all((3 - 4 * xv) > 0 and qx.subs(x, xv) > 0 for xv in (sp.Rational(0), sp.Rational(1, 2))) and
          sp.diff(qx, x, 2) < 0 and
          sp.expand(al**2 * qx.subs(x, be * cc_ / al) - (al**2 + 2 * al * be * cc_ - 2 * be**2 * cc_**2)) == 0)
    rates = {al: 17, be: 400, ga: 70000, de: 10**6}
    pt = {a: 1, b: 1, c: sp.Rational(1, 100)}
    flux = [17 * 1, 400 * 1 * sp.Rational(1, 100), 70000 * sp.Rational(1, 10**4), 10**6 * sp.Rational(1, 10**6)]
    claim("O.eq  (1,1,1/100) is an equilibrium at rates (17,400,70000,10^6), T = 201/100, fluxes (17,4,7,1)",
          [e.subs(rates).subs(pt) for e in f] == [0, 0, 0] and flux == [17, 4, 7, 1] and
          sum(pt.values()) == sp.Rational(201, 100))
    TT = sp.Rational(201, 100)
    fb = f[1].subs(rates).subs(a, TT - b - c)
    fc = f[2].subs(rates).subs(a, TT - b - c)
    J = sp.Matrix([[sp.diff(fb, b), sp.diff(fb, c)], [sp.diff(fc, b), sp.diff(fc, c)]]).subs({b: 1, c: sp.Rational(1, 100)})
    claim("O.J  Jacobian in (b,c) = [[-59,3658],[-25,79]], trace 20, det 86789, eigenvalues 10 +- i sqrt(86689)",
          J == sp.Matrix([[-59, 3658], [-25, 79]]) and J.trace() == 20 and J.det() == 86789 and
          set(J.eigenvals()) == {10 + sp.I * sp.sqrt(86689), 10 - sp.I * sp.sqrt(86689)})
    H = sp.Matrix([[93655, -290457], [-290457, 13471234]])
    claim("O.H  H positive definite (93655 > 0, det = %d > 0), J^T H + H J = 3471560 I" % H.det(),
          H[0, 0] > 0 and H.det() > 0 and J.T * H + H * J == 3471560 * sp.eye(2))
    # uniqueness at these rates: c-range (0, 17/800); the equation a(c)+b(c)+c = T has the root c = 1/100
    claim("O.root  c = 1/100 lies in (0, alpha/(2 beta)) = (0, 17/800) and gives b = 1, a^3 = 1",
          sp.Rational(1, 100) < sp.Rational(17, 800) and bfun.subs(rates).subs(cc_, sp.Rational(1, 100)) == 1 and
          a3fun.subs(rates).subs(cc_, sp.Rational(1, 100)) == 1)
    # ------------------------------------------------------------------ L
    okL = True
    tested = 0
    for _ in range(4000):
        n = rnd.randint(3, 5)
        k = rnd.randint(3, 5)
        pts = []
        while len(pts) < k:
            y = tuple(rnd.randint(0, 3) for _ in range(n))
            if sum(y) > 0 and y not in pts:
                pts.append(y)
        # one linkage class: random connected digraph (spanning path + extras), not necessarily WR
        order = list(range(k))
        rnd.shuffle(order)
        E = set()
        for i in range(k - 1):
            u_, v_ = order[i], order[i + 1]
            E.add((u_, v_) if rnd.random() < 0.5 else (v_, u_))
        for _e in range(rnd.randint(0, 3)):
            u_, v_ = rnd.sample(range(k), 2)
            E.add((u_, v_))
        Rn = [(pts[u_], pts[v_]) for u_, v_ in sorted(E)]
        Gm = L.Gamma(n, Rn)
        if L.rank(Gm) != 2 or len(L.linkage_classes(Rn)) != 1:
            continue
        ndd, Rdd = L.dual(n, Rn)
        tested += 1
        for Z, t in L.facet_types_any_rank(n, Rn):
            if not L.is_siphon(ndd, Rdd, Z):
                continue
            act = [i for i in range(n) if all(Rn[j][0][i] == 0 for j in Z) and any(Gm[i])]
            for i, h in combinations(act, 2):
                if sum(Gm[i][j] * Gm[h][j] for j in range(len(Rn))) < 0:
                    okL = False
    claim("L  no siphon edge with two opposite active dual reaction vectors in %d random single-linkage rank-two "
          "networks" % tested, okL and tested > 500)
    # ------------------------------------------------------------------ S7
    Gs = [(-2, -1, 0, 0, 2, 1), (2, 1, -2, -1, 0, 0), (0, 0, 2, 1, -2, -1)]
    Ss = [(2, 1, 0, 0, 1, 2), (1, 2, 2, 1, 0, 0), (0, 0, 1, 2, 2, 1)]
    N7 = []
    for j, order in enumerate([(0, 2, 1), (1, 0, 2), (2, 1, 0)]):
        z = Ss[j]
        for d in order:
            zz = tuple(x + y for x, y in zip(z, Gs[d]))
            N7.append((z, zz))
            z = zz
    G7 = L.Gamma(6, N7)
    claim("S7.net  six species, nine reactions (three 3-cycles), WR, every complex has total 6 ((H1) with ones), "
          "(H2), (H3), rank 2", len(N7) == 9 and L.is_wr(N7) and all(sum(y) == 6 for s_, p_ in N7 for y in (s_, p_))
          and L.consistent(6, N7) and L.H3(6, N7) and L.rank(G7) == 2)
    n7, R7 = L.dual(6, N7)
    p = [1] * 9
    p[0] = p[2] = 0
    zv = [0, -1, 0, -2, 0, 0]
    direction = [sum(G7[i][j] * zv[i] for i in range(6)) for j in range(9)]
    pos = [Fr(x) + Fr(dd, 4) for x, dd in zip(p, direction)]
    kap = sp.symbols('kap1:7', positive=True)
    ks = sp.symbols('k1:10', positive=True)
    f7 = dual_field(9, R7, kap, ks)
    psub = {ks[j]: p[j] for j in range(9)}
    eq7 = all(sp.simplify(e.subs(psub)) == 0 for e in f7)
    # chart: S* basis from rows of Gamma (two independent rows), coordinates (k1, k3)
    rows = [list(r) for r in G7]
    i1 = next(i for i in range(6) if any(rows[i]))
    i2 = next(i for i in range(6) if L.rank([rows[i1], rows[i]]) == 2)
    B = sp.Matrix([[rows[i1][j], rows[i2][j]] for j in range(9)])
    Df = sp.Matrix([[sp.diff(e, kk) for kk in ks] for e in f7]).subs(psub)
    P = sp.zeros(2, 9)
    P[0, 0] = 1
    P[1, 2] = 1
    PB = P * B
    J7 = sp.simplify(P * Df * B * PB.inv())
    claim("S7.vertex  p = (0,1,0,1,1,1,1,1,1) is a boundary point of a positive class (p + Gamma^T z/4 > 0), an "
          "equilibrium for every kappa, with Jacobian in (k1,k3) = [[-kappa2, 0],[kappa2, -kappa4]] (Hurwitz)",
          all(x > 0 for x in pos) and eq7 and PB.det() != 0 and
          J7 == sp.Matrix([[-kap[1], 0], [kap[1], -kap[3]]]), str(J7.tolist()))
    print("# runtime %.1f s" % (time.time() - t0))
    print("ALL CLAIMS PASS" if FAILS == 0 else "%d FAIL" % FAILS)
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
