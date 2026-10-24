#!/usr/bin/env python3
"""check_formulas.py -- independent checks of every formula used by allm_core / l1reg / h2_lib.

Environment: Python 3.13, sympy, python-flint 0.9.  Runtime: about 1-2 minutes.  1 thread.
Claims (PASS/FAIL):
  F1  (sympy, exact) Pego-Velazquez sorted representation F = -P1 + phi^M P2 + A^-M R1 + (A phi)^-M R2 with the
      closed forms of P1, P2, R1, R2 used here equals the determinant definition (f1 g2 - g1 f2, paper (76)),
      for symbolic phi, K and M = 3..12 (exact polynomial identity after multiplying by phi^M (A phi)^M).
  F2  (sympy, exact) regular forms: P2/K^3 = P2t, P1/K^3 = P1t, (R1 + phi^-M R2)/K^3 = W/(eps E) ... with
      K = kappa eps, M = eps^-2, phi = 1 + eps^2 z, E = phi^M treated as a symbol.
  F3  (sympy, exact) lambda = (A - 1/phi)(phi - 1) = w (K + A w)/(1 + w) and Re lambda >= 0  <=>
      a^2 (1 - A u) <= u (K + 2 A u)  (w = u + i v, a = |w|);  F(1) = F'(1) = 0 for symbolic M (double root).
  F4  (sympy, exact) (w - s2)(w + s2 - A - 1 - mu) = -(mu - lambda) w for the roots w = A phi, s2 of the two
      characteristic quadratics; the small root zeta formula; Xi formula.
  F5  (ball arithmetic, dense, no use of result 8's code) at M = 40, 100, 300 (eps = M^-1/2 exact) and the Hopf
      point enclosed by allm_core: the atoms of q, p S, r(0), r(2 i omega) produced by l1reg agree entrywise with the
      dense solutions of the linear systems (B - i omega) q = 0 (q_1 = 1), p B = i omega p (<p,q> = 1),
      (mu - B) r = S q (r in H for mu = 0), and l1 (dense Kuznetsov formula, M = 40, 100) = M Lhat (l1reg).
  F6  closed forms of the pair sums / mass sums / values (l1reg.psum, mass, val) agree with direct summation
      for random atoms at M = 37, 200 (exponentially small parts included in the enclosures).
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from flint import arb, acb, acb_mat, ctx
import allm_core as C
import l1reg as R

fails = 0
def report(name, ok, extra=""):
    global fails
    fails += (not ok)
    print("%-4s %s %s" % (name, "PASS" if ok else "FAIL", extra), flush=True)

# ------------------------------------------------------------------------------------------------ F1
phi, K = sp.symbols("phi K")
A = 1 + K
okF1 = True
for Mv in range(3, 13):
    f1 = (phi - 1) + phi**Mv * A
    f2 = (1 - A*phi) + (A*phi)**(-Mv) * A**2 * phi
    G1 = phi*(1 - A*phi)*(K*phi + (K*Mv + 1)*(phi - 1))
    G2 = (1 - A*phi)*(A*Mv*(phi - 1) + K*phi + A*phi*(phi - 1))
    H1 = (phi - 1)*(K + (K*Mv + 1)*(1 - A*phi))
    H2 = (phi - 1)*((A + Mv*A**2*phi)*(1 - A*phi) + K*A*phi)
    g1 = G1 - phi**Mv*G2; g2 = H1 - (A*phi)**(-Mv)*H2
    Fdet = f1*g2 - g1*f2
    # closed forms used (paper (e:P1-1), (e:P2-1), (e:r1), (e:r2))
    S = (phi - 1)*(A*phi - 1)*(A*phi**2 - 1)
    P1 = K*((A*phi - 1)**2*phi**2 - (phi - 1)**2) + (K*Mv + 1)*S
    P2 = phi*(A*phi - 1)**3 - A**2*(phi - 1)**2 + Mv*(A*phi - 1)*A**2*(phi - 1)**2
    R1 = -A**2*((A*phi - 1)**2*phi**2 - (phi - 1)**2)
    R2 = A*(phi - 1)**3 + A**2*phi*(A*phi - 1)**2*(Mv*(phi - 1) + phi)
    Fsort = -P1 + phi**Mv*P2 + A**(-Mv)*R1 + (A*phi)**(-Mv)*R2
    d = sp.expand((Fdet - Fsort)*phi**Mv*(A*phi)**Mv*A**Mv)
    okF1 &= (d == 0)
report("F1", okF1, "sorted representation of F (paper (76), (77), (e:r1), (e:r2)) exact for M = 3..12")

# ------------------------------------------------------------------------------------------------ F2
eps, kap, z, E, Mv = sp.symbols("epsilon kappa z E M", positive=True)
Ksym = kap*eps; Asym = 1 + Ksym; ph = 1 + eps**2*z; Msym = eps**-2
g = kap + Asym*eps*z; h = g*ph + eps*z
P2s = ph*(Asym*ph - 1)**3 - Asym**2*(ph - 1)**2 + Msym*(Asym*ph - 1)*Asym**2*(ph - 1)**2
P1s = Ksym*((Asym*ph - 1)**2*ph**2 - (ph - 1)**2) + (Ksym*Msym + 1)*(ph - 1)*(Asym*ph - 1)*(Asym*ph**2 - 1)
R1s = -Asym**2*((Asym*ph - 1)**2*ph**2 - (ph - 1)**2)
R2s = Asym*(ph - 1)**3 + Asym**2*ph*(Asym*ph - 1)**2*(Msym*(ph - 1) + ph)
P2t = (ph*g**3 - eps*Asym**2*z**2 + Asym**2*z**2*g)/kap**3
P1t = (kap*(g**2*ph**2 - eps**2*z**2) + (kap + eps)*z*g*h)/kap**3
W = -Asym**2*(g**2*ph**2 - eps**2*z**2) + (Asym*eps**4*z**3 + Asym**2*ph*g**2*(z + ph))/E
ok2 = sp.simplify(P2s/Ksym**3 - P2t) == 0 and sp.simplify(P1s/Ksym**3 - P1t) == 0 \
      and sp.simplify((R1s + R2s/E)/Ksym**3 - W/(eps*kap**3)) == 0
ok2 &= sp.simplify((Asym*ph - 1) - eps*g) == 0 and sp.simplify((Asym*ph**2 - 1) - eps*h) == 0
lam = (Asym - 1/ph)*(ph - 1)
ok2 &= sp.simplify(lam - eps**3*g*z/ph) == 0
ok2 &= sp.simplify(Msym*lam/Ksym - (z + eps*z**2/(kap*ph))) == 0
report("F2", ok2, "eps-regular forms P2t, P1t, W, g, h, Lambda (exact)")

# ------------------------------------------------------------------------------------------------ F3
u, v, Kp = sp.symbols("u v K", real=True)
w = u + sp.I*v; Ap = 1 + Kp
lamw = (Ap - 1/(1 + w))*w
ok3 = sp.simplify(lamw - w*(Kp + Ap*w)/(1 + w)) == 0
rel = sp.simplify(sp.re(sp.expand_complex(lamw))*((1 + u)**2 + v**2) - (Ap*u*((1 + u)**2 + v**2) - (u + u**2 + v**2)))
ok3 &= rel == 0
a2 = u**2 + v**2
ok3 &= sp.expand((Ap*u*(1 + 2*u + a2) - (u + a2)) - (u*(Kp + 2*Ap*u) - a2*(1 - Ap*u))) == 0
# double root at phi = 1 for symbolic M (uses phi^M, (A phi)^-M as functions)
phs, Ms = sp.symbols("phi M", positive=True)
Ksb = sp.symbols("K", positive=True); Ab = 1 + Ksb
Sb = (phs - 1)*(Ab*phs - 1)*(Ab*phs**2 - 1)
P1b = Ksb*((Ab*phs - 1)**2*phs**2 - (phs - 1)**2) + (Ksb*Ms + 1)*Sb
P2b = phs*(Ab*phs - 1)**3 - Ab**2*(phs - 1)**2 + Ms*(Ab*phs - 1)*Ab**2*(phs - 1)**2
R1b = -Ab**2*((Ab*phs - 1)**2*phs**2 - (phs - 1)**2)
R2b = Ab*(phs - 1)**3 + Ab**2*phs*(Ab*phs - 1)**2*(Ms*(phs - 1) + phs)
Fb = -P1b + phs**Ms*P2b + Ab**(-Ms)*R1b + (Ab*phs)**(-Ms)*R2b
ok3 &= sp.simplify(Fb.subs(phs, 1)) == 0 and sp.simplify(sp.diff(Fb, phs).subs(phs, 1)) == 0
report("F3", ok3, "lambda(w), region Re lambda >= 0 <=> a^2(1-Au) <= u(K+2Au), F(1) = F'(1) = 0 for real M")

# ------------------------------------------------------------------------------------------------ F4
wv, s2, mu, lm, Av = sp.symbols("w s2 mu lambda A")
# w^2 - (A+1+lambda) w + A = 0 and s2^2 - (A+1+mu) s2 + A = 0  => identity
expr = (wv - s2)*(wv + s2 - Av - 1 - mu) + (mu - lm)*wv - ((wv**2 - (Av + 1 + lm)*wv + Av) - (s2**2 - (Av + 1 + mu)*s2 + Av))
ok4 = sp.expand(expr) == 0
zeta, m = sp.symbols("zeta m")
s1 = 1 + eps**2*zeta
q4 = sp.expand(((s1 - 1)*(s1 - Asym) - eps**3*m*s1)/eps**3)
ok4 &= sp.simplify(q4 - (eps*zeta**2 - (kap + eps**2*m)*zeta - m)) == 0
cm = kap + eps**2*m
zsmall = -2*m/(cm + sp.sqrt(cm**2 + 4*eps*m))
ok4 &= sp.simplify(sp.expand(eps*zsmall**2 - cm*zsmall - m)) == 0
Xi = (kap - eps*zeta)/(1 + eps**2*zeta) + Asym*eps*z - eps**2*m
ok4 &= sp.simplify((Asym*ph + Asym/s1 - Asym - 1 - eps**3*m)/eps - Xi) == 0
report("F4", ok4, "divided-difference identity, quadratic for zeta and its small root, Xi")

# ------------------------------------------------------------------------------------------------ F5 dense
ctx.prec = 256
def Bmat(K, M):
    A = 1 + K; N = M - 1
    rows = [[acb(0)]*M for _ in range(M)]
    rows[0][0] = -A*(M + 2); rows[0][1] = 1 - K
    for c in range(2, N): rows[0][c] = -K
    rows[0][M - 1] = M*K + 1
    for l in range(1, N):
        rows[l][l - 1] += A; rows[l][l] += -A - 1; rows[l][l + 1] += 1
    rows[M - 1][0] += A; rows[M - 1][N - 1] += A; rows[M - 1][M - 1] += -A
    return acb_mat(rows)
def Smat(M):
    N = M - 1
    rows = [[acb(0)]*M for _ in range(M)]
    rows[0][0] = acb(-2)
    for c in range(1, N): rows[0][c] = acb(-1)
    for l in range(1, N): rows[l][l - 1] = acb(1); rows[l][l] = acb(-1)
    rows[M - 1][N - 1] = acb(1)
    return acb_mat(rows)

def atom_at(cx, Y, l, M):
    e = acb(cx.eps_box); A = acb(cx.A)
    c = Y.c.v * e**Y.c.p
    if isinstance(Y, R.Slow):
        return c*(1 + e*e*Y.bt)**l
    if isinstance(Y, R.BM):
        return c*(A*(1 + e*e*Y.gt))**(l - M)
    if isinstance(Y, R.Dd):
        g1 = A*(1 + e*e*Y.g1); g2 = A*(1 + e*e*Y.g2)
        return c*(g1**(l - M) - g2**(l - M))/(g1 - g2)
    if isinstance(Y, R.B1):
        return c*((1 + e*e*Y.pt)/A)**l
def seq_at(cx, seq, l, M):
    return sum((atom_at(cx, Y, l, M) for Y in seq), acb(0))

def dense_solve(Amat, rhs):
    return Amat.solve(rhs)

okF5 = True
X0 = [arb("1.8982"), arb("0.0"), arb("4.4934")]
for Mv in [40, 100, 300]:
    eps_b = 1/arb(Mv).sqrt()
    # continuation from the limit
    X = X0
    for t in [0.2, 0.4, 0.6, 0.8, 1.0]:
        X = C.newton(X, eps_b*t)
    enc, box = C.certify_hopf(X, eps_b)
    kapb = enc[0]; zb = acb(enc[1], enc[2])
    out = R.blocks(eps_b, kapb, zb)
    cx = out["cx"]
    K = kapb*eps_b; Bm = Bmat(K, Mv); Sm = Smat(Mv)
    lam = (1 + K - 1/(1 + zb/Mv))*(zb/Mv)
    iom = acb(0, lam.imag)
    # q: (B - i omega) q = 0, q_1 = 1 : replace first equation by q_1 = 1
    Mq = Bm - iom*acb_mat([[1 if i == k else 0 for k in range(Mv)] for i in range(Mv)])
    rows = [[Mq[i, k] for k in range(Mv)] for i in range(Mv)]
    rows[0] = [acb(1)] + [acb(0)]*(Mv - 1)
    rhs = acb_mat([[1]] + [[0]]*(Mv - 1))
    qd = acb_mat(rows).solve(rhs)
    errq = max(float(abs(qd[l - 1, 0] - seq_at(cx, out["q"], l, Mv)).mid()) for l in range(1, Mv + 1))
    # p: p (B - i omega) = 0, normalized <p,q> = 1: transpose system with one equation replaced by <p,q> = 1
    rowsT = [[Mq[k, i] for k in range(Mv)] for i in range(Mv)]
    rowsT[1] = [qd[k, 0] for k in range(Mv)]
    rhsT = acb_mat([[0], [1]] + [[0]]*(Mv - 2))
    pd = acb_mat(rowsT).solve(rhsT)
    pS = acb_mat([[pd[k, 0] for k in range(Mv)]])*Sm
    errpS = max(float(abs(pS[0, k - 1] - (seq_at(cx, out["pS"], k, Mv) if k < Mv else 0)).mid()) for k in range(1, Mv + 1))
    Sq = Sm*qd
    # r(2 i omega) = (2 i omega - B)^{-1} S q
    I = acb_mat([[1 if i == k else 0 for k in range(Mv)] for i in range(Mv)])
    r2d = (2*iom*I - Bm).solve(Sq)
    errr2 = max(float(abs(r2d[l - 1, 0] - seq_at(cx, out["r2"], l, Mv)).mid()) for l in range(1, Mv + 1))
    # r(0): (-B) r = S q with sum l r_l = 0 replacing the first row
    rows0 = [[-Bm[i, k] for k in range(Mv)] for i in range(Mv)]
    rows0[0] = [acb(k + 1) for k in range(Mv)]
    rhs0 = acb_mat([[0]] + [[Sq[i, 0]] for i in range(1, Mv)])
    r0d = acb_mat(rows0).solve(rhs0)
    errr0 = max(float(abs(r0d[l - 1, 0] - seq_at(cx, out["r0"], l, Mv)).mid()) for l in range(1, Mv + 1))
    scale = max(float(abs(r0d[l, 0]).mid()) for l in range(Mv))
    # dense Kuznetsov l1 (normalization q_1 = 1)
    line = ""
    if Mv <= 100:
        def Bil(x, y):
            return Sm*y*x[0, 0] + Sm*x*y[0, 0]
        qbar = acb_mat([[qd[i, 0].conjugate()] for i in range(Mv)])
        w11 = -(r0d + acb_mat([[r0d[i, 0].conjugate()] for i in range(Mv)]))
        w20 = 2*r2d
        pv = acb_mat([[pd[k, 0] for k in range(Mv)]])
        c = -2*(pv*Bil(qd, w11))[0, 0] + (pv*Bil(qbar, w20))[0, 0]
        l1d = c.real/(2*lam.imag)
        rel = float(abs(l1d/Mv - out["Lhat"].real).mid()/abs(l1d/Mv).mid())
        okF5 &= rel < 1e-30
        line = " l1(dense)/M = %s, Lhat = %s (rel %.1e)" % (arb(l1d/Mv).str(12), arb(out["Lhat"].real).str(12), rel)
    ok = errq < 1e-40 and errpS < 1e-40 and errr2 < 1e-35*max(1, scale) and errr0 < 1e-35*max(1, scale)
    okF5 &= ok
    print("   M=%d: max entry errors q %.1e, pS %.1e, r(2iw) %.1e, r(0) %.1e%s" % (Mv, errq, errpS, errr2, errr0, line))
report("F5", okF5, "atoms of q, pS, r(0), r(2 i omega) = dense solutions entrywise; Lhat = dense Kuznetsov l1/M")

# ------------------------------------------------------------------------------------------------ F6
random.seed(20261024)
okF6 = True
ctx.prec = 192
for Mv in [37, 200]:
    eps_b = 1/arb(Mv).sqrt(); kapb = arb("2.7")
    cx = R.Ctx(eps_b, kapb)
    def rnd():
        return acb(random.uniform(-3, 3), random.uniform(-6, 6))
    for trial in range(6):
        Xs = [R.Slow(R.SN(rnd(), random.randint(-2, 2)), rnd()), R.B1(R.SN(rnd(), 1), rnd())]
        g1 = rnd(); g2 = g1 + acb(random.uniform(-1, 1), random.uniform(-1, 1))
        Ys = [R.Slow(R.SN(rnd(), -1), rnd()), R.BM(R.SN(rnd(), 0), rnd()), R.Dd(R.SN(rnd(), -1), g1, g2)]
        for X in Xs:
            for Y in Ys:
                cf = R.psum(cx, X, Y)
                direct = sum((atom_at(cx, X, k, Mv)*atom_at(cx, Y, k, Mv) for k in range(1, Mv)), acb(0))
                cv = cf.v*acb(eps_b)**cf.p
                okF6 &= cv.overlaps(direct)
        for Y in Ys:
            cf = R.mass(cx, Y); cv = cf.v*acb(eps_b)**cf.p
            direct = sum((l*atom_at(cx, Y, l, Mv) for l in range(1, Mv + 1)), acb(0))
            okF6 &= cv.overlaps(direct)
            for where, l in (("1", 1), ("N", Mv - 1), ("M", Mv)):
                cf = R.val(cx, Y, where); cv = cf.v*acb(eps_b)**cf.p
                okF6 &= cv.overlaps(atom_at(cx, Y, l, Mv))
report("F6", okF6, "pair sums, mass sums and boundary values = direct summation (M = 37, 200, random atoms)")

print("SUMMARY:", "PASS" if fails == 0 else "FAIL")
sys.exit(0 if fails == 0 else 1)
