#!/usr/bin/env python3
"""Symbolic checks of every algebraic identity used in the manuscript (exact, SymPy).

Model: the n-reactor generalisation of Bisbis (1987), eqs. (5)-(8):
  dy_k/dt   = (y_{k-1}-y_k)/alpha_k + Da0*tau*(1-y_k)*E(theta_k)
  dth_k/dt  = (th_{k-1}-th_k)/alpha_k + Da0*tau*B*(1-y_k)*E(theta_k) - (beta_k*tau/alpha_k)*(th_k-eta_k)
  E(th) = exp(th/(1+eps*th)),  alpha_1 = 1,  y_0 = th_0 = 0.

Claims checked (each prints PASS/FAIL; exit code 1 on any FAIL):
  A1  reduced variables c = 1-y, z = eps+eps^2*theta turn E into exp(1/eps-1/z) and the steady-state
      equations of reactor k into the reduced recursion of Lemma 1.
  A2  one-reactor scalar reduction: on y = 1-c_in/(1+D E), alpha*f_y = 0 and alpha*f_th = (1+b)*phi(theta).
  A3  Jacobian block identities: det(alpha J) = -(1+DE)(1+b) phi',
      tr(alpha J) = (1+DE)(1+b)(phi' - lambda + DE/(1+DE)).
  A4  g' = g(1-g)/(1+eps th)^2 and g'' = g(1-g)(1+eps th)^-4 [1 - 2g - 2 eps (1+eps th)]; eps = 0 version.
  A5  F(z) = eps^2 phi(theta), F'(z) = phi'(theta), F'(z) = X c psi(z)/(1+k)^2 - 1 with X = lambda*Delta*K.
  A6  G''(z) = k z^-4 (1+k)^-3 [1 - 2z - k(1+2z)],  psi'(z) = e^{-1/z} z^-4 (1-2z).
  A7  the full Jacobian of the n = 3 cascade is block lower triangular, with sub-diagonal blocks I/alpha_k.
Environment: Python 3.12+, SymPy 1.14 (WSL Fedora 44, ~/envs/sci). Runtime: about 20 s.
"""
import sys
import sympy as sp

fails = 0


def report(name, ok):
    global fails
    print(("PASS " if ok else "FAIL ") + name)
    if not ok:
        fails += 1


def is_zero(expr):
    e = sp.simplify(expr)
    if e == 0:
        return True
    return sp.simplify(sp.expand(sp.powsimp(sp.expand_power_exp(e)))) == 0


eps, B, Da, tau, al, be = sp.symbols('epsilon B Da0 tau alpha beta', positive=True)
eta = sp.symbols('eta', real=True)
y, th, yin, thin = sp.symbols('y theta y_in theta_in', real=True)


def E(t):
    return sp.exp(t/(1 + eps*t))


def f_y(yk, tk, yprev, tprev, a=al):
    return (yprev - yk)/a + Da*tau*(1 - yk)*E(tk)


def f_th(yk, tk, yprev, tprev, a=al, b_=be, et=eta):
    return (tprev - tk)/a + Da*tau*B*(1 - yk)*E(tk) - (b_*tau/a)*(tk - et)


# ---------- A1: reduced variables ----------
c, z, cin, zin, zc = sp.symbols('c z c_in z_in z_c', positive=True)
sub = {y: 1 - c, th: (z - eps)/eps**2, yin: 1 - cin, thin: (zin - eps)/eps**2}
etaval = (zc - eps)/eps**2
report("A1a  E((z-eps)/eps^2) = exp(1/eps - 1/z)",
       is_zero(E((z - eps)/eps**2) - sp.exp(1/eps - 1/z)))
K = al*Da*tau*sp.exp(1/eps)          # infinite-temperature Damkohler number of the reactor
Dl = eps**2*B                         # Delta = eps^2 B
bb = be*tau                           # b = beta*tau
kk = K*sp.exp(-1/z)
mass = sp.simplify(al*f_y(y, th, yin, thin).subs(sub))
energy = sp.simplify(al*eps**2*f_th(y, th, yin, thin, et=etaval).subs(sub))
report("A1b  alpha*f_y = -(c_in - c - k c),  k = K exp(-1/z)", is_zero(mass + (cin - c - kk*c)))
report("A1c  alpha*eps^2*f_th = z_in - z + Delta k c - b (z - z_c)",
       is_zero(energy - (zin - z + Dl*kk*c - bb*(z - zc))))
lam = 1/(1 + bb)
cout = cin/(1 + kk)
Fz = lam*zin + (1 - lam)*zc + lam*Dl*cin*kk/(1 + kk) - z
report("A1d  mass balance <=> c = c_in/(1+k)", is_zero((cin - c - kk*c).subs(c, cout)))
report("A1e  on c = c_in/(1+k): energy equation = (1+b) F(z)",
       is_zero((zin - z + Dl*kk*c - bb*(z - zc)).subs(c, cout) - (1 + bb)*Fz))

# ---------- A2/A3: one reactor in thesis variables ----------
D = al*Da*tau
Eth = E(th)
g = D*Eth/(1 + D*Eth)
phi = (thin + bb*eta)/(1 + bb) + B*(1 - yin)*g/(1 + bb) - th
yss = 1 - (1 - yin)/(1 + D*Eth)
report("A2a  alpha*f_y vanishes on y = 1 - c_in/(1+D E)",
       is_zero((al*f_y(y, th, yin, thin)).subs(y, yss)))
report("A2b  alpha*f_th = (1+b) phi(theta) on that curve",
       is_zero((al*f_th(y, th, yin, thin)).subs(y, yss) - (1 + bb)*phi))
J = sp.Matrix([[sp.diff(al*f_y(y, th, yin, thin), v) for v in (y, th)],
               [sp.diff(al*f_th(y, th, yin, thin), v) for v in (y, th)]])
phip = sp.diff(phi, th)
DE = D*Eth
report("A3a  det(alpha J) = -(1+DE)(1+b) phi'(theta) at a steady state",
       is_zero(J.det().subs(y, yss) + (1 + DE)*(1 + bb)*phip))
report("A3b  tr(alpha J) = (1+DE)(1+b)(phi' - lambda + DE/(1+DE)) at a steady state",
       is_zero(J.trace().subs(y, yss) - (1 + DE)*(1 + bb)*(phip - lam + DE/(1 + DE))))

# ---------- A4: sign factor of g'' ----------
gp = sp.diff(g, th)
gpp = sp.diff(g, th, 2)
report("A4a  g' = g(1-g)/(1+eps theta)^2", is_zero(gp - g*(1 - g)/(1 + eps*th)**2))
report("A4b  g'' = g(1-g)(1+eps theta)^-4 [1 - 2g - 2 eps (1+eps theta)]",
       is_zero(gpp - g*(1 - g)/(1 + eps*th)**4*(1 - 2*g - 2*eps*(1 + eps*th))))
g0 = D*sp.exp(th)/(1 + D*sp.exp(th))
report("A4c  eps = 0: g'' = g(1-g)(1-2g)", is_zero(sp.diff(g0, th, 2) - g0*(1 - g0)*(1 - 2*g0)))

# ---------- A5: F(z) = eps^2 phi(theta) ----------
phi_c = phi.subs({yin: 1 - cin, thin: (zin - eps)/eps**2, eta: etaval})
report("A5a  F(z) = eps^2 phi(theta) with z = eps + eps^2 theta",
       is_zero(Fz.subs(z, eps + eps**2*th) - eps**2*phi_c))
report("A5b  F'(z) = phi'(theta)",
       is_zero(sp.diff(Fz, z).subs(z, eps + eps**2*th) - sp.diff(phi_c, th)))
X = lam*Dl*K
psi = sp.exp(-1/z)/z**2
report("A5c  F'(z) = X c_in psi(z)/(1+k)^2 - 1,  X = lambda Delta K",
       is_zero(sp.diff(Fz, z) - (X*cin*psi/(1 + kk)**2 - 1)))

# ---------- A6: G'' and psi' ----------
Kp = sp.symbols('K', positive=True)
k2 = Kp*sp.exp(-1/z)
G = k2/(1 + k2)
report("A6a  G''(z) = k z^-4 (1+k)^-3 [1 - 2z - k(1+2z)]",
       is_zero(sp.diff(G, z, 2) - k2/z**4/(1 + k2)**3*(1 - 2*z - k2*(1 + 2*z))))
report("A6b  psi'(z) = exp(-1/z) z^-4 (1 - 2z)",
       is_zero(sp.diff(psi, z) - sp.exp(-1/z)/z**4*(1 - 2*z)))
report("A6c  e^{-1/z}/(1+K e^{-1/z}) = 1/(e^{1/z}+K)",
       is_zero(1/(sp.exp(1/z) + Kp) - sp.exp(-1/z)/(1 + Kp*sp.exp(-1/z))))

# ---------- A7: block lower triangular Jacobian, n = 3 ----------
n = 3
ys = sp.symbols('y1:%d' % (n + 1), real=True)
ts = sp.symbols('t1:%d' % (n + 1), real=True)
als = (sp.Integer(1),) + sp.symbols('alpha2:%d' % (n + 1), positive=True)
bes = sp.symbols('beta1:%d' % (n + 1), positive=True)
ets = sp.symbols('eta1:%d' % (n + 1), real=True)
rhs = []
for k in range(n):
    yp = ys[k - 1] if k > 0 else 0
    tp = ts[k - 1] if k > 0 else 0
    rhs.append(f_y(ys[k], ts[k], yp, tp, a=als[k]))
    rhs.append(f_th(ys[k], ts[k], yp, tp, a=als[k], b_=bes[k], et=ets[k]))
vars_ = [v for pair in zip(ys, ts) for v in pair]
Jfull = sp.Matrix(rhs).jacobian(vars_)
upper_zero = all(Jfull[i, j] == 0 for i in range(2*n) for j in range(2*n) if j//2 > i//2)
report("A7a  n=3 Jacobian: every entry above the 2x2 diagonal blocks is identically 0", upper_zero)
sub_ok = all(sp.simplify(Jfull[2*k + i, 2*(k - 1) + j] - (sp.Integer(1)/als[k] if i == j else 0)) == 0
             for k in range(1, n) for i in range(2) for j in range(2))
report("A7b  n=3 Jacobian: sub-diagonal blocks are (1/alpha_k) I", sub_ok)
blk_ok = all(sp.simplify(Jfull[2*k:2*k + 2, 2*k:2*k + 2]
                         - (J/als[k]).subs({al: als[k], be: bes[k], eta: ets[k], y: ys[k], th: ts[k]})) == sp.zeros(2, 2)
             for k in range(n))
report("A7c  n=3 Jacobian: diagonal block k equals the one-reactor block J/alpha_k", blk_ok)

print("\nSUMMARY: %s (%d failure(s))" % ("ALL PASS" if fails == 0 else "FAILURES", fails))
sys.exit(1 if fails else 0)
