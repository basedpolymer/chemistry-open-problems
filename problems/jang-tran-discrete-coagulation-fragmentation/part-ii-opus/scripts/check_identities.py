#!/usr/bin/env python3
"""Exact symbolic checks of every algebraic identity used in part-ii.md.

Environment: Python 3 with sympy (tested with the project's WSL 'sci' environment).
Runtime: about 1-2 minutes on one core. Deterministic. Prints PASS/FAIL per claim and
exits with a non-zero status if any check fails.

Notation (see part-ii.md): z = exp(-x), tau(x) = 2 tanh(x/2), kappa = 1/2 + 1/(e^x - 1),
H_m(x,F,p) = (p-m)(p-m-1)/2 + kappa F - m.
"""
import sys
import sympy as sp

FAIL = 0


def report(name, ok):
    global FAIL
    print(("PASS " if ok else "FAIL ") + name)
    if not ok:
        FAIL += 1


def is_zero(e):
    return sp.simplify(sp.together(sp.expand(e))) == 0


x, z, m = sp.symbols('x z m', positive=True)
X, P, U = sp.symbols('X P U', real=True)

# ---------------------------------------------------------------- A1: kappa = 1/tau
tau_x = 2 * sp.tanh(x / 2)
kap_x = sp.Rational(1, 2) + 1 / (sp.exp(x) - 1)
report("A1 kappa(x) = 1/2 + 1/(e^x-1) equals 1/tau(x), tau = 2 tanh(x/2)",
       sp.simplify((kap_x - 1 / tau_x).rewrite(sp.exp)) == 0)
# in the z variable: tau = 2(1-z)/(1+z), kappa = (1+z)/(2(1-z)), tau' = 1 - tau^2/4
tau_z = 2 * (1 - z) / (1 + z)
report("A1b tau(x) = 2(1-z)/(1+z) with z = e^{-x}",
       sp.simplify((tau_x - tau_z.subs(z, sp.exp(-x))).rewrite(sp.exp)) == 0)
report("A1c tau' = 1 - tau^2/4",
       sp.simplify((sp.diff(tau_x, x) - (1 - tau_x**2 / 4)).rewrite(sp.exp)) == 0)

# ------------------------------------------- A2: Hamilton-Jacobi equation of the full system
# For a finitely supported rho (sizes 1..K), compute the right-hand side of (1.1) (support up
# to 2K) and check F_t = -H_m(F) with m = sum j rho_j, F = sum (1 - z^j) rho_j.
K = 5
r = sp.symbols('r1:%d' % (K + 1), nonnegative=True)
rho = {j: r[j - 1] for j in range(1, K + 1)}
L = 2 * K


def rhs_full(rho, L):
    """Right-hand side of (1.1) with a(j,k)=jk, b=1 for finitely supported rho."""
    get = lambda j: rho.get(j, 0)
    m1 = sum(j * get(j) for j in range(1, L + 1))
    out = {}
    for l in range(1, L + 1):
        gain = sp.Rational(1, 2) * sum(i * (l - i) * get(i) * get(l - i) for i in range(1, l))
        loss = l * get(l) * m1
        fl = -sp.Rational(l - 1, 2) * get(l) + sum(get(j) for j in range(l + 1, L + 1))
        out[l] = sp.expand(gain - loss + fl)
    return out


d = rhs_full(rho, L)
mm = sum(j * rho[j] for j in rho)
F = sum((1 - z**j) * rho[j] for j in rho)
Fx = sum(j * z**j * rho[j] for j in rho)          # d/dx of F, since dz/dx = -z
Ft = sum((1 - z**l) * d[l] for l in d)
kap_z = (1 + z) / (2 * (1 - z))
Hm = sp.Rational(1, 2) * (Fx - mm) * (Fx - mm - 1) + kap_z * F - mm
report("A2 transform of (1.1) solves F_t + H_m(F) = 0 with m = m_1(rho) (K=5 generic sizes)",
       is_zero((Ft + Hm) * (1 - z)))
# mass and zeroth moment
report("A2b d/dt sum l rho_l = 0 (finitely supported rho)",
       is_zero(sum(l * d[l] for l in d)))
m0 = sum(rho.values())
report("A2c d/dt m_0 = (m - m^2 - m_0)/2  (Jang-Tran (1.6))",
       is_zero(sum(d[l] for l in d) - (mm - mm**2 - m0) / 2))

# --------------------------------------------------- A3: truncated system (sizes 1..N)
N = 6
s = sp.symbols('s1:%d' % (N + 1), nonnegative=True)
rhoN = {j: s[j - 1] for j in range(1, N + 1)}


def rhs_trunc(rho, N):
    get = lambda j: rho.get(j, 0)
    out = {}
    for l in range(1, N + 1):
        gain = sp.Rational(1, 2) * sum(i * (l - i) * get(i) * get(l - i) for i in range(1, l))
        loss = l * get(l) * sum(k * get(k) for k in range(1, N - l + 1))
        fl = -sp.Rational(l - 1, 2) * get(l) + sum(get(j) for j in range(l + 1, N + 1))
        out[l] = sp.expand(gain - loss + fl)
    return out


dN = rhs_trunc(rhoN, N)
mN = sum(j * rhoN[j] for j in rhoN)
FN = sum((1 - z**j) * rhoN[j] for j in rhoN)
FNx = sum(j * z**j * rhoN[j] for j in rhoN)
FNt = sum((1 - z**l) * dN[l] for l in dN)
EN = sp.Rational(1, 2) * sum((1 - z**j) * (1 - z**k) * j * k * rhoN[j] * rhoN[k]
                             for j in range(1, N + 1) for k in range(1, N + 1) if j + k > N)
HmN = sp.Rational(1, 2) * (FNx - mN) * (FNx - mN - 1) + kap_z * FN - mN
report("A3 truncated system: F^N_t + H_m(F^N) = E_N, E_N = (1/2) sum_{j+k>N} (1-z^j)(1-z^k) jk rho_j rho_k",
       is_zero((FNt + HmN - EN) * (1 - z)))
report("A3b truncated system conserves sum l rho_l",
       is_zero(sum(l * dN[l] for l in dN)))

# ------------------------------------------------ A4: identities along characteristics
tau = 2 * sp.tanh(X / 2)
tp = sp.diff(tau, X)
kap = 1 / tau
H = sp.Rational(1, 2) * (P - m) * (P - m - 1) + kap * U - m
Xd = sp.diff(H, P)
Pd = -sp.diff(H, X) - P * sp.diff(H, U)
Ud = P * sp.diff(H, P) - H


def Dt(e):
    return sp.diff(e, X) * Xd + sp.diff(e, P) * Pd + sp.diff(e, U) * Ud


def zero_exp(e):
    return sp.simplify(sp.expand(e).rewrite(sp.exp)) == 0


q = m - P
report("A4a dX/dt = -(q + 1/2)", zero_exp(Xd + (q + sp.Rational(1, 2))))
report("A4b dP/dt = tau' U/tau^2 - P/tau", zero_exp(Pd - (tp * U / tau**2 - P / tau)))
report("A4c dU/dt = (P^2 + m - m^2)/2 - U/tau",
       zero_exp(Ud - ((P**2 + m - m**2) / 2 - U / tau)))
W = tau * (m + P) - 2 * U
report("A4d dW/dt + (2+tau')/(2 tau) W = q (2m - q) (1 - tau'),  1 - tau' = tau^2/4",
       zero_exp(Dt(W) + (2 + tp) / (2 * tau) * W - q * (2 * m - q) * (1 - tp))
       and zero_exp(1 - tp - tau**2 / 4))
# Riccati for r = F_xx via the variational system (J, Pi) with strip condition U_y = P J
J, Pi, rr = sp.symbols('J Pi r', real=True)


def lin(f):
    return sp.diff(f, X) * J + sp.diff(f, P) * Pi + sp.diff(f, U) * P * J


Jd, Pid, Uyd = lin(Xd), lin(Pd), lin(Ud)
report("A4e strip condition A = U_y - P J satisfies dA/dt = 0 when A = 0 (in fact dA/dt = -H_U A)",
       zero_exp(Uyd - (Pd * J + P * Jd)))
rd = sp.simplify(((Pid * J - Pi * Jd) / J**2).subs(Pi, rr * J))
report("A4f Riccati: dr/dt = -r^2 - kappa r - 2 kappa' P - kappa'' U",
       zero_exp(rd - (-rr**2 - kap * rr - 2 * sp.diff(kap, X) * P - sp.diff(kap, X, 2) * U)))
G = tau * rr
Gd = tp * Xd * rr + tau * rd
Dq = U / tau - P
report("A4g barrier: tau dGamma/dt = -Gamma^2 - Gamma (1 + tau'(q+1/2)) - 2 tau' D,  D = U/tau - P",
       zero_exp(tau * Gd - (-G**2 - G * (1 + tp * (q + sp.Rational(1, 2))) - 2 * tp * Dq)))
report("A4h W = tau (q - 2D)  (so W >= 0 iff D <= q/2)", zero_exp(W - tau * (q - 2 * Dq)))
eta = m - U / tau
eps = tau**2 / 4 * (q + sp.Rational(1, 2))
report("A4i d eta/dt = -(1/tau)[(1/2 - q + eps) eta + q^2/2] + m eps/tau,  eta = m - U/tau",
       zero_exp(Dt(eta) - (-(1 / tau) * ((sp.Rational(1, 2) - q + eps) * eta + q**2 / 2) + m * eps / tau)))
lam = U / tau
report("A4k d lambda/dt = [(P-lambda)^2 + (m-lambda)(1-m+lambda)]/(2 tau) - eps lambda/tau,  lambda = U/tau",
       zero_exp(Dt(lam) - (((P - lam)**2 + (m - lam) * (1 - m + lam)) / (2 * tau) - eps * lam / tau)))
report("A4l (1/2-q)eta + q^2/2 = [(eta-q)^2 + eta(1-eta)]/2",
       sp.expand((sp.Rational(1, 2) - q) * eta + q**2 / 2 - ((eta - q)**2 + eta * (1 - eta)) / 2) == 0)
Gam1 = -1
val = (-G**2 - G * (1 + tp * (q + sp.Rational(1, 2))) - 2 * tp * Dq).subs(rr, Gam1 / tau)
report("A4j at Gamma = -1 the barrier right-hand side equals tau'(q + 1/2 - 2D)",
       zero_exp(val - tp * (q + sp.Rational(1, 2) - 2 * Dq)))

# --------------------------------------------- A5: initial half-slope inequality, per size j
ok = True
for j in range(1, 61):
    lhs = (tau_z * j * (1 + z**j) - 2 * (1 - z**j)) * (1 + z)
    rhs = 2 * (1 - z) * sum((1 - z**i) * (1 - z**(j - i)) for i in range(1, j))
    if sp.cancel(lhs - rhs) != 0:
        ok = False
        print("   mismatch at j =", j)
report("A5 tau(x) j (1+z^j) - 2(1-z^j) = 2(1-z)/(1+z) sum_{i=1}^{j-1} (1-z^i)(1-z^{j-i}), j = 1..60",
       ok)
# general j: the bracket identity j(1+z^j) - (1+z)(1+...+z^{j-1}) = sum (1-z^i)(1-z^{j-i}) is the
# telescoping (1+z)(1+...+z^{j-1}) = 1 + 2(z+...+z^{j-1}) + z^j; checked symbolically for j<=60.

# ------------------------------------------------------- A6: elementary inequalities
qq, ee = sp.symbols('q eta', real=True)
report("A6a (m-2q) eta + q^2 - (m-q)^2 = (m-2q)(eta-m)  [so (m-2q)eta+q^2 >= 0 on 0<=q,eta<=m]",
       sp.expand((m - 2 * qq) * ee + qq**2 - (m - qq)**2 - (m - 2 * qq) * (ee - m)) == 0)
report("A6b (1/2-q) eta + q^2/2 - (1-m)/2 eta = ((m-2q) eta + q^2)/2",
       sp.expand((sp.Rational(1, 2) - qq) * ee + qq**2 / 2 - (1 - m) / 2 * ee
                 - ((m - 2 * qq) * ee + qq**2) / 2) == 0)
# U >= 0 face: dU/dt at U=0 equals (P^2 + m(1-m))/2 >= 0 for m <= 1 (A4c). P-faces: A4b.

# ----------------------------------------------- A7: closed system and recursion (1.9)
l = sp.symbols('l', positive=True)
report("A7 2 c_l = (2m+1) l + 1 with c_l = m l + (l+1)/2",
       sp.expand(2 * (m * l + (l + 1) / 2) - ((2 * m + 1) * l + 1)) == 0)

print("\nSUMMARY:", "ALL PASS" if FAIL == 0 else "%d FAIL" % FAIL)
sys.exit(1 if FAIL else 0)
