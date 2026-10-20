#!/usr/bin/env python3
"""Early verifier: symbolic / high-precision checks of the algebraic ingredients (claims 1-2).

(A) sympy, N = 2..5: for a liquid x (sum x = 1) and the vapor y = alpha x / S(x) in equilibrium with it,
    phi_liq/l = -G_x(theta) and phi_vap/v = -(theta/S) G_x(theta) identically (sign uniformity of E units);
    and the vapor of the SAME composition x (total vaporizer outlet) has phi/v = 1 - G_x(theta) (so a unit
    releasing a liquid and a vapor of equal composition would not be sign-uniform: it must be (M) + (R)).
(B) sympy: node balance.  For a generic unit with symbolic inflow/outflow component flows obeying the component
    balances and the CMO energy balance, sum_out phi - sum_in phi - Q == 0, feed term included, any q.
(C) mpmath (50 digits): feed roots, phi_F(theta_j) = 0 for q in {-1.5,-0.3,0,0.4,1,1.7}; Corollary 5
    (R >= H_j - c_j eps) on 4000 random impure product sets (liquid and vapor products, several streams per
    component, impurities moved between products so that the overall balances hold exactly), including
    adversarial impurity placements; ratio (H_j - B_j)/(c_j eps) for the worst binary placement.
Prints PASS/FAIL per claim, exits 1 on failure.
Environment: Python 3, sympy, mpmath.  Runtime: about 1 minute.
"""
import sys
import random
import sympy as sp
import mpmath as mp

ok_all = True


def report(name, ok):
    global ok_all
    ok_all &= bool(ok)
    print(('PASS ' if ok else 'FAIL ') + name)


# ---------------------------------------------------------------- (A)
th = sp.Symbol('theta', positive=True)
okA = True
for N in range(2, 6):
    al = sp.symbols('a1:%d' % (N + 1), positive=True)
    xs = sp.symbols('x1:%d' % N, positive=True)
    x = list(xs) + [1 - sum(xs)]
    S = sum(a * xi for a, xi in zip(al, x))
    G = sum(a * xi / (a - th) for a, xi in zip(al, x))
    y = [a * xi / S for a, xi in zip(al, x)]
    phi_v = 1 - sum(a * yi / (a - th) for a, yi in zip(al, y))
    phi_l = -sum(a * xi / (a - th) for a, xi in zip(al, x))
    phi_v_same = 1 - sum(a * xi / (a - th) for a, xi in zip(al, x))
    e1 = sp.simplify(sp.together(phi_v + th * G / S))
    e2 = sp.simplify(phi_l + G)
    e3 = sp.simplify(phi_v_same - (1 - G))
    okA &= (e1 == 0 and e2 == 0 and e3 == 0)
report('(A) phi_vap/phi_liq per unit flow = theta/S(x) > 0 for equilibrium phases, N=2..5', okA)

# ---------------------------------------------------------------- (B)
N = 3
al = sp.symbols('b1:4', positive=True)
q, f = sp.symbols('q f', real=True)
z = sp.symbols('z1:4', positive=True)
# two liquid inflows, one vapor inflow, outflows: liquid (2 streams) + vapor (1 stream); component flows symbolic
Lin = [sp.symbols('l%d_1:4' % k) for k in range(2)]
Vin = [sp.symbols('v0_1:4')]
Lout = [sp.symbols('m%d_1:4' % k) for k in range(2)]
Vout_free = sp.symbols('w_1:4')
# component balance fixes the vapor outflow: out = in + feed
Vout = [[sum(L[i] for L in Lin) + sum(V[i] for V in Vin) + f * z[i] - sum(L[i] for L in Lout) for i in range(N)]]
Q = sum(sum(V) for V in Vout) - sum(sum(V) for V in Vin) - (1 - q) * f      # CMO energy balance defines Q


def phi(n, vap):
    return (sum(n) if vap else 0) - sum(a * ni / (a - th) for a, ni in zip(al, n))


phiF = f * ((1 - q) - sum(a * zi / (a - th) for a, zi in zip(al, z)))
bal = sum(phi(n, False) for n in Lout) + sum(phi(n, True) for n in Vout) - sum(phi(n, False) for n in Lin) \
    - sum(phi(n, True) for n in Vin) - phiF - Q
report('(B) node balance sum_out phi - sum_in phi - phi_F = Q (symbolic, any q)', sp.simplify(bal) == 0)

# ---------------------------------------------------------------- (C)
mp.mp.dps = 50
random.seed(7)


def roots(a, z, q):
    g = lambda t: sum(ai * zi / (ai - t) for ai, zi in zip(a, z)) - (1 - q)
    out = []
    for j in range(len(a) - 1):
        lo, hi = a[j + 1], a[j]          # g increases from -inf to +inf on (lo, hi): bisection
        for _ in range(200):
            mid = (lo + hi) / 2
            if g(mid) > 0:
                hi = mid
            else:
                lo = mid
        out.append((lo + hi) / 2)
    return out


def phiP(a, n, vap, t):
    return (sum(n) if vap else 0) - sum(ai * ni / (ai - t) for ai, ni in zip(a, n))


okF = True
okC = True
worst = mp.inf
tight = []
for trial in range(4000):
    N = random.randint(2, 5)
    a = sorted([mp.mpf(random.randint(11, 80)) / 10 for _ in range(N)], reverse=True)
    if len(set(a)) < N:
        continue
    zr = [mp.mpf(random.random()) + mp.mpf('0.05') for _ in range(N)]
    z = [v / mp.fsum(zr) for v in zr]
    q = mp.mpf(random.choice([-1.5, -0.3, 0, 0.4, 1, 1.7]))
    ths = roots(a, z, q)
    for t in ths:
        okF &= abs((1 - q) - sum(ai * zi / (ai - t) for ai, zi in zip(a, z))) < mp.mpf('1e-35')
    # products: each component p split into 1-2 streams of random phase
    prods = []
    for p in range(N):
        k = random.randint(1, 2)
        fr = [mp.mpf(random.random()) for _ in range(k)]
        for fk in fr:
            n = [mp.mpf(0)] * N
            n[p] = z[p] * fk / mp.fsum(fr)
            prods.append([p, n, random.random() < 0.3])
    # move impurities: amounts of component i moved from a product assigned to i to another product
    eps_target = mp.mpf(10) ** (-random.uniform(1, 4))
    for _ in range(random.randint(1, 3 * N)):
        src = random.choice([P for P in prods if sum(P[1]) > 0])
        i = src[0]
        dst = random.choice(prods)
        if dst is src:
            continue
        amt = src[1][i] * eps_target * mp.mpf(random.random())
        src[1][i] -= amt
        dst[1][i] += amt
    eps = max(1 - P[1][P[0]] / sum(P[1]) for P in prods)
    for j, t in enumerate(ths):
        H = sum(a[p] * z[p] / (t - a[p]) for p in range(j + 1, N))
        c = max(a[j] / (a[j] - t), a[j + 1] / (t - a[j + 1]))
        B = sum(max(mp.mpf(0), phiP(a, P[1], P[2], t)) for P in prods)
        slack = B - (H - c * eps)
        worst = min(worst, slack)
        if slack < -mp.mpf('1e-30'):
            okC = False
            print('   violation', N, [mp.nstr(v, 6) for v in a], mp.nstr(q, 3), j, mp.nstr(slack, 5), mp.nstr(eps, 5),
                  'sum|P| =', mp.nstr(mp.fsum(sum(P[1]) for P in prods) - 1, 5))
# tightness: binary, light component moved into the heavy liquid product at purity 1 - eps
a = [mp.mpf(3), mp.mpf(1)]; z = [mp.mpf('0.4'), mp.mpf('0.6')]; q = mp.mpf(1)
t = roots(a, z, q)[0]
H = a[1] * z[1] / (t - a[1])
c = max(a[0] / (a[0] - t), a[1] / (t - a[1]))
for e in (mp.mpf('1e-3'), mp.mpf('1e-5')):
    # worst placement: light component 1 into the heavy product, heavy component 2 into the light product
    m, r = e * z[1], e * z[0]
    light = [z[0] - m, r]
    heavy = [m, z[1] - r]
    eps = max(1 - light[0] / sum(light), 1 - heavy[1] / sum(heavy))
    B = max(0, phiP(a, light, False, t)) + max(0, phiP(a, heavy, False, t))
    tight.append((H - B) / (c * eps))
report('(C1) phi_F(theta_j) = 0 at the feed roots for q in {-1.5,...,1.7}', okF)
report('(C2) Corollary 5: sum_P max(0,phi_P(theta_j)) >= H_j - c_j eps on 4000 random impure product sets '
       '(min slack %s)' % mp.nstr(worst, 5), okC)
print('     worst-placement ratio (H_j - B_j)/(c_j eps) for a binary example (must be <= 1; c_j is of the right order):',
      [mp.nstr(v, 5) for v in tight])
print('ALL PASS' if ok_all else 'FAIL')
sys.exit(0 if ok_all else 1)
