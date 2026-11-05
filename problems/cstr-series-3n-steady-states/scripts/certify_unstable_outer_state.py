#!/usr/bin/env python3
"""Certificate for Proposition 3: an outer steady state of one reactor can be unstable.

Single reactor (n = 1) of Bisbis (1987), eqs. (5)-(6), feed y_0 = theta_0 = 0, thesis variables (y, theta):
  f_y  = -y + Da0 tau (1-y) E(theta)
  f_th = -theta + Da0 tau B (1-y) E(theta) - beta tau (theta - eta),   E(th) = exp(th/(1+eps th)).
For each example below (exact rational parameters) the script
 1. finds three sign changes of the scalar function phi(theta) (Lemma 1/2: at most three zeros, so exactly three
    steady states), encloses the zeros by bisection,
 2. verifies each steady state by a Krawczyk test in (y, theta),
 3. evaluates trace, determinant and discriminant of the Jacobian on each box and prints the classification:
    lower state, middle state (det < 0: saddle), upper state.
Claims (PASS/FAIL, exit code 1 on failure):
  Example 1 (eps=1/20, B=20, Da0=1/10, tau=1, beta=2, eta=-1): exactly three steady states; the lower one is a
  hyperbolic sink, the middle one a saddle, the upper one an unstable focus (det > 0, trace > 0, trace^2 < 4 det).
  Example 2 (eps=1/20, B=30, Da0=69/500, tau=1, beta=5, eta=-7/20): exactly three steady states, all unstable:
  unstable node, saddle, unstable node (no stable steady state at all).
Environment: python-flint 0.9 (Arb), WSL Fedora 44 ~/envs/sci. Runtime: a few seconds.
"""
import json
import os
import sys

from flint import arb, arb_mat, fmpq, ctx

ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, '..', 'certificates')
fails = 0


def report(name, ok):
    global fails
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok:
        fails += 1


EXAMPLES = [
    dict(name='Example 1', eps=fmpq(1, 20), B=fmpq(20), Da0=fmpq(1, 10), tau=fmpq(1), beta=fmpq(2), eta=fmpq(-1),
         expect=['sink', 'saddle', 'unstable focus']),
    dict(name='Example 2', eps=fmpq(1, 20), B=fmpq(30), Da0=fmpq(69, 500), tau=fmpq(1), beta=fmpq(5),
         eta=fmpq(-7, 20), expect=['unstable node', 'saddle', 'unstable node']),
]


def make(P):
    eps, B, Da0, tau, beta, eta = (arb(P[k]) for k in ('eps', 'B', 'Da0', 'tau', 'beta', 'eta'))
    D, b = Da0 * tau, beta * tau

    def E(th):
        return (th / (1 + eps * th)).exp()

    def phi(th):
        e = E(th)
        return (b * eta + B * D * e / (1 + D * e)) / (1 + b) - th

    def rhs(x):
        y, th = x
        e = E(th)
        return [-y + D * (1 - y) * e, -th + D * B * (1 - y) * e - b * (th - eta)]

    def jac(x):
        y, th = x
        e = E(th)
        ep = e / (1 + eps * th)**2
        return [[-1 - D * e, D * (1 - y) * ep], [-B * D * e, -1 - b + B * D * (1 - y) * ep]]
    return E, phi, rhs, jac, D


def classify(tr, det):
    if det < 0:
        return 'saddle'
    if det > 0 and tr < 0:
        return 'sink'
    if det > 0 and tr > 0:
        disc = tr * tr - 4 * det
        if disc < 0:
            return 'unstable focus'
        if disc > 0:
            return 'unstable node'
        return 'unstable (undecided type)'
    return 'undecided'


results = []
for P in EXAMPLES:
    E, phi, rhs, jac, D = make(P)
    eps = P['eps']
    lo = -1 / eps + fmpq(1, 100)
    grid = [arb(lo + fmpq(i, 32)) for i in range(0, 32 * 150)]
    vals = [phi(t) for t in grid]
    decided = all((v > 0) or (v < 0) for v in vals)
    brackets = [(grid[i], grid[i + 1]) for i in range(len(grid) - 1) if bool(vals[i] > 0) != bool(vals[i + 1] > 0)]
    tail = phi(arb(10**6)) < 0
    report("%s: phi has decided signs on the grid, %d sign changes, phi(10^6) < 0" % (P['name'], len(brackets)),
           decided and len(brackets) == 3 and tail)
    states = []
    for (l, r) in brackets:
        sl = phi(l) > 0
        for _ in range(300):
            m = (l + r) / 2
            v = phi(m)
            if not ((v > 0) or (v < 0)):
                break
            if bool(v > 0) == sl:
                l = m
            else:
                r = m
        th = l.union(r).mid()
        e = E(arb(th))
        y = (D * e / (1 + D * e)).mid()
        xt = [y, th]
        rad = [(arb(fmpq(1, 10**40)) * (1 + abs(v))).mid() for v in xt]
        X = [arb(xt[i], rad[i]) for i in range(2)]
        fx = arb_mat([[v] for v in rhs([arb(v) for v in xt])])
        Jm = jac([arb(v) for v in xt])
        Y = arb_mat([[Jm[i][j].mid() for j in range(2)] for i in range(2)]).inv()
        Y = arb_mat([[Y[i, j].mid() for j in range(2)] for i in range(2)])
        JX = arb_mat(jac(X))
        I = arb_mat([[1, 0], [0, 1]])
        Kv = arb_mat([[xt[0]], [xt[1]]]) - Y * fx + (I - Y * JX) * arb_mat([[X[0] - xt[0]], [X[1] - xt[1]]])
        kok = all(bool(abs(Kv[i, 0].mid() - X[i].mid()) + Kv[i, 0].rad() < X[i].rad()) for i in range(2))
        tr = JX[0, 0] + JX[1, 1]
        det = JX[0, 0] * JX[1, 1] - JX[0, 1] * JX[1, 0]
        states.append(dict(y=xt[0].str(25, radius=False), theta=xt[1].str(25, radius=False), krawczyk=kok,
                           trace=tr.str(12), det=det.str(12), discriminant=(tr * tr - 4 * det).str(12),
                           type=classify(tr, det)))
    kr = all(s['krawczyk'] for s in states)
    types = [s['type'] for s in states]
    report("%s: three Krawczyk-verified steady states, types %s (expected %s)" % (P['name'], types, P['expect']),
           kr and types == P['expect'])
    for s in states:
        print("     theta = %s  y = %s  tr = %s  det = %s  -> %s"
              % (s['theta'][:14], s['y'][:14], s['trace'], s['det'], s['type']))
    results.append(dict(example=P['name'], parameters={k: str(P[k]) for k in ('eps', 'B', 'Da0', 'tau', 'beta', 'eta')},
                        states=states))

with open(os.path.join(CERT, 'unstable_outer_state.json'), 'w') as fh:
    json.dump(results, fh, indent=1)
print("\nSUMMARY: %s (%d failure(s))" % ("ALL PASS" if fails == 0 else "FAILURES", fails))
sys.exit(1 if fails else 0)
