#!/usr/bin/env python3
"""Interval-arithmetic certificate: the explicit family of Theorem 2 has exactly 3^n steady states,
exactly 2^n of them asymptotically stable, for n = 1..6 (eps = 1/20) and n = 1..4 (eps = 1 and 5).

The certificate works directly on the steady-state equations of the n-reactor version of Bisbis (1987),
eqs. (5)-(8), in the thesis variables (y_k, theta_k):
  f_{y,k}  = (y_{k-1}-y_k)/alpha_k + Da0 tau (1-y_k) E(theta_k)
  f_{th,k} = (th_{k-1}-th_k)/alpha_k + Da0 tau B (1-y_k) E(theta_k) - (beta_k tau/alpha_k)(th_k - eta_k)
with E(th) = exp(th/(1+eps th)), y_0 = th_0 = 0, and the parameters of Theorem 2:
  B = 6400 n m/eps^2, Da0 = e^{-1/eps}/(4n), tau = 1, alpha_k = 1, beta_k = 200 m - 1,
  eta_k = (1-20 eps)/(20 eps^2), m = max(1, eps).

Steps (all in Arb ball arithmetic, 256-bit):
 1. Reactor by reactor, for every inlet ball, the scalar function F (Lemma 1) is shown to have the sign pattern
    (+,-,+,-) at z = 1/1000, 3/20, 3, 100; by Lemma 2 (at most three zeros) it has exactly three zeros, which are
    enclosed by bisection. This yields 3^n candidate states, labelled by words in {L,M,H}^n.
 2. Each candidate is verified independently by the Krawczyk test on the full 2n-dimensional system in the
    thesis variables: K(X) contained in the interior of X proves a unique zero of f in the box X.
 3. The 3^n boxes are pairwise disjoint (so there are at least 3^n steady states; Theorem 1 gives at most 3^n).
 4. Stability: the Jacobian is block lower triangular (symbolic check A7 in check_model_algebra.py, and the
    entries above the diagonal blocks are structurally zero here); on each box the 2x2 diagonal blocks are
    evaluated: labels L/H must give det > 0 and trace < 0 (Hurwitz), label M must give det < 0 (saddle).
    Hence the states with words in {L,H}^n (2^n of them) are hyperbolic sinks and all others have a positive
    real eigenvalue.
Each (n, eps) case prints PASS/FAIL; exit code 1 on any FAIL. Certificates (box midpoints/radii, block
traces/determinants) are written to certificates/cascade_n{n}_eps{eps}.json for eps = 1/20 and n = 3, 4, and a
summary of all cases to certificates/cascade_summary.json.
Environment: python-flint 0.9 (Arb), WSL Fedora 44 ~/envs/sci. Runtime: about 1 minute (n = 6: 729 states).
"""
import itertools
import json
import os
import sys
import time

from flint import arb, arb_mat, fmpq, ctx

ctx.prec = 256
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, '..', 'certificates')
fails = 0


def report(name, ok):
    global fails
    print(("PASS " if ok else "FAIL ") + name, flush=True)
    if not ok:
        fails += 1


def family(n, eps):
    m = max(eps, fmpq(1))
    return dict(n=n, eps=eps, B=fmpq(6400 * n) * m / eps**2,
                Da0=(-1 / arb(eps)).exp() / (4 * n), tau=fmpq(1),
                alpha=[fmpq(1)] * n, beta=[200 * m - 1] * n,
                eta=[(1 - 20 * eps) / (20 * eps**2)] * n)


# ---------------- thesis-variable system ----------------
def rhs(P, x):
    n, eps = P['n'], arb(P['eps'])
    B, Da0, tau = arb(P['B']), P['Da0'], arb(P['tau'])
    f = []
    for k in range(n):
        y, th = x[2 * k], x[2 * k + 1]
        yp, tp = (x[2 * k - 2], x[2 * k - 1]) if k > 0 else (arb(0), arb(0))
        al, be, et = arb(P['alpha'][k]), arb(P['beta'][k]), arb(P['eta'][k])
        E = (th / (1 + eps * th)).exp()
        f.append((yp - y) / al + Da0 * tau * (1 - y) * E)
        f.append((tp - th) / al + Da0 * tau * B * (1 - y) * E - (be * tau / al) * (th - et))
    return f


def jac(P, x):
    n, eps = P['n'], arb(P['eps'])
    B, Da0, tau = arb(P['B']), P['Da0'], arb(P['tau'])
    J = [[arb(0)] * (2 * n) for _ in range(2 * n)]
    for k in range(n):
        y, th = x[2 * k], x[2 * k + 1]
        al, be = arb(P['alpha'][k]), arb(P['beta'][k])
        E = (th / (1 + eps * th)).exp()
        Ep = E / (1 + eps * th)**2
        J[2 * k][2 * k] = -1 / al - Da0 * tau * E
        J[2 * k][2 * k + 1] = Da0 * tau * (1 - y) * Ep
        J[2 * k + 1][2 * k] = -Da0 * tau * B * E
        J[2 * k + 1][2 * k + 1] = -1 / al + Da0 * tau * B * (1 - y) * Ep - be * tau / al
        if k > 0:
            J[2 * k][2 * k - 2] = 1 / al
            J[2 * k + 1][2 * k - 1] = 1 / al
    return J


# ---------------- step 1: reactor-by-reactor isolation (reduced variables) ----------------
def reduced(P, k):
    eps = arb(P['eps'])
    K = arb(P['alpha'][k]) * P['Da0'] * arb(P['tau']) * (1 / eps).exp()
    lam = 1 / (1 + arb(P['beta'][k]) * arb(P['tau']))
    zc = eps * (1 + eps * arb(P['eta'][k]))
    Delta = eps**2 * arb(P['B'])
    return K, lam, zc, Delta


def Ffun(z, cin, zin, red):
    K, lam, zc, Delta = red
    k = K * (-1 / z).exp()
    return lam * zin + (1 - lam) * zc + lam * Delta * cin * k / (1 + k) - z


TEST = [arb(fmpq(1, 1000)), arb(fmpq(3, 20)), arb(3), arb(100)]


def isolate(cin, zin, red):
    vals = [Ffun(t, cin, zin, red) for t in TEST]
    pattern = (vals[0] > 0) and (vals[1] < 0) and (vals[2] > 0) and (vals[3] < 0)
    if not pattern:
        return None
    roots = []
    for i in range(3):
        l, r = TEST[i], TEST[i + 1]
        sl = 1 if i % 2 == 0 else -1          # sign of F at l
        for _ in range(400):
            mid = (l + r) / 2                 # exact dyadic at this precision
            v = Ffun(mid, cin, zin, red)
            if v > 0:
                s = 1
            elif v < 0:
                s = -1
            else:
                break
            if s == sl:
                l = mid
            else:
                r = mid
        roots.append(l.union(r))
    K = red[0]
    out = []
    for z in roots:
        c = cin / (1 + K * (-1 / z).exp())
        out.append((c, z))
    return out


def krawczyk(P, xt, rad):
    N = len(xt)
    X = [arb(xt[i], rad[i]) for i in range(N)]
    fx = arb_mat([[v] for v in rhs(P, xt)])
    Jx = jac(P, xt)
    Jmid = arb_mat([[Jx[i][j].mid() for j in range(N)] for i in range(N)])
    Yb = Jmid.inv()
    Y = arb_mat([[Yb[i, j].mid() for j in range(N)] for i in range(N)])
    JX = arb_mat(jac(P, X))
    I = arb_mat([[arb(1) if i == j else arb(0) for j in range(N)] for i in range(N)])
    D = arb_mat([[X[i] - xt[i]] for i in range(N)])
    Kv = arb_mat([[xt[i]] for i in range(N)]) - Y * fx + (I - Y * JX) * D
    ok = all(bool(abs(Kv[i, 0].mid() - X[i].mid()) + Kv[i, 0].rad() < X[i].rad()) for i in range(N))
    return ok, X, JX


def run_case(n, eps, save):
    t0 = time.time()
    P = family(n, eps)
    reds = [reduced(P, k) for k in range(n)]
    # step 1
    states = [((), arb(1), arb(P['eps']), [])]       # (word, c_in, z_in, list of (c,z))
    for k in range(n):
        new = []
        for word, cin, zin, hist in states:
            rts = isolate(cin, zin, reds[k])
            if rts is None:
                return dict(ok=False, reason='sign pattern failed at reactor %d, word %s' % (k + 1, word))
            for lab, (c, z) in zip('LMH', rts):
                new.append((word + (lab,), c, z, hist + [(c, z)]))
        states = new
    words_ok = len(states) == 3**n and len(set(s[0] for s in states)) == 3**n
    window_ok = True
    for word, _, _, hist in states:
        for lab, (c, z) in zip(word, hist):
            if lab == 'L':
                window_ok &= bool(z < arb(fmpq(3, 20)))
            elif lab == 'M':
                window_ok &= bool(z > arb(fmpq(3, 20))) and bool(z < 3)
            else:
                window_ok &= bool(z > 3)
    # step 2: Krawczyk in thesis variables
    epsb = arb(P['eps'])
    records = []
    kr_ok = True
    boxes = []
    stab_ok = True
    nstable = 0
    for word, _, _, hist in states:
        xt, rad = [], []
        for (c, z) in hist:
            yv = (1 - c).mid()
            thv = ((z - epsb) / epsb**2).mid()
            xt += [yv, thv]
        for v in xt:
            rad.append(arb(fmpq(1, 10**40)) * (1 + abs(v)))
        rad = [r.mid() for r in rad]
        ok, X, JX = krawczyk(P, xt, rad)
        kr_ok &= ok
        boxes.append(X)
        # step 4: stability from diagonal blocks evaluated on the box
        blocks = []
        all_hurwitz = True
        for k, lab in enumerate(word):
            a11, a12 = JX[2 * k, 2 * k], JX[2 * k, 2 * k + 1]
            a21, a22 = JX[2 * k + 1, 2 * k], JX[2 * k + 1, 2 * k + 1]
            tr, det = a11 + a22, a11 * a22 - a12 * a21
            if lab in 'LH':
                good = bool(det > 0) and bool(tr < 0)
            else:
                good = bool(det < 0)
                all_hurwitz = False
            stab_ok &= good
            blocks.append(dict(reactor=k + 1, label=lab, trace=tr.str(12, radius=True),
                               det=det.str(12, radius=True), certified=good))
        if all_hurwitz:
            nstable += 1
        records.append(dict(word=''.join(word), stable=all_hurwitz,
                            y=[X[2 * k].mid().str(30, radius=False) for k in range(n)],
                            theta=[X[2 * k + 1].mid().str(30, radius=False) for k in range(n)],
                            box_radius=[X[i].rad().str(5, radius=False) for i in range(2 * n)],
                            krawczyk=ok, blocks=blocks))
    # step 3: disjointness
    disj = True
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if all(boxes[i][q].overlaps(boxes[j][q]) for q in range(2 * n)):
                disj = False
    res = dict(n=n, eps=str(eps), states=len(states), words_ok=words_ok, window_ok=window_ok,
               krawczyk_all=kr_ok, disjoint=disj, stability_certified=stab_ok, stable=nstable,
               expected_states=3**n, expected_stable=2**n, seconds=round(time.time() - t0, 2))
    res['ok'] = (words_ok and window_ok and kr_ok and disj and stab_ok and nstable == 2**n)
    if save:
        out = dict(model='Bisbis (1987) eqs. (5)-(8), n reactors in series, thesis variables (y_k, theta_k)',
                   parameters=dict(eps=str(eps), B=str(P['B']), Da0='exp(-1/eps)/%d' % (4 * n), tau='1',
                                   alpha=[str(v) for v in P['alpha']], beta=[str(v) for v in P['beta']],
                                   eta=[str(v) for v in P['eta']]),
                   summary=res, states=records)
        fn = os.path.join(CERT, 'cascade_n%d_eps%s.json' % (n, str(eps).replace('/', '_')))
        with open(fn, 'w') as fh:
            json.dump(out, fh, indent=1)
    return res


if __name__ == '__main__':
    summary = []
    cases = [(fmpq(1, 20), n) for n in (1, 2, 3, 4, 5, 6)] + [(e, n) for e in (fmpq(1), fmpq(5)) for n in (1, 2, 3, 4)]
    for eps, n in cases:
        res = run_case(n, eps, save=(eps == fmpq(1, 20) and n in (3, 4)))
        summary.append(res)
        report("n=%d eps=%s: %s steady states (expected %d), all Krawczyk-verified=%s, pairwise disjoint=%s, "
               "stability certified=%s, stable=%s (expected %d)  [%.1f s]"
               % (n, eps, res.get('states'), 3**n, res.get('krawczyk_all'), res.get('disjoint'),
                  res.get('stability_certified'), res.get('stable'), 2**n, res.get('seconds', 0)),
               res['ok'])
    with open(os.path.join(CERT, 'cascade_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)
    print("\nSUMMARY: %s (%d failure(s))" % ("ALL PASS" if fails == 0 else "FAILURES", fails))
    sys.exit(1 if fails else 0)
