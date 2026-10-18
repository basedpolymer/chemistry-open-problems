#!/usr/bin/env python3
"""Independent verifier's checks for the DGH charged di-block (written from scratch by the verifier;
does not import or reuse any solver script).

Environment: WSL Fedora, Python >= 3.10, numpy, mpmath, python-flint. One thread, < 500 MB RAM.
Runtime: about 10 seconds. Prints PASS/FAIL per claim, exits non-zero on failure.

What is checked (finite / numerical evidence, not a proof):
  1. Original pair Hamiltonian by brute-force path enumeration (N <= 14) == an independent
     site-charge transfer-matrix computation (exact integer polynomials), for the diblock and for
     random words; conversion to Poisat's alternative convention (2.13).
  2. The exact mixed-block identity (16) of the manuscript, as an identity of integer polynomials,
     for 2n <= 40, using the transfer-matrix polynomials (not the composition recursion).
  3. The threshold: beta_* from mpmath (60 digits) and from the Jacobi theta function
     theta_3(0, e^{-beta/2}) = 3; a rigorous Arb bracket with a verifier-derived tail bound;
     the Poisson-summation explanation of beta_* - 2*pi/9 ~ (8*pi/9) e^{-9*pi}.
  4. Cut probabilities by a full-chain forward/backward DP (log-space, float) for random and adversarial
     binary words at beta < beta_*, compared with the claimed uniform bound delta_beta of (15).
  5. Di-block at large n by full-chain DP vs. the formula p_{n,2n} = A_n^2/D_n; 3 n p_{n,2n} -> 1 at
     beta_*; p_{n,2n} -> 1/K_beta below; exponential decay above; free energy formula (4).
"""
import itertools
import math
import random
import sys
import time
from collections import defaultdict
from fractions import Fraction

import mpmath as mp
import numpy as np

FAIL = False


def check(ok, label):
    global FAIL
    print(('PASS ' if ok else 'FAIL ') + label, flush=True)
    if not ok:
        FAIL = True


# ---------------------------------------------------------------- 1. brute force & transfer matrix
def brute_force_H_poly(charges):
    """dict H -> number of paths; H = sum_{i<j} w_i w_j 1{s_i=s_j}; s_1=1, increments in {0,1}."""
    n = len(charges)
    out = defaultdict(int)
    for steps in itertools.product((0, 1), repeat=n - 1):
        pos = [1]
        for st in steps:
            pos.append(pos[-1] + st)
        H = 0
        for i in range(n):
            for j in range(i + 1, n):
                if pos[i] == pos[j]:
                    H += charges[i] * charges[j]
        out[H] += 1
    return dict(out)


def transfer_H_poly(charges):
    """Site-charge transfer matrix: state = charge of the current (last) site; exact integer counts."""
    state = {charges[0]: {0: 1}}  # current-site charge -> {H: count}
    for w in charges[1:]:
        new = defaultdict(lambda: defaultdict(int))
        for c, poly in state.items():
            for H, cnt in poly.items():
                new[c + w][H + w * c] += cnt   # stay on the same site: new pairs w*c
                new[w][H] += cnt               # step to a new site
        state = {c: dict(p) for c, p in new.items()}
    out = defaultdict(int)
    for poly in state.values():
        for H, cnt in poly.items():
            out[H] += cnt
    return dict(out)


def to_reduced(Hpoly, n):
    """reduced weight y^{sum_r Q_r^2} = y^{n + 2H} with y = e^{-beta/2}  (D = e^{-beta n/2} Z)."""
    return {n + 2 * H: c for H, c in Hpoly.items()}


def polymul(a, b, shift=0):
    out = defaultdict(int)
    for i, x in a.items():
        for j, y in b.items():
            out[i + j + shift] += x * y
    return dict(out)


def polyadd(a, b):
    out = defaultdict(int, a)
    for k, v in b.items():
        out[k] += v
    return {k: v for k, v in out.items() if v}


def part1_2():
    rng = random.Random(12345)
    for N in range(2, 15):
        words = [tuple([1] * (N // 2) + [-1] * (N - N // 2))] + [tuple(rng.choice((-1, 1)) for _ in range(N)) for _ in range(3)]
        ok = all(brute_force_H_poly(w) == transfer_H_poly(w) for w in words)
        check(ok, f'N={N}: brute-force pair Hamiltonian == site-charge transfer matrix (diblock + 3 random words)')
    # Poisat (2.13): Zbar_n(beta) = 2^{1-n} e^{-beta n} Z_n(2 beta), with Zbar = E exp(-beta sum_{i,j} w_i w_j 1{S_i=S_j})
    # computed directly from the ordered double sum over all (i,j) (diagonal included), path by path.
    mp.mp.dps = 50
    for w in [(1, 1, -1, -1), (1, -1, -1, 1, 1), (1, 1, 1, -1, -1, -1), (1, -1, 1, 1, -1, -1, 1)]:
        n = len(w)
        beta = mp.mpf('0.37')
        Z2b = mp.mpf(0)
        Zbar = mp.mpf(0)
        for steps in itertools.product((0, 1), repeat=n - 1):
            pos = [1]
            for st in steps:
                pos.append(pos[-1] + st)
            H = sum(w[i] * w[j] for i in range(n) for j in range(i + 1, n) if pos[i] == pos[j])
            Halt = sum(w[i] * w[j] for i in range(n) for j in range(n) if pos[i] == pos[j])
            Z2b += mp.e ** (-2 * beta * H)
            Zbar += mp.e ** (-beta * Halt) / 2 ** (n - 1)
        check(abs(Zbar - 2 ** (1 - n) * mp.e ** (-beta * n) * Z2b) < mp.mpf(10) ** -40,
              f'Poisat (2.13) convention identity (direct ordered double sum) for word {w}')
    # (16) as integer-polynomial identity in y, 2n <= 40
    homo = [{0: 1}] + [to_reduced(transfer_H_poly([1] * m), m) for m in range(1, 21)]
    for n in range(1, 21):
        rhs = polymul(homo[n], homo[n])
        for i in range(n):
            for j in range(n):
                rhs = polyadd(rhs, polymul(homo[i], homo[j], (i - j) ** 2))
        lhs = to_reduced(transfer_H_poly([1] * n + [-1] * n), 2 * n)
        if lhs != rhs:
            check(False, f'(16) polynomial identity, 2n={2 * n}')
            return homo
        # central cut count: paths with a cut at the interface = A_n^2
    check(True, '(16) exact integer-polynomial identity for all 2n <= 40 (transfer-matrix side)')
    return homo


# ---------------------------------------------------------------- 3. threshold
def part3():
    mp.mp.dps = 60
    W1 = lambda b: mp.nsum(lambda l: mp.e ** (-b * l * l / 2), [1, mp.inf]) - 1
    bstar = mp.findroot(W1, mp.mpf('0.7'))
    claimed = mp.mpf('0.698131700799199439091341370180112797')
    check(abs(bstar - claimed) < mp.mpf(10) ** -35, f'mpmath root of sum e^(-b l^2/2)=1: {mp.nstr(bstar, 40)}')
    th = mp.findroot(lambda b: mp.jtheta(3, 0, mp.e ** (-b / 2)) - 3, mp.mpf('0.7'))
    check(abs(th - bstar) < mp.mpf(10) ** -50, 'theta_3(0,e^(-b/2)) = 3 gives the same root')
    diff = bstar - 2 * mp.pi / 9
    approx = 8 * mp.pi / 9 * mp.e ** (-9 * mp.pi)
    print(f'INFO beta_* - 2pi/9 = {mp.nstr(diff, 20)};  (8pi/9) e^(-9pi) = {mp.nstr(approx, 20)}')
    check(diff > 0 and abs(diff / approx - 1) < 1e-2, 'beta_* - 2pi/9 = (8pi/9)e^(-9pi)(1+O(e^(-9pi)))-type agreement (<1%)')
    check(abs(diff - mp.mpf('1.4676083e-12')) < mp.mpf('1e-18'), 'difference ~ 1.4676083e-12 as claimed')
    # Poisson summation: sum_{l in Z} e^{-b l^2/2} = sqrt(2 pi/b) sum_{k in Z} e^{-2 pi^2 k^2/b}; at b=2pi/9 => 3(1+2e^{-9pi}+...) > 3
    b = 2 * mp.pi / 9
    lhs = 1 + 2 * mp.nsum(lambda l: mp.e ** (-b * l * l / 2), [1, mp.inf])
    rhs = mp.sqrt(2 * mp.pi / b) * (1 + 2 * mp.nsum(lambda k: mp.e ** (-2 * mp.pi ** 2 * k * k / b), [1, mp.inf]))
    check(abs(lhs - rhs) < mp.mpf(10) ** -50, 'Poisson summation identity at 2pi/9 (Theta(2pi/9) = 3(1+2e^{-9pi}+...))')
    # Rigorous Arb bracket with a tail bound derived here: for l >= L+1 the ratio of consecutive
    # terms is e^{-b(2l+1)/2} <= e^{-b(2L+3)/2}, so tail <= first/(1-ratio).
    from flint import arb, ctx
    ctx.prec = 300
    def bounds(bb, L=40):
        s = sum(((-bb * k * k / 2).exp() for k in range(1, L + 1)), arb(0))
        first = (-bb * (L + 1) ** 2 / 2).exp()
        ratio = (-bb * (2 * L + 3) / 2).exp()
        return s, s + first / (1 - ratio)
    lo = arb('0.698131700799199439091341370180112797008487712347')
    hi = arb('0.698131700799199439091341370180112797008487712348')
    check(bounds(lo)[0] > 1 and bounds(hi)[1] < 1, 'Arb: 0.698131700799199439091341370180112797008487712347 < beta_* < ...348')
    pi9 = 2 * arb.pi() / 9
    check(bounds(pi9)[0] > 1, 'Arb: sum at 2pi/9 strictly > 1, hence beta_* > 2pi/9 (strict monotonicity)')
    return float(bstar)


# ---------------------------------------------------------------- helpers (float)
def r_of_beta(beta, L=200):
    w = [math.exp(-beta * l * l / 2) for l in range(1, L + 1)]
    f = lambda z: sum(wl * z ** l for l, wl in enumerate(w, 1)) - 1
    lo, hi = 1e-9, 1.0
    while f(hi) < 0:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def scaled_A(beta, M, L=120):
    """u_m = A_m r^m via the renewal recurrence with probabilities p_l = w_l r^l (float)."""
    r = r_of_beta(beta)
    p = np.array([0.0] + [math.exp(-beta * l * l / 2) * r ** l for l in range(1, L + 1)])
    u = np.zeros(M + 1)
    u[0] = 1.0
    for m in range(1, M + 1):
        k = min(L, m)
        u[m] = np.dot(p[1:k + 1], u[m - 1::-1][:k])
    return r, u


def cut_probs(charges, beta):
    """Full-chain forward/backward DP in log space for the reduced weights; returns p_i, i=1..N-1."""
    w = np.asarray(charges, dtype=float)
    N = len(w)
    pre = np.concatenate([[0.0], np.cumsum(w)])
    logF = np.full(N + 1, -np.inf)
    logF[0] = 0.0
    for j in range(1, N + 1):
        k = np.arange(j)
        terms = logF[:j] - beta * (pre[j] - pre[k]) ** 2 / 2
        mx = terms.max()
        logF[j] = mx + math.log(np.exp(terms - mx).sum())
    logB = np.full(N + 1, -np.inf)  # logB[i] = log Q(i, N]
    logB[N] = 0.0
    for i in range(N - 1, -1, -1):
        k = np.arange(i + 1, N + 1)
        terms = logB[i + 1:] - beta * (pre[k] - pre[i]) ** 2 / 2
        mx = terms.max()
        logB[i] = mx + math.log(np.exp(terms - mx).sum())
    p = np.exp(logF[1:N] + logB[1:N] - logF[N])
    return p, logF[N]


# ---------------------------------------------------------------- 4. uniform lower bound (2)
def part4(bstar):
    rng = np.random.default_rng(7)
    for beta in [0.3, 0.5, 0.65, 0.69]:
        r, u = scaled_A(beta, 4000)
        c = u.min()
        S = r / (c * (1 - r))
        delta = 1 / (1 + S * S)
        worst = 1.0
        words = []
        for N in (50, 400):
            words.append(np.ones(N))
            words.append(np.array([1] * (N // 2) + [-1] * (N // 2)))
            words.append(np.array(([1] * 5 + [-1] * 5) * (N // 10)))
            for _ in range(4):
                words.append(rng.choice([-1, 1], size=N))
        for wd in words:
            p, _ = cut_probs(wd, beta)
            worst = min(worst, p.min())
        print(f'INFO beta={beta}: r={r:.6f}, c={c:.6f}, delta_beta(15)={delta:.3e}, min p over test words={worst:.4e}')
        check(worst >= delta, f'beta={beta}: every computed cut probability >= delta_beta of (15)')
    # exhaustive: all 2^12 words of length 12, odds bound (14) word by word (float, relative slack 1e-12)
    for beta in [0.69, 0.75]:
        A = [1.0]
        for m in range(1, 13):
            A.append(sum(math.exp(-beta * l * l / 2) * A[m - l] for l in range(1, m + 1)))
        ok = True
        for word in itertools.product((-1, 1), repeat=12):
            p, _ = cut_probs(word, beta)
            for i in range(1, 12):
                S1 = sum(1 / A[a] for a in range(1, i + 1))
                S2 = sum(1 / A[b] for b in range(1, 12 - i + 1))
                odds = (1 - p[i - 1]) / p[i - 1]
                if odds > S1 * S2 * (1 + 1e-12):
                    ok = False
        check(ok, f'beta={beta}: no-cut/cut odds <= (sum 1/A_a)(sum 1/A_b) for all 4096 words of length 12, all cuts')
    return


# ---------------------------------------------------------------- 5. diblock asymptotics
def diblock_formula(beta, n, u, r):
    A = u[: n + 1] * r ** (-np.arange(n + 1, dtype=float)) if r >= 1 else None
    return A


def part5(bstar):
    # (a) full-chain DP vs closed formula at moderate n
    for beta in [0.5, bstar, 0.9]:
        n = 300
        r, u = scaled_A(beta, n)
        # log A_m = log u_m - m log r
        logA = np.log(u) - np.arange(n + 1) * math.log(r)
        i = np.arange(n)
        M = logA[i][:, None] + logA[i][None, :] - beta * (i[:, None] - i[None, :]) ** 2 / 2
        terms = np.concatenate([[2 * logA[n]], M.ravel()])
        mx = terms.max()
        logD = mx + math.log(np.exp(terms - mx).sum())
        p_formula = math.exp(2 * logA[n] - logD)
        p, logZred = cut_probs([1] * n + [-1] * n, beta)
        check(abs(p[n - 1] / p_formula - 1) < 1e-9 and abs(logZred - logD) < 1e-9,
              f'beta={beta:.6f}, 2n={2 * n}: full-chain DP == formula (16)/(17) (p_mid={p_formula:.6e})')
    # (b) critical asymptotics 3 n p -> 1 and D_n / n -> 3 a^2 (formula, float, scaled)
    r, u = scaled_A(bstar, 200000)
    a = u[-1]
    Lk = 40
    ker = np.array([math.exp(-bstar * d * d / 2) for d in range(Lk + 1)])
    for n in [1000, 10000, 100000]:
        A = u[: n + 1]  # r = 1 at criticality (to float accuracy)
        s = np.dot(A[:n], A[:n])
        for d in range(1, Lk + 1):
            s += 2 * ker[d] * np.dot(A[: n - d], A[d:n])
        D = A[n] ** 2 + s
        print(f'INFO beta_*: n={n}, 3n p_mid = {3 * n * A[n] ** 2 / D:.8f}, D_n/(3 a^2 n) = {D / (3 * a * a * n):.8f}')
    check(abs(3 * 100000 * A[100000] ** 2 / D - 1) < 1e-3, 'critical: 3 n p_{n,2n} -> 1 (n=1e5 within 1e-3)')
    # (c) below: p_mid -> 1/K_beta; above: exponential decay; free energy
    for beta in [0.6, 0.69]:
        r = r_of_beta(beta)
        K = 1 + sum(r ** (x + y) * math.exp(-beta * (x - y) ** 2 / 2) for x in range(1, 2000) for y in range(max(1, x - 60), min(2000, x + 61)))
        p, logZred = cut_probs([1] * 1500 + [-1] * 1500, beta)
        F = beta / 2 + logZred / 3000
        Flim = beta / 2 + max(0.0, -math.log(r))
        print(f'INFO beta={beta}: p_mid(2n=3000)={p[1499]:.8f}, 1/K={1 / K:.8f}; F_3000={F:.6f}, F_lim={Flim:.6f}')
        check(abs(p[1499] * K - 1) < 1e-6, f'beta={beta}: p_mid -> 1/K_beta (23)')
        check(abs(F - Flim) < 2e-3, f'beta={beta}: finite-volume free energy close to (4)')
    for beta in [0.75, 1.0]:
        p, logZred = cut_probs([1] * 1500 + [-1] * 1500, beta)
        F = beta / 2 + logZred / 3000
        r = r_of_beta(beta)
        print(f'INFO beta={beta}: r={r:.5f}, p_mid={p[1499]:.3e}, p_10={p[9]:.3e}, p_100={p[99]:.3e}; F_3000={F:.6f} vs beta/2={beta / 2}')
        check(p[1499] < 1e-30 and p[99] < p[9] * 1e-3, f'beta={beta}: cut probabilities decay exponentially away from the ends (3)')
        check(abs(F - beta / 2) < 2e-3, f'beta={beta}: free energy close to beta/2 (4)')
    # (d) free energy has a kink at beta_*: left derivative of -log r_beta is non-zero
    h = 1e-5
    dleft = (-math.log(r_of_beta(bstar - h)) + math.log(r_of_beta(bstar - 2 * h))) / h
    print(f'INFO d/dbeta(-log r_beta) at beta_*^- ~ {dleft:.6f} (non-zero => first-order kink of F_db)')
    check(dleft < -0.1, 'F_db has a kink (first-order non-analyticity) at beta_*')


def main():
    t = time.time()
    part1_2()
    bstar = part3()
    part4(bstar)
    part5(bstar)
    print(f'INFO elapsed {time.time() - t:.1f}s')
    if FAIL:
        print('FAIL some checks failed')
        sys.exit(1)
    print('PASS all verifier checks')


if __name__ == '__main__':
    main()
