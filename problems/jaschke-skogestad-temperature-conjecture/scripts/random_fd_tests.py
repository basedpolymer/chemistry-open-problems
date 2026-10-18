#!/usr/bin/env python3
"""
random_fd_tests.py -- independent NUMERICAL check (high-precision finite differences) of
  (A) the corrected AMTD marginal-cost formula (Theorem A of notes.md) for M = 1..8, random data, random prices;
  (B) the failure of the paper's conjectured recursion for every M >= 3;
  (C) the general effectiveness form (Proposition C of notes.md) for the LOG-mean (LMTD) counter-current
      exchanger, where the marginal cost needs the hot OUTLET temperatures as well; and the size of the error
      made by the (corrected) AMTD formula when the plant obeys LMTD.
This is numerical evidence that is independent of the symbolic code; it is not the proof.

Environment : Python 3 + mpmath (WSL Fedora sci env).  Runtime: ~10-30 s.
Usage       : python -I random_fd_tests.py [seed]       (default seed 20261018)
Exit code 1 if any claim fails.
"""
import random
import sys
import mpmath as mp

mp.mp.dps = 60
H = mp.mpf("1e-22")          # central-difference step; truncation error O(H^2) ~ 1e-44
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261018
rng = random.Random(seed)
FAILED = []


def report(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        FAILED.append(name)


# ---------------- AMTD model ----------------
def amtd_branch(w, UA, wh, thh):
    th, tho, Q = [mp.mpf(0)], [], []
    for U, W, h in zip(UA, wh, thh):
        tp = th[-1]
        k = 2 / U + (0 if W is None else 1 / W)
        c = 2 / (1 + w * k)                     # solves the three model equations (see notes.md, Lemma 1)
        t = tp + c * (h - tp)
        q = w * (t - tp)
        th.append(t)
        Q.append(q)
        tho.append(h if W is None else h - q / W)
    return th, tho, Q


def amtd_residual(w, UA, wh, thh):
    """max |residual| of the raw model equations (guards against an error in the closed-form solve)."""
    th, tho, Q = amtd_branch(w, UA, wh, thh)
    r = mp.mpf(0)
    for i in range(len(UA)):
        r = max(r, abs(Q[i] - w * (th[i + 1] - th[i])))
        r = max(r, abs(Q[i] - UA[i] * ((thh[i] - th[i + 1]) + (tho[i] - th[i])) / 2))
        if wh[i] is not None:
            r = max(r, abs(Q[i] - wh[i] * (thh[i] - tho[i])))
    return r


def jaschke_a(th, thh, mode):
    a = []
    for i in range(1, len(th)):
        c = (th[i] - th[i - 1]) / (thh[i - 1] - th[i - 1])
        sub = mp.fsum(a) if mode == "corrected" else (a[-1] if a else 0)
        a.append(c * (th[i] + th[i - 1] - sub))
    return a


def random_branch(M, allow_inf=True):
    UA = [mp.mpf(rng.uniform(5, 200)) for _ in range(M)]
    wh = [None if (allow_inf and rng.random() < 0.2) else mp.mpf(rng.uniform(20, 300)) for _ in range(M)]
    thh = sorted(mp.mpf(rng.uniform(20, 150)) for _ in range(M))
    w = mp.mpf(rng.uniform(10, 100))
    p = [mp.mpf(rng.uniform(0.5, 2)) for _ in range(M)]
    return w, UA, wh, thh, p


def part_A_B(trials=25):
    for M in range(1, 9):
        worst_cor, best_pap, worst_res = mp.mpf(0), mp.inf, mp.mpf(0)
        for _ in range(trials):
            w, UA, wh, thh, p = random_branch(M)
            worst_res = max(worst_res, amtd_residual(w, UA, wh, thh))
            th, _, _ = amtd_branch(w, UA, wh, thh)
            _, _, Qp = amtd_branch(w + H, UA, wh, thh)
            _, _, Qm = amtd_branch(w - H, UA, wh, thh)
            dQ = [(a - b) / (2 * H) for a, b in zip(Qp, Qm)]
            ac = jaschke_a(th, thh, "corrected")
            ap = jaschke_a(th, thh, "paper")
            scale = max(abs(x) for x in dQ)
            # economic marginal value with random prices (Theorem 3 form)
            econ_fd = mp.fsum(pi * d for pi, d in zip(p, dQ))
            econ_c = mp.fsum(pi * x for pi, x in zip(p, ac)) / 2
            econ_p = mp.fsum(pi * x for pi, x in zip(p, ap)) / 2
            err_c = max(max(abs(d - x / 2) for d, x in zip(dQ, ac)), abs(econ_fd - econ_c)) / scale
            worst_cor = max(worst_cor, err_c)
            if M >= 3:
                err_p = abs(mp.fsum(dQ) - mp.fsum(ap) / 2) / abs(mp.fsum(dQ))
                best_pap = min(best_pap, err_p)
        report(f"A  M={M}: corrected formula vs finite differences, worst rel. error "
               f"{mp.nstr(worst_cor, 3)} (< 1e-30); model residual {mp.nstr(worst_res, 3)}",
               worst_cor < mp.mpf("1e-30") and worst_res < mp.mpf("1e-45"))
        if M >= 3:
            report(f"B  M={M}: paper recursion wrong in all {trials} draws; smallest rel. error of 2*d(sum Q)/dw vs "
                   f"paper T_J = {mp.nstr(best_pap, 3)} (> 1e-6)", best_pap > mp.mpf("1e-6"))


# ---------------- LMTD (counter-current) model ----------------
def eps_cold(N, R):
    """cold-side effectiveness of a counter-current exchanger, N = UA/w, R = w/wh (R != 1)."""
    if R == 0:
        return 1 - mp.exp(-N)
    E = mp.exp(-N * (1 - R))
    return (1 - E) / (1 - R * E)


def lmtd_branch(w, UA, wh, thh):
    th, tho, Q = [mp.mpf(0)], [], []
    for U, W, h in zip(UA, wh, thh):
        tp = th[-1]
        R = 0 if W is None else w / W
        e = eps_cold(U / w, R)
        t = tp + e * (h - tp)
        q = w * (t - tp)
        th.append(t)
        Q.append(q)
        tho.append(h if W is None else h - q / W)
    return th, tho, Q


def lmtd_residual(w, UA, wh, thh):
    th, tho, Q = lmtd_branch(w, UA, wh, thh)
    r = mp.mpf(0)
    for i in range(len(UA)):
        d1, d2 = thh[i] - th[i + 1], tho[i] - th[i]
        lm = (d1 - d2) / mp.log(d1 / d2) if abs(d1 - d2) > mp.mpf("1e-40") else d1
        r = max(r, abs(Q[i] - UA[i] * lm) / abs(Q[i]))
    return r


def marginal_from_temperatures(th, tho, thh):
    """Proposition C: dQ_i/dw from the four temperatures of each exchanger (no UA, no heat capacities)."""
    out = []
    for i in range(1, len(th)):
        dlt = th[i] - th[i - 1]
        e = dlt / (thh[i - 1] - th[i - 1])
        R = (thh[i - 1] - tho[i - 1]) / dlt                     # = w/wh_i (0 for constant hot temperature)
        N = mp.log((1 - R * e) / (1 - e)) / (1 - R) if R != 0 else -mp.log(1 - e)
        phi = -N * mp.diff(lambda n: eps_cold(n, R), N) + (R * mp.diff(lambda r: eps_cold(N, r), R) if R != 0 else 0)
        out.append((1 + phi / e) * dlt + e * (th[i - 1] - mp.fsum(out)))
    return out


def part_C(trials=10):
    worst, worst_res, amtd_err = mp.mpf(0), mp.mpf(0), []
    for M in (1, 2, 3, 4, 5):
        for _ in range(trials):
            w, UA, wh, thh, p = random_branch(M)
            # keep R away from 1 (removable singularity of the closed form)
            wh = [W if (W is None or abs(w / W - 1) > 0.05) else W * 1.3 for W in wh]
            worst_res = max(worst_res, lmtd_residual(w, UA, wh, thh))
            th, tho, _ = lmtd_branch(w, UA, wh, thh)
            _, _, Qp = lmtd_branch(w + H, UA, wh, thh)
            _, _, Qm = lmtd_branch(w - H, UA, wh, thh)
            dQ = [(a - b) / (2 * H) for a, b in zip(Qp, Qm)]
            mt = marginal_from_temperatures(th, tho, thh)
            scale = max(abs(x) for x in dQ)
            worst = max(worst, max(abs(a - b) for a, b in zip(dQ, mt)) / scale)
            ac = jaschke_a(th, thh, "corrected")
            amtd_err.append(abs(mp.fsum(dQ) - mp.fsum(ac) / 2) / abs(mp.fsum(dQ)))
    report(f"C  LMTD: marginal cost from the 4 temperatures of each exchanger vs finite differences, worst rel. "
           f"error {mp.nstr(worst, 3)} (< 1e-25); LMTD model residual {mp.nstr(worst_res, 3)}",
           worst < mp.mpf("1e-25") and worst_res < mp.mpf("1e-40"))
    amtd_err.sort()
    print(f"      (info) corrected AMTD formula applied to an LMTD plant: rel. error of the branch marginal cost "
          f"median {mp.nstr(amtd_err[len(amtd_err)//2], 3)}, max {mp.nstr(amtd_err[-1], 3)} -- an approximation only")


if __name__ == "__main__":
    print(f"seed = {seed}, mp.dps = {mp.mp.dps}")
    part_A_B()
    part_C()
    if FAILED:
        print(f"FAIL: {len(FAILED)} claim(s) failed")
        sys.exit(1)
    print("ALL PASS")
