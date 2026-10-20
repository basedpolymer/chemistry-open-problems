#!/usr/bin/env python3
"""check_limit_jacobian.py -- numerical checks of Theorem L (Agent L; ../limit.md, Sections 3-5).

For random data (N = 3..6, random volatilities, feed, q in [-1, 2], random pinning pattern of the binary
columns, several slack/forcing sizes) the script builds a limit configuration (limitsys.design, Agent L's own test
recipe), the limit system G_inf (limitsys.System with n = None) and its Jacobian at v_inf by central differences in
60-digit arithmetic, and checks:
  C1  G_inf(v_inf) = 0, all flows positive, root locations (L-c), (L-d) of limit.md (incl. middle roots common);
  C2  block structure: d(trace rows)/d(major variables) = 0, and d(trace rows of k)/d(trace variables of k') = 0
      for k != k';
  C3  every trace row agrees with the predicted linearization (x_{c,k} + w_k G^(k)_{x_c}(alpha_k)/Psi^(k)(alpha_k))
      divided by the pivot (Lemma 3.4 of limit.md);
  C4  every trace mode is closed at a node where its component is a trace, kappa_sigma > 0, and for every component
      the matrix A_k of Proposition 4.3 is a Z-matrix with nonnegative column sums, nonsingular, with A_k^{-1} >= 0;
  C5  the major block is block lower triangular over the columns (ordered parents first) with nonsingular diagonal
      blocks;
  C6  the full Jacobian is nonsingular (smallest singular value reported).
Prints PASS/FAIL per instance and per check; exits non-zero on any FAIL.  Evidence for the proof in limit.md, not
part of it.

Usage:  python check_limit_jacobian.py [n_instances_per_N] [Nmax] [seed]
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env), 1 thread.
Runtime: about 1 min for N = 3..5 with 4 instances each; N = 6 adds a few minutes.
"""
import sys
import random
import time
import numpy as np
import mpmath as mp
import limitsys as ls

mp.mp.dps = 60
TOL0 = mp.mpf('1e-45')      # residual of G_inf at v_inf
TOLFD = mp.mpf('1e-28')     # finite-difference noise level (h = 1e-20, 60 digits)


def node_of(s, end):
    kind, i, m = s
    if kind == 'T':
        if end == 'bot':
            return ('J', i, m)
        return ('J', i, m - 1) if m - 1 > i else ('P', i)
    if end == 'top':
        return ('J', i, m)
    return ('J', i + 1, m) if i + 1 < m else ('P', m)


def other(end):
    return 'top' if end == 'bot' else 'bot'


def run_instance(N, rng, tag, verbose=True):
    al = [mp.mpf(1)]
    for _ in range(N - 1):
        al.append(al[-1] * mp.mpf(1 + rng.uniform(0.25, 1.2)))
    al = sorted(al, reverse=True)
    zz = [rng.uniform(0.2, 1.0) for _ in range(N)]
    s = sum(zz)
    z = [mp.mpf(x) / s for x in zz]
    q = mp.mpf(rng.uniform(-1.0, 2.0))
    pattern = {p: rng.choice('BT') for p in range(N - 1)}
    eta_col = mp.mpf(10) ** (-rng.uniform(2.5, 4))
    eta_bin = mp.mpf(10) ** (-rng.uniform(1.8, 3))
    delta = mp.mpf(rng.uniform(0.01, 0.08))
    D = ls.design(al, z, q, eta_col=eta_col, eta_bin=eta_bin, delta=delta, pattern=pattern)
    head = (f'[{tag}] N={N} q={mp.nstr(q, 3)} pattern={"".join(pattern[p] for p in range(N - 1))} '
            f'eta_col={mp.nstr(eta_col, 2)} eta_bin={mp.nstr(eta_bin, 2)} delta={mp.nstr(delta, 2)} '
            f'R-maxH={mp.nstr(D.R - max(D.H), 3)}')
    return theorem_L_checks(D, head, verbose)


def theorem_L_checks(D, head, verbose=True):
    """Checks C1-C6 on a limit configuration D (limitsys.Design-like object).  Returns True iff all pass."""
    N = D.N
    pattern = D.pattern
    S = ls.System(D, None)
    ok = {}
    msgs = []
    # ---------------- C1
    c1 = True
    for c in D.cols:
        for x in (D.VT[c], D.VB[c], D.LT[c], D.LB[c]):
            if not x > 0:
                c1 = False
                msgs.append(f'nonpositive flow at {c}')
    if not all(P > 0 for P in D.P):
        c1 = False
    g = S.G(S.v_inf)
    res = max(abs(x) for x in g)
    if res > TOL0:
        c1 = False
        msgs.append(f'G_inf(v_inf) = {mp.nstr(res, 3)}')
    margins = []
    for c in D.cols:
        i, m = c
        T, B = S.sec_info[('T',) + c], S.sec_info[('B',) + c]
        rT = [t[1] for t in T['modes'] if t[0] == 'maj']
        rB = [t[1] for t in B['modes'] if t[0] == 'maj']
        if m - i >= 2:
            r, u, v, s_ = rT[0], rT[-1], rB[0], rB[-1]
            # (L-c): u > al_m, v < al_i, middle roots common, r != v, u != s (we also report r - v, u - s)
            if not (u > D.al[m] and v < D.al[i]):
                c1 = False
                msgs.append(f'(L-c) root location fails at {c}')
            for a_, b_ in zip(rT[1:-1], rB[1:-1]):
                if abs(a_ - b_) > mp.mpf('1e-40'):
                    c1 = False
                    msgs.append(f'middle roots not common at {c}')
            margins.append(min(r - v, u - s_))
        else:
            p = i
            if D.pattern[p] == 'B':
                phT = rT[0]
                php, ph2 = rB[0], rB[1]
                cond = phT > D.al[p + 1] and (p == 0 or php < D.al[p - 1]) and php > D.al[p]
                margins.append(phT - ph2)
            else:
                phB = rB[0]
                ph2, php = rT[0], rT[1]
                cond = phB < D.al[p] and (p + 2 > N - 1 or php > D.al[p + 2]) and php < D.al[p + 1]
                margins.append(ph2 - phB)
            if not cond:
                c1 = False
                msgs.append(f'(L-d) root location fails at {c}')
    if not all(mg > 0 for mg in margins):
        c1 = False
        msgs.append('a column margin is not positive')
    ok['C1'] = c1
    # ---------------- Jacobian
    J = S.jac(S.v_inf)
    rc = S.row_class()
    vc = S.var_class()
    nv = S.nv
    maj_r = [a for a in range(nv) if rc[a][0] == 'maj']
    maj_v = [b for b in range(nv) if vc[b][0] == 'maj']
    # ---------------- C2
    off = mp.mpf(0)
    for a in range(nv):
        if rc[a][0] != 'tr':
            continue
        for b_ in range(nv):
            if vc[b_][0] == 'maj' or (vc[b_][0] == 'tr' and vc[b_][1] != rc[a][1]):
                off = max(off, abs(J[a, b_]))
    ok['C2'] = off < TOLFD
    # ---------------- C3: predicted trace rows
    v0 = S.v_inf
    h = mp.mpf('1e-20')

    def end_comp(vv, s, end):
        d, b, f, X, prod = S.unpack(vv)
        xb, xt = S.ends(s, X, prod)
        return (xb if end == 'bot' else xt), d, b
    d0, b0, f0, X0, prod0 = S.unpack(v0)
    worst = mp.mpf(0)
    kappas = {}
    for a, row in enumerate(S.rows):
        s, j, typ, k, end = row
        if typ != 'tr':
            continue
        kind, i, m = s
        c = (i, m)
        xc0, _, _ = end_comp(v0, s, end)
        w0 = d0[c] if kind == 'T' else [-x for x in b0[c]]
        V = S.lim['V'][s]
        Gk = sum(D.al[ii] * xc0[ii] / (D.al[ii] - D.al[k]) for ii in range(N) if ii != k)
        Pk = V - sum(D.al[ii] * w0[ii] / (D.al[ii] - D.al[k]) for ii in range(N) if ii != k)
        info = S.sec_info[s]
        phip = info['modes'][info['ipin']][1]
        lp = [D.al[ii] / (D.al[ii] - phip) for ii in range(N)]
        piv = sum(x * y for x, y in zip(lp, xc0))
        kap = (Gk / Pk) if end == 'top' else (-Gk / Pk)
        kappas[(s, k)] = (kap, end, Gk, Pk)
        for b_ in range(nv):
            vp = list(v0)
            vm = list(v0)
            vp[b_] += h
            vm[b_] -= h
            xp, dp, bp = end_comp(vp, s, end)
            xm, dm, bm = end_comp(vm, s, end)
            dx = (xp[k] - xm[k]) / (2 * h)
            wp = dp[c][k] if kind == 'T' else -bp[c][k]
            wm = dm[c][k] if kind == 'T' else -bm[c][k]
            dw = (wp - wm) / (2 * h)
            pred = (dx + dw * Gk / Pk) / piv
            scale = 1 + max(abs(J[a, bb]) for bb in range(nv))
            worst = max(worst, abs(pred - J[a, b_]) / scale)
    ok['C3'] = worst < TOLFD
    # ---------------- C4: kappa signs, closed ends at trace nodes, M-matrix per component
    def content(node, k):
        if node[0] == 'J':
            return S.X_inf[(node[1], node[2])][k]
        return prod0[node[1]][k]
    c4 = True
    for (s, k), (kap, end, Gk, Pk) in kappas.items():
        if not kap > 0:
            c4 = False
            msgs.append(f'kappa <= 0 at {s}, comp {k}')
        if content(node_of(s, end), k) != 0:
            c4 = False
            msgs.append(f'trace mode {k} of {s} closed at a node where {k} is not a trace')
    minfo = []
    for k in range(N):
        secs = [(s, kappas[(s, k)]) for (s, kk) in kappas if kk == k]
        nodes = set()
        for s, (kap, end, Gk, Pk) in secs:
            for e in ('bot', 'top'):
                nd = node_of(s, e)
                if content(nd, k) == 0:
                    nodes.add(nd)
        nodes = sorted(nodes)
        if not nodes:
            continue
        idx = {nd: t for t, nd in enumerate(nodes)}
        A = mp.matrix(len(nodes), len(nodes))
        for nd in nodes:
            if nd[0] == 'P':
                A[idx[nd], idx[nd]] += D.P[nd[1]]
        for s, (kap, end, Gk, Pk) in secs:
            cn, on = node_of(s, end), node_of(s, other(end))
            A[idx[cn], idx[cn]] += 1 / kap
            if on in idx:
                A[idx[on], idx[cn]] -= 1 / kap
        nn = len(nodes)
        zmat = all(A[a, b_] <= 0 for a in range(nn) for b_ in range(nn) if a != b_)
        colsum = [sum(A[a, b_] for a in range(nn)) for b_ in range(nn)]
        try:
            Ai = mp.inverse(A)
            invpos = all(Ai[a, b_] >= -mp.mpf('1e-40') for a in range(nn) for b_ in range(nn))
            nonsing = True
        except ZeroDivisionError:
            invpos, nonsing = False, False
        good = zmat and all(cs >= -mp.mpf('1e-40') for cs in colsum) and nonsing and invpos
        minfo.append((k, nn, good, mp.nstr(min(colsum), 3)))
        if not good:
            c4 = False
            msgs.append(f'M-matrix property fails for component {k}')
    ok['C4'] = c4
    # ---------------- C5: major block lower triangular over columns, diagonal blocks nonsingular
    rowcol = {}
    for a, row in enumerate(S.rows):
        s = row[0]
        rowcol[a] = (s[1], s[2])
    varcol = {b_: S.vars[b_][0] for b_ in range(nv)}

    def ancestor(c1_, c2_):  # c1 strictly contains c2
        return c1_ != c2_ and c1_[0] <= c2_[0] and c1_[1] >= c2_[1]
    viol = mp.mpf(0)
    for a in maj_r:
        for b_ in maj_v:
            cr, cv = rowcol[a], varcol[b_]
            if cr != cv and not ancestor(cv, cr):
                viol = max(viol, abs(J[a, b_]))
    diag_ok = True
    for c in D.cols:
        ra = [a for a in maj_r if rowcol[a] == c]
        vb = [b_ for b_ in maj_v if varcol[b_] == c]
        if len(ra) != len(vb):
            diag_ok = False
            continue
        if not ra:
            continue
        Mb = np.array([[float(J[a, b_]) for b_ in vb] for a in ra])
        sv = np.linalg.svd(Mb, compute_uv=False)
        if sv[-1] < 1e-10 * max(1.0, sv[0]):
            diag_ok = False
    ok['C5'] = viol < TOLFD and diag_ok and len(maj_r) == len(maj_v)
    # ---------------- C6
    Jf = np.array([[float(J[a, b_]) for b_ in range(nv)] for a in range(nv)])
    sv = np.linalg.svd(Jf, compute_uv=False)
    ok['C6'] = sv[-1] > 1e-12 * sv[0]
    allok = all(ok.values())
    if verbose:
        print(f'{head} nv={nv} (maj {len(maj_v)})')
        print(f'      |G_inf(v_inf)|={mp.nstr(res, 2)} min margin={mp.nstr(min(margins), 3)} '
              f'offblock={mp.nstr(off, 2)} tracepred={mp.nstr(worst, 2)} majtri={mp.nstr(viol, 2)} '
              f'sigma_min={sv[-1]:.3e} sigma_max={sv[0]:.3e}')
        print('      M-matrices (k, size, ok, min column sum): ' + ' '.join(str(t) for t in minfo))
        print('      ' + ' '.join(f'{key}:{"PASS" if val else "FAIL"}' for key, val in ok.items())
              + ('' if not msgs else '  ' + '; '.join(msgs)))
    return allok


def main():
    nper = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    Nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = random.Random(seed)
    t0 = time.time()
    allok = True
    count = 0
    for N in range(3, Nmax + 1):
        for t in range(nper):
            okk = run_instance(N, rng, f'N{N}#{t}')
            allok &= okk
            count += 1
            sys.stdout.flush()
    print(f'{count} instances, {time.time() - t0:.1f} s')
    print('ALL PASS' if allok else 'SOME FAIL')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
