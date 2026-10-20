#!/usr/bin/env python3
"""check_unforced.py -- the role of the forcing (Agent L; ../limit.md Section 7.1).

Claim checked (proved in limit.md, Remark 7.1): if a binary column is NOT forced (eta_p = 0, sharp balanced split),
the limit system G_inf has a one-parameter family of zeros (its junction ratio is free) and DG_inf(v_inf) is singular,
with a null vector along the junction-ratio coordinate X_{c,p+1}; with forcing eta_p > 0 it is invertible
(Theorem L).  The script also reports how sigma_min(DG_inf) behaves as eta_p -> 0 (continuity: it must tend to 0,
which is why the order of limits n -> infinity first, then eta -> 0, matters).

For N = 3 and N = 4 (pattern all B), R is fixed (taken from the forced design at the largest eta), the Astra slack of
the non-binary columns is fixed, and the forcing eta of every binary column is scanned down to 0.
PASS iff sigma_min > 0 for every eta > 0 and sigma_min / sigma_max < 1e-30 at eta = 0 with the null vector having
its largest entries on junction-ratio coordinates of binary columns.  Evidence only.

Usage: python check_unforced.py
Environment: Python 3.12, mpmath 1.3, numpy (WSL Fedora sci env), 1 thread.  Runtime: about 1 minute.
"""
import sys
import numpy as np
import mpmath as mp
import limitsys as ls

mp.mp.dps = 60


def case(alpha, z, q, etas):
    D0 = ls.design(alpha, z, q, eta_col=mp.mpf('1e-3'), eta_bin=mp.mpf(etas[0]), delta=mp.mpf('0.05'))
    R = D0.R
    ok = True
    print(f'alpha={alpha} q={q} R={mp.nstr(R, 8)} (fixed), pattern all B')
    for e in etas:
        D = ls.design(alpha, z, q, eta_col=mp.mpf('1e-3'), eta_bin=mp.mpf(e), R=R)
        S = ls.System(D, None)
        g = max(abs(x) for x in S.G(S.v_inf))
        J = S.jac(S.v_inf)
        A = mp.matrix(J)
        # high-precision smallest singular value via the eigenvalues of J^T J
        Jf = np.array([[float(J[i, j]) for j in range(J.cols)] for i in range(J.rows)])
        U, sv, Vt = np.linalg.svd(Jf)
        smin = mp.sqrt(min(abs(x) for x in mp.eigsy(A.T * A, eigvals_only=True)))
        nv = Vt[-1]
        top = np.argsort(-np.abs(nv))[:2]
        names = [S.vars[t][:3] for t in top]
        rel = smin / mp.mpf(sv[0])
        # same Jacobian with the REGULARIZED normalization of the forced pinch (l~ = ((alpha_p - phi')/alpha_p) l):
        # every row of a forced section is divided by its pinch-vector scale s = (alpha_p - phi')/alpha_p
        sreg = None
        if mp.mpf(e) > 0:
            Jr = mp.matrix(J)
            for a_, row in enumerate(S.rows):
                s_ = row[0]
                kind, i, m = s_
                if m - i == 1 and kind == 'B':        # forced bottom section (pattern B)
                    info = S.sec_info[s_]
                    phi = info['modes'][info['ipin']][1]
                    scale = (D.al[i] - phi) / D.al[i]
                    for b_ in range(J.cols):
                        Jr[a_, b_] = Jr[a_, b_] / scale
            sreg = mp.sqrt(min(abs(x) for x in mp.eigsy(Jr.T * Jr, eigvals_only=True)))
        tag = ''
        if mp.mpf(e) == 0:
            good = rel < mp.mpf('1e-30') and all(nm[1] == 'X' and nm[0][1] - nm[0][0] == 1 for nm in names[:1])
            tag = 'singular (expected)' if good else 'NOT singular?'
        else:
            good = rel > mp.mpf('1e-14')
            tag = 'invertible' if good else 'near-singular?'
        ok &= good
        print(f'   eta_bin={mp.nstr(mp.mpf(e), 3):>8}  |G_inf(v_inf)|={mp.nstr(g, 2):>8}  sigma_min={mp.nstr(smin, 4):>11}  '
              f'sigma_min/sigma_max={mp.nstr(rel, 3):>9}  sigma_min(regularized forced pinch)='
              f'{mp.nstr(sreg, 4) if sreg is not None else "-":>10}  null-vector top coords {names}  {tag}')
    return ok


def main():
    etas = ['1e-1', '3e-2', '1e-2', '3e-3', '1e-3', '1e-4', '1e-5', '1e-6', '0']
    ok = case([4, 2, 1], [mp.mpf(1) / 3] * 3, 1, etas)
    ok &= case([8, 4, 2, 1], [mp.mpf(1) / 4] * 4, mp.mpf('0.5'), etas)
    print('ALL PASS' if ok else 'SOME FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
