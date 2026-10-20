"""Exploration check of Lemma G1 (recovery ordering in one TG column), Agent G.

For a single column with feed (U, f) on components i..m and split [i,k]/[l,m] (l <= k: distributing components
l..k), the TG constraints are LINEAR in (d_l..d_k, Y_T) once the feed roots theta_r (r = l-1..k) are fixed:
   Y_T = g_T(theta_r) (l <= r <= k-1),  Y_T >= g_T(theta_{l-1}),  Y_T >= g_T(theta_k),
   Y_T >= 0,  Y_B = Y_T - U >= 0,  0 <= d_p <= f_p,
with g_T(theta) = sum_{p<=k} a_p d_p/(a_p - theta).  Lemma G1 claims d_p/f_p >= d_{p+1}/f_{p+1} for l <= p < k.
We minimise d_p f_{p+1} - d_{p+1} f_p (scaled) by LP (HiGHS) on random instances.
Floating point exploration (not a certificate).  Environment: Python 3, numpy, scipy.  Runtime: ~1 min.
Usage: python recovery_order_lp.py ntrials seed
"""
import sys
import random
import numpy as np
from scipy.optimize import linprog, brentq


def roots(al, f, U):
    out = {}
    n = len(al)
    for r in range(n - 1):
        lo, hi = al[r + 1], al[r]
        w = hi - lo
        g = lambda t: sum(a * x / (a - t) for a, x in zip(al, f)) - U
        out[r] = brentq(g, lo + 1e-13 * w, hi - 1e-13 * w, xtol=1e-15, maxiter=500)
    return out


def main():
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    worst = 1e9
    nlp = 0
    for tr in range(ntr):
        n = rng.randint(3, 7)
        rat = [rng.uniform(1.02, 4.0) for _ in range(n - 1)]
        al = [1.0]
        for r in reversed(rat):
            al.insert(0, al[0] * r)
        f = [rng.uniform(0.01, 1.0) for _ in range(n)]
        U = rng.uniform(-1.5, 1.5) * sum(f)
        th = roots(al, f, U)
        # split: distillate 0..k, residue l..n-1 (0-based), l <= k
        k = rng.randint(1, n - 2)
        l = rng.randint(1, k)
        dist = list(range(l, k + 1))
        nv = len(dist) + 1  # d_l..d_k, Y
        def gT_row(t):
            row = np.zeros(nv)
            const = sum(al[p] * f[p] / (al[p] - t) for p in range(0, l))
            for idx, p in enumerate(dist):
                row[idx] = al[p] / (al[p] - t)
            return row, const
        A_eq, b_eq, A_ub, b_ub = [], [], [], []
        for r in range(l, k):          # interior roots r (0-based interval r between comps r, r+1)
            row, c = gT_row(th[r])
            row[-1] = -1.0
            A_eq.append(row); b_eq.append(-c)
        for r in (l - 1, k):
            row, c = gT_row(th[r])
            row[-1] = -1.0
            A_ub.append(row); b_ub.append(-c)
        bounds = [(0, f[p]) for p in dist] + [(max(0.0, U), None)]
        for idx in range(len(dist) - 1):
            p, pp = dist[idx], dist[idx + 1]
            obj = np.zeros(nv)
            obj[idx] = f[pp]; obj[idx + 1] = -f[p]
            obj /= (f[p] * f[pp])     # = r_p - r_{p+1}
            res = linprog(obj, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                          A_eq=np.array(A_eq) if A_eq else None, b_eq=np.array(b_eq) if b_eq else None,
                          bounds=bounds, method='highs')
            if res.status == 0:
                nlp += 1
                if res.fun < worst:
                    worst = res.fun
                if res.fun < -1e-7:
                    print('VIOLATION', n, al, f, U, k, l, p, res.fun, res.x)
    print('LPs solved', nlp, 'min (r_p - r_{p+1}) =', worst)
    print('PASS' if worst > -1e-7 else 'FAIL')


if __name__ == '__main__':
    main()
