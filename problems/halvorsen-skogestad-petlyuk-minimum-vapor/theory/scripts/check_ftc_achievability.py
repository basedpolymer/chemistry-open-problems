"""Check of the attainment theorem (paper, Section "Attainment by the FTC"): the generalized Petlyuk / FTC configuration has an operating point that
satisfies every constraint of the STRICT model (hence of the permissive TG model) with total reboiler vapor
R = max_r H_r.

Construction (paper): every column [i,j] with j - i >= 2 at its preferred split (all carried-over
feed roots theta_i..theta_{j-1} active, V = Underwood value, no excess); the binary columns [r,r+1] carry the
excess e_r = max_s H_s - H_r. The script builds the point column by column (longest streams first) in
mpmath (80 digits) and checks, with tolerance 1e-50:
  - feed equations of every column at its roots (carry-over), all flows >= 0 (strictly > 0 where nonzero),
  - strict Underwood conditions (g_T(theta_r) = V_T for all roots of non-binary columns,
    g_T(theta_r) <= V_T and g_B(theta_r) <= V_B for binary ones), liquids L_T, L_B >= 0,
  - liquid side-draw balances V_T[p,p+1] = V_B[p-1,p] and product purities/amounts,
  - R = V_B[N-1,N] = max_r H_r.
Random instances: N = 3..8, 200 instances per N, alpha ratios in [1.01, 20], z in [1e-3, 1] normalized,
q in [-1, 2] (outside [0,1] included on purpose).  Prints PASS/FAIL per N; exit code 1 on any failure.
Environment: WSL Fedora, Python 3 + mpmath (outils/wsl).  Runtime: about 2 minutes.
"""
import random
import sys

import mpmath as mp

mp.mp.dps = 80
TOL = mp.mpf(10) ** -50


def feed_roots(al, z, q):
    rhs = 1 - q
    roots = []
    for r in range(len(al) - 1):
        lo, hi = al[r + 1], al[r]
        f = lambda t: sum(a * x / (a - t) for a, x in zip(al, z)) - rhs
        a, b = lo, hi
        for _ in range(400):
            m = (a + b) / 2
            if f(m) > 0:
                b = m
            else:
                a = m
        roots.append((a + b) / 2)
    return roots


def build_and_check(al, z, q):
    N = len(al)
    th = feed_roots(al, z, q)                       # th[r-1] in (al[r], al[r-1]) for r = 1..N-1 (0-based below)
    PT = [sum(al[p] * z[p] / (al[p] - th[r]) for p in range(r + 1)) for r in range(N - 1)]
    H = [P - (1 - q) for P in PT]
    C = max(H)
    ok = True
    msgs = []

    def chk(cond, msg):
        nonlocal ok
        if not cond:
            ok = False
            msgs.append(msg)

    d, b, VT, VB, U, feed = {}, {}, {}, {}, {}, {}
    cols = sorted([(i, j) for i in range(N) for j in range(i + 1, N)], key=lambda c: -(c[1] - c[0]))
    for (i, j) in cols:
        if (i, j) == (0, N - 1):
            f = {p: z[p] for p in range(N)}
            u = 1 - q
        else:
            f = {p: mp.mpf(0) for p in range(i, j + 1)}
            u = mp.mpf(0)
            if j + 1 <= N - 1:                       # distillate of [i, j+1]
                for p in range(i, j + 1):
                    f[p] += d[(i, j + 1)].get(p, 0)
                u += VT[(i, j + 1)]
            if i - 1 >= 0:                           # residue of [i-1, j]
                for p in range(i, j + 1):
                    f[p] += b[(i - 1, j)].get(p, 0)
                u -= VB[(i - 1, j)]
        feed[(i, j)], U[(i, j)] = f, u
        for p in range(i, j + 1):
            chk(f[p] > 0, f'feed {i, j} comp {p} not positive')
        # carry-over: feed equation at the column's roots theta_i..theta_{j-1}
        for r in range(i, j):
            val = sum(al[p] * f[p] / (al[p] - th[r]) for p in range(i, j + 1)) - u
            chk(abs(val) < TOL, f'feed equation {i, j} root {r}: {val}')
        if j - i >= 2:
            # preferred split: alpha_p d_p = V_T prod_r (alpha_p - th_r) / prod_{p' != p} (alpha_p - alpha_p')
            def coef(p):
                num = mp.fprod([al[p] - th[r] for r in range(i, j)])
                den = mp.fprod([al[p] - al[pp] for pp in range(i, j) if pp != p])
                return num / den / al[p]
            vt = f[i] / coef(i)
            dd = {p: vt * coef(p) for p in range(i, j)}
            bb = {p: f[p] - dd.get(p, 0) for p in range(i + 1, j + 1)}
            chk(abs(dd[i] - f[i]) < TOL, 'd_i != f_i')
            vb = vt - u
            for p in range(i, j):
                chk(dd[p] > 0, f'd {i, j} comp {p} <= 0')
            for p in range(i + 1, j + 1):
                chk(bb[p] > 0, f'b {i, j} comp {p} <= 0')
            chk(vt > 0 and vb > 0, 'vapor <= 0')
            chk(vt - sum(dd.values()) > 0 and vb + sum(bb.values()) > 0, 'liquid <= 0')
            for r in range(i, j):
                gT = sum(al[p] * dd[p] / (al[p] - th[r]) for p in dd)
                gB = -sum(al[p] * bb[p] / (al[p] - th[r]) for p in bb)
                chk(abs(gT - vt) < TOL and abs(gB - vb) < TOL, f'strict Underwood {i, j} root {r}')
        else:
            r = i
            dd = {i: f[i]}
            bb = {j: f[j]}
            yT = al[i] * f[i] / (al[i] - th[r])
            yB = al[j] * f[j] / (th[r] - al[j])
            chk(abs(yT - yB - u) < TOL, 'binary feed equation')
            e = C - H[r]
            chk(e >= 0, 'negative excess')
            vt, vb = yT + e, yB + e
            chk(vt >= yT and vb >= yB and vt - f[i] > 0 and vb > 0, 'binary column constraints')
        d[(i, j)], b[(i, j)], VT[(i, j)], VB[(i, j)] = dd, bb, vt, vb
    # side draws and products
    for p in range(1, N - 1):
        chk(abs(VT[(p, p + 1)] - VB[(p - 1, p)]) < TOL, f'side-draw vapor balance at product {p}')
        chk(abs(d[(p, p + 1)][p] + b[(p - 1, p)][p] - z[p]) < TOL, f'product amount {p}')
    chk(abs(d[(0, 1)][0] - z[0]) < TOL and abs(b[(N - 2, N - 1)][N - 1] - z[N - 1]) < TOL, 'end products')
    R = VB[(N - 2, N - 1)]
    chk(abs(R - C) < TOL, f'R = {R} differs from max H = {C}')
    return ok, msgs


def main():
    rng = random.Random(20261020)
    allok = True
    for N in range(3, 9):
        nfail = 0
        for _ in range(200):
            ratios = [mp.mpf(rng.uniform(1.01, 20.0)) for _ in range(N - 1)]
            al = [mp.mpf(1)]
            for rt in ratios:
                al.append(al[-1] * rt)
            al = al[::-1]
            z = [mp.mpf(rng.uniform(1e-3, 1)) for _ in range(N)]
            s = sum(z)
            z = [x / s for x in z]
            q = mp.mpf(rng.uniform(-1, 2))
            ok, msgs = build_and_check(al, z, q)
            if not ok:
                nfail += 1
                if nfail <= 3:
                    print('  failure:', msgs[:3])
        print(f'N = {N}: 200 instances, {nfail} failures:', 'PASS' if nfail == 0 else 'FAIL')
        allok &= nfail == 0
    print('ALL PASS' if allok else 'SOME FAIL')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
