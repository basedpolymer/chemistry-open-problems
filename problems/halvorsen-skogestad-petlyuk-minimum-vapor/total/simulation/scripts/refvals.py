"""
refvals.py -- reference (Underwood / Halvorsen-Skogestad) values, used ONLY to compare with the stage-by-stage
results (never inside the physical model).  Agent P (Claude Opus 5.5, effort xhigh).
"""
import mpmath as mp

mp.mp.dps = 40


def feed_roots(alpha, z, q):
    """Roots theta_j in (alpha_{j+1}, alpha_j) of sum alpha z/(alpha - theta) = 1 - q."""
    a = [mp.mpf(v) for v in alpha]; zz = [mp.mpf(v) for v in z]
    N = len(a)
    f = lambda t: sum(a[p] * zz[p] / (a[p] - t) for p in range(N)) - (1 - mp.mpf(q))
    roots = []
    for j in range(N - 1):
        lo, hi = a[j + 1], a[j]   # f -> -inf at lo+, +inf at hi- (f increasing on the interval)
        for _ in range(200):
            mid = (lo + hi) / 2
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        roots.append((lo + hi) / 2)
    return roots


def peaks(alpha, z, q):
    """H_j = sum_{p>j} alpha_p z_p/(theta_j - alpha_p) (reboiler-referenced peaks), list and max."""
    th = feed_roots(alpha, z, q)
    N = len(alpha)
    H = [sum(mp.mpf(alpha[p]) * mp.mpf(z[p]) / (th[j] - mp.mpf(alpha[p])) for p in range(j + 1, N))
         for j in range(N - 1)]
    return [float(h) for h in H], float(max(H)), [float(t) for t in th]


def binary_vmin(alpha, z, q, eps):
    """Infinite-stage minimum reboiler vapor of a binary column with symmetric purity 1 - eps (Underwood)."""
    th = feed_roots(alpha, z, q)[0]
    a1, a2 = map(mp.mpf, alpha)
    z1 = mp.mpf(z[0]); e = mp.mpf(eps)
    D = (z1 - e) / (1 - 2 * e)
    VT = a1 * D * (1 - e) / (a1 - th) + a2 * D * e / (a2 - th)
    return float(VT - (1 - mp.mpf(q)))
