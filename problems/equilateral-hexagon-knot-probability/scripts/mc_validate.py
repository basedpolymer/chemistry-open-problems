#!/usr/bin/env python3
"""NUMERICAL validation (not part of the proof) of the knotting criteria used in result.md.

Samples equilateral hexagons from the symplectic (standard) measure through T135 action-angle
coordinates (d uniform in the moment polytope P, theta uniform in [0, 2pi)^3), and for every sample
computes, in floating point:
  * Calvo's algebraic intersection numbers Delta_1..Delta_6 of the triangles T_i = (v_{i-1}, v_i, v_{i+1});
  * the geometric piercing counts of all six triangles;
  * the pierce patterns A12, A23, A13, B12, B23, B13 of result.md (Section 3);
  * the number of Fox 3-colourings of the projection on a fixed random plane (an exact knot invariant:
    3 for the unknot, 9 for the trefoil), computed from the Gauss diagram with integer arithmetic mod 3.
Checks (each prints PASS/FAIL):
  C1  Fox-3 count in {3, 9} for every sample;
  C2  Fox-3 == 9  <=>  Delta_2, Delta_4, Delta_6 all nonzero (Calvo's criterion);
  C3  knotted => every triangle T_1..T_6 is pierced exactly once (Lemma 2.3 of result.md);
  C4  never (A12 and B12), (A23 and B23), (A13 and B13): the two conditions on the same angle pair exclude
      each other (Lemma 2.1 / Proposition 3.2);
  C5  for theta in (0,pi)^3: knotted <=> (A12 & A23 & A13) or (B12 & B23 & B13); R-knots have Delta = +1;
  C6  knotted => theta in (0,pi)^3 or (pi,2pi)^3;
  C7  the Monte Carlo knotting frequency is compatible (|z| < 4) with the certified interval given on the
      command line (optional).
Environment: Python 3.13, numpy 2.5, numba 0.68.  Runtime: ~3 min for 1e8 samples (one thread).
Usage: python3 mc_validate.py N SEED [LOWER UPPER]
"""
import sys, time
import numpy as np
import numba as nb

@nb.njit(cache=True)
def hexagon(d1, d2, d3, t1, t2, t3, V):
    D2 = 2*d1*d1*d2*d2 + 2*d1*d1*d3*d3 + 2*d2*d2*d3*d3 - d1**4 - d2**4 - d3**4
    D = np.sqrt(max(D2, 0.0))
    v5x = (d1*d1 - d2*d2 + d3*d3) / (2*d1); v5y = D / (2*d1)
    V[0, 0] = 0.; V[0, 1] = 0.; V[0, 2] = 0.
    V[2, 0] = d1; V[2, 1] = 0.; V[2, 2] = 0.
    V[4, 0] = v5x; V[4, 1] = v5y; V[4, 2] = 0.
    for k in range(3):
        if k == 0:
            ax, ay, bx, by, dd, t = 0., 0., d1, 0., d1, t1
        elif k == 1:
            ax, ay, bx, by, dd, t = d1, 0., v5x, v5y, d2, t2
        else:
            ax, ay, bx, by, dd, t = v5x, v5y, 0., 0., d3, t3
        mx = (ax + bx)/2; my = (ay + by)/2
        nx = -(by - ay)/dd; ny = (bx - ax)/dd
        h = np.sqrt(max(1 - dd*dd/4, 0.0))
        V[2*k+1, 0] = mx + h*np.cos(t)*nx
        V[2*k+1, 1] = my + h*np.cos(t)*ny
        V[2*k+1, 2] = h*np.sin(t)

@nb.njit(cache=True)
def orient(V, a, b, c, d):
    ux = V[b, 0]-V[a, 0]; uy = V[b, 1]-V[a, 1]; uz = V[b, 2]-V[a, 2]
    vx = V[c, 0]-V[a, 0]; vy = V[c, 1]-V[a, 1]; vz = V[c, 2]-V[a, 2]
    wx = V[d, 0]-V[a, 0]; wy = V[d, 1]-V[a, 1]; wz = V[d, 2]-V[a, 2]
    return ux*(vy*wz-vz*wy) - uy*(vx*wz-vz*wx) + uz*(vx*wy-vy*wx)

@nb.njit(cache=True)
def pierce(V, p, q, a, b, c):
    """signed crossing of the directed segment p->q through the oriented triangle (a,b,c); 0 if none"""
    sp = orient(V, a, b, c, p); sq = orient(V, a, b, c, q)
    if sp*sq >= 0:
        return 0
    s1 = orient(V, p, q, a, b); s2 = orient(V, p, q, b, c); s3 = orient(V, p, q, c, a)
    if (s1 > 0 and s2 > 0 and s3 > 0) or (s1 < 0 and s2 < 0 and s3 < 0):
        return 1 if sq > 0 else -1
    return 0

@nb.njit(cache=True)
def tri_stats(V, Dl, G):
    for c in range(6):
        a = (c-1) % 6; b = (c+1) % 6
        s = 0; g = 0
        for e in range(6):
            p = e; q = (e+1) % 6
            if p == a or p == b or p == c or q == a or q == b or q == c:
                continue
            x = pierce(V, p, q, a, c, b)
            s += x
            if x != 0:
                g += 1
        Dl[c] = s; G[c] = g

@nb.njit(cache=True)
def fox3(V, R):
    W = np.zeros((6, 3))
    for i in range(6):
        for j in range(3):
            W[i, j] = R[j, 0]*V[i, 0] + R[j, 1]*V[i, 1] + R[j, 2]*V[i, 2]
    ncr = 0
    cr_tu = np.zeros(9); cr_to = np.zeros(9)
    for i in range(6):
        for j in range(i+2, 6):
            if i == 0 and j == 5:
                continue
            p = W[i]; q = W[(i+1) % 6]; r = W[j]; s = W[(j+1) % 6]
            ex = q[0]-p[0]; ey = q[1]-p[1]; fx = s[0]-r[0]; fy = s[1]-r[1]
            den = ex*fy - ey*fx
            if den == 0:
                continue
            gx = r[0]-p[0]; gy = r[1]-p[1]
            t = (gx*fy - gy*fx)/den; u = (gx*ey - gy*ex)/den
            if t > 0 and t < 1 and u > 0 and u < 1:
                zi = p[2] + t*(q[2]-p[2]); zj = r[2] + u*(s[2]-r[2])
                if zi > zj:
                    cr_tu[ncr] = j + u; cr_to[ncr] = i + t
                else:
                    cr_tu[ncr] = i + t; cr_to[ncr] = j + u
                ncr += 1
    if ncr == 0:
        return 3
    us = np.sort(cr_tu[:ncr])
    M = np.zeros((ncr, ncr), np.int64)
    for k in range(ncr):
        x = cr_tu[k]
        idx = 0
        for m in range(ncr):
            if us[m] < x:
                idx += 1
        a_out = idx % ncr; a_in = (idx-1) % ncr
        y = cr_to[k]
        idy = 0
        for m in range(ncr):
            if us[m] < y:
                idy += 1
        a_ov = (idy-1) % ncr
        M[k, a_ov] = (M[k, a_ov]+2) % 3
        M[k, a_in] = (M[k, a_in]+2) % 3
        M[k, a_out] = (M[k, a_out]+2) % 3
    rank = 0
    for col in range(ncr):
        piv = -1
        for r in range(rank, ncr):
            if M[r, col] % 3 != 0:
                piv = r; break
        if piv < 0:
            continue
        for c2 in range(ncr):
            tmp = M[rank, c2]; M[rank, c2] = M[piv, c2]; M[piv, c2] = tmp
        inv = 1 if M[rank, col] % 3 == 1 else 2
        for c2 in range(ncr):
            M[rank, c2] = (M[rank, c2]*inv) % 3
        for r in range(ncr):
            if r != rank and M[r, col] % 3 != 0:
                f = M[r, col]
                for c2 in range(ncr):
                    M[r, c2] = (M[r, c2] - f*M[rank, c2]) % 3
        rank += 1
    return 3**(ncr-rank)

@nb.njit(cache=True)
def run(n, seed, R):
    np.random.seed(seed)
    V = np.zeros((6, 3)); Dl = np.zeros(6, np.int64); G = np.zeros(6, np.int64)
    # counters
    c = np.zeros(16, np.int64)
    # 0 samples, 1 knotted(fox9), 2 fox not in {3,9}, 3 fox9 != calvo, 4 knotted & some triangle count != 1,
    # 5 incompatibility violated, 6 upper-region knots, 7 upper knots not in R u L, 8 R u L (upper) not knotted,
    # 9 knots outside (0,pi)^3 u (pi,2pi)^3, 10 R-knots with Delta != +1, 11 L-knots with Delta != -1,
    # 12 R count (upper), 13 L count (upper)
    for it in range(n):
        while True:
            d1 = 2*np.random.random(); d2 = 2*np.random.random(); d3 = 2*np.random.random()
            if d1 <= d2 + d3 and d2 <= d1 + d3 and d3 <= d1 + d2:
                break
        t1 = 2*np.pi*np.random.random(); t2 = 2*np.pi*np.random.random(); t3 = 2*np.pi*np.random.random()
        hexagon(d1, d2, d3, t1, t2, t3, V)
        tri_stats(V, Dl, G)
        f = fox3(V, R)
        kn = f == 9
        calvo = Dl[1] != 0 and Dl[3] != 0 and Dl[5] != 0
        c[0] += 1
        if kn:
            c[1] += 1
        if f != 3 and f != 9:
            c[2] += 1
        if kn != calvo:
            c[3] += 1
        if kn:
            for i in range(6):
                if G[i] != 1:
                    c[4] += 1
                    break
        A12 = pierce(V, 3, 4, 0, 1, 2) != 0; B12 = pierce(V, 0, 1, 2, 3, 4) != 0
        A23 = pierce(V, 5, 0, 2, 3, 4) != 0; B23 = pierce(V, 2, 3, 4, 5, 0) != 0
        A13 = pierce(V, 1, 2, 4, 5, 0) != 0; B13 = pierce(V, 4, 5, 0, 1, 2) != 0
        if (A12 and B12) or (A23 and B23) or (A13 and B13):
            c[5] += 1
        upper = t1 < np.pi and t2 < np.pi and t3 < np.pi
        lower = t1 > np.pi and t2 > np.pi and t3 > np.pi
        if kn and not (upper or lower):
            c[9] += 1
        if upper:
            inR = A12 and A23 and A13
            inL = B12 and B23 and B13
            if inR:
                c[12] += 1
            if inL:
                c[13] += 1
            if kn:
                c[6] += 1
                if not (inR or inL):
                    c[7] += 1
            if (inR or inL) and not kn:
                c[8] += 1
            if inR and kn and not (Dl[1] == 1 and Dl[3] == 1 and Dl[5] == 1):
                c[10] += 1
            if inL and kn and not (Dl[1] == -1 and Dl[3] == -1 and Dl[5] == -1):
                c[11] += 1
    return c

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10**8
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 2026
    lo = float(sys.argv[3]) if len(sys.argv) > 3 else None
    hi = float(sys.argv[4]) if len(sys.argv) > 4 else None
    rng = np.random.default_rng(12345)
    R, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    t0 = time.time()
    tot = np.zeros(16, np.int64)
    chunk = 10**7
    k = 0
    while k < n:
        m = min(chunk, n - k)
        tot += run(m, seed + k // chunk, R)
        k += m
    c = tot
    p = c[1]/c[0]; se = np.sqrt(p*(1-p)/c[0])
    print("samples %d  knotted %d  p_MC = %.6e +- %.2e (1 s.e.)  time %.0f s" % (c[0], c[1], p, se, time.time()-t0))
    print("upper-region knots %d, R-pattern %d, L-pattern %d" % (c[6], c[12], c[13]))
    ok = True
    def chk(name, cond, info=""):
        nonlocal ok
        print(("PASS" if cond else "FAIL"), name, info)
        ok = ok and cond
    chk("C1 Fox-3 count in {3,9}", c[2] == 0, "(violations %d)" % c[2])
    chk("C2 Fox-3 knotted <=> Calvo Delta_2,4,6 nonzero", c[3] == 0, "(disagreements %d)" % c[3])
    chk("C3 knotted => each T_i pierced exactly once", c[4] == 0, "(violations %d)" % c[4])
    chk("C4 incompatibility of A/B on a shared angle pair", c[5] == 0, "(violations %d)" % c[5])
    chk("C5 upper region: knotted <=> R or L", c[7] == 0 and c[8] == 0, "(knot not in RuL %d, RuL not knot %d)" % (c[7], c[8]))
    chk("C5' R-knots right-handed, L-knots left-handed (Delta = +1 / -1)", c[10] == 0 and c[11] == 0, "(%d, %d)" % (c[10], c[11]))
    chk("C6 knotted => theta in (0,pi)^3 or (pi,2pi)^3", c[9] == 0, "(violations %d)" % c[9])
    if lo is not None:
        z_lo = (p - lo)/se; z_hi = (hi - p)/se
        chk("C7 MC frequency compatible with certified interval [%.4e, %.4e]" % (lo, hi), z_lo > -4 and z_hi > -4,
            "(z_lo %.2f, z_hi %.2f)" % (z_lo, z_hi))
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
