"""
v2_focal_numeric.py -- VERIFIER's generic floating-point focal values (mpmath, 60 digits), real-coordinate
Lyapunov-function method, for planar systems  x' = sum c_i x^al_i y^be_i,  y' = sum d_i x^al_i y^be_i  with an
equilibrium at (1,1) and tr J = 0.  Purpose: check that the normalisation L_k = pi * V_k (coordinates
xi = x-1, eta = -(J11 u + J12 v)/omega, time omega t, F = xi^2+eta^2+...) reproduces the focal values printed by
Boros-Hofbauer for OTHER systems (Sec. 3.2 quadrangle, Sec. 4.1 chain), hence that our L_k are exactly theirs,
and evaluate L_1..L_4 of (BH) at the certified point p* and at the paper's two checkpoints.
Not part of any proof (floating point).  Runtime: ~10 s.
"""
import sys
import mpmath as mp

mp.mp.dps = 60
results = []


def report(name, ok, info=""):
    results.append(bool(ok))
    print(("PASS" if ok else "FAIL") + "  " + name + (("  [" + info + "]") if info else ""), flush=True)


def binom(al, i):
    r = mp.mpf(1)
    for k in range(i):
        r = r * (al - k) / (k + 1)
    return r


def smul(P, Q, N):
    out = {}
    for (i1, j1), c1 in P.items():
        for (i2, j2), c2 in Q.items():
            if i1 + j1 + i2 + j2 <= N:
                k = (i1 + i2, j1 + j2)
                out[k] = out.get(k, 0) + c1 * c2
    return out


def expand(terms, N):
    """sum coef*(1+u)^al (1+v)^be as series in (u, v) up to degree N."""
    S = {}
    for cf, al, be in terms:
        for i in range(N + 1):
            bi = binom(al, i)
            for j in range(N + 1 - i):
                S[(i, j)] = S.get((i, j), 0) + cf * bi * binom(be, j)
    return S


def compose_linear(S, lu, lv, N):
    """S(u, v) with u = lu(xi,eta), v = lv(xi,eta) (linear dicts) -> series in (xi, eta)."""
    up = [{(0, 0): mp.mpf(1)}]
    vp = [{(0, 0): mp.mpf(1)}]
    for _ in range(N):
        up.append(smul(up[-1], lu, N))
        vp.append(smul(vp[-1], lv, N))
    out = {}
    for (i, j), cf in S.items():
        if cf == 0:
            continue
        for k, v in smul(up[i], vp[j], N).items():
            out[k] = out.get(k, 0) + cf * v
    return out


def focal(fterms, gterms, kmax):
    N = 2 * kmax + 1
    Sf, Sg = expand(fterms, N), expand(gterms, N)
    j11, j12, j21, j22 = Sf.get((1, 0), 0), Sf.get((0, 1), 0), Sg.get((1, 0), 0), Sg.get((0, 1), 0)
    tr, det = j11 + j22, j11 * j22 - j12 * j21
    om = mp.sqrt(det)
    lu = {(1, 0): mp.mpf(1)}
    lv = {(1, 0): -j11 / j12, (0, 1): -om / j12}
    F1 = compose_linear(Sf, lu, lv, N)
    G1 = compose_linear(Sg, lu, lv, N)
    P = {k: v / om for k, v in F1.items()}
    Q = {k: -(j11 * F1.get(k, 0) + j12 * G1.get(k, 0)) / om**2 for k in set(F1) | set(G1)}
    # linear part check
    lin_ok = (abs(P.get((1, 0), 0)) < mp.mpf(10)**-50 and abs(P.get((0, 1), 0) + 1) < mp.mpf(10)**-50
              and abs(Q.get((1, 0), 0) - 1) < mp.mpf(10)**-50 and abs(Q.get((0, 1), 0)) < mp.mpf(10)**-50)
    if abs(tr) > mp.mpf(10)**-40 or not lin_ok:
        raise RuntimeError("not a weak focus in normal form (tr=%s)" % tr)
    Ph, Qh = {}, {}
    for (i, j), v in P.items():
        if i + j >= 2:
            Ph.setdefault(i + j, [0] * (i + j + 1))[i] += v
    for (i, j), v in Q.items():
        if i + j >= 2:
            Qh.setdefault(i + j, [0] * (i + j + 1))[i] += v

    def hmul(p, q):
        out = [0] * (len(p) + len(q) - 1)
        for i1, c1 in enumerate(p):
            for i2, c2 in enumerate(q):
                out[i1 + i2] += c1 * c2
        return out

    def dxi(p):
        return [p[i] * i for i in range(1, len(p))]

    def deta(p):
        n = len(p) - 1
        return [p[i] * (n - i) for i in range(n)]

    def cmean(i, j):
        if i % 2 or j % 2:
            return mp.mpf(0)
        p, q = i // 2, j // 2
        return mp.factorial(2 * p) * mp.factorial(2 * q) / (4 ** (p + q) * mp.factorial(p) * mp.factorial(q)
                                                              * mp.factorial(p + q))
    F = {2: [mp.mpf(1), mp.mpf(0), mp.mpf(1)]}
    V = []
    for n in range(3, 2 * kmax + 3):
        T = [mp.mpf(0)] * (n + 1)
        for m in range(2, n):
            k = n - m + 1
            t1 = hmul(dxi(F[m]), Ph[k])
            t2 = hmul(deta(F[m]), Qh[k])
            for i in range(n + 1):
                T[i] += t1[i] + t2[i]
        dim = n + 1 + (1 if n % 2 == 0 else 0)
        M = mp.zeros(dim, dim)
        for i in range(n + 1):
            if i >= 1:
                M[i - 1, i] += -i
            if i + 1 <= n:
                M[i + 1, i] += (n - i)
        rhs = mp.matrix([-t for t in T] + ([0] if n % 2 == 0 else []))
        if n % 2 == 0:
            for i in range(0, n + 1, 2):
                M[i, n + 1] = -mp.binomial(n // 2, i // 2)
            for i in range(n + 1):
                M[n + 1, i] = cmean(i, n - i)
        sol = mp.lu_solve(M, rhs)
        F[n] = [sol[i] for i in range(n + 1)]
        if n % 2 == 0:
            V.append(sol[n + 1])
    return [mp.pi * v for v in V], om


def bh(a, b, d):
    K = a / (1 + b * d)
    return [(1, 0, -1), (-1, a, b)], [(-K * (d - 1), 0, 0), (-K, 0, -1), (K * d, a, b)]


def main():
    # ---- quadrangle (Sec. 3.2): x' = 1 - x y^5, y' = K[-g y + 2 + (g+2) x y^2 - 4 x y^5], g = 16 + 1/K
    def quad(K):
        g = 16 + 1 / K
        return [(1, 0, 0), (-1, 1, 5)], [(-K * g, 0, 1), (2 * K, 0, 0), (K * (g + 2), 1, 2), (-4 * K, 1, 5)]
    for Kv in (mp.mpf('0.3'), mp.mpf('1.7')):
        L, _ = focal(*quad(Kv), 1)
        ref = mp.pi * (3416 * Kv**3 + 1250 * Kv**2 - 29 * Kv - 5) / (20 * mp.sqrt(2 * (2 + 35 * Kv)**3))
        report("quadrangle L_1 at K=%s matches the paper's closed form" % Kv, abs(L[0] / ref - 1) < 1e-40,
               "ours %s, paper %s" % (mp.nstr(L[0], 12), mp.nstr(ref, 12)))
    K0 = mp.findroot(lambda K: 3416 * K**3 + 1250 * K**2 - 29 * K - 5, 0.0686)
    L, _ = focal(*quad(K0), 2)
    report("quadrangle L_2 at K_0=%s is ~0.01293 (paper)" % mp.nstr(K0, 6), abs(L[1] - mp.mpf('0.01293')) < 5e-6,
           "ours L_1 = %s, L_2 = %s" % (mp.nstr(L[0], 3), mp.nstr(L[1], 8)))

    # ---- chain (Sec. 4.1): x' = y^-q - x y^(1/2), y' = K[-(q+r+1/2) + (q+1/2) y^-q + r x y^(1/2)], K = 2/(r-q(2q+1))
    def chain(q, r):
        K = 2 / (r - q * (2 * q + 1))
        return ([(1, 0, -q), (-1, 1, mp.mpf(1) / 2)],
                [(-K * (q + r + mp.mpf(1) / 2), 0, 0), (K * (q + mp.mpf(1) / 2), 0, -q), (K * r, 1, mp.mpf(1) / 2)])
    q, r = mp.mpf('0.3'), mp.mpf('2.2')
    L, _ = focal(*chain(q, r), 1)
    ref = mp.pi * r * (3 * r * (1 - 2 * q) - q * (4 * q**2 + 16 * q + 7)) / (
        8 * (2 * q + 1) * (r - q * (2 * q + 1))**mp.mpf(1.5) * mp.sqrt(2 * q * (q + r + mp.mpf(1) / 2)))
    report("chain L_1 at (q,r)=(0.3,2.2) matches the paper's closed form", abs(L[0] / ref - 1) < 1e-40,
           "ours %s, paper %s" % (mp.nstr(L[0], 12), mp.nstr(ref, 12)))
    for q in (mp.mpf('0.1'), mp.mpf('0.4')):
        r = q * (4 * q**2 + 16 * q + 7) / (3 * (1 - 2 * q))
        L, _ = focal(*chain(q, r), 2)
        ref = mp.pi * (2 * q + 7)**2 * (3 - 2 * q) * (4 * q - 1) * mp.sqrt(2 * q + 3) / (1536 * (2 * q + 1)**4)
        report("chain L_2 on {L_1=0} at q=%s matches the paper's closed form" % q,
               abs(L[0]) < 1e-45 and abs(L[1] / ref - 1) < 1e-40,
               "ours %s, paper %s" % (mp.nstr(L[1], 12), mp.nstr(ref, 12)))
    q = mp.mpf(1) / 4
    r = q * (4 * q**2 + 16 * q + 7) / (3 * (1 - 2 * q))
    L, _ = focal(*chain(q, r), 3)
    ref = -625 * mp.pi / 110592 * mp.sqrt(mp.mpf(7) / 2)
    report("chain L_3 at q=1/4, r=15/8 equals -625 pi/110592 sqrt(7/2) (paper)",
           abs(L[0]) < 1e-45 and abs(L[1]) < 1e-45 and abs(L[2] / ref - 1) < 1e-40,
           "ours %s, paper %s" % (mp.nstr(L[2], 15), mp.nstr(ref, 15)))

    # ---- BH: the paper's checkpoints on {L1 = L2 = 0}
    def b_of(a, d):   # L_1 = 0 branch given in the paper
        return (-2 + a * (1 + mp.sqrt(1 + 8 * d))) / (2 * d)
    a, d = mp.mpf(1), mp.mpf(165) / 49
    L, _ = focal(*bh(a, b_of(a, d), d), 3)
    report("BH checkpoint (a,d)=(1,165/49): L1 = L2 = 0, L3 < 0", abs(L[0]) < 1e-45 and abs(L[1]) < 1e-45 and L[2] < 0,
           "b = %s, L = %s" % (mp.nstr(b_of(a, d), 12), [mp.nstr(x, 6) for x in L]))
    a, d = (1 + mp.sqrt(3961)) / 60, mp.mpf(3)
    L, _ = focal(*bh(a, b_of(a, d), d), 3)
    report("BH checkpoint (a,d)=((1+sqrt3961)/60,3): L1 = L2 = 0, L3 > 0",
           abs(L[0]) < 1e-45 and abs(L[1]) < 1e-45 and L[2] > 0, "L = %s" % [mp.nstr(x, 6) for x in L])
    # ---- BH at the certified point p* (centre of the solver's enclosures, 40 digits)
    a = mp.mpf('1.012821998680314600068208268899051820974')
    b = mp.mpf('0.654625621104238210727259033435915709692')
    d = mp.mpf('3.28861043930993578772416135623261176832')
    L, om = focal(*bh(a, b, d), 4)
    report("BH at p*: |L1|,|L2|,|L3| < 1e-35 and L4 = -2.580126499494e-6",
           max(abs(x) for x in L[:3]) < 1e-35 and abs(L[3] / mp.mpf("-2.5801264994944244863565347627e-6") - 1) < 1e-27,
           "L = %s, omega = %s" % ([mp.nstr(x, 8) for x in L], mp.nstr(om, 15)))
    # ---- the paper's own approximate point (5 digits): L1..L3 small, L4 < 0
    L, _ = focal(*bh(mp.mpf('1.01282'), mp.mpf('0.65463'), mp.mpf('3.28862')), 4)
    print("     at the paper's (1.01282, 0.65463, 3.28862): L =", [mp.nstr(x, 6) for x in L])
    print("ALL PASS" if all(results) else "SOME CHECKS FAILED")
    return all(results)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
