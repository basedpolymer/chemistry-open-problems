#!/usr/bin/env python3
"""v02_dense_certify.py -- verifier's INDEPENDENT certification of Hopf points and of the sign of l1 (dense, small M).

Written by the verifier without using bdlib.py.  Different method from the solver's:
  * coordinates: the dynamics is written in reduced coordinates on the mass hyperplane H,
        y = (v_2, ..., v_M),   v_1 = -sum_{l>=2} l v_l      (v = n - (1+K)),
    so the Jacobian A_hat (size M-1) is invertible and no "w11 in H" choice is needed;
  * the Jacobian and the bilinear form are obtained from the vector field of the paper (eqs. (1)-(4)) itself:
        B e_m = (f(A + e_m) - f(A - e_m))/2   (exact for a quadratic field, computed in exact integers at K=0,1),
        Bil(X, Y) = f(A+X+Y) - f(A+X) - f(A+Y) + f(A)   (exact second derivative of a quadratic field);
    no use of the solver's matrices B, S;
  * Hopf point: Krawczyk operator on the full eigen-system in the 2(M-1)+2 real unknowns (Re q, Im q, omega, K):
        A_hat(K) q = i omega q,  q_k = 1;
  * l1: generic dense Kuznetsov formula (3.20) in ball arithmetic (acb_mat.solve), with the adjoint vector from a
    bordered linear system; trilinear term computed as well (third difference of f; it must be 0);
  * spectral hypothesis (optional, --eig): acb_mat.eig on the interval matrix A_hat(K box).
All statements printed as PASS are rigorous (ball arithmetic, Arb via python-flint), except the float label check.

Usage:  python v02_dense_certify.py OUT.csv [--eig | --eigmax MMAX] [--cert DIR] M:j [M:j ...]
        python v02_dense_certify.py OUT.csv [--eig | --eigmax MMAX] [--cert DIR] --range J MMIN MMAX [STEP]
Initial guesses: kappa from the solver's certificate files (DIR, default ../certificates) -- used only as a
starting point; the enclosure is re-proved here from scratch and compared with the solver's values afterwards.
Environment: WSL Fedora, Python 3.12, python-flint 0.9, numpy.  Runtime: ~0.2 s (M=15) to ~25 s (M=300) per point,
plus eig (~7 s at M=120).
"""
import sys, os, csv, gzip, math, time
import numpy as np
from flint import arb, acb, arb_mat, acb_mat, ctx

PREC = 192

# ------------------------------------------------------------------ the model (paper eqs. (1)-(4)), generic arithmetic
def field(n, K):
    """n: list n_1..n_M (index 0..M-1).  Returns f(n) as a list."""
    M = len(n); N = M - 1
    J = [n[l]*n[0] - n[l + 1] for l in range(N)]          # J[l-1] = J_l = n_l n_1 - n_{l+1}, l = 1..N
    f = [None]*M
    sJ = J[0]
    for l in range(1, N):
        sJ = sJ + J[l]
    f[0] = -J[0] - sJ + M*K*n[M - 1]                      # dn_1/dt = -J_1 - sum_{l=1}^N J_l + M K n_M
    for l in range(1, N):
        f[l] = J[l - 1] - J[l]                            # dn_l/dt = J_{l-1} - J_l, 2 <= l <= N
    f[M - 1] = J[N - 1] - K*n[M - 1]                      # dn_M/dt = J_N - K n_M
    return f

def jac_int(M, K):
    """exact integer Jacobian at n = 1+K for integer K (central difference is exact for a quadratic field)"""
    A = 1 + K
    cols = []
    for m in range(M):
        vp = [A]*M; vm = [A]*M; vp[m] += 1; vm[m] -= 1
        fp = field(vp, K); fm = field(vm, K)
        cols.append([(a - b)//2 for a, b in zip(fp, fm)])
        assert all((a - b) % 2 == 0 for a, b in zip(fp, fm))
    return [[cols[m][i] for m in range(M)] for i in range(M)]

def reduced_int(M):
    """A_hat(K) = A0 + K A1 on y = (v_2..v_M), v_1 = -sum l v_l; exact integers. Checks affinity in K."""
    B0 = jac_int(M, 0); B1 = [[a - b for a, b in zip(r1, r0)] for r1, r0 in zip(jac_int(M, 1), B0)]
    B2 = jac_int(M, 2)
    assert all(B2[i][m] == B0[i][m] + 2*B1[i][m] for i in range(M) for m in range(M)), "B not affine in K"
    def red(Bm):
        # (P B E)_{i,m} for i, m = 2..M  (0-based 1..M-1):  B_{i,m} + B_{i,1} * (-(m))   (l-index of column m is m+1)
        return [[Bm[i][m] - (m + 1)*Bm[i][0] for m in range(1, M)] for i in range(1, M)]
    return red(B0), red(B1), B0, B1

def E(y):
    """C^{M-1} -> H subset C^M"""
    v1 = 0
    for l, yl in enumerate(y):
        v1 = v1 - (l + 2)*yl
    return [v1] + list(y)

def bil_full(X, Y, K):
    """exact second derivative D^2 f(A)[X, Y] of the quadratic field (polarisation), ball arithmetic"""
    M = len(X); A = 1 + K
    a = [A]*M
    f_xy = field([A + x + y for x, y in zip(X, Y)], K)
    f_x = field([A + x for x in X], K)
    f_y = field([A + y for y in Y], K)
    f_0 = field(a, K)
    return [p - q - r + s for p, q, r, s in zip(f_xy, f_x, f_y, f_0)]

def tri_full(X, Y, Z, K):
    """third derivative by the third-order difference; must vanish for a quadratic field"""
    M = len(X); A = 1 + K
    def g(s):
        return field([A + sum(c*w[i] for c, w in zip(s, (X, Y, Z))) for i in range(M)], K)
    out = None
    for s in [(1,1,1),(1,1,0),(1,0,1),(0,1,1),(1,0,0),(0,1,0),(0,0,1),(0,0,0)]:
        sign = (-1)**(3 - sum(s))
        v = g(s)
        out = [sign*x for x in v] if out is None else [o + sign*x for o, x in zip(out, v)]
    return out

def bil_red(x, y, K):
    return bil_full(E(x), E(y), K)[1:]

# ------------------------------------------------------------------ float helpers
def Ahat_float(A0, A1, K):
    return np.array(A0, float) + K*np.array(A1, float)

def unstable_pairs_float(A0, A1, K):
    w = np.linalg.eigvals(Ahat_float(A0, A1, K))
    sc = max(1.0, np.abs(w).max())
    return int(np.sum((w.real > 1e-13*sc) & (w.imag > 0)))

# ------------------------------------------------------------------ Hopf system (realified) and Krawczyk
def amat(A0, A1, K):
    return arb_mat(A0) + arb(K)*arb_mat(A1) if not isinstance(K, arb) else arb_mat(A0) + K*arb_mat(A1)

def F_and_J(A0m, A1m, a, b, om, K, k):
    """A0m, A1m: arb_mat; a, b: lists of arb; om, K: arb.  Returns F (list) and J (arb_mat)."""
    n = len(a)
    Ah = A0m + K*A1m
    av = arb_mat([[x] for x in a]); bv = arb_mat([[x] for x in b])
    Aa = Ah*av; Ab = Ah*bv; A1a = A1m*av; A1b = A1m*bv
    F = [Aa[i, 0] + om*b[i] for i in range(n)] + [Ab[i, 0] - om*a[i] for i in range(n)] + [a[k] - 1, b[k]]
    N2 = 2*n + 2
    rows = []
    z = arb(0)
    for i in range(n):
        r = [Ah[i, j] for j in range(n)] + [z]*n + [b[i], A1a[i, 0]]
        r[n + i] = om
        rows.append(r)
    for i in range(n):
        r = [z]*n + [Ah[i, j] for j in range(n)] + [-a[i], A1b[i, 0]]
        r[i] = -om
        rows.append(r)
    r = [z]*N2; r[k] = arb(1); rows.append(r)
    r = [z]*N2; r[n + k] = arb(1); rows.append(r)
    return F, arb_mat(rows)

def certify_hopf(M, A0, A1, K0, om0, q0, rel=1e-36):
    """Newton (float then arb midpoint) + Krawczyk. Returns dict with boxes, or None."""
    n = M - 1
    k = int(np.argmax(np.abs(q0)))
    q0 = q0/q0[k]
    # float Newton
    Af0 = np.array(A0, float); Af1 = np.array(A1, float)
    x = np.concatenate([q0.real, q0.imag, [om0, K0]])
    for _ in range(30):
        a, b, om, K = x[:n], x[n:2*n], x[2*n], x[2*n + 1]
        Ah = Af0 + K*Af1
        F = np.concatenate([Ah@a + om*b, Ah@b - om*a, [a[k] - 1, b[k]]])
        J = np.zeros((2*n + 2, 2*n + 2))
        J[:n, :n] = Ah; J[:n, n:2*n] = om*np.eye(n); J[:n, 2*n] = b; J[:n, 2*n + 1] = Af1@a
        J[n:2*n, :n] = -om*np.eye(n); J[n:2*n, n:2*n] = Ah; J[n:2*n, 2*n] = -a; J[n:2*n, 2*n + 1] = Af1@b
        J[2*n, k] = 1; J[2*n + 1, n + k] = 1
        dx = np.linalg.solve(J, F)
        x = x - dx
        if np.abs(dx).max() < 1e-13*max(1, np.abs(x).max()):
            break
    # arb Newton (midpoints)
    ctx.prec = PREC
    A0m = arb_mat(A0); A1m = arb_mat(A1)
    X = [arb(float(v)) for v in x]
    for it in range(6):
        a, b, om, K = X[:n], X[n:2*n], X[2*n], X[2*n + 1]
        F, J = F_and_J(A0m, A1m, a, b, om, K, k)
        dX = J.solve(arb_mat([[f] for f in F]), algorithm="approx")
        X = [arb((X[i] - dX[i, 0]).mid()) for i in range(2*n + 2)]
        if max(abs(float(dX[i, 0].mid())) for i in range(2*n + 2)) < 1e-50:
            break
    # Krawczyk
    scale_q = max(abs(float(v.mid())) for v in X[:2*n])
    rad = [scale_q*rel]*(2*n) + [abs(float(X[2*n].mid()))*rel, abs(float(X[2*n + 1].mid()))*rel]
    box = [arb(X[i].mid(), rad[i]) for i in range(2*n + 2)]
    a, b, om, K = X[:n], X[n:2*n], X[2*n], X[2*n + 1]
    Fm, Jm = F_and_J(A0m, A1m, a, b, om, K, k)
    Y = Jm.solve(arb_mat([[1 if i == j else 0 for j in range(2*n + 2)] for i in range(2*n + 2)]), algorithm="approx")
    ab, bb, omb, Kb = box[:n], box[n:2*n], box[2*n], box[2*n + 1]
    _, JX = F_and_J(A0m, A1m, ab, bb, omb, Kb, k)
    I = arb_mat([[1 if i == j else 0 for j in range(2*n + 2)] for i in range(2*n + 2)])
    dev = arb_mat([[box[i] - X[i]] for i in range(2*n + 2)])
    Kx = arb_mat([[X[i]] for i in range(2*n + 2)]) - Y*arb_mat([[f] for f in Fm]) + (I - Y*JX)*dev
    ok = all(box[i].contains_interior(Kx[i, 0]) for i in range(2*n + 2))
    if not ok:
        return None
    enc = [Kx[i, 0] for i in range(2*n + 2)]
    return dict(a=enc[:n], b=enc[n:2*n], om=enc[2*n], K=enc[2*n + 1], k=k, A0m=A0m, A1m=A1m)

# ------------------------------------------------------------------ l1 (Kuznetsov), reduced coordinates, ball arithmetic
def lyapunov(M, h, with_tri=False):
    n = M - 1
    K, om = h["K"], h["om"]
    q = [acb(x, y) for x, y in zip(h["a"], h["b"])]
    qb = [x.conjugate() for x in q]
    Ah = h["A0m"] + K*h["A1m"]
    I = acb(0, 1)
    # adjoint row vector u: u A_hat = i om u  <=>  (A_hat^T - i om) u^T = 0 ; row r replaced by u_k' = 1
    qf = np.array([complex(float(z.real.mid()), float(z.imag.mid())) for z in q])
    Af = np.array([[float(Ah[i, j].mid()) for j in range(n)] for i in range(n)])
    wl, vl = np.linalg.eig(Af.T)
    iu = int(np.argmin(np.abs(wl - 1j*float(om.mid()))))
    kk = int(np.argmax(np.abs(vl[:, iu]))); r = int(np.argmax(np.abs(qf)))
    rows = [[acb(Ah[j, i]) - (I*om if i == j else 0) for j in range(n)] for i in range(n)]
    rows[r] = [acb(1) if j == kk else acb(0) for j in range(n)]
    rhs = acb_mat([[acb(1) if i == r else acb(0)] for i in range(n)])
    U = acb_mat(rows).solve(rhs)
    u = [U[i, 0] for i in range(n)]
    uq = sum((x*y for x, y in zip(u, q)), acb(0))
    simple = not uq.contains(0)     # geometric simplicity: bordered matrix nonsingular; algebraic: u q != 0
    u = [x/uq for x in u]
    A1q = acb_mat(h["A1m"])*acb_mat([[x] for x in q])
    dlam = sum((u[i]*A1q[i, 0] for i in range(n)), acb(0))      # d lambda / dK = u A1 q  (u q = 1)
    # residual checks (ball): must contain 0
    Acm = acb_mat(Ah)
    res_q = Acm*acb_mat([[x] for x in q]) - acb_mat([[I*om*x] for x in q])
    resq_ok = all(res_q[i, 0].contains(0) for i in range(n))
    # w11 = A_hat^{-1} Bil(q, qbar) ; w20 = (2 i om - A_hat)^{-1} Bil(q, q)
    b11 = bil_red(q, qb, K)
    w11 = Acm.solve(acb_mat([[x] for x in b11]))
    w11 = [w11[i, 0] for i in range(n)]
    b20 = bil_red(q, q, K)
    C20 = acb_mat([[(2*I*om if i == j else 0) - Acm[i, j] for j in range(n)] for i in range(n)])
    w20 = C20.solve(acb_mat([[x] for x in b20]))
    w20 = [w20[i, 0] for i in range(n)]
    t1 = sum((x*y for x, y in zip(u, bil_red(q, w11, K))), acb(0))
    t2 = sum((x*y for x, y in zip(u, bil_red(qb, w20, K))), acb(0))
    c = -2*t1 + t2
    tri = None
    if with_tri:
        T = tri_full(E(q), E(q), E(qb), K)[1:]
        tri = sum((x*y for x, y in zip(u, T)), acb(0))
        c = c + tri
    l1 = c.real/(2*om)
    Eq = E(q)
    nrm2 = sum((abs(x)**2 for x in Eq), arb(0))
    return dict(l1_qk=l1, l1_unit=l1/nrm2, l1_q1=l1/abs(Eq[0])**2, resq_ok=resq_ok, tri=tri, simple=simple,
                transv=dlam.real)

def spectrum(M, h, j):
    """rigorous H2 + label: eigenvalues of A_hat(K box)"""
    Ah = h["A0m"] + h["K"]*h["A1m"]
    ctx.prec = 128
    Ev = acb_mat(Ah).eig()
    ctx.prec = PREC
    on = [e for e in Ev if e.real.contains(0)]
    pos = [e for e in Ev if e.real > 0]
    neg = [e for e in Ev if e.real < 0]
    om = h["om"]
    ok = (len(on) == 2 and all(abs(e.imag).overlaps(om) for e in on) and len(pos) == 2*(j - 1)
          and len(on) + len(pos) + len(neg) == M - 1)
    return ok, (len(on), len(pos), len(neg))

# ------------------------------------------------------------------ driver
def load_cert(dirn, j):
    out = {}
    for fn in ("sweep_j%d.csv.gz" % j, "sweep_j%d_b.csv.gz" % j):
        p = os.path.join(dirn, fn)
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(gzip.open(p, "rt")):
            if r["status"] in ("PASS", "FAIL"):
                out[int(r["M"])] = r
    return out

def run_point(M, j, cert, do_eig):
    t0 = time.time()
    ctx.prec = PREC
    A0, A1, _, _ = reduced_int(M)
    r = cert[M]
    K0 = float(r["kappa_mid"])/math.sqrt(M); om0 = float(r["omega_lo"])
    Af = Ahat_float(A0, A1, K0)
    w, V = np.linalg.eig(Af)
    i = int(np.argmin(np.abs(w - 1j*om0)))
    h = certify_hopf(M, A0, A1, K0, om0, V[:, i])
    if h is None:
        return dict(M=M, j=j, status="KRAWCZYK_FAIL", t=time.time() - t0)
    L = lyapunov(M, h, with_tri=(M <= 40))
    # float label check: j-1 unstable pairs just below, j just above
    Kc = float(h["K"].mid())
    lab = (unstable_pairs_float(A0, A1, Kc*(1 - 1e-6)), unstable_pairs_float(A0, A1, Kc*(1 + 1e-6)))
    eigok, eiginfo = (None, None)
    if do_eig:
        eigok, eiginfo = spectrum(M, h, j)
    sub = (M, j) == (15, 1)
    l1 = L["l1_unit"]
    sign_ok = (l1 > 0) if sub else (l1 < 0)
    # comparison with the solver's certificate (float64 stored values)
    their_l1 = float(r["l1_unit_lo"]); their_q1 = float(r["l1_q1_lo"])
    mine = float(l1.mid()); mine_q1 = float(L["l1_q1"].mid())
    agree = abs(mine - their_l1) <= 1e-9*abs(their_l1) and abs(mine_q1 - their_q1) <= 1e-9*abs(their_q1)
    kap = h["K"]*arb(M).sqrt()
    kap_agree = abs(float(kap.mid()) - float(r["kappa_mid"])) <= 1e-15*float(r["kappa_mid"])
    om_pos = bool(h["om"] > 0)
    ok = (sign_ok and om_pos and L["resq_ok"] and L["simple"] and bool(L["transv"] > 0) and lab == (j - 1, j)
          and (eigok in (None, True)))
    tri = L["tri"]
    return dict(M=M, j=j, kappa=float(kap.mid()), kappa_rad=float(kap.rad()), omega=float(h["om"].mid()),
                l1_unit=mine, l1_unit_rad=float(l1.rad()), l1_q1=mine_q1, their_l1=their_l1, agree=agree,
                kap_agree=kap_agree, label=lab, simple=L["simple"], transv_lo=float(L["transv"].lower()), eig=eiginfo, eigok=eigok,
                tri=("0" if tri is None or tri.contains(0) else "NONZERO") if tri is not None else "",
                status="PASS" if ok else "FAIL", t=time.time() - t0)

def main():
    args = sys.argv[1:]
    out = args.pop(0)
    do_eig = "--eig" in args
    if do_eig:
        args.remove("--eig")
    eigmax = 10**9
    if "--eigmax" in args:
        i = args.index("--eigmax"); eigmax = int(args[i + 1]); del args[i:i + 2]; do_eig = True
    certdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificates")
    if "--cert" in args:
        i = args.index("--cert"); certdir = args[i + 1]; del args[i:i + 2]
    pts = []
    if args and args[0] == "--range":
        J, a, b = int(args[1]), int(args[2]), int(args[3]); st = int(args[4]) if len(args) > 4 else 1
        pts = [(M, J) for M in range(a, b + 1, st)]
    else:
        pts = [tuple(int(s) for s in x.split(":")) for x in args]
    certs = {}
    nfail = 0; t0 = time.time()
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh)
        keys = ["M", "j", "kappa", "kappa_rad", "omega", "l1_unit", "l1_unit_rad", "l1_q1", "their_l1", "agree",
                "kap_agree", "label", "simple", "transv_lo", "eig", "eigok", "tri", "status", "t"]
        wr.writerow(keys)
        for (M, j) in pts:
            if j not in certs:
                certs[j] = load_cert(certdir, j)
            res = run_point(M, j, certs[j], do_eig and M <= eigmax)
            ok = res["status"] == "PASS" and res.get("agree", False) and res.get("kap_agree", False)
            nfail += not ok
            wr.writerow([res.get(k_, "") for k_ in keys]); fh.flush()
            print("M=%d j=%d kappa=%.12g l1_unit=%.10e (+-%.1e) solver=%.10e agree=%s label=%s simple=%s dRe/dK>%.3g eig=%s tri=%s %s %.1fs"
                  % (M, j, res.get("kappa", float("nan")), res.get("l1_unit", float("nan")),
                     res.get("l1_unit_rad", float("nan")), res.get("their_l1", float("nan")), res.get("agree"),
                     res.get("label"), res.get("simple"), res.get("transv_lo", float("nan")), res.get("eig"), res.get("tri"), res["status"] if ok else "FAIL", res["t"]),
                  flush=True)
    print("SUMMARY: %d points, %s (%d failures), %.0fs" % (len(pts), "PASS" if nfail == 0 else "FAIL", nfail,
                                                            time.time() - t0))
    sys.exit(1 if nfail else 0)

if __name__ == "__main__":
    main()
