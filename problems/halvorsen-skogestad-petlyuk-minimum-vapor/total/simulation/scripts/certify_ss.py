"""
certify_ss.py -- computer-assisted proof (interval Krawczyk test) that a finite-stage FTC (generalized Petlyuk)
arrangement, operated at exactly given rational flows, has a steady state in which every product has impurity at
most eps_cert, with total vaporization R exactly known.

Agent P (Claude Opus 5.5, effort xhigh).  Environment: WSL Fedora, Python 3.12, numpy/scipy; runtime ~10 s - few
min depending on size.  Usage:
    python certify_ss.py DATASET DELTA N_STAGES [r]
(DATASET as in attain.py; the flows are result 12's Underwood attainment flows + a circulation DELTA through
every column (ftc_fixed_flows.py), rounded to binary64 and then corrected EXACTLY (rational arithmetic) to satisfy
every node balance; the Underwood construction only chooses the flows.)

Model (stagenet.py): ideal VLE y_p = alpha_p x_p / S(x), S = sum alpha_k x_k; CMO component balances at every
stage  F_{s,p}(x) = feed_{s,p} + sum_c f_c c_{src,p} - Lout_s x_{s,p} - Vout_s y_{s,p} = 0.  With exactly balanced
flows, sum_p x_{s,p} = 1 at any solution (sum of the component balances).
Proof: change of variables x = m o u (m = numerical solution, rounded scalings are powers of two), equations
scaled by powers of two w; Krawczyk operator K(X) = 1 - Y Ft(1) + (I - Y Jt(X))(X - 1) on X = [1-r, 1+r]^M,
Y ~ Jt(1)^{-1}.  If |Y Ft(1)| + |I - Y Jt(X)| r < r componentwise, F has a unique zero in X (Krawczyk/Moore).
Interval enclosures: ivl.py (outward rounding by one ulp per operation); dense products bounded with the classical
error bound |fl(AB) - AB| <= gamma_k |A||B|, gamma_k = k u/(1 - k u), u = 2^-53.
Prints PASS/FAIL; writes ../runs/cert_DATASET_delta_n.json; exit code 1 on failure.
"""
import sys, json, time
from fractions import Fraction
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from ivl import I, isum, frac_enclosure, up, dn
from arrangements import regular_network, ftc_columns
from attain import DATA, cex5
from ftc_fixed_flows import flows_for
from refvals import peaks

U = 2.0 ** -53


def gamma(k):
    return k * U / (1 - k * U)


def exact_data(ds):
    if ds == 'CEX5':
        d = json.load(open('../../../verification/permissive_counterexample_N5.json'))
        return [Fraction(s) for s in d['alpha']], [Fraction(s) for s in d['z']], Fraction(d['q'])
    a, z, q = DATA[ds]
    if ds.startswith('T61'):
        zz = [Fraction(1, 3)] * 3
    else:
        zz = [Fraction(str(v)) for v in z]
    return [Fraction(v) for v in a], zz, Fraction(str(q))


def exact_balanced_flows(net, f_float, F, q):
    """Rationalize the float flows and correct them exactly so that every node balance holds; adiabatic duties
    stay exactly 0.  Returns a list of Fractions."""
    A, _ = net.balance_matrix()
    K = len(net.nodes)
    b = [Fraction(0)] * (2 * K)
    for fd in net.feeds:
        b[2 * fd['node']] += F * (1 - q)
        b[2 * fd['node'] + 1] += F * q
    fr = [Fraction(float(v)) for v in f_float]
    fixed = set()
    for k, nm in enumerate(net.nodes):
        if net.node_kind[k] == 'adiabatic':
            i = net.fidx['Q:' + nm]; fixed.add(i); fr[i] = Fraction(0)
    cols = [j for j in range(net.nf) if j not in fixed]
    Ar = [[Fraction(int(A[r, j])) for j in cols] for r in range(2 * K)]
    e = [b[r] - sum(Fraction(int(A[r, j])) * fr[j] for j in range(net.nf) if A[r, j] != 0) for r in range(2 * K)]
    # exact RREF of [Ar | e]
    M = [row[:] + [e[r]] for r, row in enumerate(Ar)]
    nr, nc = len(M), len(cols)
    piv = []
    r = 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [v * inv for v in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                fct = M[i][c]
                M[i] = [vi - fct * vr for vi, vr in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == nr:
            break
    for i in range(r, nr):
        assert M[i][-1] == 0, 'inconsistent balances'
    for i, c in enumerate(piv):
        fr[cols[c]] += M[i][-1]
    # verify
    for rr in range(2 * K):
        assert sum(Fraction(int(A[rr, j])) * fr[j] for j in range(net.nf) if A[rr, j] != 0) == b[rr]
    return fr


def pow2(v):
    return 2.0 ** np.round(np.log2(v))


def main():
    ds = sys.argv[1]; delta = float(sys.argv[2]); n = int(sys.argv[3])
    r = float(sys.argv[4]) if len(sys.argv) > 4 else None
    t0 = time.time()
    aF, zF, qF = exact_data(ds)
    alpha = [float(v) for v in aF]; z = [float(v) for v in zF]; q = float(qF)
    N = len(alpha)
    net = regular_network(alpha, z, q, ftc_columns(N), {}, n=n)
    f_float, Hmax, ncol = flows_for(net, alpha, z, q, delta)
    fr = exact_balanced_flows(net, f_float, Fraction(1), qF)
    f = np.array([float(v) for v in fr])
    R_exact = sum(fr[net.fidx['Q:' + nm]] for k, nm in enumerate(net.nodes) if net.node_kind[k] == 'heater')
    # numerical solution
    x, info = net.solve(f)
    assert info['ok'], 'no numerical steady state'
    nst, M = net.nst, net.nst * N
    # interval data
    flo = np.empty(net.nf); fhi = np.empty(net.nf)
    for i, v in enumerate(fr):
        flo[i], fhi[i] = frac_enclosure(v)
    fI = I(flo, fhi)
    alo = np.array([frac_enclosure(v)[0] for v in aF]); ahi = np.array([frac_enclosure(v)[1] for v in aF])
    assert np.all(alo == ahi), 'alpha must be binary64-exact'
    a = alo
    feedI = I(np.zeros((nst, N)))
    for fd in net.feeds:
        k = fd['node']
        for p in range(N):
            lo, hi = frac_enclosure(zF[p])        # F = 1
            feedI.lo[k, p] += lo; feedI.hi[k, p] += hi
    # outflow totals (rigorous sums of nonnegative flow intervals)
    Lo_lo = np.zeros(nst); Lo_hi = np.zeros(nst); Vo_lo = np.zeros(nst); Vo_hi = np.zeros(nst)
    for (s, i) in net.Lo:
        Lo_lo[s] = dn(Lo_lo[s] + flo[i]); Lo_hi[s] = up(Lo_hi[s] + fhi[i])
    for (s, i) in net.Vo:
        Vo_lo[s] = dn(Vo_lo[s] + flo[i]); Vo_hi[s] = up(Vo_hi[s] + fhi[i])
    LoutI = I(Lo_lo, Lo_hi); VoutI = I(Vo_lo, Vo_hi)
    # connection layers (k-th incoming connection of each stage)
    incoming = [[] for _ in range(nst)]
    for c in range(len(net.c_dst)):
        incoming[net.c_dst[c]].append(c)
    layers = []
    for l in range(max(len(v) for v in incoming)):
        st = np.array([s for s in range(nst) if len(incoming[s]) > l], int)
        cc = np.array([incoming[s][l] for s in st], int)
        layers.append((st, cc))

    def yS(X):
        aX = X * I.point(np.broadcast_to(a, X.shape))
        S = isum(aX, 1)
        Sb = I(np.repeat(S.lo[:, None], N, 1), np.repeat(S.hi[:, None], N, 1))
        return aX / Sb, S

    def yS_signed(X):
        aX = X * I.point(np.broadcast_to(a, X.shape))
        S = isum(aX, 1)
        assert np.all(S.lo > 0), 'S not positive on the box'
        Sb = I(np.repeat(S.lo[:, None], N, 1), np.repeat(S.hi[:, None], N, 1))
        return aX / Sb, S

    def residual(X):
        Y, S = yS(X)
        Lb = I(np.repeat(LoutI.lo[:, None], N, 1), np.repeat(LoutI.hi[:, None], N, 1))
        Vb = I(np.repeat(VoutI.lo[:, None], N, 1), np.repeat(VoutI.hi[:, None], N, 1))
        Rr = feedI - Lb * X - Vb * Y
        for st, cc in layers:
            src = net.c_src[cc]; ph = net.c_ph[cc]
            C = I(np.where(ph[:, None] == 0, X.lo[src], Y.lo[src]), np.where(ph[:, None] == 0, X.hi[src], Y.hi[src]))
            fc = fI[net.c_fi[cc]]
            fcb = I(np.repeat(fc.lo[:, None], N, 1), np.repeat(fc.hi[:, None], N, 1))
            add = fcb * C
            part = Rr[st] + add
            Rr.lo[st] = part.lo; Rr.hi[st] = part.hi
        return Rr

    # scalings (powers of two): x = mh * u, F scaled by 1/wh
    Lout, Vout = net.outflows(f)
    w = net.weights(x, f, Lout, Vout)
    TAU = 1e-30                     # compositions below TAU are measured in absolute terms
    mh = pow2(np.maximum(x, TAU)); wh = pow2(np.maximum(w, TAU))
    # ---- high-precision center c (python-flint arb, 256 bits) by iterative refinement with the double LU
    from flint import arb, fmpq, ctx
    import scipy.sparse.linalg as spla
    ctx.prec = 256
    fA = [arb(fmpq(v.numerator, v.denominator)) for v in fr]
    aA = [arb(fmpq(v.numerator, v.denominator)) for v in aF]
    feedA = [[arb(0)] * N for _ in range(nst)]
    for fd in net.feeds:
        for p in range(N):
            feedA[fd['node']][p] = feedA[fd['node']][p] + arb(fmpq(zF[p].numerator, zF[p].denominator))
    LoutA = [arb(0)] * nst; VoutA = [arb(0)] * nst
    for (s_, i) in net.Lo:
        LoutA[s_] = LoutA[s_] + fA[i]
    for (s_, i) in net.Vo:
        VoutA[s_] = VoutA[s_] + fA[i]
    conns = [[] for _ in range(nst)]
    for cdx in range(len(net.c_dst)):
        conns[net.c_dst[cdx]].append((int(net.c_src[cdx]), int(net.c_ph[cdx]), fA[net.c_fi[cdx]]))

    def residual_hp(C):
        Yh = []
        for s_ in range(nst):
            ax = [aA[p] * C[s_][p] for p in range(N)]
            S_ = ax[0]
            for p in range(1, N):
                S_ = S_ + ax[p]
            Yh.append([ax[p] / S_ for p in range(N)])
        Rh = []
        for s_ in range(nst):
            row = []
            for p in range(N):
                v = feedA[s_][p] - LoutA[s_] * C[s_][p] - VoutA[s_] * Yh[s_][p]
                for (src, ph, fl) in conns[s_]:
                    v = v + fl * (C[src][p] if ph == 0 else Yh[src][p])
                row.append(v)
            Rh.append(row)
        return Rh

    def mag_up(v):
        return up(up(abs(float(v.mid())) + float(v.rad())))

    C = [[arb(float(x[s_, p])) for p in range(N)] for s_ in range(nst)]
    J = net.jacobian(x, f)
    Jt = sp.csc_matrix(sp.diags(1.0 / wh.ravel()) @ J @ sp.diags(mh.ravel()))
    luJ = spla.splu(Jt)
    for it in range(12):
        Rh = residual_hp(C)
        Ft = np.array([[float(Rh[s_][p].mid()) / wh[s_, p] for p in range(N)] for s_ in range(nst)])
        if np.max(np.abs(Ft)) < 1e-40:
            break
        du = luJ.solve(-Ft.ravel()).reshape(x.shape)
        dx = du * mh
        C = [[C[s_][p] + arb(float(dx[s_, p])) for p in range(N)] for s_ in range(nst)]
    # make the center an exact (zero-radius) high-precision point
    C = [[arb(C[s_][p].mid()) for p in range(N)] for s_ in range(nst)]
    Rh = residual_hp(C)
    Ft_mag = np.array([[mag_up(Rh[s_][p]) / wh[s_, p] for p in range(N)] for s_ in range(nst)])
    cf = np.array([[float(C[s_][p].mid()) for p in range(N)] for s_ in range(nst)])
    x = cf                                    # float image of the center (used for the box enclosure)
    uc = np.maximum(cf, TAU) / mh      # box half-width in u-space = uc * r
    print('center residual (scaled) max', float(np.max(Ft_mag)), 'refinement iterations', it, flush=True)
    # approximate inverse of the scaled Jacobian at the center
    J = net.jacobian(x, f)
    Jt = sp.diags(1.0 / wh.ravel()) @ J @ sp.diags(mh.ravel())
    Jc_dense = Jt.toarray()
    Y = sla.inv(Jc_dense)
    absY = np.abs(Y)
    b_vec = (absY @ Ft_mag.ravel()) * (1 + 2 * gamma(M) + 4 * U)

    def jac_enclosure(rad):
        """Mid/rad sparse enclosure of Jt over the box u in uc*[1-rad,1+rad] (relative), i.e. x in x(1 -+ rad)."""
        sp_ = np.spacing(x)
        rho_x = up(np.maximum(x, TAU) * rad)
        X = I(dn(dn(x - rho_x) - sp_), up(up(x + rho_x) + sp_))
        Yv, S = yS_signed(X)
        rows, cols, lo, hi = [], [], [], []
        Ib = np.eye(N)
        # Jy[s,p,q] = (alpha_p delta_pq - y_p alpha_q)/S
        Sb = I(np.repeat(np.repeat(S.lo[:, None, None], N, 1), N, 2), np.repeat(np.repeat(S.hi[:, None, None], N, 1), N, 2))
        yb = I(np.repeat(Yv.lo[:, :, None], N, 2), np.repeat(Yv.hi[:, :, None], N, 2))
        aq = I.point(np.broadcast_to(a[None, None, :], (nst, N, N)))
        ad = I.point(np.broadcast_to(np.diag(a)[None], (nst, N, N)))
        Jy = (ad - yb * aq) / Sb
        Lb = I(np.repeat(np.repeat(LoutI.lo[:, None, None], N, 1), N, 2), np.repeat(np.repeat(LoutI.hi[:, None, None], N, 1), N, 2))
        Vb = I(np.repeat(np.repeat(VoutI.lo[:, None, None], N, 1), N, 2), np.repeat(np.repeat(VoutI.hi[:, None, None], N, 1), N, 2))
        D = (-(Lb * I.point(np.broadcast_to(Ib[None], (nst, N, N))))) - Vb * Jy
        st = np.arange(nst)
        ii, jj = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
        rows.append((st[:, None, None] * N + ii[None]).ravel()); cols.append((st[:, None, None] * N + jj[None]).ravel())
        lo.append(D.lo.ravel()); hi.append(D.hi.ravel())
        nc = len(net.c_dst)
        fc = fI[net.c_fi]
        fcb = I(np.repeat(np.repeat(fc.lo[:, None, None], N, 1), N, 2), np.repeat(np.repeat(fc.hi[:, None, None], N, 1), N, 2))
        src = net.c_src
        JyS = Jy[src]
        Bm = I(np.where(net.c_ph[:, None, None] == 0, np.broadcast_to(Ib, (nc, N, N)), JyS.lo),
               np.where(net.c_ph[:, None, None] == 0, np.broadcast_to(Ib, (nc, N, N)), JyS.hi))
        B = fcb * Bm
        rows.append((net.c_dst[:, None, None] * N + ii[None]).ravel())
        cols.append((src[:, None, None] * N + jj[None]).ravel())
        lo.append(B.lo.ravel()); hi.append(B.hi.ravel())
        rows = np.concatenate(rows); cols = np.concatenate(cols); lo = np.concatenate(lo); hi = np.concatenate(hi)
        # scale: entry (row i, col j) times mh_j / wh_i (powers of two: exact)
        sc = mh.ravel()[cols] / wh.ravel()[rows]
        lo = lo * sc; hi = hi * sc
        key = rows.astype(np.int64) * M + cols
        assert len(np.unique(key)) == len(key), 'duplicate Jacobian entries'
        mid = 0.5 * (lo + hi)
        rad_ = np.maximum(up(hi - mid), up(mid - lo))
        Jm = sp.csc_matrix((mid, (rows, cols)), shape=(M, M))
        Jr = sp.csc_matrix((rad_, (rows, cols)), shape=(M, M))
        return Jm, Jr

    def test(rad):
        Jm, Jr = jac_enclosure(rad)
        k = int(np.max(np.diff(Jm.indptr)))
        absJm = abs(Jm).tocsc()
        ucr = uc.ravel() * rad
        Er = np.empty(M)
        B = 1024
        JmT = Jm.T.tocsr(); absJmT = absJm.T.tocsr(); JrT = Jr.T.tocsr()
        for i0 in range(0, M, B):
            i1 = min(M, i0 + B)
            Yb = Y[i0:i1]; aYb = absY[i0:i1]
            P = (JmT @ Yb.T).T                                  # rows i0:i1 of Y @ Jm
            Eye = np.zeros_like(P); Eye[np.arange(i1 - i0), np.arange(i0, i1)] = 1.0
            E = np.abs(Eye - P) * (1 + 2 * U)
            E += gamma(k) * (1 + 2 * gamma(k)) * (absJmT @ aYb.T).T
            E += (JrT @ aYb.T).T * (1 + 2 * gamma(k))
            Er[i0:i1] = (E @ ucr) * (1 + 2 * gamma(M))
        ok = bool(np.all(up(b_vec + Er) < ucr * (1 - 1e-12)))
        rho = float(np.max(Er / ucr))
        return ok, rho, float(np.max(b_vec / uc.ravel()))

    if r is None:
        okt, rho0, bmax = test(1e-12)
        r = max(1e-12, 4 * bmax / max(1e-300, 1 - min(rho0, 0.5)))
        print('trial radius', r, 'rho0', rho0, 'bmax', bmax, flush=True)
    ok, rho, bmax = test(r)
    # impurities over the box
    imps = {}
    for pr in net.products:
        k, p = pr['node'], pr['comp']
        hiv = up(up(x[k] + up(np.maximum(x[k], TAU) * r)) + np.spacing(x[k]))
        tot = 0.0
        for qq in range(N):
            if qq != p:
                tot = up(tot + hiv[qq])
        imps[pr['name']] = tot
    H, VF, th = peaks(alpha, z, q)
    out = dict(dataset=ds, delta=delta, n=n, M=M, R_exact=str(R_exact), R=float(R_exact), V_FTC=VF,
               R_minus_VFTC=float(R_exact) - VF, box_relative_radius=r, max_row_sum=rho, max_Y_F=bmax,
               krawczyk_pass=ok, impurity_upper_bounds=imps, eps_cert=max(imps.values()),
               flows={nm: str(fr[i]) for i, nm in enumerate(net.fnames)}, seconds=round(time.time() - t0, 1))
    json.dump(out, open('../runs/cert_%s_%g_%d.json' % (ds, delta, n), 'w'), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'flows'}))
    print('Krawczyk inclusion (existence of a steady state in the box):', 'PASS' if ok else 'FAIL')
    print('every product impurity <= %.3e, R = %.10f (V_FTC = %.10f):' % (out['eps_cert'], float(R_exact), VF),
          'PASS' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
