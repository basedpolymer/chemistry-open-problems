"""deletion_formula.py -- the local vertex-deletion formula (Proposition 6 of result.md).

For a connected graph G on n vertices with Gamma = L^+ and a non-cut vertex v with neighbourhood N (|N| = d):
  M  = grounded Green function at v:  M_ab = G_ab - G_av - G_vb + G_vv         (G = Gamma)
  M2 = (M^2)_ab = (G^2)_ab - (G^2)_av - (G^2)_vb + (G^2)_vv + n (G_av - G_vv)(G_bv - G_vv)
  s_a   = (M 1)_a   = n (G_vv - G_av)
  t_a   = (M^2 1)_a = -n [(G^2)_av - (G^2)_vv] - n^2 G_vv (G_av - G_vv)
  1^T M 1 = n^2 G_vv,   1^T M^2 1 = n^2 (G^2)_vv + n^3 G_vv^2,   tr M = Rt(v) = tr G + n G_vv
  Z = [[M_NN - I, s_N], [s_N^T, n - 1 + n^2 G_vv]],  W = [[M2_NN, t_N], [t_N^T, 1^T M^2 1]]
  Kf(G - v) = (n - 1) [ Rt(v) - 1 - tr(Z^{-1} W) ],      Kf(G) = n tr G.
Only entries of Gamma and Gamma^2 on the closed neighbourhood of v enter (besides tr Gamma, G_vv).
Functions take a callable gam(x, y), gam2(x, y) so that the same code serves exact (fmpq), interval
(arb) and floating-point evaluations.
"""


def deletion_quantities(n, v, N, gam, gam2, trG, one, zero, matinv_trace):
    """Return (Kf(G), Kf(G - v)) using the local formula. `one`, `zero` are the scalars of the field,
    matinv_trace(Z, W) must return tr(Z^{-1} W) for (d+1)x(d+1) nested lists Z, W."""
    gvv = gam(v, v)
    g2vv = gam2(v, v)
    gav = {a: gam(a, v) for a in N}
    g2av = {a: gam2(a, v) for a in N}
    d = len(N)
    Z = [[zero] * (d + 1) for _ in range(d + 1)]
    W = [[zero] * (d + 1) for _ in range(d + 1)]
    for i, a in enumerate(N):
        for j, b in enumerate(N):
            m = gam(a, b) - gav[a] - gav[b] + gvv
            m2 = gam2(a, b) - g2av[a] - g2av[b] + g2vv + n * (gav[a] - gvv) * (gav[b] - gvv)
            Z[i][j] = m - (one if i == j else zero)
            W[i][j] = m2
        s = n * (gvv - gav[a])
        t = -n * (g2av[a] - g2vv) - n * n * gvv * (gav[a] - gvv)
        Z[i][d] = Z[d][i] = s
        W[i][d] = W[d][i] = t
    Z[d][d] = (n - 1) * one + n * n * gvv
    W[d][d] = n * n * g2vv + n ** 3 * gvv * gvv
    Rt = trG + n * gvv
    kfG = n * trG
    kfH = (n - 1) * (Rt - one - matinv_trace(Z, W))
    return kfG, kfH
