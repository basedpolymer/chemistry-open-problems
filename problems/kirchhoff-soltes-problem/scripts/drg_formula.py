"""drg_formula.py -- exact Kf(G) and Kf(G - v) for a distance-regular graph from its intersection array
and the spectrum of its local graph (Proposition 7 of result.md).  Exact rational arithmetic (fractions).

Intersection array {b_0,...,b_{D-1}; c_1,...,c_D}, k = b_0, k_{i+1} = k_i b_i / c_{i+1}, n = sum k_i.
Radial Green functions (layer-current argument):  if L psi = f with f radial, sum_l k_l f_l = 0, then
   k_i b_i (psi_i - psi_{i+1}) = sum_{l<=i} k_l f_l,   normalised by sum_l k_l psi_l = 0.
Gamma e_v: f = e_v - 1/n.  Gamma^2 e_v: f = Gamma e_v.  Gamma_i, Gamma2_i = values at distance i.
Local graph (neighbours of v): lambda-regular, lambda = a_1 = k - b_1 - 1, eigenvalue lambda on 1 and
`theta` = the other k-1 eigenvalues (all 0 when a_1 = 0).  Two neighbours of v are at distance 1 or 2.
Then Z_NN = (m0-1) I + m1 A + m2 (J-I-A), W_NN = w0 I + w1 A + w2 (J-I-A) and the (k+1)x(k+1) trace
tr(Z^{-1} W) splits into k-1 scalar terms and one 2x2 block (Z_NN 1 = 0 because M_NN 1 = 1).
"""
from fractions import Fraction as Fr


def drg_data(b, c):
    D = len(b)
    assert len(c) == D
    k = b[0]
    ks = [1]
    for i in range(D):
        ks.append(Fr(ks[-1] * b[i], c[i]))
    for x in ks:
        assert x.denominator == 1, 'infeasible array (non-integral k_i)'
    ks = [int(x) for x in ks]
    n = sum(ks)
    bb = list(b) + [0]
    return k, ks, n, bb, D


def radial_solve(ks, bb, n, f):
    D = len(ks) - 1
    psi = [Fr(0)]
    acc = Fr(0)
    for i in range(D):
        acc += ks[i] * f[i]
        psi.append(psi[-1] - acc / (ks[i] * bb[i]))
    assert ks[D] * f[D] + acc == 0 or abs(sum(ks[l] * f[l] for l in range(D + 1))) == 0
    mean = sum(ks[i] * psi[i] for i in range(D + 1)) / n
    return [p - mean for p in psi]


def drg_kf_pair(b, c, theta=None):
    """Return (n, k, Kf(G), Kf(G - v)) exactly.  theta: the k-1 local-graph eigenvalues other than a_1
    (None means a_1 = 0, local graph empty, theta = [0]*(k-1))."""
    k, ks, n, bb, D = drg_data(b, c)
    a1 = k - b[1] - 1 if D >= 2 else k - 1   # D = 1: complete graph
    if theta is None:
        assert a1 == 0, 'need local spectrum when a_1 > 0'
        theta = [0] * (k - 1)
    assert len(theta) == k - 1
    f = [Fr(1) - Fr(1, n)] + [Fr(-1, n)] * D
    G = radial_solve(ks, bb, n, f)
    G2 = radial_solve(ks, bb, n, G)
    g = lambda i: G[i] if i <= D else None
    g2 = lambda i: G2[i] if i <= D else None
    G0, G1 = G[0], G[1]
    G2_ = G[2] if D >= 2 else None
    H0, H1 = G2[0], G2[1]
    H2_ = G2[2] if D >= 2 else None
    # M entries among neighbours a,b of v (dist(a,v)=dist(b,v)=1)
    m0 = 2 * (G0 - G1)                 # a = b
    m1 = G1 - 2 * G1 + G0              # a ~ b
    m2 = (G2_ - 2 * G1 + G0) if D >= 2 else Fr(0)   # dist(a,b) = 2
    cross = n * (G1 - G0) ** 2
    w0 = H0 - 2 * H1 + H0 + cross
    w1 = H1 - 2 * H1 + H0 + cross
    w2 = ((H2_ - 2 * H1 + H0) if D >= 2 else Fr(0)) + cross
    s = n * (G0 - G1)
    t = -n * (H1 - H0) - n * n * G0 * (G1 - G0)
    z22 = (n - 1) + n * n * G0
    w22 = n * n * H0 + n ** 3 * G0 * G0
    # row sums on the all-ones vector of N
    zsum = (m0 - 1) + m1 * a1 + m2 * (k - 1 - a1)
    wsum = w0 + w1 * a1 + w2 * (k - 1 - a1)
    assert zsum == 0, 'M_NN 1 = 1 violated'
    tr = Fr(0)
    for th in theta:
        z = (m0 - 1 - m2) + (m1 - m2) * th
        w = (w0 - w2) + (w1 - w2) * th
        tr += w / z
    # 2x2 block in basis (1_N/sqrt(k), e_border): Z = [[0, s sqrt k],[s sqrt k, z22]], W = [[wsum, t sqrt k],[., w22]]
    tr += (2 * s * k * t - z22 * wsum) / (s * s * k)
    Rt = 2 * n * G0          # = tr Gamma + n Gamma_vv with tr Gamma = n Gamma_0
    kfG = n * n * G0
    kfH = (n - 1) * (Rt - 1 - tr)
    return n, k, kfG, kfH
