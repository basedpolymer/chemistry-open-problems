"""Exact geometry of a triangulated periodic partition and the DFM energy.

For each cell i the boundary is the union of the oriented triangles of its
faces; with a reference point r_i the cell is the signed union of the tetrahedra
(r_i, a, b, c), so that
    V_i = sum det(a,b,c)/6
    M_i = sum V_t (a+b+c)/4            (first moment about r_i)
    S_i = sum V_t (|a|^2+|b|^2+|c|^2+a.b+b.c+c.a)/10   (second moment about r_i)
exactly for the polyhedral cell; I_i = S_i - |M_i|^2 / V_i is the second moment
about the centroid.  A(X) = (2 Atot/n)/(4 pi R0^2), I(X) = (sum I_i/n)/(4 pi R0^5/5),
R0 = (3 Vbox/(4 pi n))^(1/3), F = (A^2 I)^(1/3).
"""
import numpy as np
import numba as nb

FOURPI = 4.0 * np.pi


@nb.njit(cache=True)
def _cell_sums(X, L, tri, tsh, tface, fi, fj, fsj, ref, ncell):
    T = tri.shape[0]
    V = np.zeros(ncell)
    M = np.zeros((ncell, 3))
    S = np.zeros(ncell)
    Atot = 0.0
    P = np.zeros((3, 3))
    for t in range(T):
        for k in range(3):
            v = tri[t, k]
            for c in range(3):
                P[k, c] = X[v, c] + L[c, 0] * tsh[t, k, 0] + L[c, 1] * tsh[t, k, 1] + L[c, 2] * tsh[t, k, 2]
        e1x = P[1, 0] - P[0, 0]; e1y = P[1, 1] - P[0, 1]; e1z = P[1, 2] - P[0, 2]
        e2x = P[2, 0] - P[0, 0]; e2y = P[2, 1] - P[0, 1]; e2z = P[2, 2] - P[0, 2]
        nx = e1y * e2z - e1z * e2y
        ny = e1z * e2x - e1x * e2z
        nz = e1x * e2y - e1y * e2x
        Atot += 0.5 * np.sqrt(nx * nx + ny * ny + nz * nz)
        f = tface[t]
        for side in range(2):
            if side == 0:
                ci = fi[f]
                ox = ref[ci, 0]; oy = ref[ci, 1]; oz = ref[ci, 2]
                sg = 1.0
            else:
                ci = fj[f]
                ox = ref[ci, 0] + L[0, 0] * fsj[f, 0] + L[0, 1] * fsj[f, 1] + L[0, 2] * fsj[f, 2]
                oy = ref[ci, 1] + L[1, 0] * fsj[f, 0] + L[1, 1] * fsj[f, 1] + L[1, 2] * fsj[f, 2]
                oz = ref[ci, 2] + L[2, 0] * fsj[f, 0] + L[2, 1] * fsj[f, 1] + L[2, 2] * fsj[f, 2]
                sg = -1.0
            ax = P[0, 0] - ox; ay = P[0, 1] - oy; az = P[0, 2] - oz
            bx = P[1, 0] - ox; by = P[1, 1] - oy; bz = P[1, 2] - oz
            cx = P[2, 0] - ox; cy = P[2, 1] - oy; cz = P[2, 2] - oz
            det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
            vt = sg * det / 6.0
            V[ci] += vt
            M[ci, 0] += vt * (ax + bx + cx) / 4.0
            M[ci, 1] += vt * (ay + by + cy) / 4.0
            M[ci, 2] += vt * (az + bz + cz) / 4.0
            q = (ax * ax + ay * ay + az * az + bx * bx + by * by + bz * bz + cx * cx + cy * cy + cz * cz
                 + ax * bx + ay * by + az * bz + bx * cx + by * cy + bz * cz + cx * ax + cy * ay + cz * az)
            S[ci] += vt * q / 10.0
    return Atot, V, M, S


@nb.njit(cache=True)
def _grad(X, L, tri, tsh, tface, fi, fj, fsj, ref, ncell, gA, gV, gM, gS):
    """Gradient of  gA*Atot + sum_i (gV_i V_i + gM_i . M_i + gS_i S_i)  w.r.t. X."""
    T = tri.shape[0]
    G = np.zeros(X.shape)
    P = np.zeros((3, 3))
    R = np.zeros((3, 3))
    D = np.zeros((3, 3))
    for t in range(T):
        for k in range(3):
            v = tri[t, k]
            for c in range(3):
                P[k, c] = X[v, c] + L[c, 0] * tsh[t, k, 0] + L[c, 1] * tsh[t, k, 1] + L[c, 2] * tsh[t, k, 2]
        for k in range(3):
            for c in range(3):
                D[k, c] = 0.0
        # area gradient
        e1x = P[1, 0] - P[0, 0]; e1y = P[1, 1] - P[0, 1]; e1z = P[1, 2] - P[0, 2]
        e2x = P[2, 0] - P[0, 0]; e2y = P[2, 1] - P[0, 1]; e2z = P[2, 2] - P[0, 2]
        nx = e1y * e2z - e1z * e2y
        ny = e1z * e2x - e1x * e2z
        nz = e1x * e2y - e1y * e2x
        nn = np.sqrt(nx * nx + ny * ny + nz * nz)
        if nn > 0.0:
            ux = nx / nn; uy = ny / nn; uz = nz / nn
            for k in range(3):
                kp = (k + 1) % 3
                km = (k + 2) % 3
                # dA/dP_k = 0.5 * u x (P_km - P_kp)
                wx = P[km, 0] - P[kp, 0]; wy = P[km, 1] - P[kp, 1]; wz = P[km, 2] - P[kp, 2]
                D[k, 0] += gA * 0.5 * (uy * wz - uz * wy)
                D[k, 1] += gA * 0.5 * (uz * wx - ux * wz)
                D[k, 2] += gA * 0.5 * (ux * wy - uy * wx)
        f = tface[t]
        for side in range(2):
            if side == 0:
                ci = fi[f]
                ox = ref[ci, 0]; oy = ref[ci, 1]; oz = ref[ci, 2]
                sg = 1.0
            else:
                ci = fj[f]
                ox = ref[ci, 0] + L[0, 0] * fsj[f, 0] + L[0, 1] * fsj[f, 1] + L[0, 2] * fsj[f, 2]
                oy = ref[ci, 1] + L[1, 0] * fsj[f, 0] + L[1, 1] * fsj[f, 1] + L[1, 2] * fsj[f, 2]
                oz = ref[ci, 2] + L[2, 0] * fsj[f, 0] + L[2, 1] * fsj[f, 1] + L[2, 2] * fsj[f, 2]
                sg = -1.0
            for k in range(3):
                R[k, 0] = P[k, 0] - ox; R[k, 1] = P[k, 1] - oy; R[k, 2] = P[k, 2] - oz
            ax = R[0, 0]; ay = R[0, 1]; az = R[0, 2]
            bx = R[1, 0]; by = R[1, 1]; bz = R[1, 2]
            cx = R[2, 0]; cy = R[2, 1]; cz = R[2, 2]
            det = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
            vt = sg * det / 6.0
            sx = ax + bx + cx; sy = ay + by + cy; sz = az + bz + cz
            q = (ax * ax + ay * ay + az * az + bx * bx + by * by + bz * bz + cx * cx + cy * cy + cz * cz
                 + ax * bx + ay * by + az * bz + bx * cx + by * cy + bz * cz + cx * ax + cy * ay + cz * az)
            gv = gV[ci]; gmx = gM[ci, 0]; gmy = gM[ci, 1]; gmz = gM[ci, 2]; gs = gS[ci]
            gms = gmx * sx + gmy * sy + gmz * sz
            # coefficient multiplying dV_t
            cV = gv + gms / 4.0 + gs * q / 10.0
            for k in range(3):
                kp = (k + 1) % 3
                km = (k + 2) % 3
                # dV_t/dR_k = sg/6 * (R_kp x R_km)
                dvx = sg / 6.0 * (R[kp, 1] * R[km, 2] - R[kp, 2] * R[km, 1])
                dvy = sg / 6.0 * (R[kp, 2] * R[km, 0] - R[kp, 0] * R[km, 2])
                dvz = sg / 6.0 * (R[kp, 0] * R[km, 1] - R[kp, 1] * R[km, 0])
                # dq/dR_k = 2 R_k + R_kp + R_km = R_k + s
                dqx = R[k, 0] + sx; dqy = R[k, 1] + sy; dqz = R[k, 2] + sz
                D[k, 0] += cV * dvx + vt / 4.0 * gmx + gs * vt / 10.0 * dqx
                D[k, 1] += cV * dvy + vt / 4.0 * gmy + gs * vt / 10.0 * dqy
                D[k, 2] += cV * dvz + vt / 4.0 * gmz + gs * vt / 10.0 * dqz
        for k in range(3):
            v = tri[t, k]
            G[v, 0] += D[k, 0]; G[v, 1] += D[k, 1]; G[v, 2] += D[k, 2]
    return G


def geometry(M, X=None):
    if X is None:
        X = M.X
    return _cell_sums(X, M.L, M.tri, M.tsh, M.tface, M.fi, M.fj, M.fsj, M.ref, M.ncell)


def dfm_values(M, X=None):
    """Return dict with F, A, I and per-cell volumes etc. (fixed box)."""
    Atot, V, Mo, S = geometry(M, X)
    n = M.ncell
    Vbox = abs(np.linalg.det(M.L))
    R0 = (3.0 * Vbox / (FOURPI * n)) ** (1.0 / 3.0)
    I = S - np.sum(Mo ** 2, axis=1) / V
    Abar = (2.0 * Atot / n) / (FOURPI * R0 ** 2)
    Ibar = (np.sum(I) / n) / (FOURPI * R0 ** 5 / 5.0)
    F = (Abar ** 2 * Ibar) ** (1.0 / 3.0)
    return {"F": F, "A": Abar, "I": Ibar, "V": V, "Vsum": V.sum(), "Vbox": Vbox,
            "Icell": I, "Atot": Atot, "R0": R0}


def dfm_F_and_grad(M, X):
    """F and dF/dX for fixed box."""
    Atot, V, Mo, S = geometry(M, X)
    n = M.ncell
    Vbox = abs(np.linalg.det(M.L))
    R0 = (3.0 * Vbox / (FOURPI * n)) ** (1.0 / 3.0)
    Mo2 = np.sum(Mo ** 2, axis=1)
    I = S - Mo2 / V
    sumI = I.sum()
    Abar = (2.0 * Atot / n) / (FOURPI * R0 ** 2)
    Ibar = (sumI / n) / (FOURPI * R0 ** 5 / 5.0)
    F = (Abar ** 2 * Ibar) ** (1.0 / 3.0)
    gA = 2.0 * F / (3.0 * Atot)
    cI = F / (3.0 * sumI)
    gS = np.full(n, cI)
    gM = -cI * 2.0 * Mo / V[:, None]
    gV = cI * Mo2 / V ** 2
    G = _grad(X, M.L, M.tri, M.tsh, M.tface, M.fi, M.fj, M.fsj, M.ref, n, gA, gV, gM, gS)
    return F, G


def area_volume_grad(M, X, gA, gV):
    """Gradient of gA*Atot + sum gV_i V_i (used for constrained problems)."""
    n = M.ncell
    return _grad(X, M.L, M.tri, M.tsh, M.tface, M.fi, M.fj, M.fsj, M.ref, n, gA,
                 np.asarray(gV, dtype=float), np.zeros((n, 3)), np.zeros(n))
