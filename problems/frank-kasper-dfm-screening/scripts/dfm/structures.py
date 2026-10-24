"""Periodic point sets: crystallographic structures and space-group expansion."""
import numpy as np


def expand(spacegroup, sites, tol=1e-6):
    """sites: list of (label, (x,y,z)) fractional; returns labels, P (orbits, in [0,1))."""
    import gemmi
    sg = gemmi.find_spacegroup_by_name(spacegroup)
    ops = list(sg.operations())
    labels, pts = [], []
    for lab, xyz in sites:
        orb = []
        for op in ops:
            p = np.array(op.apply_to_xyz(list(map(float, xyz))))
            p = p - np.floor(p + 1e-12)
            p[np.abs(p - 1.0) < 1e-9] = 0.0
            if not any(np.linalg.norm(((p - q) + 0.5) % 1.0 - 0.5) < tol for q in orb):
                orb.append(p)
        for p in orb:
            labels.append(lab)
            pts.append(p)
    return labels, np.array(pts)


def lattice_from_params(a, b, c, alpha=90.0, beta=90.0, gamma=90.0):
    al, be, ga = np.radians([alpha, beta, gamma])
    ax = np.array([a, 0, 0])
    bx = np.array([b * np.cos(ga), b * np.sin(ga), 0])
    cx_ = c * np.cos(be)
    cy_ = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    cz_ = np.sqrt(max(c * c - cx_ ** 2 - cy_ ** 2, 0.0))
    return np.column_stack([ax, bx, np.array([cx_, cy_, cz_])])


def hex_to_ortho(L, P, labels):
    """Equivalent orthorhombic (C-centred) cell a, b*sqrt(3), c with twice the sites."""
    a = L[0, 0]
    c = L[2, 2]
    Lo = np.diag([a, a * np.sqrt(3.0), c])
    pts, labs = [], []
    for lab, p in zip(labels, P):
        cart = L @ p
        for sh in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (-1, 0, 0), (0, -1, 0), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (2, 1, 0), (1, 2, 0)]:
            x = cart + L @ np.array(sh, dtype=float)
            f = np.linalg.solve(Lo, x)
            f = f - np.floor(f + 1e-12)
            f[np.abs(f - 1.0) < 1e-9] = 0.0
            if not any(np.linalg.norm(((f - q) + 0.5) % 1.0 - 0.5) < 1e-6 for q in pts):
                pts.append(f)
                labs.append(lab)
    return Lo, np.array(pts), labs


def bcc():
    L = np.eye(3)
    return L, np.array([[0, 0, 0], [0.5, 0.5, 0.5]]), ["a", "b"]


def a15():
    labels, P = expand("P m -3 n", [("2a", (0, 0, 0)), ("6c", (0.25, 0, 0.5))])
    return np.eye(3), P, labels


def c15(x=0.125):
    # MgCu2, Fd-3m origin choice 2: Mg 8a (1/8,1/8,1/8), Cu 16d (1/2,1/2,1/2)
    labels, P = expand("F d -3 m:2", [("8a", (x, x, x)), ("16d", (0.5, 0.5, 0.5))])
    return np.eye(3), P, labels
