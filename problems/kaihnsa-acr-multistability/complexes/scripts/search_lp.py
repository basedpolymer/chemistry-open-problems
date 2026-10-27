"""Exploration (floating point LP, not a certificate; prints candidate networks only, no PASS/FAIL).
Runtime: about 7 minutes for 20 trials per design in one process; the runs used for report.md were
  py -3 search_lp.py 100 <10..15> 8 6 <0..5>      (three points, six parallel parts, about 25 minutes each)
  py -3 search_lp.py 60 <40..45> 6 6 <0..5> 2     (two points, six parallel parts, about 20 minutes each)

Search for 6-complex bimolecular networks on species X, Y, Z with *structural* ACR in X of the following type:
an SF pair p, q = p + e_X of complexes, a weight vector w with w_X != 0, and every reaction whose source is not p or q
keeps the weight w.y unchanged. Then w.f(x) = x^p (A + B x_X) with A, B linear in the rate constants, so every
positive steady state has x_X = -A/B.

Inverse approach: fix x_X = c and three rational points on the plane, require f = 0 there (linear in kappa),
take LP vertices with random objectives (sparse networks), keep those with two exponentially stable points.

Usage: py -3 search_lp.py [n_trials_per_design] [seed] [max_reactions] [n_parts] [part_index] [n_points (2 or 3)]
"""
import sys, itertools, random
import numpy as np
from scipy.optimize import linprog

SP = "XYZ"


def cname(c):
    parts = []
    for k in range(3):
        if c[k] == 1:
            parts.append(SP[k])
        elif c[k] == 2:
            parts.append("2" + SP[k])
    return "+".join(parts) if parts else "0"


ALL = [c for c in itertools.product(range(3), repeat=3) if sum(c) <= 2]
EX = np.array([1, 0, 0])


def affdim(cs):
    a = np.array(cs, dtype=float)
    return np.linalg.matrix_rank(a[1:] - a[0]) if len(cs) > 1 else 0


def designs():
    seen = set()
    for C in itertools.combinations(ALL, 6):
        if affdim(C) != 3:
            continue
        for p in C:
            q = tuple(np.array(p) + EX)
            if q not in C:
                continue
            others = [c for c in C if c not in (p, q)]
            for w in itertools.product(range(1, 4), range(-3, 4), range(-3, 4)):
                w = np.array(w)
                edges = []
                for s in C:
                    for t in C:
                        if s == t:
                            continue
                        if s in (p, q) or np.dot(w, np.array(t) - np.array(s)) == 0:
                            edges.append((s, t))
                # p must have an edge raising the weight or q one lowering it (else no positive x)
                key = (C, p, tuple(sorted(edges)))
                if key in seen:
                    continue
                seen.add(key)
                # every non-pair complex needs at least one allowed outgoing edge or may be non-reactant
                yield C, p, q, w, edges


def jac(V, Cr, P):
    m = np.prod(P[:, None] ** Cr, axis=0)
    return V @ np.diag(m) @ Cr.T @ np.diag(1.0 / P)


VALS = [0.25, 1 / 3, 0.5, 2 / 3, 1.0, 1.5, 2.0, 3.0, 4.0]


def trial(C, p, q, w, edges, rng, maxr, npts=3):
    Carr = {c: np.array(c) for c in C}
    cval = rng.choice([0.5, 1.0, 2.0])
    pts = []
    while len(pts) < npts:
        y, z = rng.choice(VALS), rng.choice(VALS)
        if (y, z) not in pts:
            pts.append((y, z))
    P = [np.array([cval, y, z]) for (y, z) in pts]
    E = len(edges)
    M = np.zeros((3 * npts, E))
    for j, Pj in enumerate(P):
        for e, (s, t) in enumerate(edges):
            M[3 * j:3 * j + 3, e] = np.prod(Pj ** Carr[s]) * (Carr[t] - Carr[s])
    A_eq = np.vstack([M, np.ones((1, E))])
    b_eq = np.zeros(3 * npts + 1)
    b_eq[-1] = 1.0
    obj = np.array([rng.expovariate(1.0) for _ in range(E)])
    res = linprog(obj, A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * E, method="highs")
    if res.status != 0:
        return None
    k = res.x
    supp = [e for e in range(E) if k[e] > 1e-9]
    if len(supp) > maxr:
        return None
    used = set()
    for e in supp:
        used.add(edges[e][0]); used.add(edges[e][1])
    reactants = sorted(set(edges[e][0] for e in supp))
    # check B != 0 (ACR relation nondegenerate)
    B = sum(k[e] * np.dot(w, Carr[edges[e][1]] - Carr[edges[e][0]]) for e in supp if edges[e][0] == q)
    A = sum(k[e] * np.dot(w, Carr[edges[e][1]] - Carr[edges[e][0]]) for e in supp if edges[e][0] == p)
    if abs(B) < 1e-9:
        return None
    Cr = np.array(reactants).T.astype(float)
    V = np.zeros((3, len(reactants)))
    for e in supp:
        s, t = edges[e]
        V[:, reactants.index(s)] += k[e] * (Carr[t] - Carr[s])
    nst = 0
    stab = []
    for Pj in P:
        J = jac(V, Cr, Pj)
        ev = np.linalg.eigvals(J)
        nJ = np.linalg.norm(J, 2)
        ok = max(ev.real) < -1e-6 * nJ and min(abs(ev)) > 1e-6 * nJ
        stab.append(ok)
        nst += ok
    if nst < 2:
        return None
    return dict(C=C, used=sorted(used), w=w, cval=cval, pts=pts, supp=[(edges[e], k[e]) for e in supp], stab=stab,
                A=A, B=B)


def main():
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    maxr = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    npart = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    ipart = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    npts = int(sys.argv[6]) if len(sys.argv) > 6 else 3
    rng = random.Random(seed)
    nd = 0
    hits = []
    for C, p, q, w, edges in designs():
        nd += 1
        if nd % npart != ipart:
            continue
        for _ in range(ntr):
            h = trial(C, p, q, w, edges, rng, maxr, npts)
            if h is not None and len(h["used"]) == 6:
                hits.append(h)
                print("HIT", len(h["supp"]), "reactions; complexes", [cname(c) for c in h["used"]], "w", h["w"],
                      "c", h["cval"], "pts", h["pts"], "stab", h["stab"])
                for (s, t), kk in h["supp"]:
                    print("    %s -> %s  %.6g" % (cname(s), cname(t), kk))
                sys.stdout.flush()
    print("designs", nd, "hits", len(hits))


if __name__ == "__main__":
    main()
