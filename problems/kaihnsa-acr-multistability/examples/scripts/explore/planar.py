"""planar.py -- exploration: 2 species (Y,Z), bimolecular, reactant complexes {0, Y, 2Z, Y+Z} (any number of
reactions): ydot = a0 + a1 y + a2 z^2 + a3 y z, zdot = b0 + b1 y + b2 z^2 + b3 y z, columns in the bimolecular
cones.  Looks for two exponentially stable positive steady states."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1])); NS = int(sys.argv[2])
gen = {  # generators (dY, dZ) per reactant complex
    '0': [(1, 0), (0, 1), (2, 0), (0, 2), (1, 1)],
    'Y': [(-1, 0), (1, 0), (-1, 1), (0, 1), (-1, 2)],
    '2Z': [(0, -2), (1, -2), (0, -1), (2, -2), (1, -1)],
    'YZ': [(-1, -1), (0, -1), (-1, 0), (1, -1), (-1, 1)],
}
hits = 0; three = 0
for it in range(NS):
    cols = {}
    for k, g in gen.items():
        m = rng.integers(1, 3)
        idx = rng.choice(len(g), size=m, replace=False)
        w = np.exp(rng.uniform(-4, 4, m))
        cols[k] = sum(w[i] * np.array(g[idx[i]], float) for i in range(m))
    a0, b0 = cols['0']; a1, b1 = cols['Y']; a2, b2 = cols['2Z']; a3, b3 = cols['YZ']
    # eliminate y: y (a1 + a3 z) = -(a0 + a2 z^2); y (b1 + b3 z) = -(b0 + b2 z^2)
    # (a0 + a2 z^2)(b1 + b3 z) - (b0 + b2 z^2)(a1 + a3 z) = 0
    P = np.polynomial.polynomial
    R = P.polysub(P.polymul([a0, 0, a2], [b1, b3]), P.polymul([b0, 0, b2], [a1, a3]))
    zs = [r.real for r in P.polyroots(R) if abs(r.imag) < 1e-10 and r.real > 1e-12]
    sts = []
    for z in zs:
        den = a1 + a3 * z
        if abs(den) < 1e-14: continue
        y = -(a0 + a2 * z * z) / den
        if y <= 1e-12: continue
        J = np.array([[a1 + a3 * z, 2 * a2 * z + a3 * y], [b1 + b3 * z, 2 * b2 * z + b3 * y]])
        tr, det = np.trace(J), np.linalg.det(J)
        sts.append((y, z, tr < 0 and det > 0))
    if len(sts) >= 3:
        three += 1
    if sum(s for _, _, s in sts) >= 2:
        hits += 1
        if hits <= 5:
            print("HIT", {k: np.round(v, 5).tolist() for k, v in cols.items()}, sts, flush=True)
print("samples", NS, "three", three, "hits", hits)
