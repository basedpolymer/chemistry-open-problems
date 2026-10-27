import sys, ast, itertools, numpy as np
from scipy.optimize import nnls
n = 5
COMP = [np.array(c) for c in itertools.product(range(3), repeat=n) if sum(c) <= 2]
REAC = [np.array(r) for r in [(0, 0, 0, 0, 1), (0, 0, 0, 1, 1), (0, 0, 2, 0, 0), (1, 1, 0, 0, 0)]]
NM = "ABCDE"
def cname(c):
    s = []
    for i, a in enumerate(c):
        if a == 1: s.append(NM[i])
        elif a == 2: s.append("2" + NM[i])
    return "+".join(s) if s else "0"
for line in open(sys.argv[1]):
    if not line.startswith("HIT"): continue
    body = line[4:]
    i = body.index("]] [")
    N = np.array(ast.literal_eval(body[:i + 2]))
    tot = 0; desc = []
    for a in range(4):
        col = N[:, a]
        gens = [p - REAC[a] for p in COMP if not np.array_equal(p, REAC[a])]
        names = [cname(p) for p in COMP if not np.array_equal(p, REAC[a])]
        found = None
        for k in range(1, 5):
            for S in itertools.combinations(range(len(gens)), k):
                G = np.array([gens[j] for j in S]).T
                w, r = nnls(G, col)
                if r < 1e-6 * (1 + np.linalg.norm(col)) and np.all(w > 1e-9):
                    found = (k, [(names[j], round(float(wj), 4)) for j, wj in zip(S, w)])
                    break
            if found: break
        if not found:
            found = (99, None)
        tot += found[0]
        desc.append((cname(REAC[a]), found[1]))
    print(tot, desc)
