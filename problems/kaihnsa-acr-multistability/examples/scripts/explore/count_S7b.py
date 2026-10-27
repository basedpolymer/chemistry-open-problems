import numpy as np, json, sys
from fractions import Fraction
from scipy.optimize import fsolve
rng = np.random.default_rng(int(sys.argv[1])); NS = int(sys.argv[2])
Vt = np.array([[-1, 1, 0, 0, 0], [1, 0, 0.5, 0, 1]], float).T
valid = json.load(open("valid_S7.json"))
sys.argv = ["x", "4", "1", "0", "0", "1"]
exec(open("enum_S7.py").read().split("prods = [[p")[0])
hist = {}
for v in valid:
    W = np.array([[float(Fraction(x)) for x in row] for row in v['W']]).T
    # skip degenerate: resultant identically zero
    xs = np.exp(rng.uniform(-2, 2, 5))
    if class_solutions(W, xs) is None:
        continue
    mx = 0
    for s in range(NS):
        xs = np.exp(rng.uniform(-5, 5, 5))
        F = lambda t: W.T @ (xs * np.exp(Vt @ t) - xs)
        sols = [np.zeros(2)]
        for st in rng.uniform(-8, 8, size=(30, 2)):
            t, info, ier, msg = fsolve(F, st, full_output=True)
            if ier == 1 and np.all(np.isfinite(t)) and np.max(np.abs(F(t))) < 1e-10 * (1 + np.max(xs)):
                if all(np.max(np.abs(t - q)) > 1e-6 for q in sols):
                    sols.append(t)
        mx = max(mx, len(sols))
    hist[mx] = hist.get(mx, 0) + 1
    if mx >= 3:
        print(v['W'], mx, flush=True)
print(hist)
