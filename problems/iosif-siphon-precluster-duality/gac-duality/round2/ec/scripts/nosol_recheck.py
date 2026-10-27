"""Re-check (floating point) of the 'NOSOL' lines of logs/dual_search_n4_r4.log: for each (N, kappa, k0) where the
multistart Newton solver of dual_search.py found no positive steady state of N*, integrate N* from k0 (LSODA,
t <= 1e4) and report the residual |Gamma^T omega| at the end point and its minimal coordinate.  Theorem T*' and U*
predict convergence to a positive steady state.  Agent Q3, Claude Opus 5.5 (effort xhigh).  Environment: Python 3.12,
numpy, scipy (WSL).  Usage (from round2/ec/): python scripts/nosol_recheck.py logs/dual_search_n4_r4.log
"""
import re
import sys
import numpy as np
from scipy.integrate import solve_ivp
from eclib import Network


def parse_net(desc, n):
    names = [f"X{i+1}" for i in range(n)]
    rx = []
    for part in desc.split(", "):
        lhs, rhs = part.split(" -> ")
        vecs = []
        for side in (lhs, rhs):
            v = [0] * n
            for term in side.split(" + "):
                m = re.fullmatch(r"(\d*)X(\d+)", term.strip())
                v[int(m.group(2)) - 1] += int(m.group(1) or 1)
            vecs.append(tuple(v))
        rx.append(tuple(vecs))
    return Network(n, rx)


worst = 0.0
cnt = 0
for line in open(sys.argv[1], encoding="utf-8"):
    if not line.startswith("NOSOL"):
        continue
    m = re.match(r"NOSOL  (.*)  kappa=\[(.*)\] k0=\[(.*)\]", line.strip())
    desc, ks, k0s = m.groups()
    num = lambda s: [float(x) for x in re.findall(r"[-+]?\d*\.\d+(?:e[-+]?\d+)?|\d+", s.replace("np.float64", ""))]
    kap, k0 = np.array(num(ks)), np.array(num(k0s))
    N = parse_net(desc, len(kap))
    G = np.array(N.Gamma, float)
    Ye = np.array(N.Ye, float)
    F = lambda t, q: G.T @ (kap * np.exp(Ye @ np.log(np.maximum(q, 1e-300))))
    sol = solve_ivp(F, (0, 1e4), k0, method="LSODA", rtol=1e-11, atol=1e-14)
    kend = sol.y[:, -1]
    res = np.max(np.abs(F(0, kend))) / max(1.0, np.max(np.abs(kap)))
    worst = max(worst, res)
    cnt += 1
    print(f"{desc}: end point min={kend.min():.3e} residual={res:.2e} class drift={np.max(np.abs(np.linalg.svd(G)[2][N.rank():] @ (kend - k0))):.1e}")
print(f"# {cnt} cases re-checked; worst residual {worst:.2e}")
