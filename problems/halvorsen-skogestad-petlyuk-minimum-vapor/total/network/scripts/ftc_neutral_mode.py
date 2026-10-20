#!/usr/bin/env python3
"""Smallest singular value of the steady-state Jacobian of the finite-stage FTC at fixed flows (Agent N, network.md
§10.4): evidence that sections carrying excess vapor between two pinch zones have an exponentially weakly determined
"front position" (a neutral mode of the infinite-stage limit system).  Floating point; illustration only.
Usage: python ftc_neutral_mode.py N delta n1,n2,... [seed].  Environment: Python 3, numpy, scipy, mpmath.
"""
import sys
import numpy as np
from ftc_fixed_flows import build_ftc

N = int(sys.argv[1]); delta = float(sys.argv[2]); ns = [int(v) for v in sys.argv[3].split(",")]
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
if seed == 0 and N == 3:
    al, z, q = [4.0, 2.0, 1.0], [1 / 3, 1 / 3, 1 / 3], 1.0
else:
    rng = np.random.default_rng(seed)
    al = list(np.cumprod(np.r_[1.0, rng.uniform(1.4, 3.0, N - 1)])[::-1])
    z = list(rng.dirichlet(np.ones(N) * 3)); q = float(rng.uniform(0.0, 1.0))
for n in ns:
    net, P, H, C = build_ftc(al, z, q, n, delta)
    X, res = net.solve(maxit=100)
    J = net.jacobian(X.ravel()).toarray()
    # scale rows by stage throughput so that the singular values are comparable across n
    u, s, vt = np.linalg.svd(J)
    v = vt[-1].reshape(net.nst, N)
    w = np.abs(v).sum(axis=1)
    top = [net.names[t] for t in np.argsort(-w)[:3]]
    comp = np.round(v[np.argmax(w)], 3)
    print(f"n = {n:3d}: residual {res:.1e}; smallest singular values {s[-1]:.2e} {s[-2]:.2e}; "
          f"null vector concentrated at {top}, components {comp}", flush=True)
