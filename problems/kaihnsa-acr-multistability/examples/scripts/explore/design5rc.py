#!/usr/bin/env python3
"""
design5rc.py -- numerical exploration: all integral-feedback designs with 6 reactions and EXACTLY the reactant
complexes {X+Y, Y, c0, 2Z, Y+Z}, c0 in {0, X, 2X} (configuration #4, the only one with degree-3 reduced circuit
that the design can realize), heavy kappa sampling (vectorized) for two exponentially stable steady states.
Usage: python design5rc.py samples seed
"""
import numpy as np, itertools, sys, time

NSAMP = int(sys.argv[1])
SEED = int(sys.argv[2])
sys.argv = ["x", "6", "5", "10", "10", "1", "0", "1"]
_seed = int(sys.argv[2]) if False else None
exec(open("design_fast.py").read().split("t0 = time.time()")[0])
rng = np.random.default_rng(SEED)
XY, YY, Yc, ZZ, YZ = (1, 1, 0), (0, 2, 0), (0, 1, 0), (0, 0, 2), (0, 1, 1)
free = [(0, 0, 0), (1, 0, 0), (0, 0, 1), (2, 0, 0), (0, 0, 2), (1, 0, 1)]
t0 = time.time()
tested = 0
hits = 0
for c0 in [(0, 0, 0), (1, 0, 0), (2, 0, 0)]:
    opts = {YZ: [(YZ, Yc), (YZ, XY)], ZZ: [(ZZ, p) for p in free if p != ZZ], c0: [(c0, p) for p in free if p != c0]}
    extra_pool = sorted(set(
        [(XY, Yc), (XY, YZ), (XY, (1, 0, 0)), (XY, (0, 0, 1)), (XY, (1, 0, 1)), (XY, (0, 0, 0)), (XY, (2, 0, 0)),
         (XY, (0, 0, 2))]
        + [(Yc, XY), (Yc, YZ), (Yc, YY)] + [(Yc, w) for w in free]
        + opts[YZ] + opts[ZZ] + opts[c0]))
    for W in free:
        base = [(XY, YY), (Yc, W)]
        for r3 in opts[YZ]:
            for r4 in opts[ZZ]:
                for r5 in opts[c0]:
                    for r6 in extra_pool:
                        rx = base + [r3, r4, r5, r6]
                        if len(set(rx)) < 6:
                            continue
                        if set(a for a, _ in rx) != {XY, Yc, c0, ZZ, YZ}:
                            continue
                        tested += 1
                        K = np.exp(rng.uniform(-5, 5, size=(NSAMP, 6)))
                        h, info = batch(rx, K)
                        if h:
                            hits += 1
                            comps = set([a for a, _ in rx] + [b for _, b in rx])
                            print("HIT complexes=%d" % len(comps), rx, "k=", np.round(info[0], 5).tolist(),
                                  "ss=", np.round(info[1], 5).tolist(), flush=True)
print("tested", tested, "hits", hits, "time %.0f" % (time.time() - t0))
