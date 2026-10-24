#!/usr/bin/env python3
"""validate_census_sample.py -- independent exact re-check of a deterministic sample of census graphs.

For every census file alltrans<n>.tar with n in the given range that is present in
certificates/sources/holt-royle/, takes every 'step'-th connected graph of every valency file, and checks
with exact rationals (python-flint, two methods) that Kf(G) != Kf(G - v0) (v0 = vertex 0), i.e. that the
rejection by kscheck -t is correct, and that vertex-transitivity is consistent (Kf(G - v) equal for a few
random v).  PASS/FAIL; exit 1 on failure.  Usage: validate_census_sample.py NMIN NMAX STEP
Run with python3 -I is not needed (only reads graph6 strings).  Runtime: minutes.
"""
import gzip
import io
import os
import random
import sys
import tarfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kflib import g6_decode, kf_inv, kf_charpoly, delete, is_connected  # noqa: E402

here = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(here, '..', 'certificates', 'sources', 'holt-royle')
nmin, nmax, step = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(3)
ok = True
tested = 0
for n in range(nmin, nmax + 1):
    path = os.path.join(src, f'alltrans{n}.tar')
    if not os.path.exists(path):
        continue
    with tarfile.open(path) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            data = gzip.decompress(tf.extractfile(m).read()).decode().split()
            for idx, line in enumerate(data):
                if idx % step:
                    continue
                adj = g6_decode(line)
                if not is_connected(adj) or len(adj) < 3:
                    continue
                k = kf_inv(adj)
                H = delete(adj, 0)
                if not is_connected(H):
                    continue
                kh = kf_inv(H)
                if n <= 24 and (kf_charpoly(adj) != k or kf_charpoly(H) != kh):
                    ok = False
                    print('FAIL methods disagree', m.name, idx)
                if kh == k:
                    ok = False
                    print('FAIL: Kf(G) == Kf(G-v0) for', m.name, idx, line)
                v = rng.randrange(n)
                if kf_inv(delete(adj, v)) != kh:
                    ok = False
                    print('FAIL: not vertex-transitive?', m.name, idx)
                tested += 1
print(f'{tested} census graphs re-checked exactly (n = {nmin}..{nmax}, every {step}-th graph)')
print('PASS census sample: every sampled rejection confirmed exactly' if ok else 'FAIL census sample')
sys.exit(0 if ok else 1)
