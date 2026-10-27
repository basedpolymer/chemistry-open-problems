#!/usr/bin/env python3
"""
six_complex_design.py -- numerical exploration: integral-feedback designs (ACR in X via dY/dt = y(alpha x - beta))
whose complexes are exactly {X+Y, 2Y, Y, c0, 2Z, Y+Z}, c0 in {0, X, 2X} (the only 6-complex sets for which the
reduced system on x = x* can have 3 positive solutions).  All subsets of the admissible reactions are tested with
random kappa.  Usage: python six_complex_design.py samples seed
"""
import numpy as np, itertools, sys
sys.path.insert(0, '.')
SAMP = int(sys.argv[1])
rng = np.random.default_rng(int(sys.argv[2]))
XY, YY, Yc, ZZ, YZ = (1, 1, 0), (0, 2, 0), (0, 1, 0), (0, 0, 2), (0, 1, 1)
from design_fast import batch  # noqa: E402  (design_fast parses argv at import: avoid)
