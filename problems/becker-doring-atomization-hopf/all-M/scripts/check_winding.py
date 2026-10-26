#!/usr/bin/env python3
"""Adversarial checks of the rigorous winding counter (Python/flint, <1 s).

Exact polynomial examples exercise multiplicities, reversed orientation,
large argument changes and rejection of a contour passing through a zero.
"""
from unittest.mock import patch
from flint import arb, ctx
import allm_core as C

ctx.prec = 192
square = [complex(-1, -1), complex(1, -1), complex(1, 1), complex(-1, 1), complex(-1, -1)]
cases = [
    ("multiplicity four", lambda z: z**4, square, 4),
    ("reverse orientation", lambda z: z**4, square[::-1], -4),
    ("one inside, one outside", lambda z: (z-arb("0.25"))*(z-2), square, 1),
    ("no zeros", lambda z: z-2, square, 0),
    ("zero on contour", lambda z: z-1, square, None),
    ("unclosed contour", lambda z: z, square[:-1], None),
]
for label, poly, path, expected in cases:
    with patch.object(C, "QL", lambda k, z, e: (C.DN(poly(z)), C.DN(0))):
        count, _ = C.winding_count(path, arb(0), arb(0), max_depth=12,
                                   Qfun=lambda z, k, e: poly(z))
    if count != expected:
        raise AssertionError((label, count, expected))
    print("PASS", label, count)
print("SUMMARY PASS")
