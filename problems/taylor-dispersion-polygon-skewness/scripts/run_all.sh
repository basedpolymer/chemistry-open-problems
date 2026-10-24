#!/bin/sh
# Re-run every certificate (one core, about 25 minutes in total). Each script prints PASS/FAIL and exits non-zero on failure.
# Environment: python3.13, python-flint 0.9.0, sympy 1.14, numpy 2.x.
set -e
cd "$(dirname "$0")"
python3 triangle_exact.py                 # n = 3 exact (seconds)
python3 lt_certify.py 4 40                # long-time sign, finite part (~1 min)
python3 lt_uniform.py 30                  # long-time sign, all n >= 30 (seconds)
python3 sg_certify.py 4 150               # short-time sign, finite part (~15 min)
python3 sg_uniform.py 140 6               # short-time sign, all n >= 140 (seconds)
python3 alltime_windows.py                # certified time windows (~1 min)
python3 table1.py                         # enclosures of Table 1 of the paper (~3 min)
python3 numerics_check.py                 # non-rigorous sanity checks (not a certificate)
