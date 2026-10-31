#!/bin/sh
# Re-run every certificate of all-times/ (about 1.5 hours on 3 cores).  Each script prints PASS/FAIL and exits non-zero on failure.
# Environment: Python 3.12, python-flint 0.9.0, sympy 1.14, numpy, scipy.
set -e
cd "$(dirname "$0")"
mkdir -p ../logs ../certificates
python3 run_finite.py 4 200 0.015 3   | tee ../logs/finite_4_200.log     # n = 4..200 (Prop. finite)
python3 uniform.py 200 12             | tee ../logs/uniform_200.log      # all n >= 200
python3 uniform.py 0 12               | tee ../logs/uniform_disk.log     # the disk
python3 u3_uniform.py 40              | tee ../logs/u3_uniform.log       # <u^3> > 0, n >= 40
python3 triangle.py 0.025 1.0 1.6      | tee ../logs/triangle.log         # n = 3: one zero of K
python3 ../../scripts/sg_certify.py 4 150                                # <u^3> > 0, 4 <= n <= 150 (writes ../certificates of result 14)
python3 ../../scripts/sg_uniform.py 140 6                                # <u^3> > 0, n >= 140, second proof
python3 numerics_check.py             | tee ../logs/numerics_check.log   # not a certificate
