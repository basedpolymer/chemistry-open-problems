#!/bin/bash
# Runs the PASS/FAIL checks of Agent T (theory).  Environment: WSL Fedora, Python 3 with sympy, mpmath, scipy
# (HiGHS) and the read-only modules of computation/scripts (from the project root:
# outils/wsl "bash chemistry-open-problems/problems/halvorsen-skogestad-petlyuk-minimum-vapor/theory/scripts/run_theory_checks.sh").
# Runtime: about 3 minutes.  Exit code non-zero if any check fails.
set -e
cd "$(dirname "$0")"
python -u check_algebra.py                # exact identities used in the recovery-ordering and family-K lemmas and the attainment theorem (sympy)
python -u check_ftc_achievability.py      # attainment theorem: FTC point, N = 3..8, 1200 instances, 80 digits
python -u signcheck.py certify4           # N = 4: sign-rule-hard pairs = the 16 variants of family K at j = 2
echo "ALL THEORY CHECKS PASSED"
