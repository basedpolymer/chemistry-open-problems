#!/usr/bin/env bash
# Invoke from the attainment-astra task folder in the scientific environment.
# Deterministic checks, <=2 threads and <=3 GiB virtual memory.
set -eu
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1
ulimit -v 3145728
mkdir -p runs
python scripts/check_exact.py > runs/check_exact.log 2>&1
cat runs/check_exact.log
python scripts/check_prefractionator.py > runs/check_prefractionator.log 2>&1
cat runs/check_prefractionator.log
python scripts/check_ternary_returns.py > runs/check_ternary_returns.log 2>&1
cat runs/check_ternary_returns.log
python scripts/check_ftc_numerical.py > runs/check_ftc_numerical.log 2>&1
cat runs/check_ftc_numerical.log
echo 'PASS all checks (numerical checks retain their stated non-certificate status)'
