#!/bin/bash
# Run a list of cases of run_case.py in parallel (at most P processes, one thread each); logs in ../runs/<case>.log.
# Usage (from this folder, under WSL): bash run_batch.sh P case1 case2 ...
P=$1; shift
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
printf '%s\n' "$@" | xargs -P "$P" -I{} bash -c 'python run_case.py {} > ../runs/{}.log 2>&1; echo "{}: exit $?"'
