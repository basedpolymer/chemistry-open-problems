#!/bin/bash
# Regenerate the d = 3, K = 16 certificate and check it.
# Environment: WSL Fedora 44, g++ 16 + OpenMP, 6 threads.  Runtime: ~7 min (generation) + ~10 min (check);
# memory 3.5 GB.  Run scripts/build.sh first.  Output in ~/work/cubic-saw/repro.
set -e
BIN=~/work/cubic-saw/bin
OUT=~/work/cubic-saw/repro
mkdir -p "$OUT"; cd "$OUT"
export OMP_NUM_THREADS=6
"$BIN/saw_prod" 3 16 d3_k16 5000 4 | tee d3_k16.gen.log | grep -E "RESULT|integer"
sha256sum d3_k16.cert
"$BIN/verify_cert" d3_k16.cert 6 | tee d3_k16.verify.log
