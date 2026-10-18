#!/bin/bash
# Regenerate the d = 3, K = 18 certificate (needs the K = 16 certificate from reproduce_k16.sh as proxy for the
# initial policy) and check it.
# Environment: WSL Fedora 44, g++ 16 + OpenMP, 5 threads, ~10.5 GB RAM, ~40 GB free disk (temporary candidate
# cache).  Runtime: several hours (first improvement pass), see certificates/d3_k18.gen.log.
set -e
BIN=~/work/cubic-saw/bin
OUT=~/work/cubic-saw/repro
mkdir -p "$OUT"; cd "$OUT"
[ -f d3_k16.cert ] || { echo "run scripts/reproduce_k16.sh first"; exit 1; }
export OMP_NUM_THREADS=5
"$BIN/saw_big" 3 18 d3_k18 5 d3_k16.cert | tee d3_k18.gen.log | grep -E "pass .: (impr|eval)|RESULT|integer"
sha256sum d3_k18.cert d3_k18.cert2
"$BIN/verify_cert" d3_k18.cert2 5 | tee d3_k18.verify.log
