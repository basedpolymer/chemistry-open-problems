#!/bin/bash
# Regenerate and check the record certificates in dimensions 4, 5, 6 (d=4 K=14, d=5 K=14, d=6 K=12) and the
# Ponitz-Tittmann reproductions (policy 0).  Environment: WSL Fedora 44, g++ 16.  Runtime ~45 min with 1 thread
# (d=5 K=14 is the longest, ~25 min, 3.6 GB).  Run scripts/build.sh first.
set -e
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BIN=~/work/cubic-saw/bin
OUT=~/work/cubic-saw/repro_hd
mkdir -p "$OUT"; cd "$OUT"
export OMP_NUM_THREADS=${THREADS:-1}
for job in "4 10 0" "4 12 0" "5 10 0" "6 10 0" "4 14 4" "5 14 4" "6 12 4"; do
  set -- $job
  "$BIN/saw_prod" $1 $2 d$1_k$2_p$3 5000 $3 > d$1_k$2_p$3.gen.log
  grep RESULT d$1_k$2_p$3.gen.log
  python3 -I "$HERE/scripts/check_counts.py" d$1_k$2_p$3.gen.log "$HERE/certificates/counts/saw_d$1.txt"
  if [ $3 = 4 ]; then "$BIN/verify_cert" d$1_k$2_p$3.cert "$OMP_NUM_THREADS" | tail -1; fi
done
