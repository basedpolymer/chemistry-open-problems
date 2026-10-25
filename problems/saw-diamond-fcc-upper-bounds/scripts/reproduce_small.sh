#!/bin/bash
# Regenerate the small certificates (diamond K = 6..20, fcc K = 3..9; optimised children, value iteration),
# check each with latsaw_check and print the certified bounds.  Bitwise identity with the stored files is expected
# on the same CPU/compiler only (the vector is floating point before rounding); the checker decides validity.
# Environment: WSL Fedora 44, g++ 16 + OpenMP, OMP_NUM_THREADS threads (default 6).  Runtime: about 3 min, < 1 GB.
set -e
BIN=~/work/saw-diamond-fcc/bin
OUT=~/work/saw-diamond-fcc/repro_small
mkdir -p "$OUT"; cd "$OUT"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-6}"
for K in 6 8 10 12 14 16 18 20; do
  "$BIN/latsaw_gen" diamond $K dia_k$K opt > dia_k$K.gen.log
  "$BIN/latsaw_check" dia_k$K.lcert 4 > dia_k$K.check.log && echo "PASS diamond K=$K: $(tail -1 dia_k$K.check.log)" || { echo "FAIL diamond K=$K"; exit 1; }
done
for K in 3 4 5 6 7 8 9; do
  "$BIN/latsaw_gen" fcc $K fcc_k$K opt > fcc_k$K.gen.log
  "$BIN/latsaw_check" fcc_k$K.lcert 4 > fcc_k$K.check.log && echo "PASS fcc K=$K: $(tail -1 fcc_k$K.check.log)" || { echo "FAIL fcc K=$K"; exit 1; }
done
sha256sum *.lcert
