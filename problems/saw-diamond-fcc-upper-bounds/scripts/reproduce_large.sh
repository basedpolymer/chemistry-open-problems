#!/bin/bash
# Regenerate and check the large certificates (not stored in the repository):
#   diamond K = 22, 24 (value iteration, all candidates in memory), diamond K = 26 (policy iteration, candidates
#   cached on disk, proxy K = 24), fcc K = 10 (value iteration), fcc K = 11 (policy iteration, proxy K = 10, restricted candidates LATSAW_MONO=2:
#   paths through the remembered set move to an older vertex or to the next younger one).
# Environment: WSL Fedora 44, g++ 16 + OpenMP, 5-6 threads.  Resources and runtimes measured here (shared machine):
#   diamond K = 22: 2 min, 1.3 GB;  K = 24: see certificates/dia_k24.gen.log;  K = 26: see dia_k26.gen.log
#   fcc K = 10: 15 min, 4.2 GB;  fcc K = 11: see certificates/fcc_k11.gen.log (disk cache of the candidates)
# Set STEPS to a subset of "d22 d24 d26 f10 f11" to run only some of them (default: all but d26).
set -e
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BIN=~/work/saw-diamond-fcc/bin
OUT=~/work/saw-diamond-fcc/repro_large
mkdir -p "$OUT"; cd "$OUT"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-6}"
STEPS="${STEPS:-d22 d24 f10 f11}"
check() { "$BIN/latsaw_check" "$1.lcert" 4 > "$1.check.log" && echo "PASS $1: $(tail -1 "$1.check.log")" || { echo "FAIL $1"; exit 1; }; }
for s in $STEPS; do
  case $s in
    d22) "$BIN/latsaw_gen" diamond 22 dia_k22 opt > dia_k22.gen.log; check dia_k22 ;;
    d24) "$BIN/latsaw_gen" diamond 24 dia_k24 opt > dia_k24.gen.log; check dia_k24 ;;
    d26) [ -f dia_k24.lcert ] || { echo "run d24 first"; exit 1; }
         "$BIN/latsaw_gen" diamond 26 dia_k26 pi 12 dia_k24.lcert > dia_k26.gen.log; check dia_k26 ;;
    f10) "$BIN/latsaw_gen" fcc 10 fcc_k10 opt > fcc_k10.gen.log; check fcc_k10 ;;
    f11) [ -f fcc_k10.lcert ] || { echo "run f10 first"; exit 1; }
         LATSAW_MONO=2 "$BIN/latsaw_gen" fcc 11 fcc_k11 pi 12 fcc_k10.lcert > fcc_k11.gen.log; check fcc_k11 ;;
  esac
done
sha256sum *.lcert
echo "compare with $HERE/certificates/SHA256SUMS.txt (identical on the same CPU and compiler; the checker decides validity)"
