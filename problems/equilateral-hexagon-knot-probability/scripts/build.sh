#!/bin/bash
# Builds the certifier and the trigonometric table.  Environment: g++ >= 12 (C++20), Python 3 with mpmath.
# Usage: bash scripts/build.sh [OUTDIR]   (default OUTDIR = build/ next to scripts/)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-$HERE/../build}"
mkdir -p "$OUT"
# production variant: rounding mode FE_UPWARD (requires -frounding-math); no fast-math, no FMA contraction
g++ -O2 -std=c++20 -ffp-contract=off -fno-fast-math -frounding-math -DHEX_UPWARD -Wall \
    -o "$OUT/hexcert" "$HERE/../src/hexcert.cpp"
# reference variant: nextafter outward rounding in round-to-nearest mode (slower, same enclosures up to 1 ulp)
g++ -O2 -std=c++20 -ffp-contract=off -fno-fast-math -Wall \
    -o "$OUT/hexcert_nextafter" "$HERE/../src/hexcert.cpp"
python3 "$HERE/make_trig_table.py" 14 "$OUT/trig14.bin"
echo "built $OUT/hexcert $OUT/hexcert_nextafter $OUT/trig14.bin"
