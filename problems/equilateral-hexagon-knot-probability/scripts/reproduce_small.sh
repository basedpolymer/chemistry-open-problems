#!/bin/bash
# Small end-to-end reproduction (one thread, ~15-30 min): build, run a short branch and bound (5 rounds),
# export and fully re-verify the resulting (weaker) certificate.  Prints PASS/FAIL lines.
# Usage: bash scripts/reproduce_small.sh [WORKDIR]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
W="${1:-$HERE/../build/small}"
mkdir -p "$W"
bash "$HERE/build.sh" "$W"
rm -f "$W/state.pkl"
python3 "$HERE/bnb_driver.py" --state "$W/state.pkl" --bin "$W/hexcert" --table "$W/trig14.bin" \
    --batch 1000 --init-level 3 --offset 2 --lth-min 6 --lth-max 10 --rounds 5 | tee "$W/driver.log"
python3 "$HERE/export_certificate.py" "$W/state.pkl" "$W/small"
python3 "$HERE/verify_certificate.py" "$W/small.leaves.gz" "$W/small.summary.json" "$W/hexcert" "$W/trig14.bin"
