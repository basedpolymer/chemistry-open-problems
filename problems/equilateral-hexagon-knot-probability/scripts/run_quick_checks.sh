#!/bin/bash
# Runs the fast checks (a few minutes, one thread) and writes logs to certificates/.
# Usage: bash scripts/run_quick_checks.sh [BUILDDIR]
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
B="${1:-$HERE/../build}"
C="$HERE/../certificates"
[ -x "$B/hexcert" ] || bash "$HERE/build.sh" "$B"
fail=0
python3 "$HERE/test_interval.py" 200000      2>&1 | tee "$C/test_interval.log"   || fail=1
python3 "$HERE/check_geometry.py"            2>&1 | tee "$C/check_geometry.log"  || fail=1
python3 "$HERE/hake_lemma_check.py"          2>&1 | tee "$C/hake_lemma_check.log" || fail=1
python3 "$HERE/check_boxes_mc.py" "$B/hexcert" "$B/trig14.bin" 400000 2>&1 | tee "$C/check_boxes_mc.log" || fail=1
if grep -q FAIL "$C"/test_interval.log "$C"/check_geometry.log "$C"/hake_lemma_check.log "$C"/check_boxes_mc.log; then fail=1; fi
[ $fail -eq 0 ] && echo "ALL QUICK CHECKS PASS" || { echo "SOME CHECK FAILED"; exit 1; }
