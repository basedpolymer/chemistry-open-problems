#!/bin/bash
# run_all.sh -- recompute everything for the Boros-Hofbauer four-limit-cycle certificate.
# Environment: WSL Fedora, Python 3.12 with python-flint 0.9.0, mpmath 1.3.0, numpy, scipy.
# Usage (from this problem's directory):  bash scripts/run_all.sh [--with-illustration]
# Runtime: ~15 s (plus ~3 min with --with-illustration).  Exit code != 0 if any certificate check fails.
set -u
cd "$(dirname "$0")"
status=0
run() {
  echo "=================== $* ==================="
  python "$@" || status=1
}
tmp=$(mktemp -d)
run focal_exact.py 4 "$tmp"          # exact focal values V~_1..V~_4 (written to a temp dir)
run crosscheck.py                      # validation against the paper and the float implementation
run certify.py                         # interval certificate (Krawczyk, L4 < 0, Jacobian)
run uniqueness.py                      # Theorem B: the zero is unique in the whole domain
run validated_cycles.py                # Theorem C: validated integration at explicit rational parameters
run polar_check.py                     # independent float method (non-rigorous cross-check)
if [ "${1:-}" = "--with-illustration" ]; then
  run illustrate_cycles.py             # numerical illustration of four explicit cycles (non-rigorous)
fi
rm -rf "$tmp"
echo
if [ $status -eq 0 ]; then echo "RUN_ALL: ALL PASS"; else echo "RUN_ALL: SOME CHECK FAILED"; fi
exit $status
