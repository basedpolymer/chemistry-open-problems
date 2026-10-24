#!/bin/bash
# Packages the main branch-and-bound state into certificates/ (leaves, summary, run log, checksums).
# Usage: bash scripts/package_certificate.sh STATE.pkl RUN.log
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
C="$HERE/../certificates"
python3 "$HERE/export_certificate.py" "$1" "$C/hexknot"
cp "$2" "$C/bnb_run.log"
( cd "$C" && sha256sum hexknot.leaves.gz hexknot.summary.json bnb_run.log *.log 2>/dev/null | sort -u -k2 > SHA256SUMS.txt )
( cd "$HERE/../src" && sha256sum hexcert.cpp interval.hpp test_interval.cpp >> "$C/SHA256SUMS.txt" )
ls -la "$C"
