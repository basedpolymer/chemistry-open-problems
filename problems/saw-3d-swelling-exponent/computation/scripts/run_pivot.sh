#!/usr/bin/env bash
# run_pivot.sh -- production pivot-algorithm runs (NUMERICS ONLY, a sanity check of the setup, not a result).
# Environment: Linux, g++ >= 9, bash, xargs.  Runs at most 3 processes at a time (resource cap of the lab).
# Runtime: about 15 min wall clock with 3 parallel processes on a free Xeon 2.8 GHz (much longer on a shared
# machine); memory < 100 MB.
# Deterministic: fixed seeds; every line of the output is reproducible bit for bit on the same compiler/CPU.
# Output: ../data/pivot_runs.txt, one line per run:
#   R n attempts acceptance_fraction mean_R2 err tau_int(attempts) var/mean^2 seconds
#   S n E[s_1.s_n] err E[s_1.(w_n-w_0)] err
# (err = largest batch-means error over the blocking levels; see src/pivot.cpp).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
g++ -O3 -march=native -std=c++17 "$ROOT/src/pivot.cpp" -o "$TMP/pivot"
OUT="$ROOT/data/pivot_runs.txt"
# n attempts burnin seed  (burn-in >= 10 n attempts; from a straight rod)
cat > "$TMP/jobs" <<'EOF'
32000 2500000 640000 12
16000 4000000 320000 14
8000 6000000 160000 15
4000 8000000 80000 16
2000 12000000 40000 17
1000 16000000 20000 18
500 20000000 10000 19
250 25000000 10000 20
125 25000000 10000 21
EOF
{
  echo "# pivot runs, $(date -u +%Y-%m-%dT%H:%MZ), $(g++ --version | head -1), $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
  echo "# columns: R n attempts acc_frac mean_R2 err tau_int var/mean^2 seconds ; S n E[s_1.s_n] err E[s_1.(w_n-w_0)] err"
  xargs -P 3 -L 1 "$TMP/pivot" < "$TMP/jobs" | grep '^[RS] ' | sort -k1,1 -k2,2n
} > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
cat "$OUT"
