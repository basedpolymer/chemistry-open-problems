#!/bin/bash
# run_geng_search.sh -- exhaustive search over all 2-connected graphs on n vertices (nauty geng -C),
# each graph tested by scripts/kscheck (exact modular criterion, see kscheck.c header).
# Usage: scripts/run_geng_search.sh NMIN NMAX   (from the problem folder; nauty built in certificates/sources/nauty)
# Runtime: n<=10 about 40 s; n=11 about 1-2 h (1 thread).  Output: certificates/search/geng_C_n<n>.{log,out}
set -euo pipefail
NAUTY=${NAUTY:-certificates/sources/nauty/nauty2_8_9}
mkdir -p certificates/search
for n in $(seq "$1" "$2"); do
  start=$(date +%s)
  nice -n 10 "$NAUTY/geng" -C "$n" 2> "certificates/search/geng_C_n$n.geng.log" | nice -n 10 scripts/kscheck > "certificates/search/geng_C_n$n.out" 2> "certificates/search/geng_C_n$n.log"
  end=$(date +%s)
  echo "elapsed_seconds $((end-start))" >> "certificates/search/geng_C_n$n.log"
  echo "n=$n done: $(head -1 certificates/search/geng_C_n$n.log)"
done
