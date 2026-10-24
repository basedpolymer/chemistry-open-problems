#!/bin/bash
# run_regular_search.sh -- k-regular graphs on n vertices (nauty geng -d k -D k), each tested by kscheck_fast;
# with COMPL=1 also the complements (nauty complg), i.e. all (n-1-k)-regular graphs on n vertices.
# Sequential (one heavy process at a time): geng output goes to a temporary file in $TMPDIR_KS.
# Usage: scripts/run_regular_search.sh k n [C]   (C: only 2-connected graphs, geng -C)
# Output: certificates/search/regular_k<k>_n<n>[_C].{log,out} and ..._compl.{log,out}
set -uo pipefail
NAUTY=${NAUTY:-certificates/sources/nauty/nauty2_8_9}
TMPDIR_KS=${TMPDIR_KS:-/tmp/claude-0/-home-user-chemistry-open-problems/36d79390-3b54-5b69-86a2-80d7d98b1378/scratchpad}
k=$1; n=$2; flag=${3:-}
tag="regular_k${k}_n${n}${flag:+_C}"
opt="-d$k -D$k"; [ "$flag" = "C" ] && opt="-C $opt"
start=$(date +%s)
tmp="$TMPDIR_KS/$tag.g6"
nice -n 10 "$NAUTY/geng" -q $opt "$n" > "$tmp"
nice -n 10 scripts/kscheck_fast < "$tmp" > "certificates/search/$tag.out" 2> "certificates/search/$tag.log"
echo "elapsed_seconds $(( $(date +%s)-start ))" >> "certificates/search/$tag.log"
if [ "${COMPL:-0}" = "1" ]; then
  start=$(date +%s)
  nice -n 10 "$NAUTY/complg" -q < "$tmp" | nice -n 10 scripts/kscheck_fast > "certificates/search/${tag}_compl.out" 2> "certificates/search/${tag}_compl.log"
  echo "elapsed_seconds $(( $(date +%s)-start ))" >> "certificates/search/${tag}_compl.log"
fi
rm -f "$tmp"
echo "$tag done: $(head -1 certificates/search/$tag.log)"
