#!/bin/bash
# run_regular_search.sh -- k-regular graphs on n vertices (nauty geng -d k -D k), each tested by kscheck;
# with COMPL=1 also the complements (nauty complg), i.e. all (n-1-k)-regular graphs on n vertices.
# Usage: scripts/run_regular_search.sh k n [C]   (C: only 2-connected graphs, geng -C)
# Output: certificates/search/regular_k<k>_n<n>[_C].{log,out} and ..._compl.{log,out}
set -uo pipefail
NAUTY=${NAUTY:-certificates/sources/nauty/nauty2_8_9}
k=$1; n=$2; flag=${3:-}
tag="regular_k${k}_n${n}${flag:+_C}"
opt="-d$k -D$k"; [ "$flag" = "C" ] && opt="-C $opt"
start=$(date +%s)
if [ "${COMPL:-0}" = "1" ]; then
  nice -n 10 "$NAUTY/geng" -q $opt "$n" | tee >(nice -n 10 "$NAUTY/complg" -q | nice -n 10 scripts/kscheck > "certificates/search/${tag}_compl.out" 2> "certificates/search/${tag}_compl.log") \
     | nice -n 10 scripts/kscheck > "certificates/search/$tag.out" 2> "certificates/search/$tag.log"
  sleep 5
else
  nice -n 10 "$NAUTY/geng" -q $opt "$n" | nice -n 10 scripts/kscheck > "certificates/search/$tag.out" 2> "certificates/search/$tag.log"
fi
echo "elapsed_seconds $(( $(date +%s)-start ))" >> "certificates/search/$tag.log"
echo "$tag done: $(head -1 certificates/search/$tag.log)"
