#!/bin/bash
# run_census.sh -- test every graph of the Holt-Royle census of vertex-transitive graphs (Zenodo 10.5281/zenodo.4010122)
# on n vertices for n in the given list, with kscheck_fast -t (exact modular criterion at vertex 0; valid by vertex-transitivity;
# graphs passing vertex 0 are re-tested on all vertices).
# Usage (from the problem folder): scripts/run_census.sh 10 11 ... 47 ; set KEEP=1 to keep downloaded tars.
# Output: certificates/search/census_n<n>.log (one line per member file: graphs read, disconnected, rejected, survivors)
#         certificates/search/census_n<n>.out (survivors / inconclusive graphs, expected empty)
# Runtime: about 3-4 h for n<=47 on one thread (n=44 and 46 dominate); download 11.4 GB in total.
set -uo pipefail
here=$(pwd)
export KSCHECK="${KSCHECK:-$here/scripts/kscheck_fast}"
SRC=certificates/sources/holt-royle
mkdir -p "$SRC" certificates/search
for n in "$@"; do
  tarf="$SRC/alltrans$n.tar"
  if [ ! -s "$tarf" ]; then
    curl -sS -m 7200 -L --retry 5 -o "$tarf" "https://zenodo.org/api/records/4010122/files/alltrans$n.tar/content" || { echo "download failed $n"; continue; }
  fi
  sha256sum "$tarf" >> "$SRC/sha256.txt"
  export CENSUS_LOG="$here/certificates/search/census_n$n.log" CENSUS_OUT="$here/certificates/search/census_n$n.out"
  : > "$CENSUS_LOG"; : > "$CENSUS_OUT"
  start=$(date +%s)
  nice -n 10 tar -xf "$tarf" --warning=no-unknown-keyword --to-command="$here/scripts/census_member.sh"
  echo "elapsed_seconds $(( $(date +%s)-start ))" >> "$CENSUS_LOG"
  [ "${KEEP:-0}" = "1" ] || [ "$n" -le 39 ] || rm -f "$tarf"
  echo "n=$n done"
done
