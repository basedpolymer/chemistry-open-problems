#!/bin/bash
# Rigorous modular upper bound on the number of good vertices of every connected graph with 3 <= n <= NMAX
# (nauty geng -c).  Writes certificates/search/goodcount_n<N>.out (graphs with at least 2/5 of their vertices
# possibly good) and .log (histogram of the upper bound pg).  Usage (from this folder): scripts/run_goodcount.sh 3 10
# Build first: gcc -O3 -march=native -o scripts/goodcount scripts/goodcount.c
# Runtime: n <= 10 about 4 min; n = 11 about 5 h of one core (run as 4 slices geng -c 11 r/4).
set -e
for n in $(seq "$1" "$2"); do
  geng -c -q "$n" | scripts/goodcount 2 5 > "certificates/search/goodcount_n$n.out" 2> "certificates/search/goodcount_n$n.log"
  echo "n=$n done"; cat "certificates/search/goodcount_n$n.log"
done
