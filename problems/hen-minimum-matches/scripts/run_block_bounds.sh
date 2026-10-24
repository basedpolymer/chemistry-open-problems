#!/bin/bash
# Lemma 1 lower bounds (exact enumeration), PASS/FAIL per instance; exit 1 on any FAIL.
# Environment: python3 (3.10+) with numpy. Runtime ~1-2 min total (single thread), < 1.5 GB RAM.
cd "$(dirname "$0")"
D=${HEN_DATA:-../certificates/sources/min_matches_heuristics/data/mip_instances}
fail=0
for x in furman_sahinidis/37sp-yfyv:36 furman_sahinidis/20sp1:19 furman_sahinidis/23sp1:22 furman_sahinidis/22sp1:23 \
         chen_grossmann_miller/balanced12:26 grossmann_random/balanced12_random0:26 grossmann_random/balanced12_random1:26 \
         grossmann_random/balanced12_random2:26 chen_grossmann_miller/balanced15:32 grossmann_random/balanced15_random0:32 \
         grossmann_random/balanced15_random1:32 grossmann_random/balanced15_random2:32 chen_grossmann_miller/unbalanced15:32; do
  f=${x%%:*}; c=${x##*:}
  python3 block_bound.py $D/$f.dat $c || fail=1
done
[ $fail = 0 ] && echo "ALL PASS" || { echo "SOME FAIL"; exit 1; }
