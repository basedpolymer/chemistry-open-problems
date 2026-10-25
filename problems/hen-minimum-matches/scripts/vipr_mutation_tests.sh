#!/bin/bash
# Negative tests: both VIPR checkers must REJECT corrupted copies of the unbalanced10 certificate.
# usage: vipr_mutation_tests.sh <unbalanced10_le24.vipr> <scratch-dir>
# Mutations (line numbers refer to the certificate with sha256 17967e67...; the script checks the original lines first):
#   m1  problem row maxm  "sum y <= 24"  ->  "<= 25"   (the relaxed problem is feasible: no valid proof can exist)
#   m2  first multiplier of derivation ACT_L4624 changed from -1 to -2   (wrong linear combination)
#   m3  the two branching assumptions of the first unsplit step swapped (assumptions are never discharged)
#   m4  rounding step ACT_R5354 claims y >= 2 instead of y >= 1        (rounding does not give that)
# Each mutated copy is checked by viprchk (patched) and by vipr_check_independent.py; PASS means "rejected".
# Environment: VIPRCHK = path of the patched viprchk; python3 (+ gmpy2).  Runtime: about 6 min, 0.5 GB scratch space.
set -u
cd "$(dirname "$0")"
CERT=$1; W=$2; mkdir -p "$W"
VIPRCHK=${VIPRCHK:-viprchk}
[ "$(sha256sum "$CERT" | cut -c1-16)" = "17967e67e607271c" ] || { echo "FAIL: unexpected certificate (sha256)"; exit 1; }
fail=0
mutate() {  # name line old new
  local name=$1 line=$2 old=$3 new=$4
  if ! sed -n "${line}p" "$CERT" | grep -q -F -- "$old"; then echo "FAIL $name: line $line does not contain '$old'"; fail=1; return; fi
  awk -v L="$line" -v o="$old" -v n="$new" 'NR==L{i=index($0,o); $0=substr($0,1,i-1) n substr($0,i+length(o))} {print}' "$CERT" > "$W/mut.vipr"
  r1=$($VIPRCHK "$W/mut.vipr" 2>&1 | grep -c "^Successfully verified")
  python3 vipr_check_independent.py "$W/mut.vipr" > "$W/ind.log" 2>&1; r2=$?
  if [ "$r1" = 0 ] && [ "$r2" != 0 ]; then echo "PASS $name: rejected by viprchk and by the independent checker ($(tail -1 $W/ind.log | cut -c1-110))"
  else echo "FAIL $name: viprchk accepted=$r1, independent exit=$r2"; fail=1; fi
  rm -f "$W/mut.vipr"
}
mutate m1 3666 "C2452 L 24 " "C2452 L 25 "
mutate m2 5841 "{ lin 1 4623 -1 }" "{ lin 1 4623 -2 }"
mutate m3 6675 "{ uns 5431 5427  5457 5426  }" "{ uns 5431 5426  5457 5427  }"
mutate m4 6571 "ACT_R5354 G 1 1 43 1" "ACT_R5354 G 2 1 43 1"
[ $fail = 0 ] && echo "ALL PASS: every corrupted certificate rejected by both checkers" || { echo "SOME FAIL"; exit 1; }
