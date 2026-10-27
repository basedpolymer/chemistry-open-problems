#!/bin/bash
# Runs check_wfec_families.py on the EC families of agent P (round2/ec, agent Q3). Run from round2/ec/ under WSL.
for fam in "2 3 5 0" "3 2 5 0" "3 3 6 1" "4 2 5 1" "4 2 4 0" "3 3 4 0"; do
  set -- $fam
  log="logs/wfec_n$1_m$2_r$3_ord$4.log"
  timeout 1800 python3 scripts/check_wfec_families.py $1 $2 $3 $4 12 > "$log" 2>&1
  echo "family $fam exit $?" >> logs/run_families_status.log
done
