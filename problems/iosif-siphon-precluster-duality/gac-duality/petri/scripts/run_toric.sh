#!/bin/bash
# Runs check_toric.py on EC families (agent P). Logs in ../logs/. Usage (WSL): bash run_toric.sh
cd "$(dirname "$0")"
for fam in "2 3 5 0" "3 2 5 0" "3 3 6 1" "4 2 5 1" "3 3 4 0" "4 2 4 0"; do
  set -- $fam
  log=../logs/toric_n$1_m$2_r$3_ord$4.log
  timeout 1800 python check_toric.py $1 $2 $3 $4 6 3 > $log 2>&1
  echo "family $fam exit=$?" >> ../logs/run_toric_status.log
done
