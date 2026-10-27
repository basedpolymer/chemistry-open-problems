#!/bin/bash
# Runs run_toric.sh and the general families (agent P). WSL. (No longer waits for run_enum.sh.)
cd "$(dirname "$0")"
bash run_toric.sh
for fam in "2 2 5" "3 2 4" "3 3 3" "4 1 8" "4 2 3"; do
  set -- $fam
  log=../logs/general_n$1_m$2_r$3.log
  timeout 1800 python check_general.py $1 $2 $3 6 > $log 2>&1
  echo "family $fam exit=$?" >> ../logs/run_general_status.log
done
