#!/bin/bash
# Runs the EC/FC enumeration families (agent P). Logs in ../logs/.  Usage (WSL): bash run_enum.sh
cd "$(dirname "$0")"
for fam in "3 3 6 1" "2 3 5 0" "3 2 5 0" "4 2 5 1" "3 3 4 0" "4 2 4 0" "5 2 5 1"; do
  set -- $fam
  log=../logs/enum_ec_n$1_m$2_r$3_ord$4.log
  timeout 1800 python enumerate_ec.py $1 $2 $3 $4 12 > $log 2>&1
  echo "family $fam exit=$?" >> ../logs/run_enum_status.log
done
