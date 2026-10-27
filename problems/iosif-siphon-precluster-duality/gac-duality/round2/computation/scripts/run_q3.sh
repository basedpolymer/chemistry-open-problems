#!/bin/bash
# Q3 candidate tests on EC families (q3_tests.py).  Run under WSL from round2/computation/: bash scripts/run_q3.sh
cd scripts
for fam in "2 3 5 0" "3 3 4 0" "4 2 4 0" "2 4 4 0" "3 3 5 0"; do
  set -- $fam
  timeout 1500 python q3_tests.py $1 $2 $3 $4 4 > ../logs/q3_tests_$1_$2_$3_$4.log 2>&1
  echo "family $fam exit=$?" >> ../logs/run_q3_status.log
done
