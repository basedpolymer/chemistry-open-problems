#!/bin/bash
# WF-EC multistationarity families (ms_wfec.py).  Run under WSL from round2/computation/: bash scripts/run_ms.sh
# Each family is capped at 25 minutes; the exit status of each run is appended to logs/run_ms_status.log.
for fam in "2 3 5 0" "3 3 4 0" "4 2 4 0" "2 4 5 0" "3 2 6 0" "4 2 5 1" "5 2 5 1" "3 3 5 0"; do
  set -- $fam
  log=logs/ms_wfec_$1_$2_$3_$4.log
  timeout 1500 python scripts/ms_wfec.py $1 $2 $3 $4 10 > $log 2>&1
  echo "family $fam exit=$?" >> logs/run_ms_status.log
done
