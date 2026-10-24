#!/bin/bash
cd /home/user/chemistry-open-problems/problems/kirchhoff-soltes-problem
while [ ! -f certificates/search/queue2.status ]; do sleep 60; done
export COMPL=1
for args in "3 20" "4 15" "5 14" "4 16" "6 14" "3 22"; do scripts/run_regular_search.sh $args >> certificates/search/regular_run.txt 2>&1; done
unset COMPL
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for f in "tori4 48 1000" "gp 24 1000" "tri6 48 600"; do
  set -- $f
  start=$(date +%s)
  nice -n 10 python3 scripts/fourier_families.py $1 $2 $3 1e-9 > certificates/search/fourier_$1_$2_$3.log 2>&1
  echo "elapsed_seconds $(( $(date +%s)-start ))" >> certificates/search/fourier_$1_$2_$3.log
done
echo done-queue3 > certificates/search/queue3.status
