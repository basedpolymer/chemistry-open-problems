#!/bin/bash
# post-census queue (sequential, one heavy process at a time)
cd /home/user/chemistry-open-problems/problems/kirchhoff-soltes-problem
while [ ! -f certificates/search/queue2.status ]; do sleep 60; done
export COMPL=1
for args in "3 20" "4 15" "5 14"; do scripts/run_regular_search.sh $args >> certificates/search/regular_run.txt 2>&1; done
for m in 3 4 5 6 7 8 9 10 11 12 13; do
  start=$(date +%s)
  nice -n 10 python3 scripts/bicirculant_enum.py $m 2> certificates/search/bicirculant_m$m.gen.log | nice -n 10 scripts/kscheck_fast > certificates/search/bicirculant_m$m.out 2> certificates/search/bicirculant_m$m.log
  echo "elapsed_seconds $(( $(date +%s)-start ))" >> certificates/search/bicirculant_m$m.log
done
echo done-bicirculants > certificates/search/queue5.bic.status
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for f in "tori4 48 1000" "gp 24 1000" "tri6 48 600"; do
  set -- $f
  start=$(date +%s)
  nice -n 10 python3 scripts/fourier_families.py $1 $2 $3 1e-9 > certificates/search/fourier_$1_$2_$3.log 2>&1
  echo "elapsed_seconds $(( $(date +%s)-start ))" >> certificates/search/fourier_$1_$2_$3.log
done
echo done-fourier > certificates/search/queue5.fourier.status
export COMPL=1
for args in "4 16" "6 14" "3 22"; do scripts/run_regular_search.sh $args >> certificates/search/regular_run.txt 2>&1; done
echo done-queue5 > certificates/search/queue5.status
