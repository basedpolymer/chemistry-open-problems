#!/bin/bash
cd /home/user/chemistry-open-problems/problems/kirchhoff-soltes-problem
while [ ! -f certificates/search/queue3.status ]; do sleep 60; done
for m in 3 4 5 6 7 8 9 10 11 12 13; do
  start=$(date +%s)
  nice -n 10 python3 scripts/bicirculant_enum.py $m 2> certificates/search/bicirculant_m$m.gen.log | nice -n 10 scripts/kscheck_fast > certificates/search/bicirculant_m$m.out 2> certificates/search/bicirculant_m$m.log
  echo "elapsed_seconds $(( $(date +%s)-start ))" >> certificates/search/bicirculant_m$m.log
done
echo done-queue4 > certificates/search/queue4.status
