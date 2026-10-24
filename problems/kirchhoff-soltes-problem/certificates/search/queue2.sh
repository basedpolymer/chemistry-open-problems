#!/bin/bash
cd /home/user/chemistry-open-problems/problems/kirchhoff-soltes-problem
while [ ! -f certificates/search/queue1.status ]; do sleep 30; done
scripts/run_census.sh 32 33 34 35 36 37 38 39 41 43 45 47 40 42 44 46 > certificates/search/census_run.txt 2>&1
echo done-census > certificates/search/queue2.status
