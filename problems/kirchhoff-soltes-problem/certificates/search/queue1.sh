#!/bin/bash
cd /home/user/chemistry-open-problems/problems/kirchhoff-soltes-problem
while pgrep -f "run_geng_search.sh" > /dev/null; do sleep 30; done
start=$(date +%s)
nice -n 10 python3 scripts/kernel_enum.py 5 2> certificates/search/kernel_c5.gen.log | nice -n 10 scripts/kscheck > certificates/search/kernel_c5.out 2> certificates/search/kernel_c5.log
echo "elapsed_seconds $(( $(date +%s)-start ))" >> certificates/search/kernel_c5.log
echo done-c5 > certificates/search/queue1.status
