#!/bin/bash
# Numerical dual multistationarity search (round2/ec, agent Q3). Run from round2/ec/ under WSL.
timeout 1800 python3 scripts/dual_search.py 3 5 3 3 3 25 11 6 > logs/dual_search_n3_r5.log 2>&1
echo "n3 r5 exit $?" >> logs/run_dual_status.log
timeout 1800 python3 scripts/dual_search.py 4 4 3 2 2 25 13 6 > logs/dual_search_n4_r4.log 2>&1
echo "n4 r4 exit $?" >> logs/run_dual_status.log
