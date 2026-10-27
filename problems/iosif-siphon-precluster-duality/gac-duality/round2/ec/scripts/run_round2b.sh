#!/bin/bash
# Second batch (round2/ec, agent Q3): dual certificates on the families, beyond-EC scan, K8 species-side search.
# Run from round2/ec/ under WSL.
timeout 1800 python3 scripts/beyond_ec.py ~/gacE-data/all_*.pkl.gz > logs/beyond_ec.log 2>&1
echo "beyond_ec exit $?" >> logs/run_round2b_status.log
timeout 1200 python3 scripts/k8s_search.py 4 2 3 6 4000 1 6 > logs/k8s_search_n4.log 2>&1
echo "k8s n4 exit $?" >> logs/run_round2b_status.log
for fam in "2 3 5 0" "3 2 5 0" "3 3 6 1" "4 2 5 1" "4 2 4 0" "3 3 4 0"; do
  set -- $fam
  timeout 1800 python3 scripts/check_dual_families.py P $1 $2 $3 $4 6 > "logs/dual_n$1_m$2_r$3_ord$4.log" 2>&1
  echo "dual P $fam exit $?" >> logs/run_round2b_status.log
done
timeout 1800 python3 scripts/check_dual_families.py U 3 5 3 6 > logs/dual_U_n3_r5.log 2>&1
echo "dual U n3 r5 exit $?" >> logs/run_round2b_status.log
