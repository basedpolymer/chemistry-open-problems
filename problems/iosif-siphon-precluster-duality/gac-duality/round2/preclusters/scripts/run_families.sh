#!/bin/bash
# Run check_families.py on agent E's families, one process per family, 30-minute cap each (WSL).
cd "$(dirname "$0")"
for f in all_2_2_5 all_2_2_6 all_2_3_4 all_3_2_4 all_3_3_3 all_4_1_12 all_4_2_3 fc_3_6 fc_4_5 \
         wr_3_2_6 wr_4_3_4 wr_5_2_5 wr_6_2_5_homog wr_2_2_5_zero wr_3_2_5_zero; do
  timeout 1800 python check_families.py "$f.pkl.gz" > "../logs/fam_$f.log" 2>&1 &
done
wait
echo "done" > ../logs/run_families_status.log
