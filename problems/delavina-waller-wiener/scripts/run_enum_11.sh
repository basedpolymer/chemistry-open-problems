#!/bin/bash
# All connected graphs on 11 vertices (1,006,700,565 graphs), four parallel slices.
# Runtime about 5-10 minutes on 4 cores.
source ~/envs/sci/bin/activate
cd ~/dw
for r in 0 1 2 3; do
  (geng -cq 11 $r/4 | ./enum_stream 11 out_11_$r.txt; echo "slice $r done") &
done
wait
echo ALLDONE
