#!/bin/bash
# Enumerate all connected graphs on n = 5..10 vertices; write max W per diameter.
# Runtime about 1 minute. Binary and output go to ~/dw (outside the Windows mount, for speed).
source ~/envs/sci/bin/activate
S=/mnt/c/Users/demor/chemistry-open-problems/problems/delavina-waller-wiener/scripts
mkdir -p ~/dw && gcc -O3 -march=native -o ~/dw/enum_stream $S/enum_stream.c || exit 1
cd ~/dw
for n in 5 6 7 8 9 10; do
  start=$(date +%s)
  geng -cq $n | ./enum_stream $n out_$n.txt
  echo "n=$n done in $(( $(date +%s) - start ))s"
done
