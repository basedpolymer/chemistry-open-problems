#!/bin/bash
# V1: second batch (larger n where the first batch's n was below n_0; nearly equal volatilities; forcing at 99.9 % of
# the explicit bound).  At most 2 processes here (a third thread was busy).  Usage (in WSL, from scripts/): bash run_batch2.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cat > /tmp/v1_cases2.txt <<'CASES'
rand:4:4:BBT:0.9 60 150,100
rand:5:2:TBTB 60 150,110
near:4:1 60 120,90
near:5:2 60
rand:4:5:BBB:0.999 60 120
rand:4:6:TTT:0.999 60 120
near:3:3:TB 60 160,120
CASES
run_one() {
  set -- $1
  name=$(echo "$1" | tr ':' '_')
  python v1_independent.py "$1" "$2" "$3" > ../runs/v1_$name.log 2>&1
  echo "$1: exit $? ($(tail -1 ../runs/v1_$name.log))"
}
export -f run_one
tr '\n' '\0' < /tmp/v1_cases2.txt | xargs -0 -P 2 -I{} bash -c 'run_one "{}"'
# High-precision rerun of rand:5:2 at n = 150 (resolves the sign of a 1e-321 trace content, see verification.md §8.3):
#   python v1_independent.py rand:5:2:TBTB 60 150 230 > ../runs/v1_rand_5_2_TBTB_n150_extra230.log 2>&1
#   python v1_independent.py rand:5:2:TBTB 60 150 30  > ../runs/v1_rand_5_2_TBTB_n150_extra30.log 2>&1   (default precision)
