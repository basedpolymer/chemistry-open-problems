#!/bin/bash
# V1: batch of independent checks (v1_independent.py), at most 3 processes in parallel, one thread each.
# Usage (from this folder's scripts/ directory, in WSL): bash run_batch.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p ../runs
cat > /tmp/v1_cases.txt <<'CASES'
n3:BT 60 60,40,30
n3:BB 60 60,40,30
n3:TB 60 60,40,30
n3:TT 60 60,40,30
n4:BBB 60 50,35
n4:TTT 60 50,35
n4:TBT 60 50,35
rand:4:1 60 50,35
rand:4:2 60 50,35
rand:4:3:TTB:0.9 60 50,35
rand:4:4:BBT:0.9 60 50,35
rand:5:1 60 45
rand:5:2:TBTB 60 45
rand:5:3:BBBB:0.9 60
rand:5:4:TTTT:0.05 60
rand:6:1 50
rand:6:2:TTTTT 50
rand:6:3:BTBTB:0.9 50
CASES
run_one() {
  set -- $1
  name=$(echo "$1" | tr ':' '_')
  python v1_independent.py "$1" "$2" "$3" > ../runs/v1_$name.log 2>&1
  echo "$1: exit $? ($(tail -1 ../runs/v1_$name.log))"
}
export -f run_one
tr '\n' '\0' < /tmp/v1_cases.txt | xargs -0 -P 3 -I{} bash -c 'run_one "{}"'
