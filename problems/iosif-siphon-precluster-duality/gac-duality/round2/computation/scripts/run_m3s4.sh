#!/bin/bash
# Task 4 (first-round loose end): s = 4 with molecularity 3.  Seeded samples (first-round sampler, imported
# unchanged) and the full first-round test suite on them.  Run under WSL from round2/computation/:
#   bash scripts/run_m3s4.sh
# Records go to ~/gacE2-data/ (outside the repository); logs to logs/.  Each step is capped at 25 minutes.
set -u
E=../../computation/scripts
LOG=$(pwd)/logs
D=~/gacE2-data
mkdir -p "$D"
cd "$D"
timeout 1500 python "$OLDPWD/$E/sample_wr.py" 5 3 5 8 20000 3 --s 4 --cons --nproc 10 --out sample_5_3_5_8_s3.pkl.gz > "$LOG/sample_5_3_5_8_s3.log" 2>&1
timeout 1500 python "$OLDPWD/$E/sample_wr.py" 6 3 5 8 10000 4 --s 4 --cons --nproc 10 --out sample_6_3_5_8_s4.pkl.gz > "$LOG/sample_6_3_5_8_s4.log" 2>&1
timeout 1500 python "$OLDPWD/$E/tests.py" sample_5_3_5_8_s3.pkl.gz sample_6_3_5_8_s4.pkl.gz --nproc 10 --json "$LOG/tests_m3s4.json" > "$LOG/tests_m3s4.log" 2>&1
timeout 1500 python "$OLDPWD/$E/tests.py" sample_5_3_5_8_s3.pkl.gz sample_6_3_5_8_s4.pkl.gz --only X- --nproc 10 --json "$LOG/tests_m3s4X.json" > "$LOG/tests_m3s4X.log" 2>&1
sha256sum sample_5_3_5_8_s3.pkl.gz sample_6_3_5_8_s4.pkl.gz > "$LOG/data_sha256_m3s4.txt"
echo done
