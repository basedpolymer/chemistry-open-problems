#!/bin/sh
# run_all.sh -- regenerate and re-check every certificate of the all-M extension (1 thread; about 2-3 hours).
# Run from this directory.  Each step prints PASS/FAIL lines and a SUMMARY; the script stops at the first failure.
set -e
C=../certificates; L=$C/logs; R8=../../certificates
mkdir -p $L
python3 check_formulas.py                                            | tee $L/check_formulas.log
python3 asymptotics.py                                               | tee $L/asymptotics.log
# branch 1: Hopf box, nondegeneracy, l1 < 0 and H2, for all M >= 161 (eps in [0, 161^-1/2])
python3 l1_sweep.py 1 0 0.0788111 $C/l1h2_j1.csv --h2 1e-5 2e-3   > $L/l1h2_sweep_j1.log
python3 verify_sweep.py $C/l1h2_j1.csv 0 0.0788111 --h2              | tee $L/verify_l1h2_j1.log
python3 check_vs_result8.py 1 $C/l1h2_j1.csv $R8/sweep_j1.csv.gz $R8/sweep_j1_b.csv.gz | tee $L/check_vs_result8_j1.log
# branches 2 and 3: l1 < 0 for all M >= 161, then the spectral certificate on the same boxes
python3 l1_sweep.py 2 0 0.0788111 $C/l1_j2.csv 1e-5 2e-3           > $L/l1_sweep_j2.log
python3 verify_sweep.py $C/l1_j2.csv 0 0.0788111                     | tee $L/verify_l1_j2.log
python3 check_vs_result8.py 2 $C/l1_j2.csv $R8/sweep_j2.csv.gz $R8/sweep_j2_b.csv.gz | tee $L/check_vs_result8_j2.log
python3 h2j_sweep.py 2 $C/l1_j2.csv $C/h2j_j2.csv 20 20.4 -0.3 0.25  > $L/h2j_sweep_j2.log
python3 l1_sweep.py 3 0 0.0788111 $C/l1_j3.csv 1e-5 2e-3           > $L/l1_sweep_j3.log
python3 verify_sweep.py $C/l1_j3.csv 0 0.0788111                     | tee $L/verify_l1_j3.log
python3 check_vs_result8.py 3 $C/l1_j3.csv $R8/sweep_j3.csv.gz $R8/sweep_j3_b.csv.gz | tee $L/check_vs_result8_j3.log
python3 h2j_sweep.py 3 $C/l1_j3.csv $C/h2j_j3.csv 26 26.6 -0.15 0.25 > $L/h2j_sweep_j3.log
echo "run_all: all steps PASS"
