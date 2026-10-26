#!/bin/sh
# run_all.sh -- regenerate and re-check every certificate of the all-M extension (1 thread; about 1 hour on a shared machine).
# Run from this directory (scripts/).  Each step writes its log to ../certificates/logs/, prints its SUMMARY line,
# and the script stops at the first step that exits non-zero (FAIL).
C=../certificates; L=$C/logs; R8=../../certificates
mkdir -p $L
run() {
    name=$1; shift
    if "$@" > $L/$name.log 2>&1; then
        echo "PASS  $name: $(tail -n 1 $L/$name.log)"
    else
        echo "FAIL  $name (see $L/$name.log)"; exit 1
    fi
}
run check_formulas      python3 check_formulas.py
run asymptotics         python3 asymptotics.py
# branch 1: Hopf box, nondegeneracy, l1 < 0 and H2, for all M >= 161 (eps in [0, 161^-1/2])
run l1h2_sweep_j1       python3 l1_sweep.py 1 0 0.0788111 $C/l1h2_j1.csv --h2 1e-5 2e-3
run verify_l1h2_j1      python3 verify_sweep.py $C/l1h2_j1.csv 0 0.078811099999999995 --h2
run check_vs_result8_j1 python3 check_vs_result8.py 1 $C/l1h2_j1.csv $R8/sweep_j1.csv.gz $R8/sweep_j1_b.csv.gz
# branches 2 and 3: l1 < 0 for all M >= 161, then the spectral certificate (Theorem 4.5) on the same boxes
run l1_sweep_j2         python3 l1_sweep.py 2 0 0.0788111 $C/l1_j2.csv 1e-5 2e-3
run verify_l1_j2        python3 verify_sweep.py $C/l1_j2.csv 0 0.078811099999999995
run check_vs_result8_j2 python3 check_vs_result8.py 2 $C/l1_j2.csv $R8/sweep_j2.csv.gz $R8/sweep_j2_b.csv.gz
run h2j_sweep_j2        python3 h2j_sweep.py 2 $C/l1_j2.csv $C/h2j_j2.csv 20 20.4 -0.3 0.25
run l1_sweep_j3         python3 l1_sweep.py 3 0 0.0788111 $C/l1_j3.csv 1e-5 2e-3
run verify_l1_j3        python3 verify_sweep.py $C/l1_j3.csv 0 0.078811099999999995
run check_vs_result8_j3 python3 check_vs_result8.py 3 $C/l1_j3.csv $R8/sweep_j3.csv.gz $R8/sweep_j3_b.csv.gz
run h2j_sweep_j3        python3 h2j_sweep.py 3 $C/l1_j3.csv $C/h2j_j3.csv 26 26.6 -0.15 0.25
echo "run_all: all steps PASS"
