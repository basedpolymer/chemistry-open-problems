#!/bin/sh
# run_branch.sh J EPS_C MJ MC M_HI R0 Y cL cR [cL2 cR2] -- all certificates for branch J (one thread; 0.5-4 h).
#   eps-chart (every M with M^-1/2 <= EPS_C, including M = infinity):
#     l1_sweep.py  (Hopf box, simplicity, transversality, Im z > 0, f(phi) != 0, Lhat < 0), verify_sweep.py (re-proof),
#     h2cover.py   (H2 and the unstable count on merged boxes);
#   onset charts (K = infinity to M = M_HI > EPS_C^-2): onset_sweep.py (the Hopf curve, junction with the eps-chart);
#   integer points MJ <= M <= MC - 1: pointwise.py (all checks at each M), pointwise.py --verify, link_check.py;
#   check_overlap.py (bdlib l1 inside the eps-chart enclosures for MC <= M <= 3 MC).
# MC must satisfy MC^-1/2 <= EPS_C (checked).  Run from scripts/ on a copy of this file (sh reads scripts
# incrementally).  Logs in ../certificates/logs/.
J=$1; EC=$2; MJ=$3; MC=$4; MHI=$5; R0=$6; Y=$7; CL=$8; CR=$9; CL2=${10}; CR2=${11}
C=../certificates; L=$C/logs
mkdir -p $L
python3 -c "import sys; sys.exit(0 if $MC ** -0.5 <= $EC else 1)" || { echo "FAIL: MC^-1/2 > EPS_C"; exit 1; }
run() {
    name=$1; shift
    if "$@" > $L/$name.log 2>&1; then
        echo "PASS  $name: $(tail -n 1 $L/$name.log)"
    else
        echo "FAIL  $name (see $L/$name.log)"; exit 1
    fi
}
run l1_sweep_j$J      python3 l1_sweep.py $J 0 $EC $C/l1_j$J.csv 1e-5 2e-3
run verify_l1_j$J     python3 verify_sweep.py $C/l1_j$J.csv 0 $EC
run h2cover_j$J       python3 h2cover.py $J $C/l1_j$J.csv $C/h2_j$J.csv $R0 $Y $CL $CR 2e-2 $CL2 $CR2
run onset_j$J         python3 onset_sweep.py $J $MHI $C/onset_j$J.csv $C/l1_j$J.csv
run pointwise_j$J     python3 pointwise.py $J $MJ $((MC - 1)) $C/pw_j$J.csv $R0 $Y $CL $CR 256 $CL2 $CR2
run verify_pw_j$J     python3 pointwise.py --verify $C/pw_j$J.csv
run link_j$J          python3 link_check.py $C/pw_j$J.csv $C/onset_j$J.csv
run overlap_j$J       python3 check_overlap.py $C/l1_j$J.csv $MC $((MC * 3))
echo "run_branch J=$J: all steps PASS"
