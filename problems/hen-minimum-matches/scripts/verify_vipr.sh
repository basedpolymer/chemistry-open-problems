#!/bin/bash
# Verify a VIPR certificate that "instance <name> has no feasible network with <= K matches".
# usage: verify_vipr.sh <instance-name> <K> [<certificate.vipr>]
# Steps (each prints PASS/FAIL; exit 1 on failure):
#  1. regenerate the strong LP file from the data (scripts/hen_lpfile.py) and compare its sha256 with the stored one;
#  2. re-derive every LP row from first principles (check_lp_model.py);
#  3. check that the certificate's problem section is a relaxation of the LP model, RTP infeas (check_vipr_model.py);
#  4. complete weak derivations with viprcomp (untrusted) and check the result with the patched viprchk (trusted).
# Environment: python3 + numpy; VIPR tools built from github.com/scipopt/vipr (commit 30f2951) with
#   certificates/viprchk_canonicalize.patch; set VIPRCOMP / VIPRCHK to their paths.
# Runtime: dominated by viprcomp + viprchk (minutes per GB of certificate).
set -u
cd "$(dirname "$0")"
NAME=$1; K=$2
CERT=${3:-../certificates/vipr/${NAME}_le${K}.vipr}
D=${HEN_DATA:-../certificates/sources/min_matches_heuristics/data/mip_instances}
DAT=$(ls $D/*/$NAME.dat)
LP=../certificates/vipr/${NAME}_le${K}.lp
TMP=$(mktemp -d)
fail=0
python3 - "$NAME" "$DAT" "$TMP/regen.lp" "$K" <<'PY'
import sys; sys.path.insert(0,'.')
from hen_instance import parse_dat
from hen_lpfile import write_lp_strong
from hen_models import strong_rows
inst=parse_dat(sys.argv[2]); write_lp_strong(inst, sys.argv[3], strong_rows(inst), maxmatches=int(sys.argv[4]), objective=True)
PY
if cmp -s "$TMP/regen.lp" "$LP"; then echo "PASS 1: LP file regenerated identically ($(sha256sum $LP | cut -c1-16)...)"; else echo "FAIL 1: LP file differs from regeneration"; fail=1; fi
python3 check_lp_model.py "$DAT" "$LP" "$K" | tail -1 | sed 's/^/step 2: /' ; [ ${PIPESTATUS[0]} = 0 ] || fail=1
python3 check_vipr_model.py "$LP" "$CERT" infeas | tail -1 | sed 's/^/step 3: /'; [ ${PIPESTATUS[0]} = 0 ] || fail=1
VIPRCOMP=${VIPRCOMP:-viprcomp}; VIPRCHK=${VIPRCHK:-viprchk}
ln -s "$(readlink -f $CERT)" "$TMP/c.vipr"
$VIPRCOMP --threads=1 --soplex=off "$TMP/c.vipr" > "$TMP/comp.log" 2>&1
out=$($VIPRCHK "$TMP/c_complete.vipr" 2>&1 | tail -3)
echo "$out" | sed 's/^/step 4 viprchk: /'
echo "$out" | grep -q "Successfully verified infeasibility" && echo "PASS 4: VIPR certificate verified (infeasible)" || { echo "FAIL 4"; fail=1; }
rm -rf "$TMP"
[ $fail = 0 ] && echo "ALL PASS: no feasible network of $NAME with <= $K matches" || { echo "SOME FAIL"; exit 1; }
