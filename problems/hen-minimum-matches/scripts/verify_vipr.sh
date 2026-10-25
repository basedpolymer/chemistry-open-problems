#!/bin/bash
# Verify a VIPR certificate that "instance <name> has no feasible network with <= K matches".
# usage: verify_vipr.sh <instance-name> <K> [<certificate.vipr>]
# Steps (each prints PASS/FAIL; exit 1 on failure):
#  1. regenerate the strong LP file from the data (scripts/hen_lpfile.py) and compare it byte for byte with the
#     stored one (certificates/vipr/<name>_le<K>.lp);
#  2. re-derive every LP row and bound from first principles (check_lp_model.py);
#  3. check that the certificate's problem section is a relaxation of the LP model, RTP infeas (check_vipr_model.py);
#  4. viprchk (VIPR, commit 30f2951, with certificates/viprchk_canonicalize.patch and
#     certificates/viprchk_harden.patch) on the raw certificate if it has no weak/incomplete derivation (4a), and on
#     the output of viprcomp (untrusted completion step) (4b);
#  5. a second, independently written checker (vipr_check_independent.py) on the same certificate.
# Environment: python3 (+ gmpy2 for speed); VIPR tools built from github.com/scipopt/vipr (commit 30f2951) with both
#   patches; set VIPRCOMP / VIPRCHK to their paths.
# Runtime (unbalanced10, K = 24, 0.52 GB certificate): steps 1-3 about 1 min, step 4 a few minutes, step 5 about
#   15 min (WSL2, one core of a desktop CPU).
set -u
cd "$(dirname "$0")"
NAME=$1; K=$2
CERT=${3:-../certificates/vipr/${NAME}_le${K}.vipr}
D=${HEN_DATA:-../certificates/sources/min_matches_heuristics/data/mip_instances}
DAT=$(ls $D/*/$NAME.dat)
LP=../certificates/vipr/${NAME}_le${K}.lp
TMP=$(mktemp -d)
fail=0
echo "certificate: $(stat -L -c %s "$CERT") bytes, sha256 $(sha256sum "$CERT" | cut -d' ' -f1)"
python3 - "$NAME" "$DAT" "$TMP/regen.lp" "$K" <<'PY'
import sys; sys.path.insert(0,'.')
from hen_instance import parse_dat
from hen_lpfile import write_lp_strong
from hen_models import strong_rows
inst=parse_dat(sys.argv[2]); write_lp_strong(inst, sys.argv[3], strong_rows(inst), maxmatches=int(sys.argv[4]), objective=True)
PY
if cmp -s "$TMP/regen.lp" "$LP"; then echo "PASS 1: LP file regenerated identically (sha256 $(sha256sum $LP | cut -d' ' -f1))"; else echo "FAIL 1: LP file differs from regeneration"; fail=1; fi
python3 check_lp_model.py "$DAT" "$LP" "$K" | tail -1 | sed 's/^/step 2: /' ; [ ${PIPESTATUS[0]} = 0 ] || fail=1
python3 check_vipr_model.py "$LP" "$CERT" infeas | tail -1 | sed 's/^/step 3: /'; [ ${PIPESTATUS[0]} = 0 ] || fail=1
VIPRCOMP=${VIPRCOMP:-viprcomp}; VIPRCHK=${VIPRCHK:-viprchk}
nweak=$(grep -a -c -E '\{ *lin +(weak|incomplete)' "$CERT")
if [ "$nweak" = 0 ]; then
  $VIPRCHK "$CERT" > "$TMP/chk_raw.log" 2>&1
  grep -E "Successfully verified|Failed|failed|Completed in" "$TMP/chk_raw.log" | tail -3 | sed 's/^/step 4a viprchk (raw certificate): /'
  grep -q "^Successfully verified infeasibility" "$TMP/chk_raw.log" && echo "PASS 4a: raw certificate verified (infeasible)" || { echo "FAIL 4a"; fail=1; }
else
  echo "step 4a: skipped ($nweak weak/incomplete derivations need viprcomp)"
fi
ln -s "$(readlink -f $CERT)" "$TMP/c.vipr"
$VIPRCOMP --threads=1 --soplex=off "$TMP/c.vipr" > "$TMP/comp.log" 2>&1
$VIPRCHK "$TMP/c_complete.vipr" > "$TMP/chk_comp.log" 2>&1
grep -E "Successfully verified|Failed|failed|Completed in" "$TMP/chk_comp.log" | tail -3 | sed 's/^/step 4b viprchk (after viprcomp): /'
grep -q "^Successfully verified infeasibility" "$TMP/chk_comp.log" && echo "PASS 4b: completed certificate verified (infeasible)" || { echo "FAIL 4b"; fail=1; }
rm -f "$TMP/c_complete.vipr"
python3 vipr_check_independent.py "$CERT" | sed 's/^/step 5: /'; [ ${PIPESTATUS[0]} = 0 ] || fail=1
rm -rf "$TMP"
[ $fail = 0 ] && echo "ALL PASS: no feasible network of $NAME with <= $K matches" || { echo "SOME FAIL"; exit 1; }
