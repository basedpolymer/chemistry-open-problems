#!/bin/bash
# Produce a VIPR certificate that instance <name> has no feasible network with <= K matches.
# usage: run_vipr_certificate.sh <instance-name> <K> <output-dir>
#   1. writes the strong LP model <output-dir>/<name>_le<K>.lp (scripts/hen_lpfile.py:write_lp_strong with
#      hen_models.strong_rows, objective sum y, row maxm: sum y <= K); if certificates/vipr/<name>_le<K>.lp exists,
#      the regenerated file must be byte-identical to it;
#   2. runs exact SCIP (exact/enable = TRUE, no presolving, no cutting planes, no conflict analysis:
#      certificates/vipr/scip_exact_safe.set) and writes the certificate <output-dir>/<name>_le<K>.vipr;
#   3. prints size and SHA-256 of the LP file and of the certificate.
# The certificate is NOT trusted by itself: check it with scripts/verify_vipr.sh.
# Environment: SCIP 10.0.2 built with -DEXACTSOLVE=on and SoPlex 8.0.2 (GMP 6.3, MPFR 4.2, Boost); set SCIP to the
#   binary.  python3 >= 3.10 (standard library).  Runtime: unbalanced10, K = 24: about 5 min, certificate 0.5 GB.
set -eu
cd "$(dirname "$0")"
NAME=$1; K=$2; OUT=$(readlink -f "$3"); mkdir -p "$OUT"
SCIP=${SCIP:-scip}
D=${HEN_DATA:-../certificates/sources/min_matches_heuristics/data/mip_instances}
DAT=$(ls $D/*/$NAME.dat)
LP=$OUT/${NAME}_le${K}.lp
CERT=$OUT/${NAME}_le${K}.vipr
python3 - "$DAT" "$LP" "$K" <<'PY'
import sys; sys.path.insert(0, '.')
from hen_instance import parse_dat
from hen_lpfile import write_lp_strong
from hen_models import strong_rows
inst = parse_dat(sys.argv[1])
info = write_lp_strong(inst, sys.argv[2], strong_rows(inst), maxmatches=int(sys.argv[3]), objective=True)
print(f"LP written: {len(info['pairs'])} y, {info['nq']} q, {info['nr']} r, {info['nrows']} extra y-rows, heat scaled by 10^-{info['scale_pow']}")
PY
STORED=../certificates/vipr/${NAME}_le${K}.lp
if [ -f "$STORED" ]; then
  if cmp -s "$STORED" "$LP"; then echo "PASS: regenerated LP identical to $STORED"; else echo "FAIL: regenerated LP differs from $STORED"; exit 1; fi
fi
$SCIP -c "set load ../certificates/vipr/scip_exact_safe.set" -c "set display freq 5000" \
      -c "set certificate filename $CERT" -c "read $LP" -c optimize -c quit
echo "LP          : $(stat -c %s "$LP") bytes, sha256 $(sha256sum "$LP" | cut -d' ' -f1)"
echo "certificate : $(stat -c %s "$CERT") bytes, sha256 $(sha256sum "$CERT" | cut -d' ' -f1)"
