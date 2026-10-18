#!/bin/bash
# Independent verification of the cubic-saw-upper-bound certificates (verifier's programs only).
# Environment: WSL Fedora 44, g++ 16, Python 3 (standard library).  Resources: <= 4 threads, peak 3.4 GB RAM
# (K = 18 SAWCERT2, streamed child table).  Runtime: ~35 min (K = 16 search ~6.5 min, K = 18 ~2.5 min, P-T d = 3
# k = 16 ~5 min, d = 2 k = 22 ~3 min, c_n brute force ~1.5 min).  Reads the large certificates from
# ~/work/cubic-saw (read only); writes only to $WORK (default ~/verif_saw).  Prints PASS/FAIL; exit 0 iff all pass.
set -u
HERE="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${WORK:-$HOME/verif_saw}"; SRC=~/work/cubic-saw
mkdir -p "$WORK/small" "$WORK/tamper"; cd "$WORK"
fail=0
g++ -O2 -march=native -std=c++17 -pthread "$HERE/verification/sawcheck.cpp" -o sawcheck || exit 1
g++ -O2 -march=native "$HERE/verification/count_saw.cpp" -o count_saw || exit 1
g++ -O2 -march=native -std=c++17 "$HERE/verification/pt_auto.cpp" -o pt_auto || exit 1
mkdir -p solver && g++ -O3 -march=native -fopenmp -o solver/verify_cert "$HERE/src/verify_cert.cpp" || exit 1

echo "== SHA-256 (against certificates/SHA256SUMS.txt)"
for f in "$HERE"/certificates/*.cert.gz "$HERE"/certificates/*.cert2.gz; do gzip -dc "$f" > "small/$(basename "$f" .gz)"; done
declare -A big=([d3_k16.cert]=$SRC/cert/d3_k16.cert [d4_k14.cert]=$SRC/hd/d4_k14.cert [d5_k14.cert]=$SRC/hd/d5_k14_p4.cert
                [d3_k18.cert]=$SRC/cert/d3_k18.cert [d3_k18.cert2]=$SRC/cert/d3_k18.cert2)
while read -r h name rest; do
  [[ "$h" =~ ^[0-9a-f]{64}$ ]] || continue
  f="small/$name"; [ -n "${big[$name]:-}" ] && f="${big[$name]}"
  [ -f "$f" ] || { echo "MISSING $name"; fail=1; continue; }
  [ "$(sha256sum "$f" | cut -d' ' -f1)" = "$h" ] && echo "PASS sha256 $name" || { echo "FAIL sha256 $name"; fail=1; }
done < "$HERE/certificates/SHA256SUMS.txt"

echo "== independent checker on every certificate"
for f in small/d3_k8.cert small/d3_k10.cert small/d3_k12.cert small/d3_k12_explicit.cert2 small/d3_k14.cert \
         small/d5_k10.cert small/d5_k12.cert small/d6_k10.cert small/d6_k12.cert \
         $SRC/cert/d3_k16.cert $SRC/hd/d4_k14.cert $SRC/hd/d5_k14_p4.cert $SRC/cert/d3_k18.cert2; do
  out=$(/usr/bin/time -f 'wall %e s, maxRSS %M kB' ./sawcheck "$f" 4 full 2>&1); rc=$?
  echo "$out" | grep -E '^(file|exact max ratio|failures|PASS|FAIL|wall)' ; [ $rc -eq 0 ] || fail=1
done

echo "== c_n by brute force vs Clisby-Liang-Slade tables"
for dn in "3 16" "4 12" "5 10" "6 9"; do set -- $dn
  ./count_saw $1 $2 > cn_d$1.txt
  bad=$(join <(sort -k1,1 cn_d$1.txt) <(tr -s ' \t' ' ' < "$HERE/certificates/counts/saw_d$1.txt" | sed 's/^ //' | sort -k1,1) | awk '$2!=$3' | wc -l)
  [ "$bad" -eq 0 ] && echo "PASS c_n d=$1 n<=$2" || { echo "FAIL c_n d=$1"; fail=1; }
done

echo "== a_n (checked automata) >= c_n"
for c in d3_k8 d3_k10 d3_k12 d3_k12_explicit d3_k14 d5_k10 d5_k12 d6_k10 d6_k12; do
  f=small/$c.cert; [ -f $f ] || f=small/$c.cert2; d=${c:1:1}; nmax=24; [ $d = 3 ] && nmax=30
  ./sawcheck $f 4 full $((1<<40)) $nmax > out_$c.txt
  printf '%-16s ' $c; python3 -I "$HERE/verification/compare_counts.py" out_$c.txt "$HERE/certificates/counts/saw_d$d.txt" || fail=1
done

echo "== Ponitz-Tittmann automaton (independent implementation; floating point, comparison with Table 2)"
for dk in "2 4" "2 6" "2 8" "2 10" "2 12" "2 14" "2 16" "2 18" "2 20" "2 22" "3 4" "3 6" "3 8" "3 10" "3 12" "3 14" "3 16" \
          "4 4" "4 6" "4 8" "4 10" "4 12" "4 14" "5 4" "5 6" "5 8" "5 10" "6 10"; do ./pt_auto $dk 60000; done

echo "== falsification tests (both checkers)"
C1="./sawcheck {} 4 full >/dev/null 2>&1"; C2="./solver/verify_cert {} 4 >/dev/null 2>&1"
T="$HERE/verification/tamper.py"
python3 -I $T small/d3_k12.cert tamper 1971676015178709/415809901234900 "$C1" "$C2" || fail=1
python3 -I $T small/d3_k12_explicit.cert2 tamper 55211756701651/11643695445633 "$C1" "$C2" || fail=1
python3 -I $T small/d5_k10.cert tamper 6602374288228685/745326036768363 "$C1" "$C2" || fail=1
python3 -I $T small/d6_k10.cert tamper 8873828103941960/815034755979641 "$C1" "$C2" || fail=1

[ $fail -eq 0 ] && echo "ALL PASS" || echo "SOME FAIL"
exit $fail
