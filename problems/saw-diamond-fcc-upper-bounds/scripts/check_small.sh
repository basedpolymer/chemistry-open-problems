#!/bin/bash
# Re-check the certificates stored in certificates/ (diamond K <= 18, fcc K <= 9):
#   * the C++ checker latsaw_check on every stored certificate, with the a_n >= c_n cross-check against the
#     exact counts of certificates/counts/;
#   * the pure-Python brute-force checker on diamond K <= 16 and fcc K <= 8 (both already beat the previous bounds);
#   * tamper tests (corrupted copies must be rejected) on diamond K = 12 and fcc K = 6.
# Prints PASS/FAIL; exit code 0 iff everything passes.
# Environment: WSL Fedora 44, g++ 16, Python 3 (standard library), gzip.  Runtime: about 2 min, < 1 GB RAM.
# Run scripts/build.sh first.  THREADS (default 4) = checker threads.
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BIN=~/work/saw-diamond-fcc/bin
TMP=~/work/saw-diamond-fcc/check
mkdir -p "$TMP"
TH="${THREADS:-4}"
fail=0
for f in "$HERE"/certificates/*.lcert.gz; do
  b=$(basename "$f" .gz)
  gzip -dc "$f" > "$TMP/$b" || { echo "FAIL cannot decompress $f"; fail=1; continue; }
  case "$b" in dia*) cnt="$HERE/certificates/counts/diamond_cn.txt"; nmax=32;; *) cnt="$HERE/certificates/counts/fcc_cn.txt"; nmax=24;; esac
  if "$BIN/latsaw_check" "$TMP/$b" "$TH" "$cnt" "$nmax" > "$TMP/$b.check.log"; then
    echo "PASS checker $b: $(grep '(info) a_n' "$TMP/$b.check.log" | cut -c6-) ; $(tail -1 "$TMP/$b.check.log")"
  else
    echo "FAIL checker $b"; cat "$TMP/$b.check.log"; fail=1
  fi
done
for b in dia_k6 dia_k8 dia_k10 dia_k12 dia_k14 dia_k16 fcc_k3 fcc_k4 fcc_k5 fcc_k6 fcc_k7 fcc_k8; do
  if python3 -I "$HERE/scripts/verify_small_bruteforce.py" "$TMP/$b.lcert" > "$TMP/$b.brute.log"; then
    echo "PASS brute-force Python checker $b.lcert: $(tail -1 "$TMP/$b.brute.log")"
  else
    echo "FAIL brute-force Python checker $b.lcert"; cat "$TMP/$b.brute.log"; fail=1
  fi
done
for b in dia_k12 fcc_k6; do
  if python3 -I "$HERE/scripts/tamper_test.py" "$BIN/latsaw_check" "$TMP/$b.lcert" "$TMP/tamper" > "$TMP/$b.tamper.log"; then
    echo "PASS tamper test on $b.lcert (all corrupted copies rejected)"
  else
    echo "FAIL tamper test on $b.lcert"; cat "$TMP/$b.tamper.log"; fail=1
  fi
done
[ $fail -eq 0 ] && echo "ALL PASS" || echo "SOME FAIL"
exit $fail
