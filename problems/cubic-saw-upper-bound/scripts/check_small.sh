#!/bin/bash
# Re-check the certificates stored in certificates/ with the C++ checker, the brute-force Python checker (K <= 10)
# and the tamper tests.  Prints PASS/FAIL; exit code 0 iff everything passes.
# Environment: WSL Fedora 44, g++ 16, Python 3 (standard library only), gzip.  Runtime: ~10 min (the K = 10
# brute force takes ~6 min), < 1 GB RAM.  Run scripts/build.sh first.  THREADS (default 6) = checker threads.
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BIN=~/work/cubic-saw/bin
TMP=~/work/cubic-saw/check
mkdir -p "$TMP"
TH="${THREADS:-6}"
fail=0
check() {   # $1 = certificate file (uncompressed), $2 = label
  if "$BIN/verify_cert" "$1" "$TH" > "$1.verify.log"; then
    echo "PASS C++ checker $2: $(tail -1 "$1.verify.log")"
  else
    echo "FAIL C++ checker $2"; cat "$1.verify.log"; fail=1
  fi
}
for f in "$HERE"/certificates/*.cert.gz "$HERE"/certificates/*.cert2.gz; do
  [ -f "$f" ] || continue
  b=$(basename "$f" .gz)
  gzip -dc "$f" > "$TMP/$b" || { echo "FAIL cannot decompress $f"; fail=1; continue; }
  check "$TMP/$b" "$b"
done
for K in 8 10; do
  if python3 -I "$HERE/scripts/verify_small_bruteforce.py" "$TMP/d3_k$K.cert" > "$TMP/brute_d3_k$K.log"; then
    echo "PASS brute-force Python checker d3_k$K.cert: $(tail -1 "$TMP/brute_d3_k$K.log")"
  else
    echo "FAIL brute-force Python checker d3_k$K.cert"; fail=1
  fi
done
for c in d3_k12.cert d3_k12_explicit.cert2; do
  if python3 -I "$HERE/scripts/tamper_test.py" "$BIN/verify_cert" "$TMP/$c" "$TMP" > "$TMP/tamper_$c.log"; then
    echo "PASS tamper test on $c (all corrupted copies rejected)"
  else
    echo "FAIL tamper test on $c"; cat "$TMP/tamper_$c.log"; fail=1
  fi
done
[ $fail -eq 0 ] && echo "ALL PASS" || echo "SOME FAIL"
exit $fail
