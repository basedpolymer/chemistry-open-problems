#!/bin/bash
# Build all programs into ~/work/cubic-saw/bin.  Environment: WSL Fedora 44, g++ 16 with OpenMP.  Runtime: < 1 min.
set -e
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BIN=~/work/cubic-saw/bin
mkdir -p "$BIN"
for p in saw_prod saw_big verify_cert saw_automaton count_states; do
  g++ -O3 -march=native -fopenmp -o "$BIN/$p.new" "$HERE/src/$p.cpp"
  mv "$BIN/$p.new" "$BIN/$p"
  echo "built $BIN/$p"
done
