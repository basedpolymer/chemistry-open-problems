#!/bin/bash
# Build the generator and the checker into ~/work/saw-diamond-fcc/bin.
# Environment: Linux (WSL Fedora 44), g++ 16 with OpenMP.  Runtime: < 1 min.
set -e
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BIN=~/work/saw-diamond-fcc/bin
mkdir -p "$BIN"
g++ -O3 -march=native -fopenmp -o "$BIN/latsaw_gen.new" "$HERE/src/latsaw_gen.cpp" && mv "$BIN/latsaw_gen.new" "$BIN/latsaw_gen"
g++ -O3 -march=native -std=c++17 -pthread -o "$BIN/latsaw_check.new" "$HERE/src/latsaw_check.cpp" && mv "$BIN/latsaw_check.new" "$BIN/latsaw_check"
echo "built $BIN/latsaw_gen $BIN/latsaw_check"
