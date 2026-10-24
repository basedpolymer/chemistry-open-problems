#!/usr/bin/env python3
"""export_tables.py -- write human-readable exact tables (TSV) derived from the enumeration data file.

Environment: Python 3.8+ standard library.  Runtime < 5 s.
Usage: python3 export_tables.py [data file, default ../data/saw_enum_N19.txt.gz]
Writes (all values exact integers or reduced rationals p/q):
  ../data/exact_summary.tsv     n, c_n, S_n, E_n = S_n/c_n, b_n, h_n, #{x_1(w_n)=0}, sum |w_n|^4, polygons p_{n+1},
                                bridge sum |w_n|^2, half-space sum |w_n|^2, max_x c_n(x) and its argmax classes
  ../data/endpoint_small_x.tsv  c_n(x) for every orbit class (a>=b>=c>=0) and every n (long format: n a b c c_n(x))
  ../data/intermediate_E.tsv    sum over n-step SAWs of |w_k|^2 for 1 <= k <= n (so E_n|w_k|^2 = value / c_n)
  ../data/bridge_span.tsv       number of n-step bridges with span x_1(w_n) = s
"""
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sawdata import SawData  # noqa: E402

ROOT = os.path.dirname(HERE)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "saw_enum_N19.txt.gz")
    D = SawData(path)
    src = os.path.basename(path)
    out = os.path.join(ROOT, "data")
    with open(os.path.join(out, "exact_summary.tsv"), "w") as f:
        f.write("# exact SAW data on Z^3 from %s (src/saw_enum.cpp); bridges: 0 < x1(w_i) <= x1(w_n), i >= 1;\n" % src)
        f.write("# half-space: x1(w_i) > 0, i >= 1; p_{n+1} = 3 c_n(e1)/(n+1) (n >= 3).\n")
        f.write("n\tc_n\tS_n\tE_n\tb_n\th_n\tZ0_n\tS4_n\tp_{n+1}\tSbridge_n\tShalf_n\tmax_x c_n(x)\targmax\n")
        for n in range(D.N + 1):
            E = F(D.S(n), D.c[n])
            p = (3 * D.cnx(n, (1, 0, 0))) // (n + 1) if n >= 3 else 0
            best, arg = D.max_point(n)
            f.write("%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n" % (
                n, D.c[n], D.S(n), E, D.b[n], D.h[n], D.ending_on_plane(n), D.S(n, 2), p,
                D.S_bridge(n) if n else 0, D.S_half(n) if n else 0, best,
                ";".join("(%d,%d,%d)" % a for a in arg)))
    with open(os.path.join(out, "endpoint_small_x.tsv"), "w") as f:
        f.write("# c_n(x) for each orbit class x = (a,b,c), a>=b>=c>=0 (value per lattice point), from %s\n" % src)
        f.write("n\ta\tb\tc\tc_n(x)\torbit_size\n")
        for n in range(D.N + 1):
            for (a, b, c), (Ns, orb) in sorted(D.A[n].items()):
                f.write("%d\t%d\t%d\t%d\t%d\t%d\n" % (n, a, b, c, Ns // orb, orb))
    with open(os.path.join(out, "intermediate_E.tsv"), "w") as f:
        f.write("# sum over n-step SAWs of |w_k|^2, 1 <= k <= n, from %s; E_n|w_k|^2 = value / c_n\n" % src)
        f.write("n\tk\tsum_|w_k|^2\n")
        for n in range(1, D.N + 1):
            for k in range(1, n + 1):
                f.write("%d\t%d\t%d\n" % (n, k, D.S_mid(k, n)))
    with open(os.path.join(out, "bridge_span.tsv"), "w") as f:
        f.write("# number of n-step bridges (0 < x1(w_i) <= x1(w_n)) with span s = x1(w_n), from %s\n" % src)
        f.write("n\ts\tcount\n")
        for n in range(1, D.N + 1):
            for s, v in sorted(D.bridge_span_dist(n).items()):
                f.write("%d\t%d\t%d\n" % (n, s, v))
    print("wrote exact_summary.tsv, endpoint_small_x.tsv, intermediate_E.tsv, bridge_span.tsv for n <= %d" % D.N)


if __name__ == "__main__":
    main()
