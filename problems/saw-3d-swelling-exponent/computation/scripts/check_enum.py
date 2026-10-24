#!/usr/bin/env python3
"""check_enum.py -- cross-checks of the exact enumeration data (PASS/FAIL per claim; exit code 0 iff all PASS).

Environment: Linux, g++ >= 9 (C++17, pthreads), Python 3.8+ standard library only.  Tested: Ubuntu 24.04,
g++ 13.3, Python 3.13, Intel Xeon 2.8 GHz.  Runtime about 1 min (3 threads at most).

Usage: python3 check_enum.py [data file, default ../data/saw_enum_N19.txt.gz] [--rerun N2 (default 16)]

Checks
  1. brute force: pure-Python enumerator (no symmetry, scripts/brute_force.py), n <= 8, against the data file:
     c_n, h_n, b_n and all class tables A, H, B, M.
  2. naive C++ (src/saw_naive.cpp: no symmetry, linear-scan self-avoidance, pointwise B_3-invariance check) against
     a fresh run of src/saw_enum.cpp for n <= 11, with four different prefix depths (1, 4, 9, 11): identical output.
  3. re-run: fresh src/saw_enum.cpp run with N = N2 (default 16) agrees with the data file for every n <= N2;
     src/saw_gram.cpp (Gram sums) against brute force (n <= 7) and against the data file's intermediate-point table.
  4. internal consistency of the data file: orbit divisibility, parity, totals of A/H/B/M tables.
  5. published values: c_n vs OEIS A001412 (= Schram-Barkema-Bisseling 2011, Clisby-Liang-Slade 2007);
     S_n = sum |omega_n|^2 vs SBB 2011 Table I (P_N) (= CLS rho_n); polygons p_{n+1} = 3 c_n(e1)/(n+1) vs CLS;
     bridges b_n vs Clisby-Conway-Guttmann 2016 Table 6; h_{n+1} vs OEIS A116904 (terminally attached walks).
"""
import os
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sawdata import SawData, read_reference  # noqa: E402
import brute_force  # noqa: E402

ROOT = os.path.dirname(HERE)
THREADS = os.environ.get("SAW_THREADS", "3")  # resource cap: at most 3
fails = 0


def report(ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails += 1


def strip_comments(text):
    return "\n".join(l for l in text.splitlines() if not l.startswith("#"))


def main():
    args = sys.argv[1:]
    n2 = 16
    if "--rerun" in args:
        i = args.index("--rerun")
        n2 = int(args[i + 1])
        del args[i:i + 2]
    path = args[0] if args else os.path.join(ROOT, "data", "saw_enum_N19.txt.gz")
    t0 = time.time()
    D = SawData(path)
    print("data file %s: N = %d; header: %s" % (os.path.relpath(path, ROOT), D.N, " | ".join(D.header)))

    # 1. brute force
    nb = min(8, D.N)
    r = brute_force.enumerate_all(nb)
    ok = all(r["c"][n] == D.c[n] and r["h"].get(n, 0) == D.h[n] and r["b"].get(n, 0) == D.b[n] for n in range(nb + 1))
    ok = ok and all(r["A"][n] == D.A[n] and r["H"].get(n, {}) == D.H.get(n, {}) and r["B"].get(n, {}) == D.B.get(n, {})
                    for n in range(nb + 1))
    ok = ok and all(r["M"][(k, n)] == D.M[(k, n)] for n in range(2, nb + 1) for k in range(1, n))
    report(ok, "[1] Python brute force (no symmetry) = data file for n <= %d (c, h, b, A, H, B, M tables)" % nb)

    # 2./3. C++ runs
    with tempfile.TemporaryDirectory() as tmp:
        exe_e = os.path.join(tmp, "saw_enum")
        exe_n = os.path.join(tmp, "saw_naive")
        subprocess.run(["g++", "-O3", "-std=c++17", "-pthread", os.path.join(ROOT, "src", "saw_enum.cpp"), "-o", exe_e],
                       check=True)
        subprocess.run(["g++", "-O2", "-std=c++17", os.path.join(ROOT, "src", "saw_naive.cpp"), "-o", exe_n], check=True)
        pn = subprocess.run([exe_n, "11"], capture_output=True, text=True)
        report(pn.returncode == 0 and "invariance check: PASS" in pn.stdout,
               "[2a] saw_naive 11: endpoint distributions invariant under all 48 signed permutations")
        same = True
        for K in (1, 4, 9, 11):
            pe = subprocess.run([exe_e, "11", THREADS, str(K)], capture_output=True, text=True)
            same = same and pe.returncode == 0 and strip_comments(pe.stdout) == strip_comments(pn.stdout)
        report(same, "[2b] saw_enum 11 (prefix depths 1,4,9,11) output identical to saw_naive 11 (T, A, H, B, M)")
        if n2 <= D.N:
            out = os.path.join(tmp, "rerun.txt")
            with open(out, "w") as f:
                pe = subprocess.run([exe_e, str(n2), THREADS], stdout=f, stderr=subprocess.DEVNULL)
            E = SawData(out)
            ok = pe.returncode == 0
            for n in range(n2 + 1):
                ok = ok and E.c[n] == D.c[n] and E.h[n] == D.h[n] and E.b[n] == D.b[n]
                ok = ok and E.A[n] == D.A[n] and E.H.get(n) == D.H.get(n) and E.B.get(n) == D.B.get(n)
                for k in range(1, n):
                    ok = ok and E.M[(k, n)] == D.M[(k, n)]
            report(ok, "[3] fresh saw_enum run with N = %d agrees with the data file for all n <= %d" % (n2, n2))

        # 3b. Gram sums (src/saw_gram.cpp) against brute force and against the M table of the data file
        exe_g = os.path.join(tmp, "saw_gram")
        subprocess.run(["g++", "-O3", "-std=c++17", "-pthread", os.path.join(ROOT, "src", "saw_gram.cpp"), "-o", exe_g],
                       check=True)
        Gb = brute_force.gram_all(7)
        ok = True
        for K in (1, 3, 6):
            pg = subprocess.run([exe_g, "7", THREADS, str(K)], capture_output=True, text=True)
            Gc = {}
            for line in pg.stdout.splitlines():
                t = line.split()
                if t and t[0] == "G":
                    Gc[(int(t[1]), int(t[2]), int(t[3]))] = int(t[4])
            ok = ok and pg.returncode == 0 and Gc == Gb
        report(ok, "[3b] saw_gram 7 (prefix depths 1,3,6) = Python brute-force Gram sums sum w_i.w_j, n <= 7")
        ng = min(D.N, 13)
        pg = subprocess.run([exe_g, str(ng), THREADS], capture_output=True, text=True)
        Gc = {}
        for line in pg.stdout.splitlines():
            t = line.split()
            if t and t[0] == "G":
                Gc[(int(t[1]), int(t[2]), int(t[3]))] = int(t[4])
        ok = all(Gc[(n, k, k)] == D.S_mid(k, n) for n in range(1, ng + 1) for k in range(1, n + 1))
        report(ok, "[3c] saw_gram %d diagonal sum |w_k|^2 = intermediate-point table of the data file, k <= n <= %d"
               % (ng, ng))

    # 4. internal consistency
    ok = all(Ns % orb == 0 and (a + b + c - n) % 2 == 0 for n in D.A for (a, b, c), (Ns, orb) in D.A[n].items())
    ok = ok and all(sum(Ns for Ns, _ in D.A[n].values()) == D.c[n] for n in D.A)
    ok = ok and all(sum(Ns for Ns, _ in D.H.get(n, {}).values()) == D.h[n] for n in range(D.N + 1))
    ok = ok and all(sum(Ns for Ns, _ in D.B.get(n, {}).values()) == D.b[n] for n in range(D.N + 1))
    ok = ok and all(Ns % orb == 0 for T in (D.H, D.B) for n in T for (Ns, orb) in T[n].values())
    ok = ok and all(sum(Ns for Ns, _ in D.M[(k, n)].values()) == D.c[n] for (k, n) in D.M)
    ok = ok and all(Ns % orb == 0 and (a + b + c - k) % 2 == 0 for (k, n) in D.M
                    for (a, b, c), (Ns, orb) in D.M[(k, n)].items())
    ok = ok and len(D.M) == D.N * (D.N - 1) // 2
    ok = ok and any("orbit divisibility check: PASS" in h for h in D.header)
    report(ok, "[4] data file internally consistent (orbit divisibility, parity, A/H/B/M totals) for n <= %d" % D.N)

    # 5. published values
    ref = {int(r[0]): r for r in read_reference(os.path.join(ROOT, "data", "reference_series.tsv"))}
    nmax = min(D.N, max(ref))
    report(all(D.c[n] == int(ref[n][1]) for n in range(nmax + 1)),
           "[5a] c_n = OEIS A001412 (= SBB 2011 Z_N = CLS 2007 c_n) for n <= %d" % nmax)
    report(all(D.S(n) == int(ref[n][2]) for n in range(1, nmax + 1)),
           "[5b] S_n = sum |omega_n|^2 = SBB 2011 Table I P_N (= CLS rho_n) for n <= %d" % nmax)
    okp, cnt = True, 0
    for n in range(3, D.N + 1):  # a polygon has at least 4 steps; for m = 2 the 'return' 0 -> e1 -> 0 is not one
        m = n + 1
        if m in ref and ref[m][3] not in ("-",):
            val = 3 * D.cnx(n, (1, 0, 0))
            okp = okp and val % m == 0 and val // m == int(ref[m][3])
            cnt += 1
    report(okp and cnt > 0, "[5c] polygons p_{n+1} = 3 c_n(e1)/(n+1) = CLS poly_d3 for 4 <= n+1 <= %d (%d values, odd ones 0)"
           % (D.N + 1, cnt))
    rb = {int(r[0]): r for r in read_reference(os.path.join(ROOT, "data", "reference_bridges.tsv"))}
    nb2 = min(D.N, 28)
    report(all(D.b[n] == int(rb[n][1]) for n in range(1, nb2 + 1)),
           "[5d] bridges b_n = Clisby-Conway-Guttmann 2016 Table 6 for 1 <= n <= %d" % nb2)
    nh = min(D.N, 27)
    report(all(D.h[n] == int(rb[n - 1][2]) for n in range(1, nh + 1)),
           "[5e] half-space walks h_n = OEIS A116904(n-1) (terminally attached walks) for 1 <= n <= %d" % nh)
    print("runtime %.1f s; %s" % (time.time() - t0, "ALL PASS" if fails == 0 else "%d FAIL" % fails))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
