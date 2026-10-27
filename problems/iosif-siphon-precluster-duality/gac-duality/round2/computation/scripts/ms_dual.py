"""Q3-Udual: is the dual N* of an unbalanced WF-EC network N ever multistationary?  Exact search.

Environment: Python 3.12 (WSL Fedora 44), standard library + msolve 0.10.1 (certified real root isolation; WSL).
Usage (from round2/computation/): python scripts/ms_dual.py N MAXMOL RMAX [TRIALS_INV] [TRIALS_MSOLVE] [NPROC]
Networks: agent P's generator (petri/scripts/enumerate_ec.gen_networks, imported unchanged), EC, connected, (H1),
(H2), s = m - 1 (WF-EC) and with unbalanced joins (some source has two different nonzero coefficients).

N* has species k_1..k_r and reactions i = 1..n (sources = rows Ye_i of Ye(N)); f*(x) = Gamma^T w(x),
w_i = kappa_i x^{Ye_i}.  Two methods, both exact:
  (inverse) choose x2 > 0 rational and x1 = x2 + d, d in im Gamma^T (same class of N*), rho_i = (x1/x2)^{Ye_i};
            both are steady states for one kappa iff some w > 0 has Gamma^T w = 0 and Gamma^T (rho o w) = 0 (exact
            LP); then kappa_i = w_i / x2^{Ye_i}, and f*(x1) = f*(x2) = 0 is verified exactly;
  (forward) random rational kappa and class x0; msolve isolates all complex-real solutions of the square system
            {s independent components of f*, conservation laws of N*}; positive boxes are counted (boxes meeting a
            coordinate hyperplane are refined, then reported as ambiguous if still unresolved).
Claims: I every inverse certificate verifies exactly; F msolve counts are consistent with the inverse certificates
(for each certified kappa and class, msolve finds >= 2 positive solutions).  Report: networks, certificates.
"""
import os
import re
import ast
import sys
import time
import random
import tempfile
import platform
import subprocess
from fractions import Fraction
from itertools import permutations
from multiprocessing import Pool
from math import lcm

import e2lib
from e2lib import CNet, crn

sys.path.insert(0, e2lib.PSCRIPTS)
import enumerate_ec as EE  # noqa: E402


def mono(x, e):
    out = Fraction(1)
    for xi, ei in zip(x, e):
        if ei:
            out *= Fraction(xi) ** ei
    return out


def fstar(n, rx, kappa, x):
    r = len(rx)
    Ye = [[rx[j][0][i] for j in range(r)] for i in range(n)]
    G = CNet(n, rx).Gamma
    w = [kappa[i] * mono(x, Ye[i]) for i in range(n)]
    return [sum(G[i][j] * w[i] for i in range(n)) for j in range(r)]


def inverse_trial(n, rx, rng):
    r = len(rx)
    G = CNet(n, rx).Gamma
    Ye = [[rx[j][0][i] for j in range(r)] for i in range(n)]
    x2 = [Fraction(rng.randint(1, 12), rng.randint(1, 12)) for _ in range(r)]
    lam = [rng.randint(-3, 3) for _ in range(n)]
    d = [sum(G[i][j] * lam[i] for i in range(n)) for j in range(r)]
    if not any(d):
        return None
    neg = [x2[j] / -d[j] for j in range(r) if d[j] < 0]
    tmax = min(neg) if neg else Fraction(5)
    t = tmax * Fraction(rng.randint(1, 99), 100)
    x1 = [x2[j] + t * d[j] for j in range(r)]
    rho = [mono(x1, Ye[i]) / mono(x2, Ye[i]) for i in range(n)]
    eqs = [([G[i][j] for i in range(n)], 0) for j in range(r)]
    eqs += [([G[i][j] * rho[i] for i in range(n)], 0) for j in range(r)]
    ineqs = [([-int(k == i) for k in range(n)], -1) for i in range(n)]
    w = e2lib.lp_point(eqs, ineqs, n)
    if w is None:
        return None
    kappa = [w[i] / mono(x2, Ye[i]) for i in range(n)]
    ok = all(v == 0 for v in fstar(n, rx, kappa, x1)) and all(v == 0 for v in fstar(n, rx, kappa, x2))
    return dict(x1=x1, x2=x2, kappa=kappa, ok=ok)


def polystr(terms):
    """msolve syntax: no parentheses, explicit signs ('+3*x1^2*x2-5*x3+7')."""
    out = ""
    for c, m in terms:
        if c == 0:
            continue
        sg = "-" if c < 0 else "+"
        body = str(abs(c)) if m == "1" else f"{abs(c)}*{m}"
        out += sg + body
    out = out.lstrip("+")
    return out or "0"


def msolve_count(n, rx, kappa, x0, prec=64):
    """Number of positive real solutions of N* in the class of x0 (msolve).  Returns (count, nreal, status)."""
    r = len(rx)
    N = CNet(n, rx)
    G = N.Gamma
    s = N.rank()
    Ye = [[rx[j][0][i] for j in range(r)] for i in range(n)]
    # s independent columns of Gamma (= rows of Gamma^T)
    J = []
    for j in range(r):
        col = [G[i][j] for i in range(n)]
        if crn.rank([[G[i][jj] for i in range(n)] for jj in J] + [col]) > len(J):
            J.append(j)
    var = [f"x{j+1}" for j in range(r)]

    def monostr(e):
        parts = [f"{var[j]}^{e[j]}" if e[j] > 1 else var[j] for j in range(r) if e[j]]
        return "*".join(parts) if parts else "1"
    polys = []
    for j in J:
        coeffs = [G[i][j] * kappa[i] for i in range(n)]
        den = 1
        for c in coeffs:
            den = lcm(den, Fraction(c).denominator)
        terms = {}
        for i in range(n):
            c = Fraction(coeffs[i]) * den
            if c:
                key = tuple(Ye[i])
                terms[key] = terms.get(key, 0) + c
        polys.append(polystr([(int(c), monostr(e)) for e, c in terms.items() if c]))
    K = crn.nullspace([[Fraction(G[i][j]) for j in range(r)] for i in range(n)], r)  # conservation laws of N*
    for c in K:
        den = 1
        for v in list(c) + [sum(cj * xj for cj, xj in zip(c, x0))]:
            den = lcm(den, Fraction(v).denominator)
        cc = [int(v * den) for v in c]
        rhs = int(sum(Fraction(cj) * xj for cj, xj in zip(c, x0)) * den)
        polys.append(polystr([(cc[j], var[j]) for j in range(r) if cc[j]] + [(-rhs, "1")]))
    with tempfile.TemporaryDirectory() as td:
        fin = os.path.join(td, "in.ms")
        fout = os.path.join(td, "out.ms")
        with open(fin, "w") as fh:
            fh.write(",".join(var) + "\n0\n" + ",\n".join(polys) + "\n")
        subprocess.run(["msolve", "-f", fin, "-o", fout, "-p", str(prec)], check=True, capture_output=True)
        txt = open(fout).read().strip().rstrip(":")
    if txt.startswith("[-1]"):
        return 0, 0, "no complex solution"
    if txt.startswith("[1,"):
        return None, None, "positive-dimensional"
    # endpoints: integers, 'a / 2^k', or 'a / b' (exact rationals); all turned into strings 'a/b/k' = a/(b^k)
    txt = re.sub(r"(-?\d+)\s*/\s*(\d+)\^(-?\d+)", r"'\1/\2/\3'", txt)
    txt = re.sub(r"(-?\d+)\s*/\s*(\d+)(?![\d^/'])", r"'\1/\2/1'", txt)
    data = ast.literal_eval(txt)
    boxes = data[1][1]

    def val(t):
        if isinstance(t, str):
            a, b, k = t.split("/")
            return Fraction(int(a)) / Fraction(int(b)) ** int(k)
        return Fraction(t)
    pos = 0
    amb = 0
    for box in boxes:
        lo = [val(iv[0]) for iv in box]
        hi = [val(iv[1]) for iv in box]
        if all(l > 0 for l in lo):
            pos += 1
        elif any(h < 0 for h in hi):
            continue
        else:
            amb += 1
    return pos, len(boxes), ("ok" if amb == 0 else f"{amb} ambiguous")


def worker(args):
    n, m, rmax, pidx, tinv, tms, seed = args
    rng = random.Random(seed + pidx)
    perms = list(permutations(range(n)))
    seen = set()
    out = []
    for rx in EE.gen_networks(n, m, rmax, 0, pidx):
        rx = tuple(sorted(rx))
        can = EE.canonical(rx, perms)
        if can != rx or can in seen:
            continue
        seen.add(can)
        if all(len({x for x in a if x > 0}) <= 1 for a, b in rx):
            continue  # balanced joins
        if not e2lib.petri_connected(n, rx):
            continue
        N = CNet(n, rx)
        if not N.is_consistent() or not N.is_conservative() or N.rank() != e2lib.n_sources(rx) - 1:
            continue
        res = dict(rx=rx, n=n, cert=None, msolve=[])
        for _ in range(tinv):
            c = inverse_trial(n, rx, rng)
            if c is not None and c["x1"] != c["x2"]:
                res["cert"] = c
                break
        for _ in range(tms):
            r = len(rx)
            kappa = [Fraction(rng.randint(1, 20), rng.randint(1, 5)) for _ in range(n)]
            x0 = [Fraction(rng.randint(1, 9), rng.randint(1, 9)) for _ in range(r)]
            res["msolve"].append(msolve_count(n, rx, kappa, x0))
        if res["cert"] is not None and res["cert"]["ok"]:
            c = res["cert"]
            res["cert_ms"] = msolve_count(n, rx, c["kappa"], c["x2"])
        out.append(res)
    return out


def main():
    n, m, rmax = (int(a) for a in sys.argv[1:4])
    tinv = int(sys.argv[4]) if len(sys.argv) > 4 else 200
    tms = int(sys.argv[5]) if len(sys.argv) > 5 else 5
    nproc = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    t0 = time.time()
    print(f"# ms_dual.py n={n} molecularity<={m} r<={rmax}; inverse trials {tinv}, msolve samples {tms}; Python "
          f"{platform.python_version()} on {platform.system()}; msolve 0.10.1")
    nparts = sum(1 for _ in EE.set_partitions(list(range(n))))
    res = []
    with Pool(nproc) as pool:
        for chunk in pool.imap_unordered(worker, [(n, m, rmax, p, tinv, tms, 7) for p in range(nparts)]):
            res.extend(chunk)
    print(f"# unbalanced WF-EC networks: {len(res)}")
    # V: msolve interface sanity (X1 <-> X2: dual k2 -> k1, k1 -> k2 has exactly one positive steady state per class;
    # G's M = 3A <-> 2A+B: dual 3U+2V -> 2U+3V, V -> U has two positive steady states on u + v = 2 for unit rates)
    v1 = msolve_count(2, (((1, 0), (0, 1)), ((0, 1), (1, 0))), [Fraction(1), Fraction(2)], [Fraction(1), Fraction(1)])
    v2 = msolve_count(2, (((3, 0), (2, 1)), ((2, 1), (3, 0))), [Fraction(1), Fraction(1)], [Fraction(1), Fraction(1)])
    okV = v1[0] == 1 and v2[0] == 2
    print(f"{'PASS' if okV else 'FAIL'} V: msolve interface sanity (expected 1 and 2 positive solutions): {v1}, {v2}")
    certs = [q for q in res if q["cert"] is not None]
    badI = [q for q in certs if not q["cert"]["ok"]]
    print(f"{'PASS' if not badI else 'FAIL'} I: inverse certificates verified exactly: {len(certs) - len(badI)}"
          f"/{len(certs)} (networks with two positive steady states of N* in one class for some kappa)")
    badF = [q for q in certs if q["cert"]["ok"] and not (q.get("cert_ms") and q["cert_ms"][0] is not None
                                                          and q["cert_ms"][0] >= 2)]
    print(f"{'PASS' if not badF else 'FAIL'} F: msolve finds >= 2 positive solutions at every certificate "
          f"({len(certs) - len(badF)}/{len(certs)})")
    hist = {}
    stat = {}
    for q in res:
        for c, nr, st in q["msolve"]:
            hist[c] = hist.get(c, 0) + 1
            stat[st] = stat.get(st, 0) + 1
    print(f"# forward msolve samples: positive-solution counts {dict(sorted(hist.items(), key=str))}; status {stat}")
    certs.sort(key=lambda q: (len(q["rx"]), q["rx"]))
    for q in certs[:4]:
        c = q["cert"]
        print(f"  multistationary N*: N = {e2lib.describe(n, q['rx'])}")
        print(f"     N* = {e2lib.describe(len(q['rx']), e2lib.dual_rx(n, q['rx']), [f'k{j+1}' for j in range(len(q['rx']))])}")
        print(f"     kappa = {[str(v) for v in c['kappa']]}; x1 = {[str(v) for v in c['x1']]}; x2 = {[str(v) for v in c['x2']]}"
              f"; msolve in that class: {q.get('cert_ms')}")
    many = [q for q in res if any(c is not None and c >= 2 for c, nr, st in q["msolve"])]
    print(f"# networks with a forward sample having >= 2 positive solutions: {len(many)}")
    ok = not badI and not badF and okV
    print(f"{'PASS' if ok else 'FAIL'}: all claims")
    print(f"# runtime {time.time() - t0:.1f} s")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
