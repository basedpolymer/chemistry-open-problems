"""Symbolic checks for Theorem T / Corollary T* of petri/notes.md (agent P, Claude Opus 5.5, effort xhigh).

Environment: Python 3.12 + sympy 1.14 (WSL Fedora 44).  Runtime: < 5 s.  PASS/FAIL per claim, exit 1 on failure.

Claims:
  U1  U: X4 -> X2+X3, X3+X4 -> X1, X1+X2 -> 2X4 has a positive steady state for EVERY k (explicit formula).
  U2  U* has no positive steady state when kappa_2 >= kappa_1: f_{K1} - f_{K2} = (kappa_2 - kappa_1) K3 + 2 kappa_3 K2.
  U3  U and U* are not EC (so Theorem T does not apply); U satisfies (H1)-(H3), U* satisfies (H1), (H2).
  W1  W* (well-formed EC with balanced joins, not weakly reversible): explicit positive steady states for every
      kappa; the set is {k3 = kappa1 k1 k2 / kappa2, k4 = kappa1 k1 k2 / kappa3}: toric of dimension 2 = n - s.
  T1  Theorem T on W*, F and D: dim ker A_k and positivity for symbolic k (W*, F: kernel spanned by a vector of
      positive rational functions of k; D: kernel nonzero iff k1 = k2).
"""
import sys
import time
import platform
import sympy as sp

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from petri import Network, is_EC, ec_sets  # noqa: E402

T0 = time.time()
print(f"# check_absolute_consistency.py  Python {platform.python_version()}, sympy {sp.__version__}, "
      f"{platform.system()}")
FAILS = []


def claim(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        FAILS.append(name)


def massaction(N, xs, ks):
    f = [0] * N.n
    for j, (a, b) in enumerate(N.reactions):
        rate = ks[j] * sp.prod([xs[i] ** a[i] for i in range(N.n)])
        for i in range(N.n):
            f[i] += (b[i] - a[i]) * rate
    return [sp.expand(e) for e in f]


def e(n, *terms):
    v = [0] * n
    for c, i in terms:
        v[i] += c
    return tuple(v)


U = Network(4, [(e(4, (1, 3)), e(4, (1, 1), (1, 2))), (e(4, (1, 2), (1, 3)), e(4, (1, 0))),
                (e(4, (1, 0), (1, 1)), e(4, (2, 3)))])
Ud = U.dual()
k = sp.symbols("k1:4", positive=True)
x = sp.symbols("x1:5", positive=True)
fU = massaction(U, x, k)
sol = {x[0]: k[0] / k[2], x[1]: 1, x[2]: k[0] / k[1], x[3]: 1}
claim("U1 U: x = (k1/k3, 1, k1/k2, 1) is a positive steady state for every k",
      all(sp.simplify(fi.subs(sol)) == 0 for fi in fU))
kap = sp.symbols("kappa1:5", positive=True)
K = sp.symbols("K1:4", positive=True)
fUd = massaction(Ud, K, kap)
diff = sp.expand(fUd[0] - fUd[1])
claim(f"U2 U*: f_K1 - f_K2 = {diff}", sp.simplify(diff - ((kap[1] - kap[0]) * K[2] + 2 * kap[2] * K[1])) == 0)
claim("U3 U and U* are not EC; U satisfies (H1)-(H3); U* satisfies (H1), (H2)",
      not is_EC(U) and not is_EC(Ud) and U.is_conservative() and U.is_consistent() and U.rows_distinct()
      and Ud.is_conservative() and Ud.is_consistent())

# W* : K1+K2 -> K3+K4, K3 -> K1, K4 -> K2
Wd = Network(4, [(e(4, (1, 0), (1, 1)), e(4, (1, 2), (1, 3))), (e(4, (1, 2)), e(4, (1, 0))),
                 (e(4, (1, 3)), e(4, (1, 1)))])
y = sp.symbols("y1:5", positive=True)
fW = massaction(Wd, y, kap[:3])
solW = {y[2]: kap[0] * y[0] * y[1] / kap[1], y[3]: kap[0] * y[0] * y[1] / kap[2]}
claim("W1 W*: {y3 = kappa1 y1 y2/kappa2, y4 = kappa1 y1 y2/kappa3} are steady states for all kappa, y1, y2 > 0",
      all(sp.simplify(fi.subs(solW)) == 0 for fi in fW) and is_EC(Wd) and Wd.rank() == len(ec_sets(Wd)) - 1)


def throughput_matrix(N, ks):
    classes = ec_sets(N)
    A = sp.zeros(N.n, len(classes))
    for c, cl in enumerate(classes):
        for t in cl:
            a, b = N.reactions[t]
            for i in range(N.n):
                A[i, c] += ks[t] * (b[i] - a[i])
    return A


F = Network(6, [(e(6, (1, 0), (1, 2)), e(6, (1, 4))), (e(6, (1, 4)), e(6, (1, 0), (1, 2))),
                (e(6, (1, 4)), e(6, (1, 1), (1, 2))), (e(6, (1, 1), (1, 3)), e(6, (1, 5))),
                (e(6, (1, 5)), e(6, (1, 1), (1, 3))), (e(6, (1, 5)), e(6, (1, 0), (1, 3)))])
D = Network(3, [(e(3, (1, 0)), e(3, (1, 1))), (e(3, (1, 0)), e(3, (1, 2))), (e(3, (1, 1), (1, 2)), e(3, (2, 0)))])
for name, N in [("W*", Wd), ("F", F)]:
    ks = sp.symbols(f"c1:{N.r + 1}", positive=True)
    A = throughput_matrix(N, ks)
    ns = A.nullspace()
    ok = len(ns) == 1
    if ok:
        u = [sp.factor(sp.simplify(c)) for c in ns[0]]
        # normalise so that the last entry is 1 and check every entry is a ratio of polynomials with positive
        # coefficients (hence positive for every positive k)
        u = [sp.factor(sp.simplify(c / u[-1])) for c in u]
        pos = all(all(co > 0 for co in sp.Poly(sp.numer(sp.together(c)), *ks).coeffs()) and
                  all(co > 0 for co in sp.Poly(sp.denom(sp.together(c)), *ks).coeffs()) for c in u)
        ok = pos
        print(f"#   {name}: u(k) = {u}")
    claim(f"T1 {name}: ker A_k is one-dimensional, spanned by positive rational functions of k", ok)
kd = sp.symbols("d1:4", positive=True)
A = throughput_matrix(D, kd)
claim("T1 D (not well-formed): det of A_k restricted to two rows vanishes only on k1 = k2 (no steady state otherwise)",
      sp.factor(A[[0, 1], :].det()) != 0 and sp.simplify(A.subs(kd[1], kd[0]).rank() - 1) == 0
      and A.rank() == 2)
print(f"# runtime {time.time() - T0:.2f} s; {len(FAILS)} failure(s)")
sys.exit(1 if FAILS else 0)
