"""Exact algebra checks for attainment.md; Python 3 + SymPy; normally < 15 s.

Run from any directory. No files are written. At most two computational threads.
PASS means the stated algebraic check, not finite-stage FTC attainment.
Attainment: GPT-6-Astra (effort xhigh).
"""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "2"
import sympy as s


def check(label, condition):
    print(("PASS " if condition else "FAIL ") + label, flush=True)
    if not condition:
        raise AssertionError(label)


def prod(xs):
    return s.prod(list(xs))


def main():
    t = s.symbols("t")
    for n in range(3, 9):
        a = [s.Integer(2) ** (n - 1 - i) for i in range(n)]
        th = [s.Rational(3, 4) * a[i] for i in range(n - 1)]
        beta = prod(a[1:-1]) / prod(th[1:-1])
        K = (t - beta) * prod(t - x for x in th[1:-1]) / prod(t - x for x in a[1:-1])
        h = [s.S(0)] * n
        for i in range(1, n - 1):
            h[i] = ((a[i] - beta) * prod(a[i] - x for x in th[1:-1])
                    / (a[i] * prod(a[i] - a[j] for j in range(1, n - 1) if j != i)))
        check(f"N={n}: slack partial fractions and zero total redistribution",
              s.cancel(K - (1 - sum(a[i] * h[i] / (a[i] - t) for i in range(n)))) == 0
              and sum(h) == 0)
        check(f"N={n}: middle roots retained, endpoint perturbations positive",
              all(s.cancel(K).subs(t, x) == 0 for x in th[1:-1])
              and s.cancel(K).subs(t, th[0]) > 0
              and s.cancel(K).subs(t, th[-1]) > 0)

        # V_T = V_B = 1, hence q=1 after feed normalization.
        d = [s.S(0)] * n
        b = [s.S(0)] * n
        for i in range(n - 1):
            d[i] = prod(a[i] - x for x in th) / (
                a[i] * prod(a[i] - a[j] for j in range(n - 1) if j != i))
        for i in range(1, n):
            b[i] = -prod(a[i] - x for x in th) / (
                a[i] * prod(a[i] - a[j] for j in range(1, n) if j != i))
        check(f"N={n}: positive preferred net flows",
              all(v > 0 for v in d[:-1] + b[1:]) and sum(d) < 1)
        psi_t = 1 - sum(a[i] * d[i] / (a[i] - t) for i in range(n))
        psi_b = 1 + sum(a[i] * b[i] / (a[i] - t) for i in range(n))
        check(f"N={n}: preferred roots and feed equation",
              all(s.cancel(psi_t).subs(t, x) == 0 and s.cancel(psi_b).subs(t, x) == 0 for x in th))

        eta = s.Rational(1, 1000)
        de = s.Matrix([d[i] + eta * h[i] for i in range(n)])
        be = s.Matrix([b[i] - eta * h[i] for i in range(n)])
        av = s.Matrix(a)
        V = 1 + eta
        LT, LB = V - sum(de), V + sum(be)
        MT, MB = V * s.diag(*a) - de * av.T, V * s.diag(*a) + be * av.T
        ec = [a[-1] * de[i] / (LT * (a[i] - a[-1])) for i in range(n - 1)]
        ec += [1 - sum(ec)]
        ea = [a[0] * be[i] / (LB * (a[0] - a[i])) for i in range(1, n)]
        ea = [1 - sum(ea)] + ea
        check(f"N={n}: positive missing-component pinch eigenvectors",
              min(ec + ea) > 0
              and MT * s.Matrix(ec) == V * a[-1] * s.Matrix(ec)
              and MB * s.Matrix(ea) == V * a[0] * s.Matrix(ea))
        # An exact identity valid also at a zero-net-flow component's volatility.
        G = sum(a[i] * ec[i] / (a[i] - t) for i in range(n))
        psi = V - sum(a[i] * de[i] / (a[i] - t) for i in range(n))
        S = sum(a[i] * ec[i] for i in range(n))
        check(f"N={n}: pinch flux identity including the absent-component eigenvalue",
              s.cancel(psi - (LT - V * t / S) * G) == 0)

        # Cauchy interpolation: any interlacing list determines a positive composition.
        chi = [th[0] + s.Rational(1, 20)] + th[1:-1] + [th[-1] - s.Rational(1, 20)]
        weights = [prod(a[i] - x for x in chi) /
                   (a[i] * prod(a[i] - a[j] for j in range(n) if j != i)) for i in range(n)]
        xx = s.Matrix([v / sum(weights) for v in weights])
        J = s.Matrix([[1] * n] + [[a[i] / (a[i] - r) for i in range(n)] for r in chi])
        check(f"N={n}: positive feed matching with nonsingular Cauchy matrix",
              min(xx) > 0 and J * xx == s.Matrix([1] + [0] * (n - 1)) and J.det() != 0)

    # Frozen-net-flow ternary endpoint conditions. Six nodes, two simplex
    # coordinates each: F,U,D,A,B,C. This is deliberately NOT the full MESH Jacobian.
    a = [s.S(4), s.S(2), s.S(1)]
    r, ss = s.S(3), s.Rational(6, 5)
    vals = s.symbols("f1 f2 u1 u2 d1 d2 a1 a2 b1 b2 c1 c2")
    nodes = [s.Matrix([vals[2*i], vals[2*i+1], 1-vals[2*i]-vals[2*i+1]]) for i in range(6)]
    F, U, D, A, B, C = nodes
    G = lambda x, rr: sum(a[i]*x[i]/(a[i]-rr) for i in range(3))
    rows = [G(F, r), U[2], D[0], G(F, ss), A[1], A[2],
            B[0], U[2], D[0], B[2], C[0], C[1]]
    rank = s.Matrix(rows).jacobian(vals).rank()
    check("frozen-net-flow ternary limiting endpoint matrix has rank 10, not 12", rank == 10)
    print("PASS all exact algebra checks; no claim of full-network attainment")


if __name__ == "__main__":
    main()
