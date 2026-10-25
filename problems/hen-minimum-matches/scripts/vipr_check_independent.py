#!/usr/bin/env python3
"""Independent checker for VIPR certificates of infeasibility ("RTP infeas"), written from the VIPR 1.0/1.1
specification (github.com/scipopt/vipr, cert_spec_v1_1.md) without using or translating viprchk's code.

usage: python3 vipr_check_independent.py <certificate.vipr>

Semantics checked, in exact rational arithmetic (gmpy2.mpq if available, else fractions.Fraction):
  * problem section: variables (VAR), integer variables (INT), objective (OBJ), constraints (CON); every
    constraint is a_x (sense) b with sense G (>=), L (<=) or E (=);
  * each derived constraint (DER), in order, with its set of assumptions A(.):
      asm            : A = {itself}; anything may be assumed;
      lin  p (i, l)* : the combination  sum_i l_i C_i  must DOMINATE the stated constraint, where every term has a
                       consistent direction (l_i > 0 on G or l_i < 0 on L gives >=, the opposite gives <=, E rows
                       take any sign; mixing >= and <= terms is rejected); only earlier indices may be used;
                       A = union of A(C_i);
      rnd  p (i, l)* : as lin, then Chvatal-Gomory rounding of the right-hand side (ceil for >=, floor for <=);
                       REQUIRED: every nonzero coefficient of the combination is an integer and belongs to an
                       integer variable (rounding is invalid otherwise);
      uns  i1 l1 i2 l2: C_i1 and C_i2 both dominate the stated constraint, and {C_l1, C_l2} is an integer
                       disjunction  a x <= b  /  a x >= b + 1  with b an INTEGER, a integral, support of a integer
                       variables;  A = (A(C_i1) - {l1}) u (A(C_i2) - {l2});
      sol            : rejected (no primal solutions are accepted by this checker);
      'weak' / 'incomplete' derivations: rejected (complete them first; viprcomp output is accepted).
    "dominates": either the derived constraint is a falsehood (0 >= b > 0, 0 <= b < 0, 0 = b != 0), or it has
    exactly the same coefficient vector and a right-hand side at least as strong in the stated direction
    (an equality dominates either inequality; only an equality with equal rhs dominates an equality).
  * success iff some derived constraint is a falsehood with an empty set of assumptions, and RTP is infeas.
Anything after the reference index on a derivation line (e.g. SCIP's 'global' flag) is ignored, i.e. it is never
used to drop assumptions (conservative).  Prints PASS/FAIL; exit code 1 on failure.
Environment: Python >= 3.10, gmpy2 (optional, for speed).  Runtime: about 15-30 min per GB of certificate.
"""
import sys
import time

try:
    from gmpy2 import mpq as Q
except ImportError:  # pragma: no cover
    from fractions import Fraction as Q

SENSE = {"G": 1, "L": -1, "E": 0}


class Fail(Exception):
    pass


def tokens(path):
    """yield (token, line_number); comment lines (first non-blank character '%') are skipped"""
    with open(path) as fh:
        for ln, line in enumerate(fh):
            s = line.lstrip()
            if not s or s[0] == "%":
                continue
            for t in s.split():
                yield t, ln


class Stream:
    def __init__(self, path):
        self.it = tokens(path)
        self.buf = None
        self.line = -1

    def next(self):
        if self.buf is not None:
            t, self.buf = self.buf, None
            return t
        t, self.line = next(self.it)
        return t

    def skip_rest_of_line(self):
        ln = self.line
        while True:
            try:
                t, l2 = next(self.it)
            except StopIteration:
                return
            if l2 != ln:
                self.buf, self.line = t, l2
                return


def is_int(q):
    return Q(q).denominator == 1


def floor_q(q):
    q = Q(q)
    return Q(q.numerator // q.denominator)


def ceil_q(q):
    q = Q(q)
    return Q(-((-q.numerator) // q.denominator))


def is_falsehood(sense, rhs, coef):
    if any(c != 0 for c in coef.values()):
        return False
    return (sense >= 0 and rhs > 0) or (sense <= 0 and rhs < 0)


def same_coef(a, b):
    a = {k: v for k, v in a.items() if v != 0}
    b = {k: v for k, v in b.items() if v != 0}
    return a == b


def dominates(d, s):
    """d, s = (sense, rhs, coef): does derived d imply stated s?"""
    ds, dr, dc = d
    ss, sr, sc = s
    if is_falsehood(ds, dr, dc):
        return True
    if not same_coef(dc, sc):
        return False
    if ss > 0:
        return ds >= 0 and dr >= sr
    if ss < 0:
        return ds <= 0 and dr <= sr
    return ds == 0 and dr == sr


def main(path):
    t0 = time.time()
    S = Stream(path)

    def expect(word):
        t = S.next()
        if t != word:
            raise Fail(f"expected {word}, read {t}")

    expect("VER")
    ver = S.next()
    if ver not in ("1.0", "1.1"):
        raise Fail(f"unsupported version {ver}")
    expect("VAR")
    nv = int(S.next())
    for _ in range(nv):
        S.next()
    expect("INT")
    ni = int(S.next())
    isint = [False] * nv
    for _ in range(ni):
        i = int(S.next())
        if not 0 <= i < nv:
            raise Fail("integer index out of range")
        isint[i] = True
    expect("OBJ")
    osense = S.next()
    if osense not in ("min", "max"):
        raise Fail("objective sense")
    obj = {}
    k = int(S.next())
    for _ in range(k):
        i = int(S.next()); c = Q(S.next())
        if not 0 <= i < nv:
            raise Fail("objective index out of range")
        obj[i] = obj.get(i, Q(0)) + c

    def read_coefs():
        t = S.next()
        if t == "OBJ":
            return dict(obj)
        n = int(t)
        coef = {}
        for _ in range(n):
            i = int(S.next()); c = Q(S.next())
            if not 0 <= i < nv:
                raise Fail("variable index out of range")
            coef[i] = coef.get(i, Q(0)) + c
        return coef

    def read_con():
        label = S.next()
        sense_t = S.next()
        if sense_t not in SENSE:
            raise Fail(f"bad sense {sense_t} in {label}")
        rhs = Q(S.next())
        return label, SENSE[sense_t], rhs, read_coefs()

    expect("CON")
    m = int(S.next()); S.next()  # number of bound constraints: informational
    cons = []   # list of (sense, rhs, coef, assumptions:frozenset)
    for _ in range(m):
        label, sense, rhs, coef = read_con()
        cons.append((sense, rhs, coef, frozenset()))
    expect("RTP")
    rtp = S.next()
    if rtp != "infeas":
        raise Fail("this checker only handles RTP infeas")
    expect("SOL")
    if int(S.next()) != 0:
        raise Fail("solutions present (not expected for RTP infeas)")
    expect("DER")
    nd = int(S.next())
    counts = {"asm": 0, "lin": 0, "rnd": 0, "uns": 0}
    nmult = 0
    for di in range(nd):
        idx = len(cons)
        label, sense, rhs, coef = read_con()
        expect("{")
        kind = S.next()
        if kind == "asm":
            expect("}")
            A = frozenset([idx])
        elif kind in ("lin", "rnd"):
            t = S.next()
            if t in ("weak", "incomplete"):
                raise Fail(f"{label}: incomplete derivation ({t})")
            p = int(t)
            agg, arhs, dsense, A = {}, Q(0), 0, set()
            for _ in range(p):
                i = int(S.next()); lam = Q(S.next())
                if not 0 <= i < idx:
                    raise Fail(f"{label}: reference {i} out of range")
                if lam == 0:
                    continue
                cs, cr, cc, ca = cons[i]
                d = cs * (1 if lam > 0 else -1)
                if d != 0:
                    if dsense == 0:
                        dsense = d
                    elif dsense != d:
                        raise Fail(f"{label}: inconsistent multiplier sign for {i}")
                for v, c in cc.items():
                    agg[v] = agg.get(v, Q(0)) + lam * c
                arhs += lam * cr
                A |= ca
                nmult += 1
            expect("}")
            if kind == "rnd":
                for v, c in agg.items():
                    if c != 0 and (not isint[v] or not is_int(c)):
                        raise Fail(f"{label}: rounding with non-integer coefficient or continuous variable {v}")
                if dsense > 0:
                    arhs = ceil_q(arhs)
                elif dsense < 0:
                    arhs = floor_q(arhs)
            if not dominates((dsense, arhs, agg), (sense, rhs, coef)):
                raise Fail(f"{label}: derived constraint does not dominate the stated one")
            A = frozenset(A)
        elif kind == "uns":
            i1, l1, i2, l2 = (int(S.next()) for _ in range(4))
            expect("}")
            for j in (i1, l1, i2, l2):
                if not 0 <= j < idx:
                    raise Fail(f"{label}: unsplit index out of range")
            for j in (i1, i2):
                cs, cr, cc, ca = cons[j]
                if not dominates((cs, cr, cc), (sense, rhs, coef)):
                    raise Fail(f"{label}: unsplit premise {j} does not dominate")
            s1, r1, c1, _ = cons[l1]
            s2, r2, c2, _ = cons[l2]
            if not same_coef(c1, c2):
                raise Fail(f"{label}: disjunction coefficients differ")
            if (s1, s2) == (-1, 1):
                ok = r2 == r1 + 1
            elif (s1, s2) == (1, -1):
                ok = r1 == r2 + 1
            else:
                ok = False
            if not ok or not is_int(r1) or not is_int(r2):
                raise Fail(f"{label}: not an integer disjunction a x <= b | a x >= b+1 with b integer")
            for v, c in c1.items():
                if c != 0 and (not isint[v] or not is_int(c)):
                    raise Fail(f"{label}: disjunction on continuous variable or fractional coefficient")
            A = (cons[i1][3] - {l1}) | (cons[i2][3] - {l2})
        else:
            raise Fail(f"{label}: derivation type {kind} not accepted")
        counts[kind] += 1
        S.next()                 # reference index (memory management only)
        S.skip_rest_of_line()    # ignore trailing flags such as 'global'
        cons.append((sense, rhs, coef, A))
        if not A and is_falsehood(sense, rhs, coef):
            print(f"certificate {path.split('/')[-1]}: VER {ver}, {nv} variables ({ni} integer), {m} problem "
                  f"constraints, {di + 1} of {nd} derivations checked ({counts}), {nmult} multipliers, "
                  f"{time.time() - t0:.0f} s")
            print("PASS: falsehood derived without assumptions: the certificate's problem is infeasible")
            return 0
    raise Fail("no assumption-free falsehood derived")


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1]))
    except (Fail, StopIteration, ValueError) as e:
        print(f"FAIL: {e!r}")
        sys.exit(1)
