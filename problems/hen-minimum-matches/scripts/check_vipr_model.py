#!/usr/bin/env python3
"""Check that the problem section of a VIPR certificate is (a relaxation of) the MILP in a given LP file.

usage: python3 check_vipr_model.py <model.lp> <certificate.vipr> <expected RTP: 'infeas' | 'range L U'>

Verified exactly (fractions.Fraction):
  * every VIPR variable name is 't_'+name or name of an LP variable;
  * every VIPR integer variable is declared binary in the LP file;
  * every VIPR constraint (CON section, including bound constraints) equals, up to a positive rational factor,
    an LP row (same sense; '=' rows also match 'E') or a variable bound / binary bound of the LP file;
    hence the certificate's problem is a relaxation of (or equal to) the LP-file MILP, and any infeasibility or
    lower bound it proves holds for the LP-file MILP;
  * for RTP 'range', the VIPR objective equals the LP objective;
  * the RTP line equals the expected one.
Only the header sections are read (streamed), so large certificates are fine.
Environment: Python >= 3.10, standard library only.
"""
import sys, re
from fractions import Fraction


def parse_lp(path):
    rows, bounds, binaries, obj = [], {}, set(), {}
    sec = None
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("\\"):
            continue
        if s in ("Minimize", "Subject To", "Bounds", "Binaries", "End"):
            sec = s
            continue
        if sec == "Minimize":
            name, expr = s.split(":", 1)
            obj = parse_terms(expr)
        elif sec == "Subject To":
            name, rest = s.split(":", 1)
            mo = re.match(r"(.*)\s(<=|>=|=)\s(\S+)$", rest.strip())
            expr, sense, rhs = mo.group(1), mo.group(2), Fraction(mo.group(3))
            rows.append((name.strip(), parse_terms(expr), {"<=": "L", ">=": "G", "=": "E"}[sense], rhs))
        elif sec == "Bounds":
            mo = re.match(r"(\S+)\s*<=\s*(\S+)\s*<=\s*(\S+)$", s)
            bounds[mo.group(2)] = (Fraction(mo.group(1)), Fraction(mo.group(3)))
        elif sec == "Binaries":
            binaries.add(s)
    return rows, bounds, binaries, obj


def parse_terms(expr):
    toks = expr.split()
    coef = {}
    sign, val = 1, None
    i = 0
    while i < len(toks):
        t = toks[i]
        if t in "+-":
            sign = 1 if t == "+" else -1
            i += 1
            continue
        try:
            val = Fraction(t)
            v = toks[i + 1]
            i += 2
        except (ValueError, IndexError):
            val, v = Fraction(1), t
            i += 1
        coef[v] = coef.get(v, 0) + sign * val
        sign = 1
    return {k: c for k, c in coef.items() if c != 0}


def normalize(coef, sense, rhs):
    """canonical key: divide by |first coefficient| (positive factor), sorted items"""
    items = sorted(coef.items())
    f = abs(items[0][1])
    return (tuple((k, c / f) for k, c in items), sense, rhs / f)


def main():
    lp, cert, expected = sys.argv[1], sys.argv[2], sys.argv[3:]
    rows, bounds, binaries, obj = parse_lp(lp)
    allowed = set()
    for name, coef, sense, rhs in rows:
        allowed.add(normalize(coef, sense, rhs))
        if sense == "E":  # an equality may also be split into two inequalities
            allowed.add(normalize(coef, "L", rhs))
            allowed.add(normalize(coef, "G", rhs))
    lpvars = set(bounds) | binaries | {v for r in rows for v in r[1]}
    for v in lpvars:
        lo, up = bounds.get(v, (Fraction(0), None))
        if v in binaries:
            lo, up = Fraction(0), Fraction(1)
        allowed.add(normalize({v: Fraction(1)}, "G", lo))
        if up is not None:
            allowed.add(normalize({v: Fraction(1)}, "L", up))
    errs = []
    with open(cert) as f:
        tok = iter(" ".join(f.readline() for _ in range(0)).split())  # placeholder

    # streaming tokenizer
    def tokens(fh):
        for line in fh:
            if line.lstrip().startswith("%"):
                continue
            for t in line.split():
                yield t
    fh = open(cert)
    T = tokens(fh)
    assert next(T) == "VER"; ver = next(T)
    assert next(T) == "VAR"; nv = int(next(T))
    names = [next(T) for _ in range(nv)]
    base = [nm[2:] if nm.startswith("t_") and nm[2:] in lpvars else nm for nm in names]
    for nm in base:
        if nm not in lpvars:
            errs.append(f"unknown variable {nm}")
    assert next(T) == "INT"; ni = int(next(T))
    ints = [int(next(T)) for _ in range(ni)]
    for i in ints:
        if base[i] not in binaries:
            errs.append(f"integer variable {base[i]} not binary in LP")
    assert next(T) == "OBJ"; osense = next(T)
    k = int(next(T)); vobj = {}
    for _ in range(k):
        i = int(next(T)); c = Fraction(next(T)); vobj[base[i]] = c
    assert next(T) == "CON"; nc = int(next(T)); nb = int(next(T))
    for _ in range(nc):
        label = next(T); sense = next(T); rhs = Fraction(next(T)); k = next(T)
        if k == "OBJ":
            coef = dict(vobj)
        else:
            coef = {}
            for _ in range(int(k)):
                i = int(next(T)); c = Fraction(next(T)); coef[base[i]] = coef.get(base[i], 0) + c
        coef = {a: c for a, c in coef.items() if c != 0}
        if not coef:
            # 0 (sense) rhs: must be trivially true
            ok = (sense == "L" and rhs >= 0) or (sense == "G" and rhs <= 0) or (sense == "E" and rhs == 0)
            if not ok:
                errs.append(f"{label}: empty constraint violated")
            continue
        if normalize(coef, sense, rhs) not in allowed:
            errs.append(f"{label}: constraint not found in LP model")
    assert next(T) == "RTP"
    rtp = [next(T)]
    if rtp[0] == "range":
        rtp += [next(T), next(T)]
    fh.close()
    if rtp != expected:
        errs.append(f"RTP {rtp} != expected {expected}")
    if rtp[0] == "range":
        if osense != "min" or {a: c for a, c in vobj.items() if c} != obj:
            errs.append("objective differs from LP objective")
    print(f"certificate {cert.split('/')[-1]}: VER {ver}, {nv} vars, {ni} integer, {nc} constraints ({nb} bounds), RTP {' '.join(rtp)}")
    print(f"LP model {lp.split('/')[-1]}: {len(rows)} rows, {len(lpvars)} variables, {len(binaries)} binaries")
    if errs:
        print("FAIL:", "; ".join(errs[:10]), f"({len(errs)} errors)")
        sys.exit(1)
    print("PASS: certificate problem = relaxation of the LP-file MILP (every constraint matched exactly)")


if __name__ == "__main__":
    main()
