#!/usr/bin/env python3
"""Rigorous FINITE coefficient enclosures; no all-order claim.

Environment: Python 3, python-flint 0.9.0 (Arb). Typical runtime: 2-20 seconds.
One computational thread; typically below 250 MB RAM.
Run from this problem folder. Inputs and precision are fixed and reproducible.
Every claimed inequality is checked with Arb balls, never floating point.
"""
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key]="1"
import json
from pathlib import Path
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
PRECISION=4096
ORDER=256
MASSES=((1,2),(3,5),(3,4),(9,10),(99,100),(999,1000))


def main():
    ctx.prec=PRECISION
    ctx.threads=1
    records=[]
    for numerator,denominator in MASSES:
        m=arb(numerator)/denominator
        total=m*(1-m)
        tail=total
        rho=[arb(0)]
        partial_mass=arb(0)
        all_positive=True
        all_tails=True
        for n in range(1,ORDER+1):
            conv=sum((i*(n-i)*rho[i]*rho[n-i] for i in range(1,n)),arb(0))
            c=(2*tail+conv)/((2*m+1)*n+1)
            rho.append(c)
            tail-=c
            partial_mass+=n*c
            all_positive=all_positive and bool(c>0)
            all_tails=all_tails and bool(tail>0)
            if not all_positive or not all_tails:
                print(f"FAIL mass={numerator}/{denominator}, coefficient={n}")
                raise AssertionError("Insufficient precision or failed finite inequality")
        for label,good in [("all 256 coefficients positive",all_positive),
                           ("all 256 number remainders positive",all_tails),
                           ("partial mass below m",bool(partial_mass<m))]:
            print(("PASS " if good else "FAIL ")+f"m={numerator}/{denominator}: {label}")
            if not good: raise AssertionError(label)
        records.append({"mass":[numerator,denominator],
                        "order":ORDER,"all_coefficients_positive":all_positive,
                        "all_number_remainders_positive":all_tails,
                        "last_coefficient":str(rho[-1]),
                        "number_remainder":str(tail),
                        "partial_mass":str(partial_mass),
                        "mass_remainder":str(m-partial_mass)})
    output={"arithmetic":"Arb real balls","precision_bits":PRECISION,
            "scope":"Finite enclosures only; no proof of infinite positivity or total mass",
            "records":records}
    (ROOT/"certificates"/"coefficients-arb.json").write_text(
        json.dumps(output,indent=2)+"\n",encoding="utf-8")
    print("PASS finite Arb checks; certificate: certificates/coefficients-arb.json")


if __name__=="__main__":
    main()
