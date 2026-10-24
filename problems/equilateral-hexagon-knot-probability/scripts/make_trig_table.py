#!/usr/bin/env python3
"""Rigorous outward-rounded enclosures of cos(k*pi/N), sin(k*pi/N), k = 0..N, N = 2**LEV.

Environment: Python 3.13, mpmath 1.3.  Runtime: ~10 s for LEV = 14.
Each value is computed with mpmath at 60 decimal digits (error < 1e-55) and then rounded
outward to IEEE doubles (lo = largest double <= value - 1e-40, hi = smallest double >= value + 1e-40),
so [lo, hi] contains the exact value.  Output: binary file of 4*(N+1) little-endian doubles
(coslo, coshi, sinlo, sinhi) per k.  Prints PASS if the self-check succeeds.
"""
import sys, struct, math
import mpmath as mp

def main():
    lev = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    out = sys.argv[2] if len(sys.argv) > 2 else "trig_table.bin"
    N = 2 ** lev
    mp.mp.dps = 60
    eps = mp.mpf(10) ** -40
    pi = mp.pi
    buf = bytearray()
    ok = True
    for k in range(N + 1):
        x = pi * k / N
        c = mp.cos(x); s = mp.sin(x)
        vals = []
        for v in (c, s):
            lo = float(v - eps); hi = float(v + eps)   # round to nearest, then push outward
            while mp.mpf(lo) > v - eps: lo = math.nextafter(lo, -math.inf)
            while mp.mpf(hi) < v + eps: hi = math.nextafter(hi, math.inf)
            if not (mp.mpf(lo) <= v <= mp.mpf(hi)): ok = False
            vals += [lo, hi]
        buf += struct.pack("<4d", *vals)
    with open(out, "wb") as f:
        f.write(buf)
    print(("PASS" if ok else "FAIL"), "trig table", out, "N =", N)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
