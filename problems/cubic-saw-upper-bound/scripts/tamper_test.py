#!/usr/bin/env python3
"""Tamper test: corrupt a valid certificate in several ways and check that verify_cert rejects every copy.
Usage: tamper_test.py <verify_cert binary> <valid.cert or .cert2> <tmpdir>.  Python 3 standard library.  Runtime ~10 s.
Prints PASS/FAIL per case; exit code 0 iff the original passes and every corrupted copy fails."""
import os
import struct
import subprocess
import sys

binary, cert, tmp = sys.argv[1], sys.argv[2], sys.argv[3]
data = bytearray(open(cert, "rb").read())
d, K = struct.unpack_from("<II", data, 8)
N, P, Q = struct.unpack_from("<QQQ", data, 16)
codes_off, w_off = 40, 40 + 8 * N


def run(buf, name):
    path = os.path.join(tmp, "tamper_" + name + ".cert")
    open(path, "wb").write(buf)
    r = subprocess.run([binary, path, "2"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return r.returncode == 0


ok = True
if run(bytes(data), "original"):
    print("PASS original certificate accepted")
else:
    print("FAIL original certificate rejected")
    ok = False
cases = []
b = bytearray(data); struct.pack_into("<Q", b, 24, P - 1); cases.append(("claim P-1", b))
b = bytearray(data)
i = N // 2
w = struct.unpack_from("<Q", b, w_off + 8 * i)[0]
struct.pack_into("<Q", b, w_off + 8 * i, max(1, w * 9 // 10)); cases.append(("weight of one state -10%", b))
b = bytearray(data)
# remove the one-step walk: the start state is then left with the empty walk as only child (2d V0 > (P/Q) V0)
ws = struct.unpack_from("<%dQ" % N, data, w_off)
cs = struct.unpack_from("<%dQ" % N, data, codes_off)
j = next(k for k in range(N) if (cs[k] >> 57) & 63 == 1)
nc = list(cs[:j]) + list(cs[j + 1:]); nw = list(ws[:j]) + list(ws[j + 1:])
b = bytearray(data[:40]); struct.pack_into("<Q", b, 16, N - 1)
b += struct.pack("<%dQ" % (N - 1), *nc) + struct.pack("<%dQ" % (N - 1), *nw)
cases.append(("remove the one-step state", b))
b = bytearray(data); struct.pack_into("<Q", b, w_off + 8 * 5, 0); cases.append(("zero weight", b))
if data[:8] == b"SAWCERT2":
    # explicit children: redirect the child of (state 1, direction 0) to the longest listed walk (not a sub-walk)
    kid_off = w_off + 8 * N
    longest = max(range(N), key=lambda k: (cs[k] >> 57) & 63)
    b = bytearray(data); struct.pack_into("<I", b, kid_off + 4 * (1 * 2 * d + 0), longest)
    cases.append(("SAWCERT2 child pointer redirected", b))
    b = bytearray(data); struct.pack_into("<I", b, kid_off + 4 * (1 * 2 * d + 0), 0xffffffff)
    cases.append(("SAWCERT2 child missing for an allowed step", b))
for name, buf in cases:
    if run(bytes(buf), name.replace(" ", "_").replace(">=", "ge").replace("%", "pc")):
        print("FAIL corrupted certificate accepted:", name)
        ok = False
    else:
        print("PASS corrupted certificate rejected:", name)
sys.exit(0 if ok else 1)
