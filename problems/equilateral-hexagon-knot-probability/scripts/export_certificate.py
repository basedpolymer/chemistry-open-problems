#!/usr/bin/env python3
"""Export the leaves of a branch-and-bound state (bnb_driver.py) as a certificate file.

Output: gzip text, one leaf per line
    L1 a1 L2 a2 L3 a3 Lth outR certR outL certL
sorted lexicographically, plus a summary JSON with the exact bounds (as fractions num/2^EXP) and SHA-256.
Usage: python3 export_certificate.py STATE.pkl OUT_PREFIX
"""
import sys, pickle, gzip, hashlib, json
from fractions import Fraction
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from bnb_driver import weight_num, EXP

def main():
    st = pickle.load(open(sys.argv[1], "rb"))
    pre = sys.argv[2]
    rows = []
    for l in st["leaves"]:
        rows.append(tuple(l["box"]) + (l["outR"], l["certR"], l["outL"], l["certL"]))
    rows += [tuple(r) for r in st["closed"]]
    # zero leaves are not stored individually in the driver state: they are the boxes dropped... (none: the
    # driver keeps every evaluated leaf either in 'leaves' or in 'closed')
    rows.sort()
    up = sum(weight_num(*r[:7]) * (r[7] + r[9]) for r in rows)
    lo = sum(weight_num(*r[:7]) * (r[8] + r[10]) for r in rows)
    text = "".join(" ".join(map(str, r)) + "\n" for r in rows).encode()
    sha = hashlib.sha256(text).hexdigest()
    with gzip.open(pre + ".leaves.gz", "wb", compresslevel=9) as f:
        f.write(text)
    U = Fraction(up, 1 << EXP); L = Fraction(lo, 1 << EXP)
    summ = {"n_leaves": len(rows), "EXP": EXP, "upper_num": str(up), "lower_num": str(lo),
            "upper_float": float(U), "lower_float": float(L), "sha256_leaves_text": sha,
            "rounds": st["rounds"], "nevals": st["nevals"], "anomalies": len(st["anomalies"]),
            "driver_args": st.get("args")}
    json.dump(summ, open(pre + ".summary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if not k.endswith("_num")}, indent=1))

if __name__ == "__main__":
    main()
