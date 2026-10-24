#!/usr/bin/env python3
"""Adaptive branch-and-bound driver for the certified enclosure of the knotting probability of
equilateral hexagons (see ../result.md).  Calls the C++ certifier src/hexcert on d-boxes and keeps
the exact dyadic sums of the upper and lower bounds.

Usage:
  python3 bnb_driver.py --state STATE.pkl --bin HEXCERT --table TRIG.bin [--rounds R] [--batch B]
         [--offset OFF] [--lth-max LMAX] [--hours H]
The state file is updated after every round, so the run can be interrupted and resumed with the
same command (resumable chunks).  One thread, memory < 1 GB.

Each leaf box D = [2 a1, 2 a1 + 2]/2^L1 x [a2, a2+1]/2^L2 x [a3, a3+1]/2^L3 in (d1, s, t) carries
theta level Lth and certified counts (outR, inR, certR, outL, inL, certL).  Its contribution to the
probability is (exact dyadic rational)
    w(D) * count,   w(D) = ((a1+1)^3 - a1^3) (2 a2 + 1) / (4 * 2^(3 L1 + 2 L2 + L3 + 3 Lth)),
upper = sum w (outR + outL), lower = sum w (certR + certL).
"""
import argparse, pickle, subprocess, sys, time, heapq, os
from fractions import Fraction

EXP = 400   # common exponent: values are stored as integers / 2^EXP

def weight_num(L1, a1, L2, a2, L3, a3, Lth):
    """w(D) * 2^EXP as an exact integer (w is a dyadic rational)."""
    num = ((a1 + 1) ** 3 - a1 ** 3) * (2 * a2 + 1)
    e = 2 + 3 * L1 + 2 * L2 + L3 + 3 * Lth
    assert e <= EXP
    return num << (EXP - e)

class Certifier:
    def __init__(self, binpath, table):
        self.p = subprocess.Popen([binpath, table, "1"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  text=True, bufsize=1)
    def run(self, boxes, chunk=200):
        out = []
        for i in range(0, len(boxes), chunk):
            part = boxes[i:i + chunk]
            for j, b in enumerate(part):
                self.p.stdin.write("%d %d %d %d %d %d %d %d\n" % ((i + j,) + tuple(b)))
            self.p.stdin.flush()
            for j in range(len(part)):
                line = self.p.stdout.readline()
                f = line.split()
                if len(f) != 14 or int(f[0]) != i + j:
                    raise RuntimeError("bad certifier output: %r" % line)
                vals = list(map(int, f[1:]))
                out.append(vals)
        return out
    def close(self):
        self.p.stdin.close(); self.p.wait()

def lth_for(L2, off, lmin, lmax):
    return max(lmin, min(lmax, L2 + off))

def leaf_record(box, res):
    L1, a1, L2, a2, L3, a3, Lth = box
    (Lth_out, outR, inR, certR, outL, inL, certL, ncR, ncL, failR, failL, unkR, unkL) = res
    assert Lth_out == Lth
    w = weight_num(*box)
    return {"box": box, "outR": outR, "inR": inR, "certR": certR, "outL": outL, "inL": inL, "certL": certL,
            "failR": failR, "failL": failL, "unkR": unkR, "unkL": unkL,
            "up": w * (outR + outL), "lo": w * (certR + certL)}

def children(box, off, lmin, lmax):
    L1, a1, L2, a2, L3, a3, Lth = box
    ch = []
    nL2 = L2 + 1
    nl = lth_for(nL2, off, lmin, lmax)
    nl = max(nl, Lth)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ch.append((L1 + 1, 2 * a1 + i, L2 + 1, 2 * a2 + j, L3 + 1, 2 * a3 + k, nl))
    return ch

def summary(state):
    up = state["up_closed"] + sum(l["up"] for l in state["leaves"])
    lo = state["lo_closed"] + sum(l["lo"] for l in state["leaves"])
    return Fraction(up, 1 << EXP), Fraction(lo, 1 << EXP)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--rounds", type=int, default=1000000)
    ap.add_argument("--batch", type=int, default=500)
    ap.add_argument("--offset", type=int, default=2)
    ap.add_argument("--lth-min", type=int, default=6)
    ap.add_argument("--lth-max", type=int, default=9)
    ap.add_argument("--hours", type=float, default=1e9)
    ap.add_argument("--init-level", type=int, default=3)
    args = ap.parse_args()
    t0 = time.time()
    cert = Certifier(args.bin, args.table)
    if os.path.exists(args.state):
        with open(args.state, "rb") as f:
            state = pickle.load(f)
        print("resumed: %d leaves, %d boxes evaluated" % (len(state["leaves"]), state["nevals"]), flush=True)
    else:
        L = args.init_level
        boxes = []
        nl = lth_for(L, args.offset, args.lth_min, args.lth_max)
        for a1 in range(2 ** (L + 1)):
            for a2 in range(2 ** L):
                for a3 in range(2 ** L):
                    boxes.append((L + 1, a1, L, a2, L, a3, nl))
        res = cert.run(boxes)
        leaves = [leaf_record(b, r) for b, r in zip(boxes, res)]
        state = {"leaves": [l for l in leaves if l["up"] > 0], "up_closed": 0, "lo_closed": 0,
                 "nevals": len(boxes), "rounds": 0, "anomalies": [], "log": []}
        for l in leaves:
            if l["unkR"] or l["unkL"]:
                state["anomalies"].append(l)
    for rnd in range(args.rounds):
        if time.time() - t0 > args.hours * 3600:
            break
        leaves = state["leaves"]
        gaps = sorted(range(len(leaves)), key=lambda i: -(leaves[i]["up"] - leaves[i]["lo"]))
        sel = set(gaps[:args.batch])
        newboxes = []
        keep = []
        for i, l in enumerate(leaves):
            if i in sel:
                newboxes += children(l["box"], args.offset, args.lth_min, args.lth_max)
            else:
                keep.append(l)
        res = cert.run(newboxes)
        for b, r in zip(newboxes, res):
            rec = leaf_record(b, r)
            if rec["unkR"] or rec["unkL"]:
                state["anomalies"].append(rec)
            if rec["up"] > 0:
                if rec["up"] == rec["lo"]:
                    state["up_closed"] += rec["up"]; state["lo_closed"] += rec["lo"]
                else:
                    keep.append(rec)
        state["leaves"] = keep
        state["nevals"] += len(newboxes)
        state["rounds"] += 1
        up, lo = summary(state)
        msg = "round %d evals %d leaves %d  upper %.6e lower %.6e ratio %.4f  anomalies %d  t=%.0fs" % (
            state["rounds"], state["nevals"], len(keep), float(up), float(lo), float(up / lo) if lo else float("inf"),
            len(state["anomalies"]), time.time() - t0)
        state["log"].append(msg)
        print(msg, flush=True)
        tmp = args.state + ".tmp"
        with open(tmp, "wb") as f:
            pickle.dump(state, f)
        os.replace(tmp, args.state)
    cert.close()
    up, lo = summary(state)
    print("FINAL upper %.8e lower %.8e" % (float(up), float(lo)))

if __name__ == "__main__":
    main()
