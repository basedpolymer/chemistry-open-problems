#!/usr/bin/env python3
"""Adaptive branch-and-bound driver for the certified enclosure of the knotting probability of
equilateral hexagons (see ../result.md, Section 5).  Calls the C++ certifier src/hexcert on d-boxes and
keeps the exact dyadic sums of the upper and lower bounds.

Usage:
  python3 bnb_driver.py --state STATE.pkl --bin HEXCERT --table TRIG.bin [--rounds R] [--batch B]
         [--init-level L0] [--lth-init T0] [--lth-max TMAX] [--hours H]
The state file is rewritten (atomically) after every round, so the run can be interrupted and resumed with
the same command (resumable chunks).  One thread; memory grows with the number of leaves (~300 B per leaf).

A leaf is a dyadic box D = [2 a1, 2 a1 + 2]/2^L1 x [a2, a2+1]/2^L2 x [a3, a3+1]/2^L3 of the (d1, s, t) chart
with a theta level Lth and certified counts (outR, certR, outL, certL).  Its contribution to the knotting
probability is the exact dyadic rational  w(D) * count,
    w(D) = ((a1+1)^3 - a1^3) (2 a2 + 1) / (4 * 2^(3 L1 + 2 L2 + L3 + 3 Lth)),
upper = sum w (outR + outL), lower = sum w (certR + certL).

Refinement (heuristic, does not affect validity): the leaves with the largest gap w*(out - cert) are refined;
a leaf is either re-evaluated with Lth + 1 or bisected along one d-axis, according to the certifier's
attribution of its remaining uncertainty (largest share wins).  The leaves always tile the chart.
"""
import argparse, pickle, subprocess, sys, time, os
from fractions import Fraction

EXP = 600   # common exponent: values are stored as integers / 2^EXP
OCT_LEVEL = 6   # boxes coarser than 2^-6 (in s, t; 2^-5 in d1/2) are split isotropically

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
    def run(self, boxes, chunk=100):
        out = []
        for i in range(0, len(boxes), chunk):
            part = boxes[i:i + chunk]
            for j, b in enumerate(part):
                self.p.stdin.write("%d %d %d %d %d %d %d %d\n" % ((i + j,) + tuple(b)))
            self.p.stdin.flush()
            for j in range(len(part)):
                line = self.p.stdout.readline()
                f = line.split()
                if len(f) != 18 or int(f[0]) != i + j:
                    raise RuntimeError("bad certifier output: %r" % line)
                out.append([int(x) for x in f[1:14]] + [float(x) for x in f[14:18]])
        return out
    def close(self):
        self.p.stdin.close(); self.p.wait()

def leaf_record(box, res):
    (Lth_out, outR, inR, certR, outL, inL, certL, ncR, ncL, failR, failL, unkR, unkL, a0, a1, a2, a3) = res
    assert Lth_out == box[6]
    w = weight_num(*box)
    return {"box": tuple(box), "outR": outR, "inR": inR, "certR": certR, "outL": outL, "inL": inL, "certL": certL,
            "failR": failR, "failL": failL, "unkR": unkR, "unkL": unkL, "attr": (a0, a1, a2, a3),
            "up": w * (outR + outL), "lo": w * (certR + certL)}

def compact(l):
    return tuple(l["box"]) + (l["outR"], l["certR"], l["outL"], l["certL"])

def refine(l, lth_max):
    """children boxes of leaf l (list), according to the uncertainty attribution"""
    L1, a1, L2, a2, L3, a3, Lth = l["box"]
    m = min(L1 - 1, L2, L3)
    if m < OCT_LEVEL:                     # large boxes: isotropic octree split, theta level tied to box size
        nl = min(lth_max, max(Lth, m + 1 + 2))
        return [(L1 + 1, 2 * a1 + i, L2 + 1, 2 * a2 + j, L3 + 1, 2 * a3 + k, nl)
                for i in range(2) for j in range(2) for k in range(2)]
    at = list(l["attr"])
    if Lth >= lth_max or Lth > min(L1 - 1, L2, L3) + 3:
        at[0] = -1.0                      # theta grid already much finer than the d-box
    if max(at) <= 0:                      # no usable attribution: bisect the coarsest d-axis
        lev = {1: L1 - 1, 2: L2, 3: L3}
        k = min((2, 3, 1), key=lambda i: lev[i])
    else:
        k = max(range(4), key=lambda i: at[i])
    if k == 0:
        return [(L1, a1, L2, a2, L3, a3, Lth + 1)]
    if k == 1:
        return [(L1 + 1, 2 * a1 + i, L2, a2, L3, a3, Lth) for i in range(2)]
    if k == 2:
        return [(L1, a1, L2 + 1, 2 * a2 + i, L3, a3, Lth) for i in range(2)]
    return [(L1, a1, L2, a2, L3 + 1, 2 * a3 + i, Lth) for i in range(2)]

def summary(state):
    up = state["up_closed"] + sum(l["up"] for l in state["leaves"])
    lo = state["lo_closed"] + sum(l["lo"] for l in state["leaves"])
    return Fraction(up, 1 << EXP), Fraction(lo, 1 << EXP)

def save(state, path):
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)

def add_results(state, boxes, res):
    for b, r in zip(boxes, res):
        rec = leaf_record(b, r)
        if rec["unkR"] or rec["unkL"]:
            state["anomalies"].append(rec)
        if rec["up"] == rec["lo"]:
            state["up_closed"] += rec["up"]; state["lo_closed"] += rec["lo"]
            state["closed"].append(compact(rec))
        else:
            state["leaves"].append(rec)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--rounds", type=int, default=10 ** 9)
    ap.add_argument("--batch", type=int, default=2000)
    ap.add_argument("--init-level", type=int, default=3)
    ap.add_argument("--lth-init", type=int, default=6)
    ap.add_argument("--lth-max", type=int, default=12)
    ap.add_argument("--hours", type=float, default=1e9)
    ap.add_argument("--oct-level", type=int, default=6)
    args = ap.parse_args()
    global OCT_LEVEL
    OCT_LEVEL = args.oct_level
    t0 = time.time()
    cert = Certifier(args.bin, args.table)
    if os.path.exists(args.state):
        with open(args.state, "rb") as f:
            state = pickle.load(f)
        print("resumed: %d open leaves, %d closed, %d boxes evaluated" % (len(state["leaves"]), len(state["closed"]), state["nevals"]), flush=True)
    else:
        L = args.init_level
        boxes = [(L + 1, a1, L, a2, L, a3, args.lth_init)
                 for a1 in range(2 ** (L + 1)) for a2 in range(2 ** L) for a3 in range(2 ** L)]
        state = {"leaves": [], "closed": [], "up_closed": 0, "lo_closed": 0, "nevals": len(boxes),
                 "rounds": 0, "anomalies": [], "log": [], "args": vars(args), "cpu": 0.0}
        add_results(state, boxes, cert.run(boxes))
        save(state, args.state)
    for rnd in range(args.rounds):
        if time.time() - t0 > args.hours * 3600:
            break
        tr = time.time()
        leaves = state["leaves"]
        leaves.sort(key=lambda l: -(l["up"] - l["lo"]))
        sel, keep = leaves[:args.batch], leaves[args.batch:]
        state["leaves"] = keep
        newboxes = []
        nth = 0
        for l in sel:
            ch = refine(l, args.lth_max)
            nth += (len(ch) == 1)
            newboxes += ch
        add_results(state, newboxes, cert.run(newboxes))
        state["nevals"] += len(newboxes)
        state["rounds"] += 1
        up, lo = summary(state)
        msg = ("round %d evals %d open %d closed %d  upper %.6e lower %.6e ratio %.4f  theta-refinements %d/%d  "
               "anomalies %d  dt=%.0fs t=%.0fs") % (
            state["rounds"], state["nevals"], len(state["leaves"]), len(state["closed"]), float(up), float(lo),
            float(up / lo) if lo else float("inf"), nth, len(sel), len(state["anomalies"]), time.time() - tr, time.time() - t0)
        state["log"].append(msg)
        print(msg, flush=True)
        save(state, args.state)
    cert.close()
    up, lo = summary(state)
    print("FINAL upper %.8e lower %.8e" % (float(up), float(lo)))

if __name__ == "__main__":
    main()
