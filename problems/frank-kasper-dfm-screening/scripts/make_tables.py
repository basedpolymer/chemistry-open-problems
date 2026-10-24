#!/usr/bin/env python3
"""Build the screen tables (CSV in data/, Markdown for README) from the run
results in work/runs (box-relaxed runs), with catalogue metadata.
usage: make_tables.py   (writes data/dfm_screen_boxrelaxed.csv and work/screen_table.md)"""
import os
import sys
import json
import glob
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analyze import load_runs, extrapolate  # noqa: E402

ROOT = os.path.join(HERE, "..")
PHYS = {"Phys01": "C14 (MgZn2)", "Phys02": "C15 (MgCu2)", "Phys03": "C36 (MgNi2)", "Phys04": "6-layer Laves (MgCuNi)",
        "Phys06": "9-layer Laves", "Phys12": "zra-d (K7Cs6)", "Phys13": "p-sigma (Th6Cd7)", "Phys14": "mu (Mo6Co7)",
        "Phys18": "Z (Zr4Al3)", "Phys20": "delta (MoNi)", "Phys25": "H (complex)", "Phys26": "sigma (Cr46Fe54)",
        "Phys27": "A15 (Cr3Si)"}
PAPERNAME = {"BCC": "BCC", "A15": "A15", "C15": "C15", "Z": "Z", "C14": "C14", "sigma": "sigma", "psigma": "p-sigma",
             "H": "H", "mu": "mu", "M": "M (Nb48Ni39Al13)", "P": "P (Mo42Cr18Ni40)", "delta": "delta"}


def group_of(name):
    if ":" not in name:
        return ""
    pref, cid = name.split(":")
    sub = "ListDATApaper" if pref in ("cat", "catR") else "ListDATApaperStr"
    f = os.path.join(ROOT, "certificates", "sources", "dutour2010", "SpaceFullerene", sub, cid + ".cgd")
    for line in open(f):
        w = line.split()
        if w and w[0].lower() == "group":
            return w[1]
    return ""


def main():
    runs = load_runs(os.path.join(ROOT, "work", "runs"))
    rows = []
    for (name, box), r in runs.items():
        if box != "relax":
            continue
        F = {L: rec["F"] for L, rec in r["levels"].items()}
        R, err, rho, K = extrapolate(F)
        last = r["levels"][K]
        cells = last.get("cells")
        comp = ""
        if cells:
            c = Counter(cells["nfaces"])
            comp = "%d/%d/%d/%d" % (c.get(12, 0), c.get(14, 0), c.get(15, 0), c.get(16, 0))
        if ":" in name:
            cid = name.split(":")[1]
            label = ("cat2 " if name.startswith("cat2") else "") + cid + (" re-emb." if name.split(":")[0].endswith("R") else "")
            phys = PHYS.get(cid, "") if name.startswith("cat:") else ""
        else:
            label = "paper " + name
            phys = PAPERNAME.get(name, name)
        rows.append({"name": name, "label": label, "phys": phys, "group": group_of(name), "n": r["meta"].get("ncell"),
                     "comp": comp, "K": K, "F_K": F[K], "R": R, "err": err, "rho": rho,
                     "min_edge": last["min_foam_edge_over_R0"], "Vrms": last["Vrms_over_V0"],
                     "cone": last["cone_min"],
                     "levels": [F[L] for L in sorted(F)], "L_unitvol": last["L_unitvol"]})
    # merge aliases: paper structures present in the catalogue, re-embedded restarts
    P2C = {"sigma": "Phys26", "A15": "Phys27", "C15": "Phys02", "Z": "Phys18", "C14": "Phys01", "psigma": "Phys13",
           "H": "Phys25", "mu": "Phys14", "delta": "Phys20"}
    def ident(nm):
        if ":" in nm:
            pref, cid = nm.split(":")
            return ("II:" if pref.startswith("cat2") else "") + cid
        return P2C.get(nm, nm)
    groups = {}
    for r in rows:
        r["key"] = r["R"] if r["R"] is not None else r["F_K"]
        groups.setdefault(ident(r["name"]), []).append(r)
    merged = []
    for idn, lst in groups.items():
        # prefer embedded (cone_min > 0) results, then the deepest level, then the lowest value
        lst.sort(key=lambda r: (r["cone"] <= 0, -r["K"], r["key"]))
        b = dict(lst[0])
        b["id"] = idn
        b["aliases"] = " ".join(sorted(x["name"] for x in lst))
        if not b["group"]:
            for x in lst:
                if x["group"]:
                    b["group"] = x["group"]
        if not b["phys"] or b["name"] in P2C:
            b["phys"] = PHYS.get(idn, PAPERNAME.get(b["name"], ""))
        b["label"] = idn
        merged.append(b)
    rows = merged
    for r in rows:
        if r["cone"] <= 0:      # relaxation left the class of embedded partitions: no valid value
            r["key"] = float("inf")
    rows.sort(key=lambda r: r["key"])
    sig = [r for r in rows if r["id"] == "Phys26"][0]["key"]
    hdr = "rank,id,physical,group,ncell,Z12/Z14/Z15/Z16,maxlevel,F_maxlevel_upper_bound,F_extrap,err_est,rel_diff_vs_sigma,min_edge_over_R0,Vrms_over_V0,run_used,aliases,F_by_level"
    lines = [hdr]
    md = ["| # | structure | group | n | Z12/14/15/16 | L | 𝓕_L (upper bound) | 𝓕 extrapolated | Δ vs σ | min edge/R₀ |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(rows, 1):
        if r["cone"] <= 0:
            lines.append("%d,%s,%s,%s,%s,%s,%d,,,,,,,%s,%s,INVALID (relaxation not embedded: cone_min=%.1e)" % (
                i, r["id"], r["phys"], r["group"], r["n"], r["comp"], r["K"], r["name"], r["aliases"], r["cone"]))
            md.append("| %d | %s | %s | %s | %s | %d | — | relaxation degenerates (not embedded) | — | %.3f |" % (
                i, r["label"], r["group"], r["n"], r["comp"], r["K"], r["min_edge"]))
            continue
        d = (r["key"] - sig) / sig
        lines.append("%d,%s,%s,%s,%s,%s,%d,%.10f,%s,%s,%+.3e,%.4f,%.4f,%s,%s,%s" % (
            i, r["id"], r["phys"], r["group"], r["n"], r["comp"], r["K"], r["F_K"],
            "%.9f" % r["R"] if r["R"] is not None else "", "%.1e" % r["err"] if r["err"] is not None else "",
            d, r["min_edge"], r["Vrms"], r["name"], r["aliases"], ";".join("%.10f" % x for x in r["levels"])))
        md.append("| %d | %s%s | %s | %s | %s | %d | %.8f | %s | %+.3f%% | %.3f |" % (
            i, r["label"], (" — " + r["phys"]) if r["phys"] else "", r["group"], r["n"], r["comp"], r["K"], r["F_K"],
            ("%.7f ± %.0e" % (r["R"], r["err"])) if (r["R"] is not None and r["err"] is not None) else
            ("%.7f" % r["R"] if r["R"] is not None else "—"), 100 * d, r["min_edge"]))
    open(os.path.join(ROOT, "data", "dfm_screen_boxrelaxed.csv"), "w").write("\n".join(lines) + "\n")
    open(os.path.join(ROOT, "work", "screen_table.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md[:40]))
    print("... %d rows" % len(rows))


if __name__ == "__main__":
    main()
