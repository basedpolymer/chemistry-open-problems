"""Combine the enumeration statistics (certificates/enumeration_*.json) into the table of the manuscript.

Environment: Python 3.12, standard library only.  Runtime: < 1 s.
Usage: python summarize_enumeration.py ../certificates/enumeration_*.json
Prints, for the formulations of Table 1 (thesis variant of preclusters, which coincides with the other two by
Lemma 4), the number of networks satisfying H for which the statement holds, summed over all enumerated families,
and checks (PASS/FAIL) that every formulation fails at least once in the class 'all'.
"""
import sys
import json

ROWS = [
    ("A|thesis|ray|P=M", "P(N*) = MIPS(N)"),
    ("A|thesis|ray|P<=M", "P(N*) refines MIPS(N)"),
    ("A|thesis|ray|M<=P", "MIPS(N) refines P(N*)"),
    ("A|thesis|ray|#P=#M", "#P(N*) = #MIPS(N)"),
    ("B|thesis|ray|P=M", "P(N) = MIPS(N*)"),
    ("B|thesis|ray|P<=M", "P(N) refines MIPS(N*)"),
    ("B|thesis|ray|M<=P", "MIPS(N*) refines P(N)"),
    ("B|thesis|ray|#P=#M", "#P(N) = #MIPS(N*)"),
    ("same-network|#P(N)=#M(N)", "#P(N) = #MIPS(N)"),
    ("A|thesis|semiflow|P=M", "P(N*) = semiflow supports of N"),
    ("A|thesis|set|P=M", "P(N*) = MIPS(N), set reading (q<=2)"),
    ("A|thesis|P=S", "P(N*) = minimal siphons of N"),
    ("A|thesis|S<=P", "each minimal siphon of N inside a block of P(N*)"),
    ("A|thesis|P<=S", "each block of P(N*) inside a minimal siphon of N"),
    ("A|thesis|P_blocks_siphons", "each block of P(N*) is a siphon of N"),
    ("A|thesis|P_blocks_unions_of_S", "each block of P(N*) is a union of minimal siphons"),
    ("A|thesis|#P=#S", "#P(N*) = #minimal siphons of N"),
    ("B|thesis|P=S", "P(N) = minimal siphons of N*"),
    ("B|thesis|S<=P", "each minimal siphon of N* inside a block of P(N)"),
    ("B|thesis|P<=S", "each block of P(N) inside a minimal siphon of N*"),
    ("B|thesis|P_blocks_siphons", "each block of P(N) is a siphon of N*"),
    ("B|thesis|#P=#S", "#P(N) = #minimal siphons of N*"),
    ("B|thesis|C<=P", "each internal cycle of N inside a precluster (Rem. 23)"),
    ("A|S<=MIPS", "each minimal siphon of N inside a MIPS of N"),
]
CLASSES = ["all", "strong", "wr", "wr_def0"]


def main(files):
    tot = {}
    counts = {c: 0 for c in CLASSES}
    for f in files:
        d = json.load(open(f))
        for c in CLASSES:
            counts[c] += d["counts"][c]
        for key, (h, fl) in d["stats"].items():
            st, c = key.split("||")
            a, b = tot.get((st, c), (0, 0))
            tot[(st, c)] = (a + h, b + fl)
    print("# families:", ", ".join(f.rsplit("/", 1)[-1] for f in files))
    print("# networks satisfying H per class:", counts)
    print(f"{'formulation':58s}" + "".join(f"{c:>15s}" for c in CLASSES))
    ok = True
    for st, label in ROWS:
        line = f"{label:58s}"
        for c in CLASSES:
            h, fl = tot.get((st, c), (0, 0))
            line += f"{f'{h}/{h + fl}':>15s}"
        print(line)
        if st.startswith("A|S<=MIPS"):
            continue
        h, fl = tot.get((st, "all"), (0, 0))
        if fl == 0:
            ok = False
    print(f"{'PASS' if ok else 'FAIL'}: every formulation of Conjecture 22 / Question 24 / Remark 23 in the table "
          f"fails for at least one enumerated network satisfying H")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
