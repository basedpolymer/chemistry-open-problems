"""Named networks used by the V1 checks (definitions copied from the manuscript paper.tex and from the notes,
re-typed by hand)."""
from v1lib import net_from_str, Net


def B(k):
    parts = ["X -> Y%d" % a for a in range(1, k + 1)] + ["Y%d -> Z" % a for a in range(1, k + 1)] + ["Z -> X"]
    return net_from_str(", ".join(parts), name="B%d" % k)


EX = {
    "R": net_from_str("X1 -> X2, X2 -> X1", name="R"),
    "T": net_from_str("X1 -> X2, X2 -> X3, X3 -> X1", name="T"),
    "W": net_from_str("X1 -> X2, X1 -> X3, X2 -> X1, X3 -> X1", name="W"),
    "U": net_from_str("X4 -> X2+X3, X3+X4 -> X1, X1+X2 -> 2X4", species_order=["X1", "X2", "X3", "X4"], name="U"),
    "D": net_from_str("X1 -> X2, X1 -> X3, X2+X3 -> 2X1", name="D"),
    "F": net_from_str("S0+E -> ES0, ES0 -> S0+E, ES0 -> S1+E, S1+F -> FS1, FS1 -> S1+F, FS1 -> S0+F", name="F"),
    "N3": net_from_str("2X1 -> X1+X2, X1+X2 -> 2X1, X1+X2 -> 2X2", name="N3"),
    "N412": net_from_str("2X1 -> X1+X2, X1+X2 -> 2X1", name="N412"),
    "B2": B(2), "B3": B(3), "B4": B(4),
    "K3": net_from_str("X -> Y, Y -> X, Y -> Z, Z -> Y, X -> Z, Z -> X", name="K3"),
    "XY4": net_from_str("X+Y -> 2X, X+Y -> 2Y, X+2Y -> 2X+Y, X+2Y -> 3Y", name="XY4"),
    "Z4": net_from_str("2X2 -> 2X1, 2X1 -> 2X2, X1+X2 -> 2X1, X1+X2 -> 2X2", species_order=["X1", "X2"], name="Z4"),
    "PiCF2": net_from_str("2X2 -> X1+X2, 2X2 -> 2X1, X1+X2 -> 2X2, X1+X2 -> 2X1", species_order=["X1", "X2"],
                          name="PiCF2"),
    "Y5": net_from_str("A -> B, B -> C, B -> D, B -> E, C -> A, D -> A, E -> A", name="Y5"),
}
EX["Wstar"] = EX["W"].dual()
EX["Wstar"].name = "W*"
# B2dup: B2 with Z -> X duplicated
EX["B2dup"] = net_from_str("X -> Y1, X -> Y2, Y1 -> Z, Y2 -> Z, Z -> X, Z -> X", name="B2dup")

# K8 pairs (agent P, section 5)
K8 = {
    "N1": net_from_str("X3 -> X4, X2+X3 -> X1+X2, X1+X4 -> 2X3", species_order=["X1", "X2", "X3", "X4"], name="N1"),
    "N2": net_from_str("2X3 -> 2X4, X2+X3 -> 2X2, X2+X3 -> 2X1, X1+X4 -> 2X3", species_order=["X1", "X2", "X3", "X4"],
                       name="N2"),
    "M1": net_from_str("X3 -> X1+X2, X2+X3 -> X1+X3, X1+X2 -> X3, X1+X2 -> 2X2", species_order=["X1", "X2", "X3"],
                       name="M1"),
    "M2": net_from_str("2X3 -> X2+X3, X1+X3 -> X1+X2, X1+X2 -> 2X3, X1+X2 -> 2X1", species_order=["X1", "X2", "X3"],
                       name="M2"),
    # agent E's species-side pair (computation/notes.md section 6)
    "E1": net_from_str("X2+X3 -> 2X3, X2+X3 -> 2X2, X1+X2+X3 -> X1+2X3", species_order=["X1", "X2", "X3"], name="E1"),
    "E2": net_from_str("X2+X3 -> 2X3, X1+X2+X3 -> 2X2+X3, X1+X2+X3 -> 2X1+X2", species_order=["X1", "X2", "X3"],
                       name="E2"),
}
