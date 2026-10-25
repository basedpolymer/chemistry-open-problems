"""Preclusters, MIPS and siphons of standard biochemical motifs (exact).

Environment: Python 3.12 (WSL Fedora 44), standard library only; exact rational arithmetic (crn.py).
Runtime: < 10 s.

For each motif the script prints the hypotheses of Conjecture 22 and the objects on both sides of the duality,
then checks (PASS/FAIL) the statements recorded in Table 2 of the manuscript.  Exits with status 1 on failure.
"""
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from crn import Network, refines, fmt  # noqa: E402


def parse(species, text):
    """Parse 'A + E -> AE, AE -> A + E' into a Network over the given species names."""
    idx = {s: i for i, s in enumerate(species)}

    def cplx(t):
        v = [0] * len(species)
        t = t.strip()
        if t == "0":
            return tuple(v)
        for term in t.split("+"):
            term = term.strip()
            k = 1
            j = 0
            while j < len(term) and term[j].isdigit():
                j += 1
            if j:
                k = int(term[:j])
            v[idx[term[j:].strip()]] += k
        return tuple(v)
    rx = []
    for part in text.split(","):
        a, b = part.split("->")
        rx.append((cplx(a), cplx(b)))
    return Network(len(species), rx, names=species)


MOTIFS = {
    "reversible Michaelis-Menten": (["E", "S", "ES", "P"],
                                    "E + S -> ES, ES -> E + S, ES -> E + P, E + P -> ES"),
    "one-site futile cycle": (["S0", "S1", "E", "F", "ES0", "FS1"],
                              "S0 + E -> ES0, ES0 -> S0 + E, ES0 -> S1 + E, "
                              "S1 + F -> FS1, FS1 -> S1 + F, FS1 -> S0 + F"),
    "two-site distributive phosphorylation": (
        ["S0", "S1", "S2", "E", "F", "ES0", "ES1", "FS2", "FS1"],
        "S0 + E -> ES0, ES0 -> S0 + E, ES0 -> S1 + E, S1 + E -> ES1, ES1 -> S1 + E, ES1 -> S2 + E, "
        "S2 + F -> FS2, FS2 -> S2 + F, FS2 -> S1 + F, S1 + F -> FS1, FS1 -> S1 + F, FS1 -> S0 + F"),
    "thesis network G (Example 1.17)": (["X1", "X2"],
                                        "3X1 -> 2X1 + X2, 2X1 + X2 -> 3X1, X1 + X2 -> 2X2, X2 -> X1"),
    "triangle with dimerisation": (["A", "B", "C"], "2A -> 2B, 2B -> A + C, A + C -> 2A"),
}

# expected values used in Table 2 of the manuscript:
# (H, P(N)=M(N*), P(N*)=M(N), P(N*) refines M(N), M(N) refines P(N*), P(N) refines M(N*), M(N*) refines P(N))
EXPECTED = {
    "reversible Michaelis-Menten": (True, False, False, True, False, False, True),
    "one-site futile cycle": (True, False, False, False, True, False, True),
    "two-site distributive phosphorylation": (True, False, False, True, False, False, True),
    "thesis network G (Example 1.17)": (True, True, False, True, False, True, True),
    "triangle with dimerisation": (True, False, False, True, False, True, False),
}

def main():
    fails = 0
    for name, (sp, text) in MOTIFS.items():
        N = parse(sp, text)
        D = N.dual()
        Hval = N.is_conservative() and N.is_consistent() and N.rows_distinct()
        PN, MD = N.preclusters("thesis"), D.mips_ray()
        PD, MN = D.preclusters("thesis"), N.mips_ray()
        rnames = [f"r{j+1}" for j in range(N.r)]
        print(f"== {name}: {N.describe()}")
        print(f"   n={N.n}, r={N.r}, rank={N.rank()}, deficiency={N.deficiency()}, weakly reversible={N.weakly_reversible()}")
        print(f"   conservative={N.is_conservative()}, consistent={N.is_consistent()}, rows distinct={N.rows_distinct()}")
        print(f"   reaction side: P(N) = {fmt(PN, rnames)};  MIPS(N*) (coupled reactions) = {fmt(MD, rnames)}")
        print(f"   species side:  P(N*) = {fmt(PD, N.names)};  MIPS(N) = {fmt(MN, N.names)}")
        print(f"   minimal siphons of N = {fmt(N.minimal_siphons(), N.names)};  semiflow supports = {fmt(N.mips_semiflow(), N.names)}")
        for v in ("onestep", "closure"):
            if set(N.preclusters(v)) != set(PN) or set(D.preclusters(v)) != set(PD):
                print(f"   note: precluster variant {v} differs")
        got = (Hval, set(PN) == set(MD), set(PD) == set(MN), refines(PD, MN), refines(MN, PD),
               refines(PN, MD), refines(MD, PN))
        ok = got == EXPECTED[name]
        fails += not ok
        print(f"{'PASS' if ok else 'FAIL'}: {name}: (H, P(N)=M(N*), P(N*)=M(N), P(N*)<=M(N), M(N)<=P(N*), P(N)<=M(N*), M(N*)<=P(N)) = {got}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
