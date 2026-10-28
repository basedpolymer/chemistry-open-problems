# Rigorous upper bound for the connective constant of the simple cubic lattice

<!-- plain-words:start -->
**The question.** A polymer chain cannot cross itself: on the cubic lattice, the number of such “self-avoiding walks” of n steps grows like μⁿ. Numerically μ ≈ 4.684, but the best proved upper bound (4.7387) dated from 2000.

**The answer.** μ ≤ 4.722048, a 30 % cut of the gap to the estimate, plus better bounds in dimensions 4, 5 and 6.

**How.** A finite automaton remembers recent sites and forbids only short loops, so it counts every self-avoiding walk plus a few extra paths. With a longer memory (206 million states) and a smarter choice of what to forget, its growth rate gives the bound, certified by an integer vector anyone can re-check.

**Overlaps.** Method reused by [result 24](../saw-diamond-fcc-upper-bounds/).

<p align="center"><img src="../../figures/saw-automaton-bound.svg" alt="A self-avoiding walk, the automaton bound and the new value" width="720"></p>

<p align="center"><sub>An automaton that only forbids short loops overcounts self-avoiding walks, so its growth rate bounds μ from above.</sub></p>

<!-- plain-words:end -->

**Status: record: new rigorous bounds with exact certificates, verified (`verification.md`, verdict
CONFIRMED WITH CORRECTIONS; the corrections are editorial and change no bound).**

**Attribution.** Method, certificates, checker and write-up: Claude Opus 5.5 (effort xhigh), sub-agent of the
orchestrator, 18 Oct 2026. LaTeX edition: Claude Opus 5.5 (effort xhigh). Verification: Claude Opus 5.5
(effort xhigh), fresh context, no access to the solver's reasoning (`verification.md`, own checker
`verification/sawcheck.cpp`).

## Problem

Let \(c_n\) be the number of \(n\)-step self-avoiding walks on \(\mathbb Z^3\) and \(\mu=\lim c_n^{1/n}\). The best
rigorous upper bound, unchanged since 2000, is \(\mu\le 4.7387\) (A. Pönitz, P. Tittmann, *Improved upper bounds
for self-avoiding walks in Z^d*, Electron. J. Combin. 7 (2000) R21, Table 2, loop memory k = 14). It is still
quoted as the state of the art in 2025:

> "Current state-of-the-art for rigorous bounds for connective constants for SAWs is 2.62002 ≤ µ ≤ 2.66235 on the
> square lattice (Couronné, 2022), and 4.57213 ≤ µ ≤ 4.73871 for the cubic lattice (Hara et al., 1993, Pönitz and
> Tittmann, 2000)."
> — S. Kim, F. Pinna, arXiv:2509.04568 (2025), introduction.

A circulating "µ < 4.7114" is erroneous (Clisby, Liang, Slade, J. Phys. A 40 (2007) 10973, §7.2: "in fact c26
gives the weaker bound µ ≤ 4.7626"). Novelty search: `references.md`.

## Status

**Record: new rigorous bounds with exact certificates, checked by the solver's checker and re-checked by the
verifier's own checker (PROTOCOL §4), with identical exact maximal ratios.**

| lattice | previous bound | new bound (certified) | memory K | states |
|---|---|---|---|---|
| \(\mathbb Z^3\) | 4.7387 (P–T 2000) | **4.722048** | 18 | 206 118 560 |
| \(\mathbb Z^3\) | | 4.726984 (smaller certificate) | 16 | 11 641 894 |
| \(\mathbb Z^4\) | 6.8040 (P–T 2000) | **6.801078** | 14 | 4 235 309 |
| \(\mathbb Z^5\) | 8.8602 (P–T 2000) | **8.850674** | 14 | 7 787 256 |
| \(\mathbb Z^6\) | 10.8886 (Noonan 1998) | **10.885008** | 12 | 192 224 |

(Numerical estimate \(\mu_3\approx4.68404\): the rigorous gap shrinks from 0.0547 to 0.0380, i.e. by 30 %.)

All five certificates passed `src/verify_cert.cpp` (logs in `certificates/`): K = 18 in 4 min 47 s on the
explicit-automaton file (8.2 GB), K = 16 in 13 min, d = 4..6 in 32 min or less (one thread). The verifier's checker (`verification/sawcheck.cpp`, written from the mathematical specification) passes all of them
with 4 threads: K = 18 in 151 s (3.3 GB, streamed), K = 16 in 384 s, d = 4, K = 14 in 91 s, d = 5, K = 14 in 193 s,
d = 6, K = 12 in 2 s; see `verification.md`.

## Result in one paragraph

We use finite "forbidden-set" automata whose paths dominate self-avoiding walks: a state remembers some visited
sites (up to the lattice symmetries), steps onto remembered sites are forbidden, and after each step the automaton
may forget any sites. The growth rate is bounded by the Collatz–Wielandt ratio \(\max_s(AV)_s/V_s\) of an integer
vector \(V>0\), checked exactly. The states are those of Pönitz–Tittmann (walks of loop size \(\le K\)); the gain
comes from letting the automaton choose, when it must forget, among all sub-walks through contacts of the
remembered walk, the choice being optimised against the Perron vector (min-sum value / policy iteration, after
Couronné's idea for \(\mathbb Z^2\)). At equal memory this gains about 0.005 in \(d=3\); with K = 16 (Pönitz–Tittmann
stopped at 14) the bound becomes 4.726984, and with K = 18 (2.06·10⁸ states) 4.722048. The method reproduces
the Pönitz–Tittmann table (d = 2..5) and Noonan's d = 6 value when restricted to their forgetting rule; the
verifier's own implementation of the P–T automaton reproduces every entry of their Table 2, including
d = 4, k = 14 (6.8039597 → 6.8040).

## How to reproduce

On Linux (g++, Python 3), from this directory:

```bash
bash scripts/build.sh
bash scripts/check_small.sh     # ~10 min
bash scripts/reproduce_k16.sh   # ~20 min, 3.5 GB
bash scripts/reproduce_hd.sh    # ~45 min, 3.6 GB
bash scripts/reproduce_k18.sh   # ~3.5 h, 10.5 GB, ~45 GB temporary disk
bash verification/run_all.sh    # verifier's checks on all certificates, ~35 min, 4 threads, 3.4 GB
```

The large certificates (K = 16, K = 18, d = 4 and 5 with K = 14; up to 8.2 GB) are not stored in the repository;
their SHA-256 sums are in `certificates/SHA256SUMS.txt` and the `reproduce_*` scripts regenerate them
deterministically.

`check_small.sh` re-checks the certificates stored in `certificates/` (d = 3 with K ≤ 14, d = 5, 6 with K ≤ 12)
with the C++ checker, the brute-force Python checker (K ≤ 10) and the tamper tests, printing PASS/FAIL. The
`reproduce_*` scripts regenerate the large certificates deterministically and check them; binaries and large files
go to `~/work/cubic-saw/`.

## Files

* `paper.tex`, `paper.pdf` — manuscript (generated from `notes.md`).
* `notes.md` — full write-up: lemmas and proofs, automaton, certificate formats, checks, results, what to verify.
* `references.md` — sources and novelty search.
* `verification.md` — verification report; `verification/` — the verifier's programs (separate
  checker `sawcheck.cpp`, `count_saw.cpp`, `pt_auto.cpp`, `tamper.py`, `compare_counts.py`, `patch_k18.py`,
  `run_all.sh`) and logs.
* `src/saw_prod.cpp` — generator storing all candidate children (K ≤ 16, d ≤ 6). `src/saw_big.cpp` — generator for
  K = 18 (policy iteration, candidates cached on disk, explicit-automaton certificate). `src/verify_cert.cpp` —
  checker. `src/saw_automaton.cpp`, `src/count_states.cpp` — exploration tools (Pönitz–Tittmann reproduction,
  variants that were tried).
* `scripts/` — `build.sh`, `check_small.sh`, `reproduce_k16.sh`, `reproduce_k18.sh`, `reproduce_hd.sh`,
  `verify_small_bruteforce.py` (pure-Python checker), `tamper_test.py`, `check_counts.py` (a_n ≥ c_n).
* `certificates/` — small certificates (`*.cert.gz`, `*.cert2.gz`), generation and verification logs,
  `SHA256SUMS.txt` (also for the large certificates kept in `~/work/cubic-saw/`), `counts/` (exact c_n tables of
  Clisby–Liang–Slade, d = 3..6).
