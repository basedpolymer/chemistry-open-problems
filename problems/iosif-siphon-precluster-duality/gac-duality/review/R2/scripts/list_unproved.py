"""R2: list the theorem-like statements of paper.tex that are not followed by a proof environment.

Environment: any Python 3 (standard library). Runtime: < 1 s.
Run from the repository root:
  py -3 problems/iosif-siphon-precluster-duality/gac-duality/review/R2/scripts/list_unproved.py
"""
import re
from pathlib import Path

TEX = Path(__file__).resolve().parents[4] / "paper.tex"
t = TEX.read_text(encoding="utf-8")
BS = chr(92)
n = 0
for m in re.finditer(re.escape(BS) + r"begin\{(theorem|proposition|lemma|corollary)\}", t):
    env = m.group(1)
    endtag = BS + "end{" + env + "}"
    end = t.find(endtag, m.end())
    lab = re.search(re.escape(BS) + r"label\{([^}]+)\}", t[m.end():end])
    after = t[end + len(endtag): end + len(endtag) + 60].lstrip()
    n += 1
    if not after.startswith(BS + "begin{proof}"):
        line = t.count("\n", 0, m.start()) + 1
        print(f"line {line}: {env} {lab.group(1) if lab else '?'}: no proof environment follows")
print(f"{n} theorem-like environments scanned")
