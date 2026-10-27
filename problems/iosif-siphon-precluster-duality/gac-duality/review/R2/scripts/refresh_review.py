"""R2: replace section 8.2 of review.md by the current output of apply_corrections.py --markdown.

Environment: any Python 3 (standard library). Runtime: < 1 s. Run from anywhere:
  py -3 problems/iosif-siphon-precluster-duality/gac-duality/review/R2/scripts/refresh_review.py
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
md = subprocess.run([sys.executable, str(HERE / "scripts" / "apply_corrections.py"), "--markdown"],
                    capture_output=True, text=True, encoding="utf-8", check=True).stdout
(HERE / "logs" / "replacements.md").write_text(md, encoding="utf-8")
review = (HERE / "review.md").read_text(encoding="utf-8")
start_marker = "Apply them in any order (no two overlap).\n\n"
end_marker = "\n### 8.3 `references.md`"
a = review.index(start_marker) + len(start_marker)
b = review.index(end_marker)
review = review[:a] + md + review[b:]
(HERE / "review.md").write_text(review, encoding="utf-8")
print("section 8.2 refreshed:", md.count("#### "), "items")
