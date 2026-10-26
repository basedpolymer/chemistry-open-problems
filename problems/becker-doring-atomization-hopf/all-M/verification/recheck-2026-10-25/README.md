# Result 8 completion certificates (25–26 October 2026)

`summary.json` assembles all 21 successful jobs of `all-M/scripts/recheck_all.py`.
Each entry names its source execution and log. This is an aggregate after corrections,
not a claim that one invocation passed unchanged. `earlier-attempts/` retains the
previous summaries, failed checks and explanations. The final verdict is in
`../../verification.md`.

The eight generated CSV files contain 337 direct matrix points, 1,423 higher-branch
spectral boxes and 8,614 onset boxes. The original 1,985 compact coefficient boxes
remain in `../../certificates/`. The logs re-prove all of them, together with all
branch links, exact small-size exclusions and formula/limit checks.

`manifest.json` records SHA-256 digests of the source, input certificates, manuscript,
report and archived outputs at completion. Paths are relative to the problem folder.
Text is normalized from CRLF to LF before hashing; `.gz` and `.pdf` files are hashed
byte for byte. The manifest itself is excluded. A later rerun can change generated
files and timings, so these digests identify this archived execution only.

To verify the archive from the repository root:

```python
from pathlib import Path
import hashlib, json
root = Path("problems/becker-doring-atomization-hopf")
archive = root / "all-M/verification/recheck-2026-10-25"
manifest = json.loads((archive / "manifest.json").read_text())
for name, expected in manifest["files"].items():
    data = (root / name).read_bytes()
    if Path(name).suffix not in (".gz", ".pdf"):
        data = data.replace(b"\r\n", b"\n")
    assert hashlib.sha256(data).hexdigest() == expected, name
print("PASS", len(manifest["files"]), "file hashes")
```

To repeat the mathematics, use a fresh output folder:

```bash
python problems/becker-doring-atomization-hopf/all-M/scripts/recheck_all.py /tmp/bd8-recheck
```
