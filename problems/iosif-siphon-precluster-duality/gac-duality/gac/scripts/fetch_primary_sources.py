"""Fetch open primary sources, into git-ignored local storage only.

Environment: Python 3.12+, standard library, network access; Poppler pdftotext.
Expected runtime: 5--90 s, network dependent. PASS/FAIL per download/conversion.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import platform
import subprocess
import sys
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "certificates" / "sources"
SOURCES = {
    "anderson-2011": "https://arxiv.org/pdf/1101.0761",
    "pantea-2012": "https://arxiv.org/pdf/1103.0603",
    "cnp-2013": "https://people.math.wisc.edu/~craciun/PAPERS/Craciun_Nazarov_Pantea_PERSISTENCE_PERMANENCE.pdf",
    "gms-2014": "https://arxiv.org/pdf/1305.5303",
    "anderson-shiu-2010": "https://web.math.wisc.edu/~dfanderson/papers/AndShiuFacetsSIAM.pdf",
    "xu-2026": "https://arxiv.org/pdf/2409.01598",
    "wiuf-2026": "https://arxiv.org/pdf/2609.24553v1",
    "craciun-v3-2026": "https://arxiv.org/pdf/1501.02860v3",
    "adls-2007": "https://www.sontaglab.org/FTPDIR/angeli_leenheer_sontag_math_biosciences_MBS-D-06-00188R1.pdf",
    "craciun-v2-2016": "https://arxiv.org/pdf/1501.02860v2",
    "openai-wr-2026": "https://raw.githubusercontent.com/openai/math/main/preprints/Boundedness-and-persistence-of-weakly-reversible-mass-action-systems-September-25-2026/paper.pdf",
    "openai-permanence-2026": "https://raw.githubusercontent.com/openai/math/main/preprints/Uniform-Permanence-in-Weakly-Reversible-Mass-Action-Systems-October-5-2026/permanence.pdf",
}

def fetch(item):
    key, url = item
    try:
        path = OUT / (key + ".pdf")
        if not path.exists():
            req = urllib.request.Request(url, headers={"User-Agent": "research-source-check/1.0"})
            with urllib.request.urlopen(req, timeout=55) as response:
                data = response.read()
            if not data.startswith(b"%PDF"):
                raise ValueError("Not a PDF")
            path.write_bytes(data)
        subprocess.run(["pdftotext", "-layout", str(path), str(path.with_suffix(".txt"))], check=True, capture_output=True)
        return True, f"PASS {key}: PDF fetched and text extracted"
    except Exception as exc:
        return False, f"FAIL {key}: {exc}"

if __name__ == "__main__":
    print(f"Environment: Python {platform.python_version()}, {platform.platform()}; expected runtime 5--90 s", flush=True)
    start = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(fetch, SOURCES.items()))
    for ok, message in results:
        print(message)
    print(f"Runtime: {time.perf_counter() - start:.3f} s")
    sys.exit(0 if all(ok for ok, _ in results) else 1)
