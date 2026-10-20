#!/usr/bin/env python3
"""Read-only literature audit. Python 3 stdlib; typical runtime 1-5 minutes.

At most three concurrent HTTP requests. Network results are date-dependent,
not mathematical certificates. Raw third-party material stays in the ignored
certificates/sources/ folder. Run from this problem folder.
"""
import concurrent.futures
import hashlib
import json
import subprocess
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "certificates" / "sources"
OUT.mkdir(parents=True, exist_ok=True)

URLS = {
    "jang-tran-arxiv.pdf": "https://arxiv.org/pdf/2409.17974",
    "jang-tran-journal.pdf": "https://link.springer.com/content/pdf/10.1186/s13662-025-03946-4.pdf",
    "van-2026.pdf": "https://arxiv.org/pdf/2607.02293",
    "openalex-journal.json": "https://api.openalex.org/works/https://doi.org/10.1186/s13662-025-03946-4",
    "openalex-search.json": "https://api.openalex.org/works?" + urllib.parse.urlencode({
        "search": "Discrete coagulation fragmentation multiplicative constant Jang Tran", "per-page": 100}),
    "semantic-arxiv.json": "https://api.semanticscholar.org/graph/v1/paper/ARXIV:2409.17974?fields=title,year,externalIds,citations.title,citations.year,citations.externalIds,citations.url",
    "semantic-journal.json": "https://api.semanticscholar.org/graph/v1/paper/DOI:10.1186/s13662-025-03946-4?fields=title,year,externalIds,citations.title,citations.year,citations.externalIds,citations.url",
    "zenodo-title.json": "https://zenodo.org/api/records?" + urllib.parse.urlencode({
        "q": '"coagulation" AND "fragmentation" AND ("Jang" OR "multiplicative")',
        "size": 100}),
    "zenodo-arxiv.json": "https://zenodo.org/api/records?" + urllib.parse.urlencode({
        "q": '"2409.17974" OR "s13662-025-03946-4"', "size": 100}),
    "arxiv-search.xml": "https://export.arxiv.org/api/query?" + urllib.parse.urlencode({
        "search_query": 'all:coagulation AND all:fragmentation AND submittedDate:[202401010000 TO 202610202359]',
        "start": 0, "max_results": 200, "sortBy": "submittedDate", "sortOrder": "descending"}),
}


def fetch(item):
    name, url = item
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "LiteratureAudit/1.0"})
        with urllib.request.urlopen(req, timeout=45) as r:
            data = r.read()
            status = r.status
        if name.endswith(".pdf") and not data.startswith(b"%PDF-"):
            (OUT / (name + ".html")).write_bytes(data)
            return {"file": name + ".html", "url": url, "http_status": status,
                    "ok": False, "error": "Response is not a PDF"}
        (OUT / name).write_bytes(data)
        return {"file": name, "url": url, "http_status": status,
                "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data), "ok": True}
    except Exception as exc:
        return {"file": name, "url": url, "ok": False, "error": str(exc)}


def main():
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch, URLS.items()))
    # OpenAlex citation lookup depends on the journal record's identifier.
    journal = OUT / "openalex-journal.json"
    if journal.exists():
        data = json.loads(journal.read_text(encoding="utf-8"))
        if data.get("id"):
            url = "https://api.openalex.org/works?" + urllib.parse.urlencode({
                "filter": "cites:" + data["id"].rsplit("/", 1)[-1], "per-page": 200})
            records.append(fetch(("openalex-citations.json", url)))
    commands = {
        "github-ipitchford.json": ["gh", "api", "users/ipitchford/repos?per_page=100", "--paginate", "--slurp"],
        "github-mmislan.json": ["gh", "api", "users/mmislan/repos?per_page=100", "--paginate", "--slurp"],
        "github-search-coagulation.json": ["gh", "search", "repos", "coagulation fragmentation", "--limit", "100", "--json", "fullName,description,url"],
        "github-search-jang-tran.json": ["gh", "search", "repos", '"Jang" "Tran" coagulation', "--limit", "100", "--json", "fullName,description,url"],
        "github-formal-conjectures-tree.json": ["gh", "api", "repos/google-deepmind/formal-conjectures/git/trees/main?recursive=1"],
    }
    for name, cmd in commands.items():
        try:
            r = subprocess.run(cmd, capture_output=True, timeout=60)
            if r.returncode == 0:
                (OUT / name).write_bytes(r.stdout)
            records.append({"file": name, "command": cmd, "ok": r.returncode == 0,
                            "error": r.stderr.decode(errors="replace").strip()})
        except Exception as exc:
            records.append({"file": name, "command": cmd, "ok": False, "error": str(exc)})
    (OUT / "audit-manifest.json").write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    for row in records:
        print(("PASS" if row["ok"] else "FAIL") + " retrieval: " + row["file"]
              + ("" if row["ok"] else " (" + row.get("error", "") + ")"))
    print("NOTE: retrieval PASS does not certify novelty; inspect every result.")


if __name__ == "__main__":
    main()
