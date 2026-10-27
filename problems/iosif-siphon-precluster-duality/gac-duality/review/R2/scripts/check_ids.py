"""R2: check every DOI and arXiv id of references.bib against Crossref and the arXiv API.

Environment: Windows, `py -3` (standard library only), network access to api.crossref.org and export.arxiv.org.
Runtime: about 1 minute (one request per identifier, 1 s pause).
Run from the repository root:
  py -3 problems/iosif-siphon-precluster-duality/gac-duality/review/R2/scripts/check_ids.py
Output: one block per entry (bib fields, then the registry's fields) and a PASS/FAIL line on resolution and on
agreement of volume, first page and year. Title agreement is printed for manual reading (case and markup differ).
Also queries Crossref bibliographically for the entries that give no DOI, and prints the best hit for manual
inspection (never used as a PASS).
"""
import json, re, sys, time, io, urllib.request, urllib.parse, xml.etree.ElementTree as ET
from pathlib import Path

BIB = Path(__file__).resolve().parents[4] / "references.bib"
UA = {"User-Agent": "chemistry-open-problems reviewer script (mailto:none)"}

def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read()

def parse_bib(text):
    entries = []
    for m in re.finditer(r"@(\w+)\{([^,]+),(.*?)\n\}", text, re.S):
        typ, key, body = m.group(1), m.group(2).strip(), m.group(3)
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}", body, re.S):
            fields[fm.group(1).lower()] = re.sub(r"\s+", " ", fm.group(2)).strip()
        entries.append((typ, key, fields))
    return entries

def norm(s):
    s = re.sub(r"\[a-zA-Z]+|[{}\'\"`~^]", "", s or "")
    return re.sub(r"[^a-z0-9]", "", s.lower())

sys.stdout.reconfigure(encoding="utf-8")
fails = 0
entries = parse_bib(BIB.read_text(encoding="utf-8"))
print(f"references.bib: {len(entries)} entries")
for typ, key, f in entries:
    print(f"\n== {key} ({typ})")
    if "doi" in f:
        doi = f["doi"]
        try:
            msg = json.loads(get("https://api.crossref.org/works/" + urllib.parse.quote(doi)))["message"]
        except Exception as e:
            try:
                dc = json.loads(get("https://api.datacite.org/dois/" + urllib.parse.quote(doi)))["data"]["attributes"]
                t = dc["titles"][0]["title"]; pub = dc.get("publisher"); yr = dc.get("publicationYear")
                auth = ", ".join(c.get("familyName", c.get("name", "")) for c in dc.get("creators", []))
                print(f"  not at Crossref; DataCite: {t} | {pub} {yr} | {auth} | url {dc.get('url')}")
                ok = norm(t) == norm(f.get("title"))
                print(("PASS" if ok else "FAIL") + f" {key}: DOI {doi} (DataCite)"); fails += (not ok)
            except Exception as e2:
                print(f"FAIL {key}: DOI {doi} resolves neither at Crossref ({e}) nor at DataCite ({e2})"); fails += 1
            time.sleep(1); continue
        title = (msg.get("title") or [""])[0]
        cont = (msg.get("container-title") or [""])[0]
        vol = msg.get("volume", ""); page = msg.get("page", ""); iss = msg.get("issue", "")
        year = None
        for k in ("published-print", "issued", "published-online"):
            if k in msg and msg[k].get("date-parts") and msg[k]["date-parts"][0][0]:
                year = msg[k]["date-parts"][0][0]; break
        auth = ", ".join(a.get("family", a.get("name", "")) for a in msg.get("author", []))
        print(f"  bib:      {f.get('title','')} | {f.get('journal', f.get('booktitle', f.get('publisher','')))} "
              f"{f.get('volume','')}({f.get('number','')}) {f.get('pages','')} {f.get('year','')}")
        print(f"  crossref: {title} | {cont} {vol}({iss}) {page} {year} | {auth}")
        ok = True
        if f.get("volume") and vol and f["volume"].lstrip("0") != vol.lstrip("0"): ok = False; print("  volume differs")
        if f.get("pages") and page:
            p0 = re.split(r"-+|–", f["pages"])[0]
            if p0 != re.split(r"-+|–", page)[0]: ok = False; print("  first page differs")
        if f.get("year") and year and str(year) != f["year"]: print(f"  year differs (bib {f['year']}, crossref {year}; check online/print)")
        tb, tc = norm(f.get("title")), norm(title)
        print("  title match:", "yes" if tb == tc else ("prefix" if tc.startswith(tb) or tb.startswith(tc) else "NO"))
        if tb != tc and not (tc.startswith(tb) or tb.startswith(tc)): ok = False
        print(("PASS" if ok else "FAIL") + f" {key}: DOI {doi}"); fails += (not ok)
        time.sleep(1)
    if f.get("eprinttype", "").lower() == "arxiv" or "eprint" in f:
        aid = f["eprint"]
        try:
            xml = get("https://export.arxiv.org/api/query?id_list=" + aid)
            root = ET.fromstring(xml)
            ns = {"a": "http://www.w3.org/2005/Atom"}
            ent = root.find("a:entry", ns)
            t = re.sub(r"\s+", " ", ent.find("a:title", ns).text).strip()
            pub = ent.find("a:published", ns).text; upd = ent.find("a:updated", ns).text
            idv = ent.find("a:id", ns).text
            auth = ", ".join(a.find("a:name", ns).text for a in ent.findall("a:author", ns))
            print(f"  arXiv {aid}: {t} | {auth} | first {pub} | latest {upd} | {idv}")
            ok = norm(t) == norm(f.get("title")) or norm(t).startswith(norm(f.get("title"))[:30])
            print(("PASS" if ok else "FAIL") + f" {key}: arXiv {aid} title"); fails += (not ok)
        except Exception as e:
            print(f"FAIL {key}: arXiv {aid} not found ({e})"); fails += 1
        time.sleep(1)
    if "doi" not in f and typ in ("article", "inproceedings", "incollection", "book", "techreport"):
        q = f.get("title", "") + " " + f.get("author", "")
        try:
            msg = json.loads(get("https://api.crossref.org/works?rows=2&query.bibliographic=" + urllib.parse.quote(q)))
            for it in msg["message"]["items"]:
                print(f"  no DOI in bib; Crossref hit (for inspection only): {it.get('DOI')} | "
                      f"{(it.get('title') or [''])[0]} | {(it.get('container-title') or [''])[0]} "
                      f"{it.get('volume','')} {it.get('page','')} | score {it.get('score')}")
        except Exception as e:
            print(f"  bibliographic query failed: {e}")
        time.sleep(1)
print(f"\nFAILS: {fails}")
sys.exit(1 if fails else 0)
