"""R2: Crossref / arXiv check of the three sources that agent L read and listed as prior art for the dual network
(literature/notes.md §10, §18) but that the paper does not cite: Famili-Palsson 2003, Jamshidi-Palsson 2009,
Duvall arXiv:2409.12487.  Environment: Windows `py -3`, standard library, network. Runtime: about 10 s.
Run from the repository root:
  py -3 problems/iosif-siphon-precluster-duality/gac-duality/review/R2/scripts/check_ids3.py
"""
import json, re, sys, time, urllib.request, urllib.parse, xml.etree.ElementTree as ET
sys.stdout.reconfigure(encoding="utf-8")
UA = {"User-Agent": "chemistry-open-problems reviewer script"}
def get(u):
    with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=30) as r: return r.read()
fails = 0
for key, doi in {"FamiliPalsson2003": "10.1016/S0006-3495(03)74450-6", "JamshidiPalsson2009": "10.1016/j.bpj.2009.06.049"}.items():
    try:
        m = json.loads(get("https://api.crossref.org/works/" + urllib.parse.quote(doi)))["message"]
        auth = ", ".join(a.get("given", "") + " " + a.get("family", "") for a in m.get("author", []))
        yr = (m.get("published-print") or m.get("issued"))["date-parts"][0]
        print(f"PASS {key}: {doi} -> {(m.get('title') or [''])[0]} | {(m.get('container-title') or [''])[0]} vol {m.get('volume','')} "
              f"no {m.get('issue','')} pp {m.get('page','')} | {yr} | {auth}")
    except Exception as e:
        print(f"FAIL {key}: {doi} ({e})"); fails += 1
    time.sleep(1)
try:
    root = ET.fromstring(get("https://export.arxiv.org/api/query?id_list=2409.12487"))
    ns = {"a": "http://www.w3.org/2005/Atom"}; e = root.find("a:entry", ns)
    print("PASS Duvall2024: arXiv 2409.12487 ->", re.sub(r"\s+", " ", e.find("a:title", ns).text).strip(), "|",
          ", ".join(a.find("a:name", ns).text for a in e.findall("a:author", ns)), "| first", e.find("a:published", ns).text,
          "| latest", e.find("a:updated", ns).text, "|", e.find("a:id", ns).text)
except Exception as ex:
    print("FAIL Duvall2024", ex); fails += 1
print("FAILS:", fails); sys.exit(1 if fails else 0)
