"""R2: direct Crossref lookup of the DOIs that references.bib does not give (candidates found by check_ids.py's
bibliographic query), and resolution of the URLs and handles of references.bib.

Environment: Windows, `py -3` (standard library only), network. Runtime: about 30 s.
Run from the repository root:
  py -3 problems/iosif-siphon-precluster-duality/gac-duality/review/R2/scripts/check_ids2.py
"""
import json, sys, time, urllib.request, urllib.parse
sys.stdout.reconfigure(encoding="utf-8")
UA = {"User-Agent": "chemistry-open-problems reviewer script"}
def get(url, method="GET"):
    req = urllib.request.Request(url, headers=UA, method=method)
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.status, r.geturl(), r.read() if method == "GET" else b""

CAND = {
    "DeselEsparza1995": "10.1017/CBO9780511526558",
    "CamposChiolaSilva1991": "10.1109/9.106153",
    "RecaldeJulvezSilva2002": "10.3182/20020721-6-ES-1901.00078",
    "MaedaKodamaOhta1978": "10.1109/TCS.1978.1084490",
    "Moreau2005": "10.1109/TAC.2004.841888",
    "CraciunNazarovPantea2013": "10.1137/100812355",
    "ZhangGaoDochain2023": "10.1109/TAC.2023.3275880",
}
fails = 0
for key, doi in CAND.items():
    try:
        _, _, b = get("https://api.crossref.org/works/" + urllib.parse.quote(doi))
        m = json.loads(b)["message"]
        yr = (m.get("published-print") or m.get("issued"))["date-parts"][0]
        auth = ", ".join(a.get("given", "") + " " + a.get("family", "") for a in m.get("author", []))
        print(f"PASS {key}: {doi} -> {(m.get('title') or [''])[0]} | {(m.get('container-title') or [''])[0]} "
              f"vol {m.get('volume','')} no {m.get('issue','')} pp {m.get('page','')} | {yr} | {auth} | ISBN {m.get('ISBN','')}")
    except Exception as e:
        print(f"FAIL {key}: {doi} ({e})"); fails += 1
    time.sleep(1)

URLS = {
    "Iosif2020 HAL v1": "https://hal.science/hal-03058670v1",
    "Iosif2024talk slides": "https://gasiull.webs.ull.es/beamer-2024-11-19.pdf",
    "Hack1972 handle": "https://hdl.handle.net/1721.1/149406",
    "CamposChiolaSilva1991 author copy": "https://webdiis.unizar.es/~jcampos/web/pdf/ccs-tac-91.pdf",
    "OpenAI2026persistence repository": "https://github.com/openai/math",
    "Iosif thesis landing (DataCite url)": "https://opendata.uni-halle.de//handle/1981185920/14078",
}
for key, url in URLS.items():
    try:
        st, final, _ = get(url, "HEAD")
        print(f"PASS {key}: {url} -> HTTP {st}, final {final}")
    except Exception as e:
        try:
            st, final, _ = get(url)
            print(f"PASS {key}: {url} -> HTTP {st} (GET), final {final}")
        except Exception as e2:
            print(f"FAIL {key}: {url} ({e2})"); fails += 1
    time.sleep(1)
print(f"FAILS: {fails}")
sys.exit(1 if fails else 0)
