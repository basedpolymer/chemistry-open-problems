#!/usr/bin/env python3
"""Read-only novelty audit for the DeLaVina-Waller conjecture (PROTOCOL.md section 1).

Queries OpenAlex, Semantic Scholar, arXiv, Zenodo and GitHub (through `gh` when available),
and writes every query, URL and result to certificates/novelty-api.json.
Results depend on the date of the run. Environment: Python 3 standard library, plus gh if on PATH.
"""
import json, subprocess, sys, time, urllib.parse, urllib.request, os, re, xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "certificates", "novelty-api.json")
UA = {"User-Agent": "chemistry-open-problems-novelty-audit (mailto:noreply@example.org)"}


def get(url, timeout=40, tries=4):
    """GET with exponential backoff, since OpenAlex and Semantic Scholar rate-limit."""
    req = urllib.request.Request(url, headers=UA)
    err = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return {"status": r.status, "body": r.read().decode("utf-8", "replace")}
        except Exception as e:
            err = e
            time.sleep(4 * 2 ** i)
    return {"status": "error", "body": repr(err)}


res = {"date": time.strftime("%Y-%m-%d"), "queries": []}


def add(name, url, summary):
    res["queries"].append({"name": name, "url": url, "summary": summary})
    print(name, "->", json.dumps(summary)[:400])


# OpenAlex: search for works
for q in ["On a conjecture of DeLaVina and Waller",
          "Wiener index diameter 2d+1 cycle conjecture",
          "maximum Wiener index graphs given diameter order 2d+1",
          "DeLaVina Waller conjecture Wiener index"]:
    url = "https://api.openalex.org/works?per-page=15&search=" + urllib.parse.quote(q)
    r = get(url)
    try:
        js = json.loads(r["body"])
        items = [{"title": w.get("title"), "year": w.get("publication_year"), "doi": w.get("doi"),
                  "id": w.get("id"), "cited_by": w.get("cited_by_count")} for w in js["results"]]
    except Exception as e:
        items = {"error": repr(e), "status": r["status"]}
    add("openalex:" + q, url, items)
    time.sleep(1.0)

# OpenAlex: citers of the three key sources, by DOI
for doi in ["10.37236/757", "10.26493/1855-3974.3077.63a", "10.48550/arXiv.2605.24855", "10.1002/jgt.3190080102",
            "10.1017/S0004972718000163", "10.1016/j.amc.2019.03.025", "10.5281/zenodo.22015516", "10.5281/zenodo.22015517"]:
    url = "https://api.openalex.org/works/https://doi.org/" + doi
    r = get(url)
    try:
        js = json.loads(r["body"])
        wid = js["id"].rsplit("/", 1)[-1]
        title = js.get("title")
        cited = js.get("cited_by_count")
        time.sleep(1.0)
        url2 = f"https://api.openalex.org/works?per-page=200&filter=cites:{wid},from_publication_date:2024-01-01"
        r2 = get(url2)
        js2 = json.loads(r2["body"])
        items = [{"title": w.get("title"), "year": w.get("publication_year"), "doi": w.get("doi")}
                 for w in js2["results"]]
        add("openalex-citers:" + doi, url2, {"title": title, "cited_by_total": cited, "citers_since_2024": items})
    except Exception as e:
        add("openalex-citers:" + doi, url, {"error": repr(e), "status": r["status"], "body": r["body"][:200]})
    time.sleep(0.3)

# Semantic Scholar: citers of the arXiv paper
for ident in ["arXiv:2605.24855", "DOI:10.37236/757"]:
    url = f"https://api.semanticscholar.org/graph/v1/paper/{ident}/citations?fields=title,year,externalIds&limit=200"
    r = get(url)
    try:
        js = json.loads(r["body"])
        items = [{"title": c["citingPaper"].get("title"), "year": c["citingPaper"].get("year"),
                  "ids": c["citingPaper"].get("externalIds")} for c in js.get("data", [])
                 if (c["citingPaper"].get("year") or 0) >= 2024]
    except Exception as e:
        items = {"error": repr(e), "status": r["status"], "body": r["body"][:200]}
    add("semanticscholar:" + ident, url, items)
    time.sleep(1.5)

# arXiv API
for q in ['au:"Liu_Mingchang" AND all:"Wiener"', 'all:"DeLaVina" AND all:"Wiener index"',
          'abs:"Wiener index" AND abs:"diameter" AND abs:"cycle" AND abs:"conjecture"',
          'ti:"Wiener index" AND ti:"diameter"',
          'all:"average distance" AND all:"diameter" AND all:"2d+1"',
          'all:"DeLaVina-Waller"', 'all:"DeLaVina and Waller"']:
    url = ("http://export.arxiv.org/api/query?max_results=30&sortBy=submittedDate&sortOrder=descending&search_query="
           + urllib.parse.quote(q))
    r = get(url)
    try:
        root = ET.fromstring(r["body"])
        ns = {"a": "http://www.w3.org/2005/Atom"}
        items = [{"id": e.find("a:id", ns).text, "published": e.find("a:published", ns).text[:10],
                  "title": " ".join(e.find("a:title", ns).text.split())} for e in root.findall("a:entry", ns)]
    except Exception as e:
        items = {"error": repr(e), "status": r["status"]}
    add("arxiv:" + q, url, items)
    time.sleep(3.2)

# Zenodo
for q in ['"DeLaViña–Waller"', '"DeLaVina-Waller"', '"Wiener index" DeLaViña', "DeLaVina Waller conjecture", "Wiener index diameter cycle conjecture 2d+1",
          "Wiener index DeLaVina", "\"Wiener index\" \"diameter\" conjecture proof cycle maximum",
          "2605.24855"]:
    url = "https://zenodo.org/api/records?size=25&q=" + urllib.parse.quote(q)
    r = get(url)
    try:
        js = json.loads(r["body"])
        items = {"total": js["hits"]["total"],
                 "hits": [{"title": h["metadata"].get("title"), "date": h["metadata"].get("publication_date"),
                           "doi": h.get("doi")} for h in js["hits"]["hits"]]}
    except Exception as e:
        items = {"error": repr(e), "status": r["status"], "body": r["body"][:200]}
    add("zenodo:" + q, url, items)
    time.sleep(1)


# GitHub
def gh(args):
    try:
        p = subprocess.run(["gh"] + args, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=180)
        return (p.stdout or "") if p.returncode == 0 else "ERROR: " + (p.stderr or "")[:300]
    except Exception as e:
        return "ERROR: " + repr(e)


KEY = re.compile(r"wiener|delavin|average.dist|plesn|cycle.*diam|diam.*cycle|graffiti", re.I)
for user in ["ipitchford", "mmislan"]:
    out = gh(["api", f"users/{user}/repos?per_page=100", "--paginate", "--jq", ".[] | [.name, (.description // \"\")] | @tsv"])
    if out.startswith("ERROR"):
        add("github-user:" + user, user, out)
    else:
        lines = out.strip().split("\n")
        add("github-user:" + user, user, {"repos": len(lines), "matches": [l for l in lines if KEY.search(l)]})
for q in ["Wiener index DeLaVina", "DeLaVina Waller conjecture", "Wiener index diameter conjecture cycle", "2605.24855",
          "delavina-waller"]:
    out = gh(["search", "repos", q, "--limit", "30", "--json", "fullName,description,updatedAt"])
    add("github-search:" + q, q, out[:3000])
for q in ["2605.24855", "\"DeLaVina\" Wiener", "\"Conjecture 7\" \"Wiener\" DeLaVina"]:
    out = gh(["search", "code", q, "--limit", "20", "--json", "repository,path"])
    add("github-code:" + q, q, out[:3000])
out = gh(["search", "code", "repo:google-deepmind/formal-conjectures", "Wiener", "--limit", "30", "--json", "path"])
add("formal-conjectures:Wiener", "repo:google-deepmind/formal-conjectures", out[:3000])

OUT_TMP = OUT
with open(OUT_TMP, "w", encoding="utf-8") as f:
    json.dump(res, f, indent=1, ensure_ascii=False)
print("written", OUT)
