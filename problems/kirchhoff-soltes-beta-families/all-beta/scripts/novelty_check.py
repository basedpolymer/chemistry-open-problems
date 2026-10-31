#!/usr/bin/env python3
"""Fresh novelty check (read-only API queries). Not a certificate; prints what it finds."""
import json, re, time, urllib.request, urllib.parse, sys
def get(url, t=40):
    req = urllib.request.Request(url, headers={"User-Agent": "chemistry-open-problems-novelty/1.0"})
    with urllib.request.urlopen(req, timeout=t) as r:
        return r.read().decode("utf8", "replace")
print("=== OpenAlex ===")
try:
    d = json.loads(get("https://api.openalex.org/works?filter=doi:10.1016/j.amc.2025.129694"))
    for x in d.get("results", []):
        print(x["id"], x.get("cited_by_count"), x.get("title"))
    d = json.loads(get("https://api.openalex.org/works?filter=cites:W4413802645&per-page=50"))
    print("citers meta", d.get("meta", {}).get("count"))
    for x in d.get("results", []):
        print(x["id"], x.get("publication_date"), x.get("title"))
    for q in ["Kirchhoff Šoltés", "Kirchhoff Soltes", "Kirchhoff index good vertices", "Kirchhoff index vertex removal preserving"]:
        d = json.loads(get("https://api.openalex.org/works?search=%s&per-page=8&sort=publication_date:desc" % urllib.parse.quote(q)))
        print("--", q, d.get("meta", {}).get("count"))
        for x in d.get("results", []):
            print("  ", x.get("publication_date"), (x.get("title") or "")[:100])
except Exception as e:
    print("ERR", e)
print("=== Crossref ===")
try:
    d = json.loads(get("https://api.crossref.org/works?query=Kirchhoff+Soltes&rows=15&select=DOI,title,created"))
    for x in d["message"]["items"]:
        print(x["DOI"], x.get("created", {}).get("date-time"), (x.get("title") or [""])[0][:100])
    d = json.loads(get("https://api.crossref.org/works/10.1016/j.amc.2025.129694"))
    print("is-referenced-by-count", d["message"].get("is-referenced-by-count"))
except Exception as e:
    print("ERR", e)
print("=== arXiv ===")
qs = ['all:"Kirchhoff" AND all:"Soltes"', 'all:"Kirchhoff index" AND all:"good vertices"',
      'all:"Kirchhoff index" AND all:"vertex removal"', 'all:"resistance distance" AND all:"Soltes"',
      'all:"Kirchhoff index" AND all:"vertex deletion"', 'all:"Kirchhoff index" AND all:"removal of a vertex"',
      'all:"Soltes"']
for q in qs:
    try:
        t = get("http://export.arxiv.org/api/query?" + urllib.parse.urlencode(
            {"search_query": q, "max_results": 12, "sortBy": "submittedDate", "sortOrder": "descending"}))
        print("--", q)
        for m in re.finditer(r"<entry>.*?</entry>", t, re.S):
            e = m.group(0)
            print("  ", re.search(r"<published>(.*?)</published>", e).group(1)[:10],
                  re.search(r"<id>(.*?)</id>", e).group(1),
                  re.search(r"<title>(.*?)</title>", e, re.S).group(1).replace("\n", " ")[:100])
    except Exception as e:
        print("ERR", q, e)
    time.sleep(4)
print("=== Zenodo ===")
qs = ['+Kirchhoff +Soltes', '+Kirchhoff +"good vertices"', '"Kirchhoff index" vertex deletion unchanged', 'beta-Kirchhoff',
      '"Soltes problem"', 'Soltes', 'Šoltés', '"Kirchhoff index" "rook graph"', '"Kirchhoff index" "vertex removal"',
      '"Kirchhoff index" Soltes graphs proportion', 'Gottwald Kirchhoff Soltes', '"Problem 4" Kirchhoff']
for q in qs:
    try:
        d = json.loads(get("https://zenodo.org/api/records?" + urllib.parse.urlencode({"q": q, "size": 12, "sort": "mostrecent"})))
        print("--", q, "total", d["hits"]["total"])
        for x in d["hits"]["hits"]:
            print("  ", x["metadata"].get("publication_date"), x["id"], x["metadata"].get("title", "")[:100],
                  [c.get("name") for c in x["metadata"].get("creators", [])][:2])
    except Exception as e:
        print("ERR", q, e)
    time.sleep(1)
