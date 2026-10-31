#!/usr/bin/env python3
"""Repeat the API part of the novelty search for the McClure-Shiu question.

Environment: Python 3.12, standard library only. Needs network access; `gh` for the GitHub part.
Runtime: about 1 minute. Output: certificates/novelty-api.json (results change with time).
This script is a record of the search, not a mathematical certificate.
"""
import json
import pathlib
import subprocess
import time
import urllib.parse
import urllib.request

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent / "certificates" / "novelty-api.json"
UA = {"User-Agent": "chemistry-open-problems-novelty-audit/1.0"}


def get(url, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # network errors are recorded, not fatal
            err = repr(exc)
            time.sleep(4 * (attempt + 1))
    return json.dumps({"error": err})


def openalex_citers():
    work = json.loads(get("https://api.openalex.org/works/doi:10.1016/j.jsc.2024.102323"))
    oid = work.get("id", "").rsplit("/", 1)[-1]
    out = []
    cursor = "*"
    error = None
    while cursor and oid:
        url = (f"https://api.openalex.org/works?filter=cites:{oid}&per-page=200&cursor={cursor}"
               "&select=id,doi,title,publication_year")
        page = json.loads(get(url))
        error = page.get("error", error)  # the free daily budget can run out
        out += page.get("results", [])
        cursor = page.get("meta", {}).get("next_cursor")
    return {"id": oid, "cited_by_count": work.get("cited_by_count"), "citers": out, "error": error}


def semantic_scholar_citers():
    url = ("https://api.semanticscholar.org/graph/v1/paper/arXiv:2303.03960/citations"
           "?fields=title,year,externalIds&limit=500")
    return json.loads(get(url))


def arxiv_search():
    queries = [
        'all:"multistationarity region" AND all:"one-species"',
        'all:"multistationarity regions" AND all:connected',
        'all:"multistationarity" AND all:"one species" AND all:connected',
        'all:"parameter region of multistationarity" AND all:connectivity',
        'all:"one-species networks"',
        'all:"one-species" AND all:"mass-action" AND all:multistationary',
        'ti:"one-species reaction networks"',
        'au:White_P AND abs:multistationarity',
    ]
    res = {}
    for q in queries:
        url = ("http://export.arxiv.org/api/query?search_query=" + urllib.parse.quote(q)
               + "&start=0&max_results=40&sortBy=submittedDate&sortOrder=descending")
        txt = get(url)
        titles = [t.split("</title>")[0].strip().replace("\n", " ")
                  for t in txt.split("<title>")[2:]]
        ids = [t.split("</id>")[0].strip() for t in txt.split("<id>")[2:]]
        res[q] = list(zip(ids, titles))
        time.sleep(4)
    return res


def zenodo_search():
    queries = ['"multistationarity region"', '"one-species" multistationarity',
               'McClure Shiu multistationarity connected', '"multistationarity" connectedness']
    res = {}
    for q in queries:
        url = "https://zenodo.org/api/records?size=25&q=" + urllib.parse.quote(q)
        data = json.loads(get(url))
        hits = data.get("hits", {}).get("hits", []) if isinstance(data, dict) else []
        res[q] = [(h.get("doi"), h.get("metadata", {}).get("title"),
                   h.get("metadata", {}).get("publication_date")) for h in hits]
    return res


def gh(*args):
    try:
        out = subprocess.run(["gh", *args], capture_output=True, timeout=120).stdout
        return out.decode("utf-8", "replace")
    except Exception as exc:
        return json.dumps({"error": repr(exc)})


def github_search():
    res = {}
    for user in ("mmislan", "ipitchford"):
        txt = gh("api", f"users/{user}/repos?per_page=100", "--paginate", "--jq",
                 ".[] | [.name, .description, .pushed_at] | @tsv")
        rows = [ln.split("\t") for ln in txt.splitlines() if ln.strip()]
        res[user] = {"total": len(rows), "all": rows}
    res["repo_search"] = {}
    for q in ("multistationarity region connected", "one-species multistationarity",
              "McClure Shiu multistationarity", "multistationarity connectedness"):
        res["repo_search"][q] = gh("search", "repos", q, "--limit", "20", "--json",
                                   "fullName,description,updatedAt")
    # Repository search misses it; code search is what found pw/crnt-connectivity.
    res["code_search"] = {}
    for q in ('"multistationarity region" "one-species"', 'multistationarity McClure',
              'multistationarity "one-species" connected'):
        res["code_search"][q] = gh("api", "search/code?q=" + urllib.parse.quote(q), "--jq",
                                   ".items[] | [.repository.full_name, .path] | @tsv")
    res["pw_crnt_connectivity"] = {
        "repo": gh("api", "repos/pw/crnt-connectivity", "--jq",
                   "[.full_name, .description, .created_at, .pushed_at]"),
        "commits": gh("api", "repos/pw/crnt-connectivity/commits", "--jq",
                      ".[] | [.sha, .commit.committer.date] | @tsv"),
    }
    return res


def main():
    result = {
        "date": time.strftime("%Y-%m-%d"),
        "openalex": openalex_citers(),
        "semantic_scholar": semantic_scholar_citers(),
        "arxiv": arxiv_search(),
        "zenodo": zenodo_search(),
        "github": github_search(),
    }
    OUT.write_text(json.dumps(result, indent=1, ensure_ascii=False), encoding="utf-8")
    print("written", OUT)
    print("openalex citers:", len(result["openalex"]["citers"]))


if __name__ == "__main__":
    main()
