#!/usr/bin/env python3
"""Bounded read-only novelty audit; Python 3.12+ stdlib and authenticated gh CLI.
Runtime: normally 1--3 minutes (network dependent), Windows or Linux.
PASS/FAIL denotes request success, never proof of novelty. API failures are logged.
Run from anywhere: python scripts/novelty_audit.py
"""
import concurrent.futures
import datetime
import json
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'certificates' / 'novelty-api.json'
HEADERS = {'User-Agent': 'chemistry-open-problems-literature-audit/1.0',
           'Accept': 'application/json'}
TITLE = 'Approximation algorithms for process systems engineering'
DOI = '10.1016/j.compchemeng.2019.106599'


def url_query(base, query):
    return base + urllib.parse.urlencode(query)


def request(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=35) as r:
        return json.load(r)


def run_job(job):
    name, kind, query = job
    try:
        if kind == 'url':
            data = request(query)
        else:
            result = subprocess.run(query, text=True, capture_output=True, timeout=80, encoding='utf-8')
            if result.returncode:
                raise RuntimeError(result.stderr.strip())
            data = json.loads(result.stdout)
        return name, {'kind': kind, 'query': query, 'status': 'ok', 'data': data}
    except Exception as exc:
        return name, {'kind': kind, 'query': query, 'status': 'unavailable', 'error': str(exc)}


def main():
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    jobs = [
        ('openalex_doi', 'url', 'https://api.openalex.org/works/https://doi.org/' + DOI),
        ('openalex_title', 'url', url_query('https://api.openalex.org/works?', {'search': TITLE, 'per-page': 10})),
        ('semantic_scholar', 'url', 'https://api.semanticscholar.org/graph/v1/paper/DOI:' + DOI + '?fields=title,year,citationCount,citations.title,citations.year,citations.externalIds,citations.url'),
    ]
    for k, phrase in enumerate(['"multistage" AND "utility"', '"1909.12328"', '"minimum utility cost"', '"heat exchanger" AND "complexity"']):
        jobs.append((f'zenodo_{k}', 'url', url_query('https://zenodo.org/api/records?', {'q': phrase, 'size': 100})))
    for k, phrase in enumerate(['"multistage minimum utility"', '"1909.12328"', '"minimum utility cost"', '"heat exchanger" complexity']):
        jobs.append((f'gh_code_{k}', 'gh', ['gh', 'search', 'code', phrase, '--limit', '100', '--json', 'repository,path,url']))
        jobs.append((f'gh_repos_{k}', 'gh', ['gh', 'search', 'repos', phrase, '--limit', '100', '--json', 'fullName,description,url']))
    jobs += [
        ('gh_ipitchford', 'gh', ['gh', 'api', 'users/ipitchford/repos?per_page=100', '--paginate', '--slurp']),
        ('gh_mmislan', 'gh', ['gh', 'api', 'users/mmislan/repos?per_page=100']),
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        results = dict(pool.map(run_job, jobs))
    work = results['openalex_doi'].get('data')
    if work:
        workid = work['id'].rsplit('/', 1)[-1]
        job = ('openalex_citers', 'url', url_query('https://api.openalex.org/works?', {'filter': 'cites:' + workid, 'per-page': 200}))
        name, result = run_job(job)
        results[name] = result
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({'date_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                               'warning': 'Read-only search evidence, not proof of novelty.',
                               'queries': results}, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    for name, item in results.items():
        print(('PASS' if item['status'] == 'ok' else 'FAIL'), name, item.get('error', ''))
        data = item.get('data', {})
        if name.startswith('gh_'):
            flat = [v for page in data for v in page] if name == 'gh_ipitchford' else data
            if isinstance(flat, list):
                print('  returned', len(flat))
                for row in flat:
                    screen = json.dumps({k: row.get(k) for k in ('name', 'fullName', 'description', 'path')}).lower()
                    if any(term in screen for term in ('heat', 'utility', 'pooling', 'letsios', '1909.12328')):
                        print(' ', json.dumps(row if 'path' in row else {k: str(row.get(k, ''))[:250] for k in ('fullName', 'full_name', 'description', 'html_url', 'url')}, ensure_ascii=False))
        elif name.startswith('zenodo'):
            print('  total', data.get('hits', {}).get('total'))
            for row in data.get('hits', {}).get('hits', []):
                print(' ', row['id'], row['metadata']['title'])
        elif name == 'openalex_citers':
            print('  count', data.get('meta', {}).get('count'))
            for row in data.get('results', []):
                print(' ', row['publication_year'], row['doi'], row['display_name'])
    # Search availability is not a mathematical certificate. Exit nonzero on missing sources.
    return int(any(item['status'] != 'ok' for item in results.values()))


if __name__ == '__main__':
    raise SystemExit(main())
