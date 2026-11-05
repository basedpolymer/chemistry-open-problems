"""Read-only audit retries and source downloads. Python stdlib; normally <2 min.
PASS/FAIL reports access only; this is not a novelty or mathematical certificate.
"""
import concurrent.futures
import datetime
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / 'certificates/sources'
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
URLS = {
    'semantic_arxiv': 'https://api.semanticscholar.org/graph/v1/paper/ARXIV:1909.12328?fields=title,citationCount,citations.title,citations.year,citations.externalIds',
    'rankine_review.pdf': 'https://pure.tudelft.nl/ws/portalfiles/portal/232149828/1-s2.0-S0360128524000595-main.pdf',
    'optimization_review.pdf': 'https://fepbl.com/index.php/csitrj/article/download/1027/1250',
    'gh_code_hit.json': 'https://raw.githubusercontent.com/chauhanpriyanshu001/Khudkibook-Web.github.io/webtest/data/syllabus_cache/2711606.json',
    'misener_publications.html': 'https://wp.doc.ic.ac.uk/rmisener/publications/',
}
for n, q in enumerate(['"1909.12328"', '"multistage" AND "utility"', '"minimum utility cost"', '"heat exchanger" AND "complexity"']):
    URLS[f'zenodo_{n}'] = 'https://zenodo.org/api/records?' + urllib.parse.urlencode({'q': q, 'size': 10})


def fetch(item):
    name, url = item
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=40) as r:
            content = r.read()
        if '.' in name:
            (SOURCES / name).write_bytes(content)
            data = {'local': name, 'bytes': len(content)}
        else:
            data = json.loads(content)
        return name, {'url': url, 'status': 'ok', 'data': data}
    except urllib.error.HTTPError as e:
        return name, {'url': url, 'status': 'unavailable', 'error': str(e), 'body': e.read().decode(errors='replace')[:1500]}
    except Exception as e:
        return name, {'url': url, 'status': 'unavailable', 'error': str(e)}


if __name__ == '__main__':
    SOURCES.mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        data = dict(pool.map(fetch, URLS.items()))
    (ROOT / 'certificates/novelty-followup.json').write_text(json.dumps({'date_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'queries': data}, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    for name, entry in data.items():
        print('PASS' if entry['status'] == 'ok' else 'FAIL', name, entry.get('error', ''))
        if name.startswith('zenodo'):
            print('  hits', entry.get('data', {}).get('hits', {}).get('total'), entry.get('body', ''))
        if name == 'semantic_arxiv':
            print(json.dumps(entry.get('data', {}), ensure_ascii=False))
    raise SystemExit(int(any(e['status'] != 'ok' for e in data.values())))
