"""Read-only Zenodo topic queries, Python 3.12+ stdlib, normally <2 minutes.
PASS/FAIL means request success only, not novelty. Failures exit nonzero.
"""
import datetime
import json
from pathlib import Path
import urllib.parse
import urllib.request

out = {}
for term in ['multistage utility', 'minimum utility cost', 'heat exchanger complexity', 'multistage AND utility', '"minimum utility cost"', '"heat exchanger" AND complexity', '"1909.12328"']:
    url = 'https://zenodo.org/api/records?' + urllib.parse.urlencode({'q': term, 'size': 10}, quote_via=urllib.parse.quote)
    try:
        with urllib.request.urlopen(url, timeout=25) as response:
            data = json.load(response)
        out[term] = {'url': url, 'status': 'ok', 'data': data}
        print('PASS', term, 'hits', data['hits']['total'])
        for hit in data['hits']['hits']:
            print(' ', hit['id'], hit['metadata']['title'].encode('ascii', errors='replace').decode()[:160])
    except Exception as exc:
        out[term] = {'url': url, 'status': 'unavailable', 'error': str(exc)}
        print('FAIL', term, str(exc))
target = Path(__file__).resolve().parents[1] / 'certificates/zenodo-topics.json'
target.write_text(json.dumps({'date_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'queries': out}, indent=2) + '\n', encoding='utf-8')
raise SystemExit(int(any(x['status'] != 'ok' for x in out.values())))
