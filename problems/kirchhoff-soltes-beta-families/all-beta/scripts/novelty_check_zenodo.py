#!/usr/bin/env python3
"""Second pass of the novelty check, Zenodo and Semantic Scholar (read-only API queries; prints all hits).  Not a certificate."""
import json, time, urllib.request, urllib.parse


def get(url, t=40):
    req = urllib.request.Request(url, headers={"User-Agent": "chemistry-open-problems-novelty/1.0"})
    with urllib.request.urlopen(req, timeout=t) as r:
        return r.read().decode("utf8", "replace")


KEY = ('kirchhoff', 'resistance', 'soltes', 'šoltés', 'wiener', 'good vert', 'transmission')
print('=== Zenodo, full listing, keyword filter on title + description ===')
for q in ['Soltes', 'Šoltés', 'Soltes Kirchhoff', 'Šoltés Kirchhoff', 'Soltes Wiener index', '"Wiener index" "vertex removal"',
          '"resistance distance" "vertex deletion"', '"Kirchhoff index" "good vertex"', '"Kirchhoff index" "deleting a vertex"',
          '"Kirchhoff index" "vertex deleted"', '"Kirchhoff index" conjecture graph proportion', '"Kirchhoff index" "open problem"']:
    try:
        total = None
        n_seen = 0
        for page in (1, 2, 3, 4):
            url = 'https://zenodo.org/api/records?' + urllib.parse.urlencode({'q': q, 'size': 25, 'page': page, 'sort': 'mostrecent'})
            d = json.loads(get(url))
            hits = d['hits']['hits']
            total = d['hits']['total']
            n_seen += len(hits)
            for x in hits:
                md = x['metadata']
                txt = (md.get('title', '') + ' ' + md.get('description', '')).lower()
                if ('kirchhoff' in txt and ('vertex' in txt or 'vertices' in txt)) or 'soltes' in txt.replace('š', 's').replace('é', 'e'):
                    print('  ', md.get('publication_date'), x['id'], md.get('title', '')[:110],
                          [c.get('name') for c in md.get('creators', [])][:2])
            if len(hits) < 25:
                break
            time.sleep(1.2)
        print('--', q, 'total', total, 'retrieved', n_seen)
    except Exception as e:
        print('ERR', q, e)
    time.sleep(1.5)
print('=== Semantic Scholar search ===')
for q in ['Kirchhoff index Soltes problem', 'graphs preserving Kirchhoff index vertex removal', 'good vertices Kirchhoff index proportion']:
    try:
        d = json.loads(get('https://api.semanticscholar.org/graph/v1/paper/search?' + urllib.parse.urlencode(
            {'query': q, 'limit': 15, 'fields': 'title,year,externalIds'})))
        print('--', q, d.get('total'))
        for x in d.get('data', []):
            print('  ', x.get('year'), x.get('title', '')[:110])
    except Exception as e:
        print('ERR', q, e)
    time.sleep(12)
