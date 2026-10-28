#!/usr/bin/env python3
"""Read-only, bounded public-source novelty audit. Python 3.12 standard library.
Run: python scripts/novelty_audit.py. Network failures are recorded, not hidden.
Absence of search hits is not proof of novelty. Outputs query URLs and summaries.
"""
import concurrent.futures
import datetime
import json
import urllib.request
import urllib.parse
from pathlib import Path

HEADERS = {'User-Agent':'chemistry-open-problems-research-audit/1.0',
           'Accept':'application/json'}
base = Path(__file__).resolve().parents[1]/'certificates'
base.mkdir(exist_ok=True)
queries = {
 'openalex_title': 'https://api.openalex.org/works?search='+urllib.parse.quote('Coagulation equations with particle emission')+'&per-page=10',
 'semantic_scholar': 'https://api.semanticscholar.org/graph/v1/paper/ARXIV:2604.14076?fields=title,year,citationCount,citations.title,citations.year,citations.url',
 'zenodo_exact': 'https://zenodo.org/api/records?q='+urllib.parse.quote('"2604.14076"')+'&size=10',
 'github_exact': 'https://api.github.com/search/repositories?q='+urllib.parse.quote('"2604.14076"')+'&per_page=20',
 'github_topic': 'https://api.github.com/search/repositories?q='+urllib.parse.quote('coagulation emission')+'&per_page=20',
}
for owner in ['mmislan','ipitchford']:
    for page in range(1,5):
        queries[f'github_{owner}_page{page}']=f'https://api.github.com/users/{owner}/repos?per_page=100&page={page}&sort=updated'

def fetch(item):
    key,url=item
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers=HEADERS),timeout=25) as r:
            data=json.load(r)
        if key.startswith('github_') and isinstance(data,list):
            summary={'count':len(data),'repositories':[{'name':x['full_name'],'description':x['description'],'url':x['html_url']} for x in data]}
        elif key=='openalex_title':
            summary={'count':data['meta']['count'],'results':[{'id':x['id'],'title':x['display_name'],'year':x['publication_year'],'cited_by_count':x['cited_by_count'],'doi':x['doi']} for x in data['results']]}
        elif key.startswith('github_'):
            summary={'count':data['total_count'],'results':[{'name':x['full_name'],'description':x['description'],'url':x['html_url']} for x in data['items']]}
        elif key.startswith('zenodo'):
            summary={'total':data['hits']['total'],'results':[{'id':x['id'],'title':x['metadata']['title']} for x in data['hits']['hits']]}
        else:
            summary=data
        return key,{'url':url,'status':'ok','summary':summary}
    except Exception as e:
        return key,{'url':url,'status':'unavailable','error':str(e)}

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    results=dict(pool.map(fetch,queries.items()))
result={'audited_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope':'Public API queries only; inspect results and record limitations.','queries':results}
(base/'novelty-api.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
for name,item in results.items():
    summary=item.get('summary',{})
    print(name,item['status'],summary.get('count',summary.get('citationCount',summary.get('total',''))),item.get('error',''))
    if isinstance(summary,dict):
        rows=summary.get('repositories',summary.get('results',[]))
        for row in rows:
            value=json.dumps(row).lower()
            if any(x in value for x in ['coagulat','emission','klobus','gelation','rakaus']):
                print(json.dumps(row,ensure_ascii=False))
