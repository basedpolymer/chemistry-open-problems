#!/usr/bin/env python3
"""Read-only public API audit; failures and incomplete indexing are recorded.
Python standard library, ~30 seconds. Not a proof of priority.
"""
import concurrent.futures as cf
import datetime,json,urllib.request,urllib.parse
from pathlib import Path

quote=urllib.parse.quote
queries={
 'openalex_citing':'https://api.openalex.org/works?filter=cites:W4390503227&per-page=100',
 'openalex_topic':'https://api.openalex.org/works?search='+quote('absolute concentration robustness multistability')+'&per-page=50',
 'semantic_scholar':'https://api.semanticscholar.org/graph/v1/paper/ARXIV:2307.04186?fields=title,citations.title,citations.year,citations.url',
 'zenodo':'https://zenodo.org/api/records?q='+quote('"absolute concentration robustness" multistability')+'&size=100',
 'github_topic':'https://api.github.com/search/repositories?q='+quote('"absolute concentration robustness"')+'&per_page=100',
 'github_four':'https://api.github.com/search/repositories?q='+quote('"four reactant"')+'&per_page=100',
 'formal_conjectures':'https://api.github.com/repos/google-deepmind/formal-conjectures/git/trees/main?recursive=1'}
for owner in ('mmislan','ipitchford'):
 for page in range(1,5): queries[f'{owner}_{page}']=f'https://api.github.com/users/{owner}/repos?per_page=100&page={page}'

def get(item):
 name,url=item
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'chemistry-open-problems-novelty-audit/1.0','Accept':'application/json'})
  with urllib.request.urlopen(req,timeout=25) as r: data=json.load(r)
  if name.startswith('openalex'):
   data={'count':data['meta']['count'],'items':[{'title':x['display_name'],'year':x['publication_year'],'doi':x.get('doi')} for x in data['results']]}
  elif isinstance(data,list):
   data={'count':len(data),'items':[{'name':x['full_name'],'description':x['description'],'url':x['html_url']} for x in data]}
  elif name.startswith('github'):
   data={'count':data['total_count'],'items':[{'name':x['full_name'],'description':x['description'],'url':x['html_url']} for x in data['items']]}
  elif name=='formal_conjectures':
   data={'truncated':data.get('truncated'),'matching_paths':[x['path'] for x in data.get('tree',[]) if any(w in x['path'].lower() for w in ('reaction','kaihn','multistab','robust'))]}
  elif name=='zenodo':
   data={'count':data['hits']['total'],'items':[{'id':x['id'],'title':x['metadata']['title']} for x in data['hits']['hits']]}
  return name,{'url':url,'status':'ok','result':data}
 except Exception as exc: return name,{'url':url,'status':'unavailable','error':str(exc)}

with cf.ThreadPoolExecutor(max_workers=5) as pool: result=dict(pool.map(get,queries.items()))
out=Path(__file__).resolve().parents[1]/'certificates';out.mkdir(exist_ok=True)
(out/'novelty-api.json').write_text(json.dumps({'date':datetime.datetime.now(datetime.timezone.utc).isoformat(),'queries':result},indent=2)+'\n')
for name,value in result.items():
 print(name,value['status'],value.get('result',{}).get('count',''),value.get('error',''),flush=True)
 for row in value.get('result',{}).get('items',[]):
  if name.startswith('openalex') or any(w in str(row).lower() for w in ('acr','multistab','robust','kaihn')): print(row,flush=True)
