"""Public metadata audit. Python 3 stdlib + GitHub CLI; typically 1-3 minutes.

Network responses are not mathematical certificates. Failures are recorded, not
interpreted as empty search results. Run from this problem folder.
"""
import concurrent.futures
import json
from pathlib import Path
import subprocess
import time
import urllib.request
import urllib.parse
import sys

sys.stdout.reconfigure(encoding='utf-8')

OUT = Path('certificates')
OUT.mkdir(exist_ok=True)

def get(url):
    try:
        req = urllib.request.Request(url, headers={'User-Agent':'research-literature-audit/1.0'})
        with urllib.request.urlopen(req, timeout=35) as r:
            body = r.read().decode()
        try:
            return {'url':url, 'status':'OK', 'data':json.loads(body)}
        except ValueError:
            return {'url':url, 'status':'OK', 'text':body}
    except Exception as e:
        return {'url':url, 'status':'UNAVAILABLE', 'error':str(e)}

def main():
    start = time.monotonic()
    queries = {
        'openalex_klo':'https://api.openalex.org/works/https://doi.org/10.1016/j.jde.2017.04.016',
        'openalex_dulou':'https://api.openalex.org/works?search='+urllib.parse.quote('Proof of a conjecture for the perturbed Gelfand equation from combustion theory')+'&per-page=10',
        'semantic_klo':'https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/j.jde.2017.04.016?fields=title,year,citationCount,citations.title,citations.year,citations.externalIds',
        'semantic_arxiv':'https://api.semanticscholar.org/graph/v1/paper/ARXIV:1612.07722?fields=title,year,citationCount,citations.title,citations.year,citations.externalIds',
        'zenodo_gelfand':'https://zenodo.org/api/records?q='+urllib.parse.quote('"perturbed Gelfand"'),
        'zenodo_critical':'https://zenodo.org/api/records?q='+urllib.parse.quote('Gelfand AND combustion'),
        'arxiv_search':'https://export.arxiv.org/api/query?search_query=all:'+urllib.parse.quote('"perturbed Gelfand"')+'&start=0&max_results=100&sortBy=submittedDate&sortOrder=descending',
    }
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results = dict(zip(queries, pool.map(get, queries.values())))
    for key in ['openalex_klo','openalex_dulou']:
        data = results[key].get('data',{})
        works = data.get('results',[data])
        for work in works:
            if not work.get('id'):
                continue
            if key.endswith('dulou') and not any(a.get('author',{}).get('display_name','').lower().endswith('du') for a in work.get('authorships',[])):
                continue
            wid = work['id'].rsplit('/',1)[-1]
            results['citers_'+wid] = get('https://api.openalex.org/works?filter=cites:'+wid+'&per-page=200')
    for owner in ['ipitchford','mmislan']:
        try:
            cp = subprocess.run(['gh','api',f'users/{owner}/repos?per_page=100','--paginate','--slurp'],capture_output=True,text=True,encoding='utf-8',check=True)
            items = [item for page in json.loads(cp.stdout) for item in page]
            results['github_'+owner] = {'status':'OK','count':len(items),'repositories':[{'name':x['name'],'description':x.get('description'),'url':x['html_url']} for x in items]}
        except Exception as e:
            results['github_'+owner] = {'status':'UNAVAILABLE','error':str(e)}
    for query in ['"perturbed Gelfand"','Gelfand combustion','repo:google-deepmind/formal-conjectures Gelfand']:
        try:
            command = ['gh','search','repos',query,'--limit','100','--json','fullName,description,url'] if not query.startswith('repo:') else ['gh','search','code','Gelfand','--repo','google-deepmind/formal-conjectures','--limit','100','--json','path,url']
            cp = subprocess.run(command,capture_output=True,text=True,encoding='utf-8',check=True)
            results['github_search_'+query] = {'status':'OK','data':json.loads(cp.stdout)}
        except Exception as e:
            results['github_search_'+query] = {'status':'UNAVAILABLE','error':str(e)}
    (OUT/'novelty-api.json').write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for key,r in results.items():
        print(('PASS retrieval ' if r['status']=='OK' else 'FAIL retrieval ')+key)
        d=r.get('data',{})
        if key.startswith('citers_'):
            print('  count:',d.get('meta',{}).get('count'))
            for w in d.get('results',[]):
                print(' ',w.get('publication_year'),w.get('title'),w.get('doi'))
        elif key.startswith('github_') and 'repositories' in r:
            hits=[x for x in r['repositories'] if any(t in (x['name']+' '+str(x['description'])).lower() for t in ['gelfand','combustion','thermal','elliptic','bifurcation'])]
            print('  repositories:',r['count'],'keyword matches:',hits)
    print('Runtime seconds:',round(time.monotonic()-start,2))

if __name__=='__main__':
    main()
