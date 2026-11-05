"""Read saved search responses without dumping entire API objects. Python stdlib, <1 s."""
import json
import sys
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
d = json.loads((Path(__file__).resolve().parents[1] / 'certificates/novelty-api.json').read_text(encoding='utf-8'))['queries']
for key, item in d.items():
    print(key, item['status'], item.get('error', ''))
    data = item.get('data', {})
    if key.startswith('gh_'):
        rows = [x for page in data for x in page] if key == 'gh_ipitchford' else data
        print('  count', len(rows))
        for row in rows:
            if key in ('gh_ipitchford', 'gh_mmislan'):
                s = str(row.get('name')) + ' ' + str(row.get('description'))
                if not any(t in s.lower() for t in ('heat', 'utility', 'pool', 'letsios')):
                    continue
            print(' ', {k: str(row[k])[:220] for k in ('fullName','full_name','description','path','url','html_url') if k in row})
    elif key == 'openalex_citers':
        print('  count', data.get('meta', {}).get('count'))
        for row in data.get('results', []):
            print(' ', row['publication_year'], row['display_name'], row['doi'])
            index = row.get('abstract_inverted_index') or {}
            text = sorted((pos, word) for word, places in index.items() for pos in places)
            print('  ABSTRACT:', ' '.join(word for _, word in text))
            print('  LOCATIONS:', [(loc['landing_page_url'], loc['pdf_url']) for loc in row['locations']])
