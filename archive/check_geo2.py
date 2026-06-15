import urllib.request
import json
import time

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=Laue+E[Author]+AND+(Nanog+OR+Sox2+OR+pluripotency+OR+CHD4)&retmode=json'
try:
    time.sleep(2)
    res = json.loads(urllib.request.urlopen(url).read().decode('utf-8'))
    ids = res['esearchresult']['idlist']
    if ids:
        url2 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={",".join(ids)}&retmode=json'
        time.sleep(2)
        res2 = json.loads(urllib.request.urlopen(url2).read().decode('utf-8'))
        for i in ids:
            try:
                acc = res2['result'][i]['accession']
                title = res2['result'][i]['title']
                print(f'{acc}: {title}')
            except:
                pass
except Exception as e:
    print(f'Error: {e}')
