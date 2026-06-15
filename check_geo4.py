import urllib.request
import json
import time

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=Hendrich+B[Author]+AND+(Nanog+OR+Sox2)&retmode=json'
try:
    time.sleep(1)
    res = json.loads(urllib.request.urlopen(url).read().decode('utf-8'))
    ids = res['esearchresult']['idlist']
    print(ids)
    if ids:
        url2 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={",".join(ids)}&retmode=json'
        time.sleep(1)
        res2 = json.loads(urllib.request.urlopen(url2).read().decode('utf-8'))
        for i in ids:
            try:
                print(res2['result'][i]['accession'], res2['result'][i]['title'])
            except: pass
except Exception as e:
    print(e)
