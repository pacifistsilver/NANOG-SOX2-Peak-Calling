import urllib.request
import json
import time

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSE203303&retmode=json'
try:
    time.sleep(1)
    res = json.loads(urllib.request.urlopen(url).read().decode('utf-8'))
    ids = res['esearchresult']['idlist']
    if ids:
        url2 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={ids[0]}&retmode=json'
        time.sleep(1)
        res2 = json.loads(urllib.request.urlopen(url2).read().decode('utf-8'))
        title = res2['result'][ids[0]]['title']
        print(f'GSE203303: {title}')
except Exception as e:
    print(f'Error: {e}')
