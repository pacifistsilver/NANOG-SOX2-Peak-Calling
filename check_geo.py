import urllib.request
import json
for gse in ['GSE311420', 'GSE203303']:
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={gse}&retmode=json'
    try:
        res = json.loads(urllib.request.urlopen(url).read().decode('utf-8'))
        ids = res['esearchresult']['idlist']
        if ids:
            url2 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={ids[0]}&retmode=json'
            res2 = json.loads(urllib.request.urlopen(url2).read().decode('utf-8'))
            print(f'{gse}: {res2["result"][ids[0]]["title"]}')
    except Exception as e:
        print(f"Error for {gse}: {e}")
