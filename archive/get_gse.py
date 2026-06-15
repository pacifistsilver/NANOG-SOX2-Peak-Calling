import urllib.request
import re
url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12867480/fullTextXML'
try:
    xml = urllib.request.urlopen(url).read().decode('utf-8')
    gse_matches = list(set(re.findall(r'GSE\d+', xml)))
    print('GSE matches in the paper:', gse_matches)
    for gse in gse_matches:
        url2 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={gse}&retmode=json'
        import json, time
        time.sleep(1)
        res = json.loads(urllib.request.urlopen(url2).read().decode('utf-8'))
        ids = res['esearchresult']['idlist']
        if ids:
            url3 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={ids[0]}&retmode=json'
            time.sleep(1)
            res3 = json.loads(urllib.request.urlopen(url3).read().decode('utf-8'))
            title = res3['result'][ids[0]]['title']
            print(f'{gse}: {title}')
except Exception as e:
    print(e)
