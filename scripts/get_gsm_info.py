import urllib.request
import json
url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=300288360,300288359,300288358,300288357,300288356,300288355,300288354,300288353,300288352,300288351,300288350,300288349,300288348,300288347,300288346,300288345&retmode=json'
response = urllib.request.urlopen(url)
data = json.loads(response.read().decode('utf-8'))
for uid, res in data['result'].items():
    if uid != 'uids':
        print(f"{res['accession']}: {res['title']}")
