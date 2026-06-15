import requests
import xml.etree.ElementTree as ET

for gsm in ['GSM1082341', 'GSM1082342']:
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={gsm}[Accession]'
    res = requests.get(url).text
    id_elem = ET.fromstring(res).find('.//Id')
    if id_elem is None:
        print(f"{gsm}: not found")
        continue
    url2 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={id_elem.text}'
    res2 = requests.get(url2).text
    title_elem = ET.fromstring(res2).find('.//Item[@Name="title"]')
    if title_elem is not None:
        print(f"{gsm}: {title_elem.text}")
    else:
        print(f"{gsm}: title not found")
