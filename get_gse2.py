import urllib.request
import re

url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12867480/fullTextXML'
try:
    xml = urllib.request.urlopen(url).read().decode('utf-8')
    # search for Data Availability
    match = re.search(r'<sec sec-type="data-availability".*?</sec>', xml, re.DOTALL)
    if match:
        text = re.sub(r'<[^>]+>', ' ', match.group(0))
        print("Data Availability:")
        print(text)
    else:
        print("No Data Availability section found.")
        
    print("Mentions of Nanog GSE or SRR:")
    for m in re.finditer(r'.{0,50}GSE\d+.{0,50}', xml):
        print(m.group(0))
    for m in re.finditer(r'.{0,50}PRJNA\d+.{0,50}', xml):
        print(m.group(0))
    for m in re.finditer(r'.{0,50}E-MTAB-\d+.{0,50}', xml):
        print(m.group(0))
except Exception as e:
    print(e)
