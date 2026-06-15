import requests

url = "https://www.encodeproject.org/search/?type=Experiment&target.label=H3K4me3&biosample_ontology.term_name=embryonic+stem+cell&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&format=json"
res = requests.get(url).json()

found = False
for exp in res.get('@graph', []):
    exp_res = requests.get(f"https://www.encodeproject.org{exp['@id']}?format=json").json()
    for file in exp_res.get('files', []):
        if file.get('file_format') == 'bed' and file.get('assembly') == 'mm9' and 'peak' in file.get('file_format_type', '').lower():
            print('Found H3K4me3:', file['href'])
            found = True
            break
    if found: break

if not found:
    print('Not found')
