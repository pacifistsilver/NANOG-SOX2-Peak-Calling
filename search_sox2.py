import requests

url_sox2 = "https://www.encodeproject.org/search/?type=Experiment&target.label=SOX2&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&format=json"
res = requests.get(url_sox2).json()
for exp in res.get('@graph', []):
    exp_res = requests.get(f"https://www.encodeproject.org{exp['@id']}?format=json").json()
    for file in exp_res.get('files', []):
        if file.get('file_format') == 'bed' and file.get('assembly') == 'mm10':
            print('SOX2:', file['href'])
            exit(0)
print('SOX2 not found')
