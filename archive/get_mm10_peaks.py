import requests

def find_file(target):
    url = f"https://www.encodeproject.org/search/?type=Experiment&target.label={target}&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&format=json"
    res = requests.get(url).json()
    for exp in res.get('@graph', []):
        exp_res = requests.get(f"https://www.encodeproject.org{exp['@id']}?format=json").json()
        for file in exp_res.get('files', []):
            if file.get('file_format') == 'bed' and file.get('assembly') == 'mm10' and 'peak' in file.get('file_format_type', '').lower():
                print(f"{target}: https://www.encodeproject.org{file['href']}")
                return

find_file('NANOG')
find_file('SOX2')
