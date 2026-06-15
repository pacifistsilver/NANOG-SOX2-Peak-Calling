import requests

def search_encode(target):
    url = f'https://www.encodeproject.org/search/?type=Experiment&target.investigated_as=transcription+factor&biosample_ontology.term_name=embryonic+stem+cell&assay_title=TF+ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&target.label={target}&format=json'
    res = requests.get(url)
    data = res.json()
    for exp in data['@graph']:
        exp_url = f"https://www.encodeproject.org{exp['@id']}?format=json"
        exp_res = requests.get(exp_url).json()
        for file in exp_res['files']:
            if file['file_format'] == 'bed' and file.get('file_format_type') == 'narrowPeak' and file.get('assembly') == 'mm10':
                print(f"{target} file: https://www.encodeproject.org{file['href']}")
                return
    print(f"No mm10 narrowPeak found for {target}")

search_encode('NANOG')
search_encode('SOX2')
