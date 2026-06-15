import urllib.request
import re
url = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/'
html = urllib.request.urlopen(url).read().decode('utf-8')
files = re.findall(r'href="([^"]*H3k4me3[^"]*Pk\.txt\.gz)"', html, re.IGNORECASE)
for f in set(files):
    if 'E14' in f or 'Bruce' in f or 'Es' in f:
        print(f)
