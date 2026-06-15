import os
import json
import gzip
import pandas as pd
from intervaltree import Interval, IntervalTree
import urllib.request
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Get SOX2 and NANOG gene targets
cache_file = 'mapping_cache.json'

with open(cache_file, 'r') as f:
    data = json.load(f)
    nanog_ensembl = set(data['nanog'])
    sox2_ensembl = set(data['sox2'])

target_genes = sox2_ensembl.union(nanog_ensembl)
print(f"Total target genes to process: {len(target_genes)}")

# 2. Get mm9 TSS coordinates
ensgene_url = "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/ensGene.txt.gz"
ensgene_file = "ensGene_mm9.txt.gz"
if not os.path.exists(ensgene_file):
    print("Downloading mm9 ensGene...")
    urllib.request.urlretrieve(ensgene_url, ensgene_file)

# We want intervals of +/- 2000 bp around the TSS for each gene
gene_intervals = {} # ensembl_id -> (chrom, start, end)
with gzip.open(ensgene_file, 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[2]
        strand = parts[3]
        txStart = int(parts[4])
        txEnd = int(parts[5])
        gene_id = parts[12]
        
        if gene_id in target_genes:
            tss = txStart if strand == '+' else txEnd
            start = max(0, tss - 2000)
            end = tss + 2000
            
            if gene_id not in gene_intervals:
                gene_intervals[gene_id] = []
            gene_intervals[gene_id].append((chrom, start, end))

print(f"Found mm9 coordinates for {len(gene_intervals)} genes")

# Create IntervalTrees for fast lookup
# We'll store a mutable list [gene_id, sox2_max, nanog_max] in the interval data
# so we can easily update it while iterating through WIG
trees = {}
results_map = {}

for gene_id, intervals in gene_intervals.items():
    # Initialize the results container for this gene
    res = [gene_id, 0.0, 0.0]
    results_map[gene_id] = res
    for chrom, start, end in intervals:
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        # Ensure start < end
        if start < end:
            trees[chrom].addi(start, end, res)

def parse_wig(wig_file, signal_index):
    print(f"Parsing {wig_file}...")
    current_chrom = None
    span = 1
    tree = None
    
    with gzip.open(wig_file, 'rt') as f:
        for line in f:
            if line.startswith('track'):
                continue
            if line.startswith('variableStep'):
                parts = line.strip().split()
                current_chrom = parts[1].split('=')[1]
                tree = trees.get(current_chrom)
                if len(parts) > 2 and parts[2].startswith('span='):
                    span = int(parts[2].split('=')[1])
                continue
            
            if tree is None:
                continue
                
            parts = line.split()
            if len(parts) == 2:
                pos = int(parts[0])
                val = float(parts[1])
                # Find overlaps
                overlaps = tree.overlap(pos, pos + span)
                for overlap in overlaps:
                    res = overlap.data
                    if val > res[signal_index]:
                        res[signal_index] = val

# SOX2 WIG is 1, NANOG is 2
parse_wig('GSM1082341_10192012_C1AVNACXX_3.GTGGCC.wig.gz', 1)
parse_wig('GSM1082342_10192012_C1AVNACXX_2.GTTTCG.wig.gz', 2)

# Convert to DataFrame
df = pd.DataFrame(list(results_map.values()), columns=['ensembl_gene_id', 'sox2_pileup', 'nanog_pileup'])

# Determine target category
def get_category(gene_id):
    if gene_id in sox2_ensembl and gene_id in nanog_ensembl:
        return 'Both'
    elif gene_id in sox2_ensembl:
        return 'SOX2 only'
    elif gene_id in nanog_ensembl:
        return 'NANOG only'
    return 'None'

df['category'] = df['ensembl_gene_id'].apply(get_category)
df.to_csv('promoter_binding_pileups.csv', index=False)

# Analysis: How many targets have binding sites?
# We'll use a pileup threshold of 5 as a reasonable minimum for a MACS peak
threshold = 5.0
df['has_sox2_site'] = df['sox2_pileup'] > threshold
df['has_nanog_site'] = df['nanog_pileup'] > threshold

print("\n--- Summary (Threshold > 5) ---")
for cat in ['SOX2 only', 'NANOG only', 'Both']:
    sub = df[df['category'] == cat]
    sox2_bound = sub['has_sox2_site'].sum()
    nanog_bound = sub['has_nanog_site'].sum()
    print(f"{cat} target genes ({len(sub)} total):")
    print(f"  With SOX2 site: {sox2_bound} ({(sox2_bound/len(sub)*100):.1f}%)")
    print(f"  With NANOG site: {nanog_bound} ({(nanog_bound/len(sub)*100):.1f}%)")

# Plotting Nature-style figure
plt.rcParams.update({
    'font.size': 7,
    'axes.labelsize': 7,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'axes.linewidth': 0.5,
    'font.sans-serif': 'Arial',
    'font.family': 'sans-serif',
})

fig, axes = plt.subplots(1, 2, figsize=(4.5, 2))

cats = ['SOX2 only', 'NANOG only', 'Both']
sox2_pcts = [df[df['category']==c]['has_sox2_site'].mean()*100 for c in cats]
nanog_pcts = [df[df['category']==c]['has_nanog_site'].mean()*100 for c in cats]

axes[0].bar(cats, sox2_pcts, color=['#1f77b4', '#ff7f0e', '#2ca02c'], width=0.6)
axes[0].set_title('Promoters with SOX2 Peak', fontsize=7)
axes[0].set_ylabel('% of Target Genes')

axes[1].bar(cats, nanog_pcts, color=['#1f77b4', '#ff7f0e', '#2ca02c'], width=0.6)
axes[1].set_title('Promoters with NANOG Peak', fontsize=7)

for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(0, 100)
    ax.tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('promoter_binding_summary_nature.pdf')
plt.savefig('promoter_binding_summary_nature.png', dpi=300)
print("Saved plots.")
