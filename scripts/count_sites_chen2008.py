import os
import json
import gzip
import pandas as pd
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns
from liftover import ChainFile

print("1. Loading SOX2 targets from cache...")
cache_file = '../data/raw/mapping_cache.json'
with open(cache_file, 'r') as f:
    data = json.load(f)
    # The user wants "all promoters which are known to be activated by sox2"
    # We use all SOX2 targets from the ESCAPE dataset
    sox2_targets = set(data['sox2'])

print(f"Total SOX2 targets to process: {len(sox2_targets)}")

print("2. Parsing mm9 ensGene TSS coordinates...")
gene_tss_windows = {} # gene_id -> list of (chrom, start, end)
with gzip.open("../data/references/ensGene_mm9.txt.gz", 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[2]
        strand = parts[3]
        txStart = int(parts[4])
        txEnd = int(parts[5])
        gene_id = parts[12]
        
        if gene_id in sox2_targets:
            tss = txStart if strand == '+' else txEnd
            start = max(0, tss - 2000)
            end = tss + 2000
            
            if gene_id not in gene_tss_windows:
                gene_tss_windows[gene_id] = []
            gene_tss_windows[gene_id].append((chrom, start, end))

print("3. Loading H3K4me3 peaks to define active promoters...")
h3k4me3_tree = {}
with gzip.open('../data/references/h3k4me3_peaks_mm9.txt.gz', 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[1]
        start = int(parts[2])
        end = int(parts[3])
        if chrom not in h3k4me3_tree:
            h3k4me3_tree[chrom] = IntervalTree()
        # Enforce start < end
        if start < end:
            h3k4me3_tree[chrom].addi(start, end)

# Define active promoters: Genes with H3K4me3 peaks overlapping the TSS +/- 2kb
# We will save the original broader window (e.g. TSS +/- 5kb) to count the TF peaks
active_promoters_windows = {} # gene_id -> list of (chrom, start, end) for TSS +/- 5kb
for gene_id, intervals in gene_tss_windows.items():
    is_active = False
    for chrom, start, end in intervals:
        if chrom in h3k4me3_tree:
            overlaps = h3k4me3_tree[chrom].overlap(start, end) # TSS +/- 2kb
            if overlaps:
                is_active = True
                break
                
    if is_active:
        # Save a broader window around the TSS to count TF sites (TSS +/- 5kb)
        # We recalculate from the original coordinates
        broader_windows = []
        for chrom, start, end in intervals:
            center = (start + end) // 2
            broader_windows.append((chrom, max(0, center - 5000), center + 5000))
        active_promoters_windows[gene_id] = broader_windows

genes_with_active_promoter = len(active_promoters_windows)
print(f"Found active H3K4me3 promoters for {genes_with_active_promoter} out of {len(sox2_targets)} SOX2 target genes")

print("4. Extracting discrete SOX2 and NANOG peaks from narrowPeak files...")
def extract_peaks_from_bed(bed_file):
    print(f"   Extracting peaks from {bed_file}...")
    peaks_tree = {} # chrom -> IntervalTree
    
    # Initialize liftover mm8 -> mm9
    converter = ChainFile('../data/references/mm8ToMm9.over.chain.gz')
    
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom = parts[0]
                # Ensure it starts with 'chr' to match the H3K4me3 and EnsGene references
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                
                # We need to lift over mm8 coordinates to mm9
                try:
                    start_mm8 = int(parts[1])
                    end_mm8 = int(parts[2])
                    
                    # liftover requires 0-based coordinates, which BED is.
                    start_res = converter[chrom][start_mm8]
                    end_res = converter[chrom][end_mm8]
                    
                    # Only keep peaks that successfully map to the same chromosome
                    if start_res and end_res:
                        new_chrom = start_res[0][0]
                        if new_chrom == end_res[0][0]:
                            start = min(start_res[0][1], end_res[0][1])
                            end = max(start_res[0][1], end_res[0][1])
                            
                            if new_chrom not in peaks_tree:
                                peaks_tree[new_chrom] = IntervalTree()
                            if start < end:
                                peaks_tree[new_chrom].addi(start, end)
                except Exception as e:
                    pass
                    
    for chrom, tree in peaks_tree.items():
        tree.merge_overlaps()
        
    return peaks_tree

sox2_peaks = extract_peaks_from_bed('../data/raw/Sox2_Chen2008_peaks.narrowPeak')
nanog_peaks = extract_peaks_from_bed('../data/raw/Nanog_Chen2008_peaks.narrowPeak')

print("5. Counting SOX2 and NANOG sites in active promoters (TSS +/- 5kb)...")
results = []

for gene_id, windows in active_promoters_windows.items():
    sox2_count = 0
    nanog_count = 0
    
    for chrom, start, end in windows:
        # Overlap with SOX2 peaks
        if chrom in sox2_peaks:
            overlaps = sox2_peaks[chrom].overlap(start, end)
            sox2_count += len(overlaps)
            
        # Overlap with NANOG peaks
        if chrom in nanog_peaks:
            overlaps = nanog_peaks[chrom].overlap(start, end)
            nanog_count += len(overlaps)
            
    results.append({
        'ensembl_gene_id': gene_id,
        'sox2_sites': sox2_count,
        'nanog_sites': nanog_count
    })

df = pd.DataFrame(results)
df.to_csv('../data/processed/h3k4me3_promoter_site_counts_chen2008.csv', index=False)

print("\n--- Summary ---")
print(f"Total SOX2 targets with H3K4me3 promoter: {len(df)}")
print(f"Average SOX2 sites per promoter: {df['sox2_sites'].mean():.2f}")
print(f"Average NANOG sites per promoter: {df['nanog_sites'].mean():.2f}")
print(f"Max SOX2 sites at a single promoter: {df['sox2_sites'].max()}")
print(f"Max NANOG sites at a single promoter: {df['nanog_sites'].max()}")

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

fig, ax = plt.subplots(figsize=(3, 3))
sns.violinplot(data=df[['sox2_sites', 'nanog_sites']], palette=['#1f77b4', '#ff7f0e'], inner='quartile', ax=ax)
ax.set_xticklabels(['SOX2 Sites', 'NANOG Sites'])
ax.set_ylabel('Number of discrete binding sites\nper active promoter')
ax.set_title('TF Binding Sites at SOX2 Target Genes', fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('../plots/discrete_sites_count_nature.pdf')
plt.savefig('../plots/discrete_sites_count_nature.png', dpi=300)
print("Saved plots to discrete_sites_count_nature.pdf/png")
