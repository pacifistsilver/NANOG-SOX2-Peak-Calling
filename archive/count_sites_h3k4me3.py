import os
import json
import gzip
import pandas as pd
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns

print("1. Loading SOX2 targets from cache...")
cache_file = 'mapping_cache.json'
with open(cache_file, 'r') as f:
    data = json.load(f)
    # The user wants "all promoters which are known to be activated by sox2"
    # We use all SOX2 targets from the ESCAPE dataset
    sox2_targets = set(data['sox2'])

print(f"Total SOX2 targets to process: {len(sox2_targets)}")

print("2. Parsing mm9 ensGene TSS coordinates...")
gene_tss_windows = {} # gene_id -> list of (chrom, start, end)
with gzip.open("ensGene_mm9.txt.gz", 'rt') as f:
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
with gzip.open('h3k4me3_peaks_mm9.txt.gz', 'rt') as f:
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

# Define active promoters: H3K4me3 peaks overlapping the TSS +/- 2kb
active_promoters = {} # gene_id -> IntervalTree of merged H3K4me3 peaks
for gene_id, intervals in gene_tss_windows.items():
    tree = IntervalTree()
    for chrom, start, end in intervals:
        if chrom in h3k4me3_tree:
            overlaps = h3k4me3_tree[chrom].overlap(start, end)
            for ov in overlaps:
                tree.addi(ov.begin, ov.end, chrom)
    
    if tree:
        tree.merge_overlaps()
        active_promoters[gene_id] = tree

genes_with_active_promoter = len(active_promoters)
print(f"Found active H3K4me3 promoters for {genes_with_active_promoter} out of {len(sox2_targets)} SOX2 target genes")

print("4. Extracting discrete SOX2 and NANOG peaks from WIG pileup...")
def extract_peaks_from_wig(wig_file, threshold=5.0, max_gap=100):
    print(f"   Extracting peaks from {wig_file}...")
    peaks_tree = {} # chrom -> IntervalTree
    
    current_chrom = None
    span = 1
    in_peak = False
    peak_start = 0
    peak_end = 0
    
    with gzip.open(wig_file, 'rt') as f:
        for line in f:
            if line.startswith('track'): continue
            if line.startswith('variableStep'):
                if in_peak:
                    if current_chrom not in peaks_tree: peaks_tree[current_chrom] = IntervalTree()
                    if peak_start < peak_end:
                        peaks_tree[current_chrom].addi(peak_start, peak_end)
                    in_peak = False
                    
                parts = line.strip().split()
                current_chrom = parts[1].split('=')[1]
                if len(parts) > 2 and parts[2].startswith('span='):
                    span = int(parts[2].split('=')[1])
                continue
            
            parts = line.split()
            if len(parts) == 2:
                pos = int(parts[0])
                val = float(parts[1])
                
                if val > threshold:
                    if not in_peak:
                        in_peak = True
                        peak_start = pos
                        peak_end = pos + span
                    else:
                        if pos <= peak_end + max_gap:
                            peak_end = max(peak_end, pos + span)
                        else:
                            if current_chrom not in peaks_tree: peaks_tree[current_chrom] = IntervalTree()
                            if peak_start < peak_end:
                                peaks_tree[current_chrom].addi(peak_start, peak_end)
                            peak_start = pos
                            peak_end = pos + span
                else:
                    if in_peak:
                        if current_chrom not in peaks_tree: peaks_tree[current_chrom] = IntervalTree()
                        if peak_start < peak_end:
                            peaks_tree[current_chrom].addi(peak_start, peak_end)
                        in_peak = False
                        
    if in_peak:
        if current_chrom not in peaks_tree: peaks_tree[current_chrom] = IntervalTree()
        if peak_start < peak_end:
            peaks_tree[current_chrom].addi(peak_start, peak_end)
            
    # Merge nearby peaks
    for chrom, tree in peaks_tree.items():
        tree.merge_overlaps()
        
    return peaks_tree

sox2_peaks = extract_peaks_from_wig('GSM1082341_10192012_C1AVNACXX_3.GTGGCC.wig.gz')
nanog_peaks = extract_peaks_from_wig('GSM1082342_10192012_C1AVNACXX_2.GTTTCG.wig.gz')

print("5. Counting SOX2 and NANOG sites in active promoters...")
results = []

for gene_id, prom_tree in active_promoters.items():
    sox2_count = 0
    nanog_count = 0
    
    for prom in prom_tree:
        chrom = prom.data
        start = prom.begin
        end = prom.end
        
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
df.to_csv('h3k4me3_promoter_site_counts.csv', index=False)

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
plt.savefig('discrete_sites_count_nature.pdf')
plt.savefig('discrete_sites_count_nature.png', dpi=300)
print("Saved plots to discrete_sites_count_nature.pdf/png")
