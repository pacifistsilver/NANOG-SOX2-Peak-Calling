import os
import gzip
import numpy as np
import pandas as pd
import urllib.request
import matplotlib.pyplot as plt
from collections import defaultdict

# ----------------- CONFIGURATION -----------------
# Relative paths since we run from the project root
ENSGENE_FILE = 'data/references/ensGene_mm9.txt.gz'
NANOG_BED = 'data/processed/NANOG_peaks_mm9.bed'
SOX2_BED = 'data/processed/SOX2_peaks_mm9.bed'
OUTPUT_CSV = 'data/processed/tf_gene_association_scores.csv'

# Bin boundaries as defined in the prompt
BINS = [-np.inf, -100000, -50000, -20000, -10000, -5000, -2000, -1000, 0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, np.inf]

def load_mm9_chrom_sizes():
    """Fetch mm9 chromosome sizes directly from UCSC."""
    print("Fetching mm9 chrom sizes from UCSC...")
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes"
    chrom_sizes = {}
    try:
        response = urllib.request.urlopen(url)
        for line in response:
            parts = line.decode('utf-8').strip().split('\t')
            if len(parts) == 2 and parts[0].startswith('chr') and '_' not in parts[0]:
                chrom_sizes[parts[0]] = int(parts[1])
    except Exception as e:
        print(f"Failed to fetch chrom sizes: {e}")
        # Fallback to standard mm9 sizes if offline
        chrom_sizes = {'chr1': 197195432, 'chr2': 181748087, 'chr3': 159599783, 'chr4': 155630120, 'chr5': 152537259, 'chr6': 149517037, 'chr7': 152524553, 'chr8': 131336506, 'chr9': 124076172, 'chr10': 129993255, 'chr11': 121843856, 'chr12': 121257530, 'chr13': 120284312, 'chr14': 125194864, 'chr15': 103494974, 'chr16': 98319150, 'chr17': 95272651, 'chr18': 90702639, 'chr19': 61342430, 'chrX': 166650296, 'chrY': 15902555, 'chrM': 16299}
    return chrom_sizes

def load_genes(filepath):
    """Load non-redundant TSS locations from ensGene."""
    print(f"Loading genes from {filepath}...")
    genes = {}
    
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 13:
                continue
            
            chrom = parts[2]
            strand = parts[3]
            txStart = int(parts[4])
            txEnd = int(parts[5])
            gene_id = parts[12]
            
            # Simple non-redundancy: take longest transcript
            length = txEnd - txStart
            if gene_id not in genes or genes[gene_id]['length'] < length:
                genes[gene_id] = {
                    'chrom': chrom,
                    'strand': strand,
                    'tss': txStart if strand == '+' else txEnd,
                    'length': length
                }
                
    print(f"Found {len(genes)} non-redundant genes.")
    return genes

def load_peaks(filepath):
    """Load peaks from a BED file and extract center coordinates."""
    print(f"Loading peaks from {filepath}...")
    peaks_by_chrom = defaultdict(list)
    total_peaks = 0
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                center = (start + end) // 2
                peaks_by_chrom[chrom].append(center)
                total_peaks += 1
                
    # Sort peaks for fast binary search
    for chrom in peaks_by_chrom:
        peaks_by_chrom[chrom].sort()
        peaks_by_chrom[chrom] = np.array(peaks_by_chrom[chrom])
        
    print(f"Loaded {total_peaks} peaks.")
    return peaks_by_chrom, total_peaks

def calculate_nearest_distances(genes, peaks_by_chrom):
    """Calculate the distance from TSS to the nearest peak for all genes."""
    distances = []
    
    for gene_id, g in genes.items():
        chrom = g['chrom']
        if chrom not in peaks_by_chrom or len(peaks_by_chrom[chrom]) == 0:
            distances.append(np.nan)
            continue
            
        tss = g['tss']
        peaks = peaks_by_chrom[chrom]
        
        # Find nearest peak
        idx = np.searchsorted(peaks, tss)
        
        if idx == 0:
            nearest_peak = peaks[0]
        elif idx == len(peaks):
            nearest_peak = peaks[-1]
        else:
            left = peaks[idx-1]
            right = peaks[idx]
            nearest_peak = left if (tss - left) < (right - tss) else right
            
        # Calculate signed distance relative to TSS
        if g['strand'] == '+':
            dist = nearest_peak - tss # Positive means peak is downstream
        else:
            dist = tss - nearest_peak # Positive means peak is downstream
            
        distances.append(dist)
        
    return np.array(distances)

def calculate_histogram(distances):
    """Calculate normalized histogram using defined bins."""
    valid_dists = distances[~np.isnan(distances)]
    counts, _ = np.histogram(valid_dists, bins=BINS)
    hist = counts / len(valid_dists) # Normalize so sum = 1
    return hist

def generate_random_peaks(chrom_sizes, total_peaks):
    """Generate uniformly distributed random peaks across the genome."""
    chrom_names = list(chrom_sizes.keys())
    chrom_lens = np.array(list(chrom_sizes.values()))
    chrom_probs = chrom_lens / chrom_lens.sum()
    
    random_peaks_by_chrom = defaultdict(list)
    # Pick random chromosomes
    chosen_chroms = np.random.choice(chrom_names, size=total_peaks, p=chrom_probs)
    
    for chrom in chosen_chroms:
        pos = np.random.randint(0, chrom_sizes[chrom])
        random_peaks_by_chrom[chrom].append(pos)
        
    for chrom in random_peaks_by_chrom:
        random_peaks_by_chrom[chrom].sort()
        random_peaks_by_chrom[chrom] = np.array(random_peaks_by_chrom[chrom])
        
    return random_peaks_by_chrom

def compute_association_scores(tf_name, genes, peaks_by_chrom, total_peaks, chrom_sizes):
    print(f"\nProcessing {tf_name} association scores...")
    
    # 1. Real histogram
    real_distances = calculate_nearest_distances(genes, peaks_by_chrom)
    hist_real = calculate_histogram(real_distances)
    
    # 2. Random histogram (averaged over 10 iterations for stability)
    hist_rands = []
    print("Running 10 randomizations for null model...")
    for _ in range(10):
        rand_peaks = generate_random_peaks(chrom_sizes, total_peaks)
        rand_dists = calculate_nearest_distances(genes, rand_peaks)
        hist_rands.append(calculate_histogram(rand_dists))
        
    hist_rand = np.mean(hist_rands, axis=0)
    
    # 3. Calculate scores for each gene
    scores = []
    bins = []
    
    # Pre-calculate score lookup table for each bin
    # Score = (H_real - H_rand)/H_real if H_real >= H_rand else 0
    bin_scores = np.zeros(len(hist_real))
    for k in range(len(hist_real)):
        if hist_real[k] >= hist_rand[k] and hist_real[k] > 0:
            bin_scores[k] = (hist_real[k] - hist_rand[k]) / hist_real[k]
        else:
            bin_scores[k] = 0.0
            
    # Assign scores and count valid physical peaks for each gene
    gene_ids = list(genes.keys())
    peak_counts = []
    
    for i, (gene_id, g) in enumerate(genes.items()):
        dist = real_distances[i]
        if np.isnan(dist):
            scores.append(0.0)
            bins.append(-1)
            peak_counts.append(0)
        else:
            k = np.digitize(dist, BINS) - 1
            scores.append(bin_scores[k])
            bins.append(k)
            
            # Count physical peaks using score threshold (Score > 0)
            chrom = g['chrom']
            tss = g['tss']
            peaks = peaks_by_chrom.get(chrom, np.array([]))
            
            if len(peaks) > 0:
                if g['strand'] == '+':
                    all_dists = peaks - tss
                else:
                    all_dists = tss - peaks
                    
                peak_bins = np.digitize(all_dists, BINS) - 1
                # Clip peak_bins to be within 0 to 15 (though digitize with np.inf handles this, 
                # but exactly matching np.inf could yield 16)
                peak_bins = np.clip(peak_bins, 0, 15)
                valid_count = np.sum(bin_scores[peak_bins] > 0)
                peak_counts.append(valid_count)
            else:
                peak_counts.append(0)
            
    # Setup Nature-style formatting
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['axes.titlesize'] = 16

    plt.figure(figsize=(10, 5))
    
    # Create clean x-axis labels to match the image precisely
    x_labels = []
    for i in range(16):
        if BINS[i] == -np.inf:
            x_labels.append(f"< -{int(BINS[i+1]//-1000)}k")
        elif BINS[i+1] == np.inf:
            x_labels.append(f"> {int(BINS[i]//1000)}k")
        else:
            x_labels.append(f"{int(BINS[i]//1000)}k to {int(BINS[i+1]//1000)}k")
            
    x_pos = np.arange(len(x_labels))
    
    # Plot as a bar chart matching the reference image
    plt.bar(x_pos, bin_scores, color='#9CA2E8', edgecolor='black', width=0.6)
    
    plt.title(f'{tf_name}', fontweight='bold', fontsize=18)
    plt.xlabel('location of nearest TF binding site in relative to TSS', fontweight='bold')
    plt.ylabel('TF-gene association score', fontweight='bold')
    
    plt.xticks(x_pos, x_labels, rotation=45, ha='right', fontweight='bold')
    
    # Set y limit slightly above max score or to 1.2 minimum to match the scale
    max_score = max(np.max(bin_scores) * 1.1, 1.2)
    plt.ylim(0, max_score)
    
    # Add light grid lines matching the image background style
    plt.grid(axis='y', linestyle='-', alpha=0.7)
    
    # Make background slightly gray like the Excel plot if desired, or keep white for Nature. 
    # Kept white for cleaner Nature look.
    
    # Remove top and right spines
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f'data/processed/{tf_name}_association_score_histogram.pdf', dpi=300)
    plt.savefig(f'data/processed/{tf_name}_association_score_histogram.png', dpi=300)
    plt.close()
    
    return gene_ids, real_distances, bins, scores, peak_counts

def main():
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    
    chrom_sizes = load_mm9_chrom_sizes()
    genes = load_genes(ENSGENE_FILE)
    
    nanog_peaks, num_nanog = load_peaks(NANOG_BED)
    sox2_peaks, num_sox2 = load_peaks(SOX2_BED)
    
    # Process NANOG
    gene_ids, nanog_dists, nanog_bins, nanog_scores, nanog_counts = compute_association_scores(
        "NANOG", genes, nanog_peaks, num_nanog, chrom_sizes)
        
    # Process SOX2
    _, sox2_dists, sox2_bins, sox2_scores, sox2_counts = compute_association_scores(
        "SOX2", genes, sox2_peaks, num_sox2, chrom_sizes)
        
    # Combine results
    df = pd.DataFrame({
        'GeneSymbol': gene_ids,
        'NANOG_Distance': nanog_dists,
        'NANOG_Bin': nanog_bins,
        'NANOG_Score': nanog_scores,
        'NANOG_PeakCount': nanog_counts,
        'SOX2_Distance': sox2_dists,
        'SOX2_Bin': sox2_bins,
        'SOX2_Score': sox2_scores,
        'SOX2_PeakCount': sox2_counts
    })
    
    try:
        df.to_csv(OUTPUT_CSV, index=False)
        print(f"\nSuccessfully wrote association scores and counts to {OUTPUT_CSV}")
    except PermissionError:
        print(f"\n[WARNING] Could not save to {OUTPUT_CSV}. Please close the file if it is open in another program.")
        fallback_csv = OUTPUT_CSV.replace('.csv', '_with_counts.csv')
        df.to_csv(fallback_csv, index=False)
        print(f"Saved results to fallback file: {fallback_csv}")

if __name__ == "__main__":
    main()
