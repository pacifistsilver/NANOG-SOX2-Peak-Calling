import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import time
import os
import json

# --- Nature Style Formatting ---
# Nature guidelines suggest:
# Single column width: 89 mm (~3.5 inches), Double column width: 183 mm (~7.2 inches)
# Font sizes: 5-7 pt. Sans-serif fonts like Arial or Helvetica.
# Removing top and right spines.
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 18,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})

def map_symbols_to_ensembl(symbols):
    ensembl_ids = set()
    symbols_list = list(symbols)
    chunk_size = 1000
    
    for i in range(0, len(symbols_list), chunk_size):
        chunk = symbols_list[i:i + chunk_size]
        q = ",".join(chunk)
        retries = 3
        while retries > 0:
            try:
                res = requests.post("https://mygene.info/v3/query", 
                                    data={'q': q, 'scopes': 'symbol,alias', 'fields': 'ensembl.gene', 'species': 'mouse'})
                res.raise_for_status()
                data = res.json()
                for item in data:
                    if 'ensembl' in item:
                        ens = item['ensembl']
                        if isinstance(ens, list):
                            for e in ens:
                                if 'gene' in e:
                                    ensembl_ids.add(e['gene'])
                        elif isinstance(ens, dict):
                            if 'gene' in ens:
                                ensembl_ids.add(ens['gene'])
                break
            except Exception as e:
                print(f"API request failed: {e}. Retrying...")
                retries -= 1
                time.sleep(2)
        time.sleep(0.1)
    return ensembl_ids

def get_mappings():
    cache_file = "mapping_cache.json"
    if os.path.exists(cache_file):
        print("Loading mapped Ensembl IDs from cache...")
        with open(cache_file, 'r') as f:
            data = json.load(f)
            return set(data['nanog']), set(data['sox2'])
            
    print("Cache not found. Parsing ESCAPE_tf_targets.txt...")
    nanog_symbols = set()
    sox2_symbols = set()
    nanog_ensembl = set()
    sox2_ensembl = set()
    
    with open("ESCAPE_tf_targets.txt", "r") as f:
        # Skip header if exists
        first_line = f.readline()
        if not first_line.startswith("sourceName"):
            f.seek(0)
        
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) > 3:
                source = parts[0].strip().upper()
                target = parts[2].strip()
                if source == "NANOG":
                    if target.startswith("ENSMUSG"): nanog_ensembl.add(target)
                    else: nanog_symbols.add(target)
                elif source == "SOX2":
                    if target.startswith("ENSMUSG"): sox2_ensembl.add(target)
                    else: sox2_symbols.add(target)
                    
    print("Mapping NANOG symbols...")
    nanog_ensembl.update(map_symbols_to_ensembl(nanog_symbols))
    print("Mapping SOX2 symbols...")
    sox2_ensembl.update(map_symbols_to_ensembl(sox2_symbols))
    
    with open(cache_file, 'w') as f:
        json.dump({'nanog': list(nanog_ensembl), 'sox2': list(sox2_ensembl)}, f)
        
    return nanog_ensembl, sox2_ensembl

nanog_ensembl, sox2_ensembl = get_mappings()

def get_category(ens_id):
    in_nanog = ens_id in nanog_ensembl
    in_sox2 = ens_id in sox2_ensembl
    if in_nanog and in_sox2: return "Both"
    elif in_nanog: return "NANOG only"
    elif in_sox2: return "SOX2 only"
    else: return "Background"

def create_burst_plot(df, x_col, y_col, c_col, filename, x_label, y_label):
    categories = ['NANOG only', 'SOX2 only', 'Both']
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharex=True, sharey=True, constrained_layout=True)
    
    target_mask = df['Target_Type'] != 'Background'
    vmin = df.loc[target_mask, c_col].min()
    vmax = df.loc[target_mask, c_col].max()
    
    for i, category in enumerate(categories):
        ax = axes[i]
        
        # Plot ALL points as background
        ax.scatter(df[x_col], df[y_col], color='#d3d3d3', s=0.5, alpha=0.2, rasterized=True, edgecolors='none')
        
        # Plot the specific category on top colored by intrinsic noise
        target_df = df[df['Target_Type'] == category]
        sc = ax.scatter(target_df[x_col], target_df[y_col], c=target_df[c_col], cmap='viridis', 
                   vmin=vmin, vmax=vmax, s=1.5, alpha=0.8, rasterized=True, edgecolors='none')
                   
        # Add density contours
        if len(target_df) > 10:
            try:
                sns.kdeplot(data=target_df, x=x_col, y=y_col, ax=ax, color='black', alpha=0.6, linewidths=0.5, levels=5)
            except Exception:
                pass
                
        ax.set_title(category)
        ax.set_xlabel(x_label)
        if i == 0:
            ax.set_ylabel(y_label)
            
    # Add a colorbar for intrinsic noise
    cbar = fig.colorbar(sc, ax=axes, shrink=0.8, aspect=20, pad=0.02)
    cbar.set_label('Intrinsic Noise\n(log.norm.int.noise)', rotation=270, labelpad=15)
            
    plt.savefig(filename.replace(".pdf", ".png"), dpi=300, bbox_inches='tight')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved {filename} and corresponding .png")


def create_exp_plot(df, x_col, y_col, filename, x_label, y_label):
    categories = ['NANOG only', 'SOX2 only', 'Both']
    colors = {'NANOG only': '#1f77b4', 'SOX2 only': '#d62728', 'Both': '#9467bd'}
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharex=True, sharey=True, constrained_layout=True)
    
    for i, category in enumerate(categories):
        ax = axes[i]
        
        # Plot ALL points as background
        ax.scatter(df[x_col], df[y_col], color='#d3d3d3', s=0.5, alpha=0.2, rasterized=True, edgecolors='none')
        
        # Plot the specific category on top
        target_df = df[df['Target_Type'] == category]
        ax.scatter(target_df[x_col], target_df[y_col], color=colors[category], s=1.5, alpha=0.6, 
                   rasterized=True, edgecolors='none', label=category)
                   
        # Add density contours 
        if len(target_df) > 10:
            try:
                # Supply log_scale if seaborn supports it, else fallback
                sns.kdeplot(data=target_df, x=x_col, y=y_col, ax=ax, color='black', alpha=0.6, linewidths=0.5, levels=5, log_scale=(True, False))
            except Exception:
                try:
                    sns.kdeplot(data=target_df, x=x_col, y=y_col, ax=ax, color='black', alpha=0.6, linewidths=0.5, levels=5)
                except Exception:
                    pass
                
        ax.set_title(category)
        ax.set_xlabel(x_label)
        ax.set_xscale('log')
        ax.set_xlim(0.01, 1000)
        
        if i == 0:
            ax.set_ylabel(y_label)
            
    plt.savefig(filename.replace(".pdf", ".png"), dpi=300, bbox_inches='tight')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved {filename} and corresponding .png")

# Process Burst Data
file_burst = "log.burst.freq vs log.burst.size.csv"
if os.path.exists(file_burst):
    print(f"Processing {file_burst}...")
    df_burst = pd.read_csv(file_burst)
    
    # Filter out the two extreme outliers
    df_burst = df_burst[(df_burst['log.burst.size'] > -10) & (df_burst['log.burst.freq'] < 10)]
    
    df_burst['Target_Type'] = df_burst['ensembl_gene_id'].apply(get_category)
    create_burst_plot(df_burst, 'log.burst.size', 'log.burst.freq', 'log.norm.int.noise',
                      "plot_burst_nature.pdf", "Log Burst Size", "Log Burst Freq")

# Process Exp Data
file_exp = "exp.norm vs. log.norm.int.noise.csv"
if os.path.exists(file_exp):
    print(f"Processing {file_exp}...")
    df_exp = pd.read_csv(file_exp)
    df_exp['Target_Type'] = df_exp['ensembl_gene_id'].apply(get_category)
    create_exp_plot(df_exp, 'exp.norm', 'log.norm.int.noise', 
                    "plot_exp_nature.pdf", "Exp Norm", "Log Norm Int Noise")

print("Done!")
