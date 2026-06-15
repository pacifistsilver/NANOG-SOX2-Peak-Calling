import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Nature guidelines (matching the residence time landscape plot)
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
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})

print("Loading data...")
df_bind = pd.read_csv('h3k4me3_promoter_site_counts_chen2008.csv')
df_burst = pd.read_csv('log.burst.freq vs log.burst.size.csv')

print("Merging...")
# We merge on ensembl_gene_id. df_burst may have multiple transcripts (tid) per gene, so we can group by gene or just use all records
# Grouping by gene and taking mean burst parameters for simplicity
df_burst_gene = df_burst.groupby('ensembl_gene_id').agg({
    'log.burst.freq': 'mean',
    'log.burst.size': 'mean'
}).reset_index()

df = pd.merge(df_bind, df_burst_gene, on='ensembl_gene_id')
print(f"Total merged genes: {len(df)}")

# Cap sites at "3+" to avoid thin tails
df['sox2_sites_cat'] = df['sox2_sites'].apply(lambda x: str(x) if x < 3 else '3+')
df['nanog_sites_cat'] = df['nanog_sites'].apply(lambda x: str(x) if x < 3 else '3+')

# Sort categories
cats = ['0', '1', '2', '3+']
df['sox2_sites_cat'] = pd.Categorical(df['sox2_sites_cat'], categories=cats, ordered=True)
df['nanog_sites_cat'] = pd.Categorical(df['nanog_sites_cat'], categories=cats, ordered=True)

# Remove NA burst data
df = df.dropna(subset=['log.burst.freq', 'log.burst.size'])

print("Plotting Burst Frequency...")
fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

sns.boxplot(data=df, x='sox2_sites_cat', y='log.burst.freq', ax=axes[0], color='#1f77b4', fliersize=0.5, linewidth=0.5)
axes[0].set_title('SOX2 Binding')
axes[0].set_xlabel('Number of SOX2 peaks\nat promoter')
axes[0].set_ylabel('Burst Frequency\n(log10 bursts / hour)')

sns.boxplot(data=df, x='nanog_sites_cat', y='log.burst.freq', ax=axes[1], color='#ff7f0e', fliersize=0.5, linewidth=0.5)
axes[1].set_title('NANOG Binding')
axes[1].set_xlabel('Number of NANOG peaks\nat promoter')

plt.tight_layout()
plt.savefig('burst_freq_vs_binding_nature.pdf')
plt.savefig('burst_freq_vs_binding_nature.png', dpi=300)

print("Plotting Burst Size...")
fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

sns.boxplot(data=df, x='sox2_sites_cat', y='log.burst.size', ax=axes[0], color='#1f77b4', fliersize=0.5, linewidth=0.5)
axes[0].set_title('SOX2 Binding')
axes[0].set_xlabel('Number of SOX2 peaks\nat promoter')
axes[0].set_ylabel('Burst Size\n(log10 mRNAs / burst)')

sns.boxplot(data=df, x='nanog_sites_cat', y='log.burst.size', ax=axes[1], color='#ff7f0e', fliersize=0.5, linewidth=0.5)
axes[1].set_title('NANOG Binding')
axes[1].set_xlabel('Number of NANOG peaks\nat promoter')

plt.tight_layout()
plt.savefig('burst_size_vs_binding_nature.pdf')
plt.savefig('burst_size_vs_binding_nature.png', dpi=300)

print("Done! Saved burst_freq_vs_binding_nature.pdf and burst_size_vs_binding_nature.pdf")
