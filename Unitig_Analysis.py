#!/usr/bin/env python3
import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import argparse
import os

# === Parse terminal arguments ===
parser = argparse.ArgumentParser(description="Process ggCaller unitig data and plot statistics.")
parser.add_argument('--tsv_dir', required=True, help='Directory containing *_assemblies_colors.tsv files')
parser.add_argument('--fasta_dir', required=True, help='Directory containing *_assemblies_untigs.fa files')
parser.add_argument('--out_dir', required=True, help='Directory where plots will be saved')

args = parser.parse_args()

tsv_dir = args.tsv_dir
fasta_dir = args.fasta_dir
out_dir = args.out_dir

# Create output directory if it doesn't exist
os.makedirs(out_dir, exist_ok=True)

# === Initialize containers ===
summary_list = []
all_scatter_dfs = []

# === Loop through all *_assemblies_colors.tsv files ===
for filepath in glob.glob(os.path.join(tsv_dir, "*_assemblies_colors.tsv")):
    df = pd.read_csv(filepath, sep='\t')

    # Get ISL/sample name
    basename = os.path.basename(filepath)
    isl_id = basename.replace("_assemblies_colors.tsv", "")
    df['isl_label'] = isl_id

    # === Load corresponding FASTA ===
    fasta_path = os.path.join(fasta_dir, f"{isl_id}_assemblies_untigs.fa")
    if not os.path.exists(fasta_path):
        print(f"⚠️ No FASTA found for {isl_id}, skipping.")
        continue

    unitig_length_dict = {
        record.id: len(record.seq)
        for record in SeqIO.parse(fasta_path, "fasta")
    }

    # Assign unitig lengths
    df['unitig_length'] = df['query_name'].map(unitig_length_dict).fillna(0).astype(int)

    # Compute presence count
    non_genome_cols = ['query_name', 'isl_label', 'unitig_length', 'presence_count']
    genome_columns = [col for col in df.columns if col not in non_genome_cols]
    df['presence_count'] = df[genome_columns].sum(axis=1)

    # Summary stats
    n_genomes = len(genome_columns)
    core_count = (df['presence_count'] == n_genomes).sum()
    unique_count = (df['presence_count'] == 1).sum()
    distribution = df['presence_count'].value_counts().sort_index().to_dict()

    summary_list.append({
        'file': filepath,
        'total_unitigs': len(df),
        'core_unitigs': core_count,
        'unique_unitigs': unique_count,
        'distribution': distribution
    })

    # Collect unitig-level info
    scatter_df = df[['unitig_length', 'presence_count', 'isl_label']]
    all_scatter_dfs.append(scatter_df)

# === Summary Table ===
summary_df = pd.DataFrame([
    {
        'file': s['file'],
        'total_unitigs': s['total_unitigs'],
        'core_unitigs': s['core_unitigs'],
        'unique_unitigs': s['unique_unitigs']
    }
    for s in summary_list
])

print("Summary of Core vs. Unique Unitigs per File:")
print(summary_df.to_string(index=False))

for s in summary_list:
    print(f"\nFile: {s['file']}")
    print(f"  Total unitigs:    {s['total_unitigs']}")
    print(f"  Core unitigs:     {s['core_unitigs']}")
    print(f"  Unique unitigs:   {s['unique_unitigs']}")
    print("  Presence-count distribution (unitigs → #genomes):")
    for count, freq in s['distribution'].items():
        print(f"    {count}: {freq}")
    print("-" * 60)

# === Barplot ===
summary_df['file_label'] = summary_df['file'].apply(lambda x: os.path.basename(x).replace('_assemblies_colors.tsv', ''))

melted_df = summary_df.melt(
    id_vars='file_label',
    value_vars=['total_unitigs', 'core_unitigs', 'unique_unitigs'],
    var_name='Unitig_Type',
    value_name='Count'
)

plt.figure(figsize=(14, 6))
sns.barplot(
    data=melted_df,
    x='file_label',
    y='Count',
    hue='Unitig_Type'
)
plt.xticks(rotation=45, ha='right')
plt.title('Comparison of Total, Core, and Unique Unitigs per File')
plt.xlabel('Sample')
plt.ylabel('Unitig Count')
plt.tight_layout()
plt.legend(title='Type')
plt.savefig(os.path.join(out_dir, "unitig_summary_barplot.png"))
plt.close()

# === Unified Scatterplot ===
combined_scatter_df = pd.concat(all_scatter_dfs, ignore_index=True)

plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=combined_scatter_df,
    x="unitig_length",
    y="presence_count",
    hue="isl_label",
    s=15,
    alpha=0.6
)
plt.xlabel("Unitig Length (bp)")
plt.ylabel("Number of Isolates Present")
plt.title("Unitig Length vs. Isolate Presence (Colored by ISL)")
plt.legend(title='ISL', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "unitig_length_vs_presence_scatterplot.png"))
plt.close()

# === Per-ISL Scatterplot (Filtered) ===
filtered_df = combined_scatter_df[combined_scatter_df['unitig_length'] > 100]

g = sns.FacetGrid(
    filtered_df,
    col="isl_label",
    col_wrap=4,
    sharex=False,
    sharey=False,
    height=4
)
g.map_dataframe(
    sns.scatterplot,
    x="unitig_length",
    y="presence_count",
    s=15,
    alpha=0.6
)
g.set_axis_labels("Unitig Length (bp)", "# of Isolates Present")
g.set_titles(col_template="{col_name}")
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle("Unitig Length > 100 bp vs. Isolate Presence (Per ISL)")
plt.tight_layout()
g.savefig(os.path.join(out_dir, "per_isl_scatterplots_filtered.png"))
plt.close()
