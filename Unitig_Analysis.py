#!/usr/bin/env python3
import os
import re
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# =========================
# Args
# =========================
parser = argparse.ArgumentParser(
    description="Process ggCaller/Bifrost unitig data and plot statistics."
)
parser.add_argument('--tsv_dir', required=True,
                    help='Directory containing *_assemblies_colors.tsv files')
parser.add_argument('--fasta_dir', required=True,
                    help='Directory containing *_assemblies_untigs.fa files')
parser.add_argument('--out_dir', required=True,
                    help='Directory where plots will be saved')
parser.add_argument('--assembly_dir', required=False, default=None,
                    help='Base directory with per-sample assembly FASTAs (for assembly size comparison)')
args = parser.parse_args()

tsv_dir = args.tsv_dir
fasta_dir = args.fasta_dir
out_dir = args.out_dir
assembly_dir = args.assembly_dir

os.makedirs(out_dir, exist_ok=True)

# =========================
# Helpers
# =========================
def _sum_fasta_bp(fp: str) -> int:
    if not fp or not os.path.exists(fp):
        return 0
    total = 0
    opener = open
    if fp.endswith('.gz'):
        import gzip
        opener = gzip.open
    with opener(fp, 'rt') as handle:
        for rec in SeqIO.parse(handle, 'fasta'):
            total += len(rec.seq)
    return total

def _normalize_sid(raw) -> str:
    """
    Normalize a TSV column/header to canonical 'S.######.#####' (or 'REFERENCE').
    Works for 'S.200213.01443', 'S_200213_01443', 'S200213.01443', 'S-200213-01443', etc.
    """
    if raw is None:
        return ''
    s = str(raw).strip()
    if s.lower().startswith('reference'):
        return 'REFERENCE'
    base = os.path.basename(s)
    base = re.sub(r'\.(fa|fna|fasta)(\.gz)?$', '', base, flags=re.IGNORECASE)
    m = re.search(r'S[._-]?\d+[._-]\d+', base)
    if m:
        tok = m.group(0)
        tok = tok.replace('_', '.').replace('-', '.')
        tok = re.sub(r'^S(?=\d)', 'S.', tok)
        return tok
    # last resort: normalize whole string’s separators
    base = base.replace('_', '.').replace('-', '.')
    base = re.sub(r'^S(?=\d)', 'S.', base)
    return base

# =========================
# Containers
# =========================
summary_list = []
all_scatter_dfs = []
all_genome_sizes = []    # proxy from unitigs (per-sample)
all_assembly_sizes = []  # per-sample assembly size
assembly_bp_by_sample = {}  # cache per (isl_id, sample_id)

# =========================
# Main loop over *_assemblies_colors.tsv
# =========================
for filepath in glob.glob(os.path.join(tsv_dir, "*_assemblies_colors.tsv")):
    df = pd.read_csv(filepath, sep='\t')

    isl_id = os.path.basename(filepath).replace("_assemblies_colors.tsv", "")
    df['isl_label'] = isl_id

    # Load corresponding unitig FASTA for lengths
    fasta_path = os.path.join(fasta_dir, f"{isl_id}_assemblies_untigs.fa")
    if not os.path.exists(fasta_path):
        print(f"⚠️ No unitig FASTA found for {isl_id}, skipping.")
        continue

    unitig_length_dict = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
    df['unitig_length'] = df['query_name'].map(unitig_length_dict).fillna(0).astype(int)

    # Identify sample columns (everything that isn’t a helper column)
    helper_cols = {'query_name','isl_label','unitig_length','presence_count'}
    genome_columns = [c for c in df.columns if c not in helper_cols]

    # Presence count per unitig
    presence_matrix = df[genome_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
    df['presence_count'] = (presence_matrix > 0.5).sum(axis=1)

    # Summary stats
    n_genomes = len(genome_columns)
    core_count = (df['presence_count'] == n_genomes).sum()
    unique_count = (df['presence_count'] == 1).sum()
    distribution = df['presence_count'].value_counts().sort_index().to_dict()

    summary_list.append({'file': filepath,
                         'total_unitigs': len(df),
                         'core_unitigs': core_count,
                         'unique_unitigs': unique_count,
                         'distribution': distribution})

    # For unitig scatter (length vs presence)
    all_scatter_dfs.append(df[['unitig_length','presence_count','isl_label']])

    # Per-sample sizes: proxy (sum unitigs present) + assembly if available
    for i, sample_col in enumerate(genome_columns, start=1):
        # Proxy genome size from unitigs present
        col_vals = pd.to_numeric(df[sample_col], errors='coerce').fillna(0)
        size_bp = df.loc[col_vals > 0.5, 'unitig_length'].sum()
        all_genome_sizes.append({'isl_label': isl_id,
                                 'sample_col': sample_col,
                                 'rank': i,
                                 'size_mbp': size_bp/1e6})

        # Actual assembly size (if assemblies exist and are named by normalized sample id)
        if assembly_dir:
            sample_id = _normalize_sid(sample_col)
            cache_key = (isl_id, sample_id)
            if cache_key not in assembly_bp_by_sample:
                asm_fp = os.path.join(assembly_dir, f"{isl_id}_assemblies", f"{sample_id}.fna")
                asm_bp = _sum_fasta_bp(asm_fp) if os.path.exists(asm_fp) else 0
                assembly_bp_by_sample[cache_key] = asm_bp
            asm_bp_val = assembly_bp_by_sample[cache_key]
            all_assembly_sizes.append({'isl_label': isl_id,
                                       'sample_col': sample_col,
                                       'rank': i,
                                       'assembly_mbp': asm_bp_val/1e6})

# =========================
# UNITIG ANALYSIS PLOTS
# =========================

# 1) Summary table (stdout) + barplot
summary_df = pd.DataFrame([
    {'file': s['file'],
     'total_unitigs': s['total_unitigs'],
     'core_unitigs': s['core_unitigs'],
     'unique_unitigs': s['unique_unitigs']}
    for s in summary_list
])
print("Summary of Core vs. Unique Unitigs per File:")
if not summary_df.empty:
    print(summary_df.to_string(index=False))
for s in summary_list:
    print(f"\nFile: {s['file']}")
    print(f"  Total unitigs:    {s['total_unitigs']}")
    print(f"  Core unitigs:     {s['core_unitigs']}")
    print(f"  Unique unitigs:   {s['unique_unitigs']}")
    print("  Presence-count distribution (unitigs → #genomes):")
    # Uncomment to see distribution details:
    # for count, freq in s['distribution'].items():
    #     print(f"    {count}: {freq}")
    # print("-" * 60)

if not summary_df.empty:
    summary_df['file_label'] = summary_df['file'].apply(
        lambda x: os.path.basename(x).replace('_assemblies_colors.tsv', '')
    )
    melted_df = summary_df.melt(
        id_vars='file_label',
        value_vars=['total_unitigs', 'core_unitigs', 'unique_unitigs'],
        var_name='Unitig_Type',
        value_name='Count'
    )
    plt.figure(figsize=(14, 6))
    sns.barplot(data=melted_df, x='file_label', y='Count', hue='Unitig_Type')
    plt.xticks(rotation=45, ha='right')
    plt.title('Comparison of Total, Core, and Unique Unitigs per File')
    plt.xlabel('Sample')
    plt.ylabel('Unitig Count')
    plt.tight_layout()
    plt.legend(title='Type')
    plt.savefig(os.path.join(out_dir, "unitig_summary_barplot.png"))
    plt.close()

# 2) Unitig length vs presence (combined + per-ISL)
if len(all_scatter_dfs):
    combined_scatter_df = pd.concat(all_scatter_dfs, ignore_index=True)

    # Combined
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=combined_scatter_df, x="unitig_length", y="presence_count",
        hue="isl_label", s=15, alpha=0.6
    )
    plt.xlabel("Unitig Length (bp)")
    plt.ylabel("Number of Isolates Present")
    plt.title("Unitig Length vs. Isolate Presence (Colored by ISL)")
    plt.legend(title='ISL', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "unitig_length_vs_presence_scatterplot.png"))
    plt.close()

    # Per-ISL small multiples (filtered length > 100)
    filtered_df = combined_scatter_df[combined_scatter_df['unitig_length'] > 100]
    g = sns.FacetGrid(filtered_df, col="isl_label", col_wrap=4, sharex=False, sharey=False, height=4)
    g.map_dataframe(sns.scatterplot, x="unitig_length", y="presence_count", s=15, alpha=0.6)
    g.set_axis_labels("Unitig Length (bp)", "# of Isolates Present")
    g.set_titles(col_template="{col_name}")
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Unitig Length > 100 bp vs. Isolate Presence (Per ISL)")
    g.savefig(os.path.join(out_dir, "per_isl_scatterplots_filtered.png"))
    plt.close()

# 3) Proxy vs Assembly size scatter (per sample), if assemblies provided
if len(all_genome_sizes) and len(all_assembly_sizes):
    sizes = pd.DataFrame(all_genome_sizes)
    asm_sizes = pd.DataFrame(all_assembly_sizes)
    merged = pd.merge(
        sizes[['isl_label','sample_col','rank','size_mbp']],
        asm_sizes[['isl_label','sample_col','rank','assembly_mbp']],
        on=['isl_label','sample_col','rank'], how='inner'
    )
    if not merged.empty:
        plt.figure(figsize=(6, 6))
        ax = sns.scatterplot(data=merged, x='size_mbp', y='assembly_mbp',
                             hue='isl_label', alpha=0.85)
        lo = min(merged['size_mbp'].min(), merged['assembly_mbp'].min())
        hi = max(merged['size_mbp'].max(), merged['assembly_mbp'].max())
        ax.plot([lo, hi], [lo, hi], linestyle='--', linewidth=1)
        ax.set_xlim(lo*0.98, hi*1.02)
        ax.set_ylim(lo*0.98, hi*1.02)
        plt.xlabel('Proxy size (sum of unitigs present, Mbp)')
        plt.ylabel('Per-sample assembly size (Mbp)')
        plt.title('Proxy vs Assembly Genome Size (per sample)')
        handles, labels = ax.get_legend_handles_labels()
        bylabel = dict(zip(labels, handles))
        plt.legend(bylabel.values(), bylabel.keys(),
                   title='ISL', loc='upper left', bbox_to_anchor=(1.02, 1),
                   frameon=False)
        plt.tight_layout(rect=(0, 0, 0.85, 1))
        plt.savefig(os.path.join(out_dir, 'proxy_vs_assembly_scatter.png'))
        plt.close()


####### new, need to test ##########

'''
#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# =========================
# Args
# =========================
parser = argparse.ArgumentParser(
    description="Process ggCaller/Bifrost unitig data (unitig-only) and plot statistics."
)
parser.add_argument('--tsv_dir', required=True,
                    help='Directory containing *_assemblies_colors.tsv files')
parser.add_argument('--fasta_dir', required=True,
                    help='Directory containing *_assemblies_untigs.fa files')
parser.add_argument('--out_dir', required=True,
                    help='Directory where plots will be saved')
args = parser.parse_args()

tsv_dir = args.tsv_dir
fasta_dir = args.fasta_dir
out_dir = args.out_dir
os.makedirs(out_dir, exist_ok=True)

# =========================
# Containers
# =========================
summary_list = []
all_scatter_dfs = []

# =========================
# Main loop over *_assemblies_colors.tsv
# =========================
for filepath in glob.glob(os.path.join(tsv_dir, "*_assemblies_colors.tsv")):
    df = pd.read_csv(filepath, sep='\t')

    isl_id = os.path.basename(filepath).replace("_assemblies_colors.tsv", "")
    df['isl_label'] = isl_id

    # Load corresponding unitig FASTA for lengths
    fasta_path = os.path.join(fasta_dir, f"{isl_id}_assemblies_untigs.fa")
    if not os.path.exists(fasta_path):
        print(f"⚠️ No unitig FASTA found for {isl_id}, skipping.")
        continue

    unitig_length_dict = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
    df['unitig_length'] = df['query_name'].map(unitig_length_dict).fillna(0).astype(int)

    # Identify sample columns (everything that isn’t a helper column)
    helper_cols = {'query_name', 'isl_label', 'unitig_length', 'presence_count'}
    genome_columns = [c for c in df.columns if c not in helper_cols]

    # Presence count per unitig (NA -> 0; present if value > 0.5)
    presence_matrix = df[genome_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
    df['presence_count'] = (presence_matrix > 0.5).sum(axis=1)

    # Summary stats
    n_genomes = len(genome_columns)
    core_count = (df['presence_count'] == n_genomes).sum()
    unique_count = (df['presence_count'] == 1).sum()

    summary_list.append({
        'file': filepath,
        'total_unitigs': len(df),
        'core_unitigs': core_count,
        'unique_unitigs': unique_count
    })

    # For unitig scatter (length vs presence)
    all_scatter_dfs.append(df[['unitig_length', 'presence_count', 'isl_label']])

# =========================
# UNITIG ANALYSIS PLOTS
# =========================

# 1) Summary table (stdout) + barplot
summary_df = pd.DataFrame(summary_list)
print("Summary of Core vs. Unique Unitigs per File:")
if not summary_df.empty:
    print(summary_df.to_string(index=False))

if not summary_df.empty:
    summary_df['file_label'] = summary_df['file'].apply(
        lambda x: os.path.basename(x).replace('_assemblies_colors.tsv', '')
    )
    melted_df = summary_df.melt(
        id_vars='file_label',
        value_vars=['total_unitigs', 'core_unitigs', 'unique_unitigs'],
        var_name='Unitig_Type',
        value_name='Count'
    )
    plt.figure(figsize=(14, 6))
    sns.barplot(data=melted_df, x='file_label', y='Count', hue='Unitig_Type')
    plt.xticks(rotation=45, ha='right')
    plt.title('Comparison of Total, Core, and Unique Unitigs per ISL')
    plt.xlabel('ISL')
    plt.ylabel('Unitig Count')
    plt.tight_layout()
    plt.legend(title='Type')
    plt.savefig(os.path.join(out_dir, "unitig_summary_barplot.png"))
    plt.close()

# 2) Unitig length vs presence (combined + per-ISL)
if len(all_scatter_dfs):
    combined_scatter_df = pd.concat(all_scatter_dfs, ignore_index=True)

    # Combined
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        data=combined_scatter_df, x="unitig_length", y="presence_count",
        hue="isl_label", s=15, alpha=0.6
    )
    plt.xlabel("Unitig Length (bp)")
    plt.ylabel("Number of Isolates Present")
    plt.title("Unitig Length vs. Isolate Presence (Colored by ISL)")
    plt.legend(title='ISL', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "unitig_length_vs_presence_scatterplot.png"))
    plt.close()

    # Per-ISL small multiples (filtered length > 100)
    filtered_df = combined_scatter_df[combined_scatter_df['unitig_length'] > 100]
    g = sns.FacetGrid(filtered_df, col="isl_label", col_wrap=4, sharex=False, sharey=False, height=4)
    g.map_dataframe(sns.scatterplot, x="unitig_length", y="presence_count", s=15, alpha=0.6)
    g.set_axis_labels("Unitig Length (bp)", "# of Isolates Present")
    g.set_titles(col_template="{col_name}")
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle("Unitig Length > 100 bp vs. Isolate Presence (Per ISL)")
    g.savefig(os.path.join(out_dir, "per_isl_scatterplots_filtered.png"))
    plt.close()
'''
