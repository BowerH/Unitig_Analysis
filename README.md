# 🧬 Unitig Analysis and Visualization

This script processes output from [`ggCaller`](https://github.com/phelimb/ggCaller) to generate summary statistics and visualizations for unitig-based bacterial pangenome analysis.

It supports:
- Automatic matching of `.tsv` and `.fa` files per sample
- Calculation of unitig presence across isolates
- Visualization of core/unique unitigs and unitig length vs. isolate presence
- Saving plots directly to a user-specified directory

---

## 📂 Input

You must provide:

- `--tsv_dir`: Directory containing one or more `*_assemblies_colors.tsv` files (from `ggCaller`)
- `--fasta_dir`: Directory containing corresponding `*_assemblies_untigs.fa` files (unitig sequences)
- `--out_dir`: Destination directory for all plots (will be created if it doesn’t exist)

> Each `.tsv` and `.fa` must share the same sample prefix.  
> Example pair:
> - `S0002_ISL_25_1_assemblies_colors.tsv`
> - `S0002_ISL_25_1_assemblies_untigs.fa`

---

## 🧪 Requirements

Install with pip:

```bash
pip install pandas biopython matplotlib seaborn
```

## Usage

```
python Unitig_Analysis.py \
  --tsv_dir ./GGCaller_Output/ \
  --fasta_dir ./GGCaller_Output/ \
  --out_dir ./Unitig_Plots/
```

## 📊 Outputs (saved to --out_dir)

#### unitig_summary_barplot.png:
Barplot comparing total, core, and unique unitigs per ISL
#### unitig_length_vs_presence_scatterplot.png:
Scatterplot of unitig length vs. isolate presence, colored by ISL
#### per_isl_scatterplots_filtered.png:
Faceted scatterplots for each ISL, showing only unitigs > 100bp


## 🧠 Example

Running on a directory of 12 .tsv/.fa pairs:

```
python Unitig_Analysis.py \
  --tsv_dir data/colors/ \
  --fasta_dir data/fasta/ \
  --out_dir results/plots/
```

Results:

3 figures saved to results/plots/

Printed summary table of unitig counts per file


## 🔧 Customization

Change the unitig length filter in the script (default: > 100 bp)

Add plt.xscale('log') if your unitigs span orders of magnitude

Modify to export CSVs or summary tables as needed
