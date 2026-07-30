import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def run_novelty_metrics(summary_file, outdir):
    print(f"🔬 Analyzing Novelty Metrics from {summary_file}...")
    df = pd.read_csv(summary_file, index_col=0)
    novel_rows = df.loc[['New_Alt_1', 'New_Alt_2', 'New_Alt_3']]
    valid_novel_alleles = (novel_rows != '-') & (novel_rows.notna())
    
    total_alleles = valid_novel_alleles.sum().sum()
    total_positions = valid_novel_alleles.any().sum()
    
    print(f"   ➡️ Unique Genomic Positions with Novel SNP: {total_positions}")
    print(f"   ➡️ Unique Novel Alleles Discovered: {total_alleles}")
    
    metrics_data = [
        {"Metric": "Total_Unique_Genomic_Positions_with_Novel_SNP", "Count": total_positions},
        {"Metric": "Total_Unique_Novel_Alleles_Discovered", "Count": total_alleles}
    ]
    pd.DataFrame(metrics_data).to_csv(os.path.join(outdir, "Novelty_Summary_Metrics.csv"), index=False)

def get_gene_positions(anno_file):
    if not anno_file or not os.path.exists(anno_file):
        return {}
        
    df = pd.read_csv(anno_file)
    df['Name'] = df['Name'].str.strip()
    
    gene_groups = {
        'HLA-G': ['HLA-G'],
        'HLA-A': ['HLA-A'],
        'HLA-C': ['HLA-C'],
        'HLA-B': ['HLA-B'],
        'MICA': ['MICA'],
        'MUC22': ['MUC22'],
        'C4': ['C4A', 'C4B'],
        'HLA-DRA': ['HLA-DRA'],
        'HLA-DRB1': ['HLA-DRB1'],
        'HLA-DQ': ['HLA-DQA1', 'HLA-DQB1', 'HLA-DQA2', 'HLA-DQB2'],
        'HLA-DP': ['HLA-DPA1', 'HLA-DPB1', 'HLA-DPA2', 'HLA-DPB2']
    }
    
    pos_dict = {}
    for label, names in gene_groups.items():
        subset = df[df['Name'].isin(names)]
        if not subset.empty:
            pos_dict[label] = (subset['Minimum'].min(), subset['Maximum'].max())
            
    return pos_dict

def run_snp_density(matrix_file, summary_file, anno_file, outdir, window):
    print(f"\n🧮 Generating Dual SNP Density Data (Window: {window}bp)...")
    
    # 1. Get TOTAL positions
    cols = pd.read_csv(matrix_file, nrows=0).columns.tolist()
    total_positions = sorted(list(set([int(c.split(':')[1]) for c in cols if c.startswith('chr6:')])))
    
    # 2. Get NOVEL positions
    df_summary = pd.read_csv(summary_file, index_col=0)
    novel_rows = df_summary.loc[['New_Alt_1', 'New_Alt_2', 'New_Alt_3']]
    valid_novel = (novel_rows != '-') & (novel_rows.notna())
    novel_cols = valid_novel.columns[valid_novel.any()].tolist()
    novel_positions = sorted(list(set([int(c.split(':')[1]) for c in novel_cols if c.startswith('chr6:')])))

    if not total_positions:
        print("❌ No valid chr6 positions found!")
        return

    min_pos, max_pos = min(total_positions), max(total_positions)
    idx = np.arange(min_pos, max_pos + 1)
    
    s_total = pd.Series(0, index=idx)
    s_total.loc[total_positions] = 1
    density_total = s_total.rolling(window=window, center=True).sum().dropna()
    
    s_novel = pd.Series(0, index=idx)
    s_novel.loc[novel_positions] = 1
    density_novel = s_novel.rolling(window=window, center=True).sum().dropna()

    df_tot = density_total[density_total > 0].reset_index()
    df_tot.columns = ['Position_chr6', f'Total_SNPs_in_{window}bp']
    df_tot.to_csv(os.path.join(outdir, "Total_SNP_Density_Data.csv"), index=False)

    df_nov = density_novel[density_novel > 0].reset_index()
    df_nov.columns = ['Position_chr6', f'Novel_SNPs_in_{window}bp']
    df_nov.to_csv(os.path.join(outdir, "Novel_SNP_Density_Data.csv"), index=False)

    # 3. Generate the Subplots
    print("🎨 making annotated plot...")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    
    step = max(1, len(density_total) // 20000)
    x = density_total.index[::step] / 1_000_000
    y_total = density_total.values[::step]
    y_novel = density_novel.values[::step]

    # Plot 1: Total
    ax1.plot(x, y_total, color='firebrick', linewidth=1.2)
    ax1.fill_between(x, y_total, color='lightcoral', alpha=0.4)
    ax1.set_title(f'Total SNP Density on chr6 ({window}bp Window)', fontsize=14, fontweight='bold', pad=10)
    ax1.set_ylabel(f'Total SNPs', fontsize=12)

    # Plot 2: Novel
    ax2.plot(x, y_novel, color='royalblue', linewidth=1.2)
    ax2.fill_between(x, y_novel, color='cornflowerblue', alpha=0.4)
    ax2.set_title(f'Novel SNP Density on chr6 ({window}bp Window)', fontsize=14, fontweight='bold', pad=10)
    ax2.set_ylabel(f'Novel SNPs', fontsize=12)
    ax2.set_xlabel('Genomic Position on chr6 (Mb)', fontsize=12)

    # 4. Add Gene Annotations (Updated to Solid Black Lines)
    genes = get_gene_positions(anno_file)
    min_x, max_x = min_pos / 1_000_000, max_pos / 1_000_000
    
    if genes:
        min_x = min(min_x, min([v[0] for v in genes.values()])/1_000_000)
        max_x = max(max_x, max([v[1] for v in genes.values()])/1_000_000)

    for ax in [ax1, ax2]:
        for i, (label, (start, end)) in enumerate(genes.items()):
            start_mb, end_mb = start / 1_000_000, end / 1_000_000
            mid_mb = (start_mb + end_mb) / 2
            
            # Draw solid black vertical line at the midpoint
            ax.axvline(x=mid_mb, color='gray', linewidth=1.5, linestyle='--', zorder=1, alpha=0.3)
            
            # Stagger text heights and add a white background box for readability
            y_pos = 0.96 - (i % 3) * 0.12 
            ax.text(mid_mb, y_pos, label, transform=ax.get_xaxis_transform(),
                    rotation=90, ha='center', va='top', fontsize=10, color='black', fontweight='bold',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1), zorder=2)
        
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda val, pos: f'{val:.2f}'))
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(min_x - 0.05, max_x + 0.05)
        ax.set_ylim(bottom=0)

    plt.tight_layout()
    plot_out = os.path.join(outdir, "Dual_Annotated_SNP_Density.png")
    plt.savefig(plot_out, dpi=300, bbox_inches='tight')
    print(f"🖼️  Plot image saved to {plot_out}")

def main():
    parser = argparse.ArgumentParser(description="Calculate Novelty Metrics and Dual SNP Density.")
    parser.add_argument("--matrix", required=True, help="Path to datatable_all_positions.csv")
    parser.add_argument("--summary", required=True, help="Path to summary_novelty_positions.csv")
    parser.add_argument("--anno", help="Path to chr6_Annotations.csv (Optional for plot overlay)")
    parser.add_argument("--outdir", default="Downstream_Analysis", help="Directory to save the outputs")
    parser.add_argument("--window", type=int, default=500, help="Rolling window size in bp (Default: 500)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    
    print("========================================")
    print("   SNP DOWNSTREAM ANALYSIS & PLOTTING   ")
    print("========================================\n")
    
    run_novelty_metrics(args.summary, args.outdir)
    run_snp_density(args.matrix, args.summary, args.anno, args.outdir, args.window)
    
    print("\n✅ All analyses complete!")

if __name__ == "__main__":
    main()
