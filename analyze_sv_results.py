import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def run_novelty_metrics(novel_file, outdir):
    print(f"🔬 Analyzing SV Novelty Metrics from {novel_file}...")
    
    if not os.path.exists(novel_file):
        print("⚠️ No novel SVs file found. Skipping metrics.")
        return
        
    df = pd.read_csv(novel_file)
    
    # Total unique variants vs total unique starting positions
    total_alleles = len(df)
    total_positions = df['POS'].nunique()
    
    print(f"   ➡️ Unique Genomic Positions with Novel SV: {total_positions}")
    print(f"   ➡️ Unique Novel SVs Discovered: {total_alleles}")
    
    metrics_data = [
        {"Metric": "Total_Unique_Genomic_Positions_with_Novel_SV", "Count": total_positions},
        {"Metric": "Total_Unique_Novel_SVs_Discovered", "Count": total_alleles}
    ]
    pd.DataFrame(metrics_data).to_csv(os.path.join(outdir, "Novelty_SV_Metrics.csv"), index=False)

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

def run_sv_density(matrix_file, novel_file, anno_file, outdir, window):
    print(f"\n🧮 Generating Dual SV Density Data (Window: {window}bp)...")
    
    # 1. Get TOTAL SV positions
    df_tot = pd.read_csv(matrix_file)
    df_tot = df_tot[df_tot['CHROM'].astype(str).str.contains('chr6')]
    
    # If multiple SVs start at the exact same position, count them all
    tot_pos_counts = df_tot['POS'].value_counts()
    
    # 2. Get NOVEL SV positions
    if os.path.exists(novel_file):
        df_nov = pd.read_csv(novel_file)
        df_nov = df_nov[df_nov['CHROM'].astype(str).str.contains('chr6')]
        nov_pos_counts = df_nov['POS'].value_counts()
    else:
        nov_pos_counts = pd.Series(dtype=int)

    if tot_pos_counts.empty:
        print("❌ No valid chr6 SV positions found!")
        return

    min_pos = int(tot_pos_counts.index.min())
    max_pos = int(tot_pos_counts.index.max())
    idx = np.arange(min_pos, max_pos + 1)
    
    # Calculate Total Density
    s_total = pd.Series(0, index=idx)
    s_total.update(tot_pos_counts)
    density_total = s_total.rolling(window=window, center=True).sum().dropna()
    
    # Calculate Novel Density
    s_novel = pd.Series(0, index=idx)
    s_novel.update(nov_pos_counts)
    density_novel = s_novel.rolling(window=window, center=True).sum().dropna()

    # Export Data Tables (Filtered to >0 to save space)
    df_tot_export = density_total[density_total > 0].reset_index()
    df_tot_export.columns = ['Position_chr6', f'Total_SVs_in_{window}bp']
    df_tot_export.to_csv(os.path.join(outdir, "Total_SV_Density_Data.csv"), index=False)

    df_nov_export = density_novel[density_novel > 0].reset_index()
    df_nov_export.columns = ['Position_chr6', f'Novel_SVs_in_{window}bp']
    df_nov_export.to_csv(os.path.join(outdir, "Novel_SV_Density_Data.csv"), index=False)

    # 3. Generate the Subplots
    print("🎨 Rendering dual-layer annotated SV plot...")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    
    step = max(1, len(density_total) // 20000)
    x = density_total.index[::step] / 1_000_000
    y_total = density_total.values[::step]
    y_novel = density_novel.values[::step]

    # Plot 1: Total
    ax1.plot(x, y_total, color='darkorange', linewidth=1.2)
    ax1.fill_between(x, y_total, color='moccasin', alpha=0.6)
    ax1.set_title(f'Total SV Density on chr6 ({window}bp Window)', fontsize=14, fontweight='bold', pad=10)
    ax1.set_ylabel(f'Total SVs', fontsize=12)

    # Plot 2: Novel
    ax2.plot(x, y_novel, color='purple', linewidth=1.2)
    ax2.fill_between(x, y_novel, color='plum', alpha=0.6)
    ax2.set_title(f'Novel SV Density on chr6 ({window}bp Window)', fontsize=14, fontweight='bold', pad=10)
    ax2.set_ylabel(f'Novel SVs', fontsize=12)
    ax2.set_xlabel('Genomic Position on chr6 (Mb)', fontsize=12)

    # 4. Add Gene Annotations
    genes = get_gene_positions(anno_file)
    min_x, max_x = min_pos / 1_000_000, max_pos / 1_000_000
    
    if genes:
        min_x = min(min_x, min([v[0] for v in genes.values()])/1_000_000)
        max_x = max(max_x, max([v[1] for v in genes.values()])/1_000_000)

    for ax in [ax1, ax2]:
        for i, (label, (start, end)) in enumerate(genes.items()):
            start_mb, end_mb = start / 1_000_000, end / 1_000_000
            mid_mb = (start_mb + end_mb) / 2
            
            # Solid black vertical line at the midpoint
            ax.axvline(x=mid_mb, color='gray', linewidth=1.5, linestyle='--', zorder=1, alpha=0.3)
            
            # Stagger text heights with white background box
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
    plot_out = os.path.join(outdir, "Dual_Annotated_SV_Density.png")
    plt.savefig(plot_out, dpi=300, bbox_inches='tight')
    print(f"🖼️  Plot image saved to {plot_out}")

def main():
    parser = argparse.ArgumentParser(description="Calculate Novelty Metrics and Dual SV Density.")
    parser.add_argument("--matrix", required=True, help="Path to datatable_all_SVs_presence.csv")
    parser.add_argument("--novel", required=True, help="Path to detailed_novel_SVs_list.csv")
    parser.add_argument("--anno", help="Path to chr6_Annotations.csv")
    parser.add_argument("--outdir", default="SV_Downstream_Analysis", help="Directory to save the outputs")
    # Bumped default window to 5000 since SVs are rarer than SNPs
    parser.add_argument("--window", type=int, default=5000, help="Rolling window size in bp (Default: 5000)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    
    print("========================================")
    print("    SV DOWNSTREAM ANALYSIS & PLOTTING   ")
    print("========================================\n")
    
    run_novelty_metrics(args.novel, args.outdir)
    run_sv_density(args.matrix, args.novel, args.anno, args.outdir, args.window)
    
    print("\n✅ All analyses complete!")

if __name__ == "__main__":
    main()
