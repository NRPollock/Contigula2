import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches

def read_fasta_length(filepath):
    """Returns the total sequence length of the FASTA."""
    length = 0
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                length += len(line.strip())
    return length

def parse_gff_genes(gff_file):
    """Extracts all genes and their start/end coordinates from a GFF."""
    genes = {}
    if not os.path.exists(gff_file):
        return genes
        
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) > 8 and parts[2] in ['CDS', 'gene', 'mRNA']:
                attr = parts[8]
                gene_name = None
                for a in attr.split(';'):
                    if a.startswith('Name=') or a.startswith('gene=') or a.startswith('gene_id='):
                        gene_name = a.split('=')[1]
                        break
                
                if gene_name:
                    start, end = int(parts[3]), int(parts[4])
                    if gene_name not in genes:
                        genes[gene_name] = [start, end]
                    else:
                        genes[gene_name][0] = min(genes[gene_name][0], start)
                        genes[gene_name][1] = max(genes[gene_name][1], end)
    return genes

def main():
    parser = argparse.ArgumentParser(description="Graph trimmed assemblies aligned by anchor genes.")
    parser.add_argument("--samples", required=True, help="TXT file with one sample per line")
    parser.add_argument("--gff_directory", required=True, help="Directory containing the reannotated GFF files")
    parser.add_argument("--trimmed_assemblies_dir", required=True, help="Directory containing the trimmed FASTA files")
    parser.add_argument("--out", default="Trimmed_Assemblies_Graph.png", help="Output image filename (.png)")
    parser.add_argument("--key", default="Graph_Sample_Key.csv", help="Output CSV mapping graph Y-axis positions to Sample Names")
    args = parser.parse_args()

    gff_dir = os.path.abspath(args.gff_directory)
    fasta_dir = os.path.abspath(args.trimmed_assemblies_dir)

    with open(args.samples, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]

    print(" Parsing FASTAs and new GFFs...")
    
    assembly_data = {}
    longest_sample = None
    max_len = 0

    # 1. Load Data & Identify Longest Assembly
    for sample in samples:
        fasta_file = os.path.join(fasta_dir, f"{sample}_CDS_trimmed_MHC.fasta")
        gff_file = os.path.join(gff_dir, f"{sample}_reannotated.gff")
        
        if not os.path.exists(fasta_file):
            print(f"⚠️ Missing FASTA for {sample}. Skipping.")
            continue
        if not os.path.exists(gff_file):
            print(f"⚠️ Missing GFF for {sample}. Skipping.")
            continue
            
        seq_len = read_fasta_length(fasta_file)
        genes = parse_gff_genes(gff_file)
        
        assembly_data[sample] = {
            'length': seq_len,
            'genes': genes,
            'shift': 0,
            'anchor_used': None,
            'is_full': False
        }
        
        if seq_len > max_len:
            max_len = seq_len
            longest_sample = sample

    if not longest_sample:
        print(" No valid assemblies found to plot.")
        return

    longest_data = assembly_data[longest_sample]
    print(f" Longest Assembly Selected: {longest_sample} ({longest_data['length']:,} bp)")

    # 2. Build the gene_markers from the Longest GFF
    target_genes = [
        'HLA-G', 'HLA-A', 'MUC22', 'HLA-C', 'HLA-B', 'MICA', 
        'C4', 'HLA-DRA', 'HLA-DRB1', 'HLA-DQ', 'HLA-DP'
    ]
    
    gene_markers = {}
    longest_genes = longest_data['genes']
    
    for tg in target_genes:
        if tg in longest_genes:
            gene_markers[tg] = longest_genes[tg][0]
        else:
            for mg in longest_genes.keys():
                if tg in mg:
                    gene_markers[tg] = longest_genes[mg][0]
                    break

    # 3. Calculate alignment shifts & determine if "Full"
    for sample, data in assembly_data.items():
        genes = data['genes']
        
        # --- NEW BIOLOGICAL FULL ASSEMBLY CHECK ---
        genes_upper = [str(g).upper() for g in genes.keys()]
        has_start = any('GPX5' in g for g in genes_upper)
        has_end = any('ZBTB9' in g for g in genes_upper)
        data['is_full'] = has_start and has_end
        # ------------------------------------------

        if sample == longest_sample:
            data['anchor_used'] = "LONGEST (REFERENCE)"
            continue
            
        shift = 0
        anchor = None
        
        if 'HLA-A' in genes and 'HLA-A' in longest_genes:
            shift = longest_genes['HLA-A'][0] - genes['HLA-A'][0]
            anchor = "HLA-A"
        elif 'HLA-DRB1' in genes and 'HLA-DRB1' in longest_genes:
            shift = longest_genes['HLA-DRB1'][0] - genes['HLA-DRB1'][0]
            anchor = "HLA-DRB1"
        else:
            center_pos = data['length'] / 2
            shared_genes = [g for g in genes.keys() if g in longest_genes]
            
            if shared_genes:
                closest_gene = min(shared_genes, key=lambda g: abs(genes[g][0] - center_pos))
                shift = longest_genes[closest_gene][0] - genes[closest_gene][0]
                anchor = f"Middle-Fallback ({closest_gene})"
            else:
                shift = 0
                anchor = "None (Unaligned)"

        data['shift'] = shift
        data['anchor_used'] = anchor

    # ==========================================
    # 4. PLOTTING AESTHETICS 
    # ==========================================
    print("Generating Plot...")
    
    # Sort samples: 
    # 1st Priority: Is it a "Full" assembly? (GPX5 to ZBTB9)
    # 2nd Priority: Earliest shifted start
    # 3rd Priority: Longest assembly
    sorted_samples = sorted(
        assembly_data.keys(),
        key=lambda s: (
            0 if assembly_data[s]['is_full'] else 1, 
            assembly_data[s]['shift'], 
            -assembly_data[s]['length']
        )
    )

    fig, ax = plt.subplots(figsize=(12, 8))
    y_positions = list(range(len(sorted_samples), 0, -1))
    
    starts = [assembly_data[s]['shift'] for s in sorted_samples]
    lengths = [assembly_data[s]['length'] for s in sorted_samples]
    
    # Color logic: "Full" if both anchor genes are present
    bar_colors = ['darkblue' if assembly_data[s]['is_full'] else 'royalblue' for s in sorted_samples]

    # Draw the main assembly bars
    ax.barh(y=y_positions, 
            width=lengths, 
            left=starts, 
            color=bar_colors,  
            edgecolor=bar_colors, 
            linewidth=0.5,
            height=0.6,
            zorder=2)

    ax.set_ylim(bottom=0)
    y_min, y_max = ax.get_ylim()

    # Draw the Vertical Gene Markers (Red Dotted)
    for gene, pos in gene_markers.items():
        ax.axvline(x=pos, color='red', linestyle=':', alpha=0.6, zorder=3)
        ax.text(x=pos, 
                y=y_max, 
                s=gene, 
                color='black', 
                fontsize=10, 
                rotation=90, 
                verticalalignment='bottom', 
                horizontalalignment='center')

    # Formatting & Aesthetics
    ax.set_title('Trimmed MHC Assemblies Aligned to Longest Assembly', fontsize=14, pad=60)
    ax.set_xlabel('Coordinates Relative to Longest Assembly (Mb)', fontsize=12)
    ax.set_ylabel('Sample Number', fontsize=12)

    # Create and export the Decoder Ring CSV
    mapping_df = pd.DataFrame({
        'Y_Axis_Position': y_positions,
        'Sample_Name': sorted_samples,
        'Length_bp': lengths,
        'Is_Full_Assembly': [assembly_data[s]['is_full'] for s in sorted_samples],
        'Aligned_Via': [assembly_data[s]['anchor_used'] for s in sorted_samples]
    })
    mapping_df.to_csv(args.key, index=False)
    print(f"Decoder Key saved as {args.key}")

    # Custom Legend
    full_patch = mpatches.Patch(color='darkblue', label='Full Assembly (GPX5 to ZBTB9)')
    partial_patch = mpatches.Patch(color='royalblue', label='Partial Assembly')
    ax.legend(handles=[full_patch, partial_patch], 
              loc='lower left', 
              framealpha=0.9, 
              edgecolor='black')

    # Force ticks every 500,000 base pairs (0.5 Mb)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(500000))
    formatter = ticker.FuncFormatter(lambda x, pos: f'{x / 1_000_000:.2f}')
    ax.xaxis.set_major_formatter(formatter)

    ax.grid(axis='y', linestyle='-', alpha=0.2, zorder=1)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    min_x = min(starts)
    max_x = max([s + l for s, l in zip(starts, lengths)])
    padding = 200000
    ax.set_xlim(min_x - padding, max_x + padding)

    # Save and Display
    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    print(f"Graph successfully saved as {args.out}")

if __name__ == "__main__":
    main()
