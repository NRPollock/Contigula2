import argparse
import os
import subprocess
import pysam
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==========================================
# 1. HELPER FUNCTIONS
# ==========================================
def run_command(cmd, sample_name="System"):
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"[{sample_name}] Error running command: {cmd}\n{e.stderr.decode()}")
        raise

def read_fasta(filepath):
    seqs = {}
    curr_header = None
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                curr_header = line
                seqs[curr_header] = []
            elif curr_header:
                seqs[curr_header].append(line)
    for k in seqs:
        seqs[k] = "".join(seqs[k])
    return seqs

def get_cds_coords_from_gff(gff_file, target_gene):
    coords = []
    if not target_gene or not os.path.exists(gff_file):
        return coords
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) > 8 and parts[2] in ['CDS', 'gene', 'mRNA']:
                attr = parts[8]
                if target_gene in attr:
                    coords.extend([int(parts[3]), int(parts[4])])
    return coords

def parse_gff_genes(gff_file):
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

def find_pre_trimmed_gff(sample, base_dir):
    flat_path = os.path.join(base_dir, f"{sample}_MHC_GPX5_ZBTB9.gff")
    if os.path.exists(flat_path): return flat_path
    nested_path = os.path.join(base_dir, f"{sample}_assembly", f"{sample}_mhca", f"{sample}_MHC_GPX5_ZBTB9.gff")
    if os.path.exists(nested_path): return nested_path
    matches = glob.glob(os.path.join(base_dir, "**", f"{sample}_*.gff"), recursive=True)
    return matches[0] if matches else None

# ==========================================
# 2. DYNAMIC REGION CALCULATION
# ==========================================
def calculate_trim_region(sample, pre_trim_fasta, pre_trimmed_gff_dir, loh_df, snp_df):
    sample_loh = loh_df[loh_df['Sample'] == sample]
    if sample_loh.empty or sample_loh['LOH_Start_Position'].values[0] == 'ERROR':
        print(f"[{sample}] Skipping: Not found in LOH summary.")
        return None

    hg38_start = int(sample_loh['LOH_Start_Position'].values[0])
    hg38_end = int(sample_loh['LOH_End_Position'].values[0])
    
    gff_file = find_pre_trimmed_gff(sample, pre_trimmed_gff_dir)
    if not gff_file or not os.path.exists(pre_trim_fasta):
        print(f"[{sample}] Skipping: Missing Pre-trimmed GFF or Pre-Trimmed FASTA.")
        return None

    seqs = read_fasta(pre_trim_fasta)
    header, full_seq = list(seqs.items())[0]
    total_len = len(full_seq)
    contig_name = header.lstrip(">").split()[0]

    annotated_snps = snp_df.dropna(subset=['gene/annotation']).copy()
    invalid_annotations = ['0', '0.0', 'NA', 'N/A', 'NAN', 'NONE', '']
    valid_mask = ~annotated_snps['gene/annotation'].astype(str).str.strip().str.upper().isin(invalid_annotations)
    clean_snps = annotated_snps[valid_mask]
    
    bypass_end = hg38_end > 33377424
    end_gene = None
    end_cds_coords = []
    
    if bypass_end:
        end_gene = "ZBTB9" 
    else:
        downstream_snps = clean_snps[clean_snps['hg38_position'] >= hg38_end].sort_values('hg38_position', ascending=True)
        for _, row in downstream_snps.iterrows():
            candidate_gene = row['gene/annotation']
            coords = get_cds_coords_from_gff(gff_file, candidate_gene)
            if coords:
                end_gene = candidate_gene
                end_cds_coords = coords
                break

    upstream_snps = clean_snps[clean_snps['hg38_position'] <= hg38_start].sort_values('hg38_position', ascending=False)
    start_gene = None
    start_cds_coords = []
    bypass_start = False
    
    for _, row in upstream_snps.iterrows():
        candidate_gene = row['gene/annotation']
        if 'GPX5' in candidate_gene.upper():
            start_gene = candidate_gene
            bypass_start = True
            break
        coords = get_cds_coords_from_gff(gff_file, candidate_gene)
        if coords:
            start_gene = candidate_gene
            start_cds_coords = coords
            bypass_start = False
            break

    if not start_gene or not end_gene:
        print(f"[{sample}] WARNING: Exhausted anchor candidates. Skipping.")
        return None

    cds_trim_start = 1
    cds_trim_end = total_len

    if not (bypass_start and bypass_end):
        if not bypass_start and not bypass_end:
            b1 = min(start_cds_coords)
            b2 = max(end_cds_coords)
            if b1 > b2:
                cds_trim_start = max(1, b2 - 500)
                cds_trim_end = min(total_len, b1 + 500)
            else:
                cds_trim_start = max(1, b1 - 500)
                cds_trim_end = min(total_len, b2 + 500)
        elif bypass_start: 
            b2 = max(end_cds_coords)
            cds_trim_start = 1
            cds_trim_end = min(total_len, b2 + 500)
        elif bypass_end: 
            b1 = min(start_cds_coords)
            cds_trim_start = max(1, b1 - 500)
            cds_trim_end = total_len

    region_string = f"{contig_name}:{cds_trim_start}-{cds_trim_end}"
    print(f"[{sample}] Resolved Coordinates: {region_string}")
    return region_string

# ==========================================
# 3. PARALLEL WORKER FUNCTION (TWO-PASS)
# ==========================================
def process_sample(sample, region, pre_trim_asm_dir, trim_asm_dir, reannotated_gff_dir, sr_dir, outdir, blastn_dir, blast_db, window_size, hetz_per_window, min_cov, maf, min_mapq, threads):
    print(f"[{sample}] Mapping & Extracting...")
    
    # Isolate sample into its own output subdirectory
    sample_outdir = os.path.join(outdir, sample)
    os.makedirs(sample_outdir, exist_ok=True)
    
    pre_trim_fasta = os.path.join(pre_trim_asm_dir, f"{sample}_trimmed_MHC.fasta") 
    trim_fasta = os.path.join(trim_asm_dir, f"{sample}_CDS_trimmed_MHC.fasta")
    reannotated_gff = os.path.join(reannotated_gff_dir, f"{sample}_reannotated.gff")
    
    raw_r1 = os.path.join(sr_dir, f"{sample}_mhc_sr_R1.fastq")
    raw_r2 = os.path.join(sr_dir, f"{sample}_mhc_sr_R2.fastq")
    
    pre_trim_bam = os.path.join(sample_outdir, f"{sample}_pre_trimmed_mapped.bam")
    filt_bam = os.path.join(sample_outdir, f"{sample}_region_extracted.bam")
    filt_n_bam = os.path.join(sample_outdir, f"{sample}_region_extracted_namesorted.bam")
    filt_r1 = os.path.join(sample_outdir, f"{sample}_filtered_R1.fastq")
    filt_r2 = os.path.join(sample_outdir, f"{sample}_filtered_R2.fastq")
    
    trim_bam = os.path.join(sample_outdir, f"{sample}_CDS_trimmed_mapped.bam")
    vcf = os.path.join(sample_outdir, f"{sample}_variants.vcf.gz")

    seqs = read_fasta(trim_fasta)
    seq_len = len(list(seqs.values())[0])
    trim_contig = list(seqs.keys())[0].lstrip(">").split()[0]
    
    # --- PASS 1: MAP TO PRE-TRIMMED ASSEMBLY & EXTRACT ---
    if not os.path.exists(filt_r1) or os.path.getsize(filt_r1) == 0:
        if not os.path.exists(f"{pre_trim_fasta}.bwt") or os.path.getsize(f"{pre_trim_fasta}.bwt") == 0:
            run_command(f"bwa index {pre_trim_fasta}", sample)
            
        if not os.path.exists(pre_trim_bam) or os.path.getsize(pre_trim_bam) == 0:
            run_command(f"bwa mem -t {threads} {pre_trim_fasta} {raw_r1} {raw_r2} | samtools sort -@ {threads} -o {pre_trim_bam}", sample)
            run_command(f"samtools index {pre_trim_bam}", sample)
            
        run_command(f"samtools view -q {min_mapq} -b {pre_trim_bam} {region} > {filt_bam}", sample)
        run_command(f"samtools sort -n -@ {threads} {filt_bam} -o {filt_n_bam}", sample)
        run_command(f"samtools fastq -@ {threads} -1 {filt_r1} -2 {filt_r2} -0 /dev/null -s /dev/null {filt_n_bam}", sample)
    else:
        print(f"[{sample}] Filtered reads found and valid. Skipping Pass 1.")

    # --- PASS 2: MAP FILTERED READS TO FINAL TRIMMED ASSEMBLY ---
    if not os.path.exists(trim_bam) or os.path.getsize(trim_bam) == 0:
        if not os.path.exists(f"{trim_fasta}.bwt") or os.path.getsize(f"{trim_fasta}.bwt") == 0:
            run_command(f"bwa index {trim_fasta}", sample)
            
        run_command(f"bwa mem -t {threads} {trim_fasta} {filt_r1} {filt_r2} | samtools view -q {min_mapq} -b | samtools sort -@ {threads} -o {trim_bam}", sample)
        run_command(f"samtools index {trim_bam}", sample)
    else:
        print(f"[{sample}] Trimmed BAM found and valid. Skipping Pass 2.")

    # --- VARIANT CALLING ---
    if not os.path.exists(vcf) or os.path.getsize(vcf) == 0:
        run_command(f"bcftools mpileup -q {min_mapq} -a FORMAT/AD,FORMAT/DP -Ou -f {trim_fasta} {trim_bam} | "
                    f"bcftools call -mv -Oz -o {vcf}", sample)
        run_command(f"bcftools index {vcf}", sample)
    else:
        print(f"[{sample}] VCF found and valid. Skipping Variant Calling.")

    # --- ZYGOSITY CALCULATIONS ---
    window_counts = {}
    max_idx = (seq_len // window_size) + 1
    for i in range(max_idx):
        window_counts[i] = 0

    with pysam.VariantFile(vcf) as vcf_in:
        for rec in vcf_in:
            if not rec.alts: continue
            if len(rec.ref) != 1 or any(len(alt) != 1 for alt in rec.alts): continue

            samp_data = rec.samples[0]
            dp = samp_data.get('DP')
            ad = samp_data.get('AD')
            
            if dp is None or dp < min_cov: continue
            if ad is None or len(ad) < 2: continue
                
            alt_reads = ad[1]
            if alt_reads is None: continue
                
            vaf = alt_reads / dp
            if maf <= vaf <= (1.0 - maf):
                w_idx = (rec.pos - 1) // window_size
                if w_idx in window_counts:
                    window_counts[w_idx] += 1

    blocks = []
    for i in range(max_idx):
        w_start = i * window_size
        w_end = min(w_start + window_size, seq_len)
        if w_start >= seq_len: break
        blocks.append({
            'start': w_start,
            'end': w_end,
            'width': w_end - w_start,
            'is_het': window_counts[i] >= hetz_per_window
        })

    # --- BLASTN LOGIC FOR HETEROZYGOUS BLOCKS (200bp FOCUS) ---
    if blastn_dir and blast_db:
        sample_blast_dir = os.path.join(blastn_dir, sample)
        os.makedirs(sample_blast_dir, exist_ok=True)
        
        for b in blocks:
            if b['is_het']:
                b_start = b['start']
                b_end = b['end']
                
                # Calculate a 200bp focused window in the absolute center of the block
                midpoint = b_start + (b['width'] // 2)
                focus_start = max(b_start, midpoint - 100)
                focus_end = min(b_end, midpoint + 100)
                
                blast_fasta = os.path.join(sample_blast_dir, f"{sample}_het_block_{b_start}_{b_end}_focus200bp.fasta")
                blast_out = os.path.join(sample_blast_dir, f"{sample}_het_block_{b_start}_{b_end}_focus200bp_blast_raw.tsv")
                clean_csv = os.path.join(sample_blast_dir, f"{sample}_het_block_{b_start}_{b_end}_focus200bp_BLAST_SCORES.csv")

                if not os.path.exists(clean_csv):
                    print(f"[{sample}] 🔬 Extracting 10 reads from 200bp focus window ({focus_start}-{focus_end}) inside Block {b_start}-{b_end}...")
                    
                    extract_cmd = f"samtools view {trim_bam} {trim_contig}:{focus_start}-{focus_end} | shuf -n 10 | awk '{{print \">\"$1\"\\n\"$10}}' > {blast_fasta}"
                    run_command(extract_cmd, sample)

                    if os.path.exists(blast_fasta) and os.path.getsize(blast_fasta) > 0:
                        blast_cmd = f"blastn -query {blast_fasta} -db {blast_db} -num_threads {threads} -max_target_seqs 10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -out {blast_out}"
                        run_command(blast_cmd, sample)
                        
                        if os.path.exists(blast_out) and os.path.getsize(blast_out) > 0:
                            df = pd.read_csv(blast_out, sep='\t', header=None, names=[
                                'Read_ID', 'Target_Transcript', 'Percent_Identity', 'Alignment_Length',
                                'Mismatches', 'Gap_Opens', 'Read_Start', 'Read_End',
                                'Transcript_Start', 'Transcript_End', 'E_Value', 'Bit_Score'
                            ])
                            df = df.sort_values(by=['Read_ID', 'Bit_Score'], ascending=[True, False])
                            df.to_csv(clean_csv, index=False)
                            
                            os.remove(blast_out)

    genes = parse_gff_genes(reannotated_gff)
    genes_upper = [str(g).upper() for g in genes.keys()]
    is_full = any('GPX5' in g for g in genes_upper) and any('ZBTB9' in g for g in genes_upper)

    print(f"[{sample}] Finished processing.")
    
    return {
        'sample': sample,
        'length': seq_len,
        'genes': genes,
        'blocks': blocks,
        'shift': 0,
        'anchor_used': None,
        'is_full': is_full
    }

def process_sample_wrapper(args):
    return process_sample(*args)

# ==========================================
# 4. MAIN SCRIPT
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Auto-Trimming Two-Pass Zygosity Plotter with Focused Scored BLAST.")
    parser.add_argument("--samples", required=True, help="TXT file with one sample per line")
    parser.add_argument("--loh", required=True, help="Path to the LOH Summary CSV")
    parser.add_argument("--snps", required=True, help="Path to the MEGA_MHC_SNPs.xlsx file")
    
    parser.add_argument("--pre_trimmed_dir", required=True, help="Directory with SAMPLE_trimmed_MHC.fasta")
    parser.add_argument("--assembly_dir", required=True, help="Directory with SAMPLE_CDS_trimmed_MHC.fasta")
    parser.add_argument("--pre_trimmed_gff_dir", required=True, help="Directory with pre-trimmed GFFs (for calculating slices)")
    parser.add_argument("--reannotated_gff_dir", required=True, help="Directory with SAMPLE_reannotated.gff (for plotting anchors)")
    parser.add_argument("--short_read_dir", required=True, help="Directory with SAMPLE_mhc_sr_R1/R2.fastq")
    parser.add_argument("--outdir", default="Zygosity_Output", help="Directory for all mapping outputs")
    
    parser.add_argument("--blastn_dir", default=None, help="Directory to save per-sample BLAST results for het blocks")
    parser.add_argument("--blast_db", default=None, help="Path to local BLAST database (required if --blastn_dir is used)")
    
    parser.add_argument("--window_size", type=int, default=10000, help="Base pairs per window (default: 10000)")
    parser.add_argument("--hetz_per_window", type=int, default=100, help="Min heterozygous SNPs to call block HET (default: 100)")
    parser.add_argument("--min_coverage", type=int, default=20, help="Min depth to consider a variant (default: 20)")
    parser.add_argument("--MAF", type=float, default=0.3, help="Minor allele frequency threshold for heterozygosity (default: 0.3)")
    parser.add_argument("--min_mapQ", type=int, default=10, help="Min mapping quality for BAM filter and variant calling (default: 10)")
    parser.add_argument("--jobs", type=int, default=2, help="Number of samples to process simultaneously")
    parser.add_argument("--threads", type=int, default=4, help="Threads per job for BWA, Samtools, and BLAST")
    
    parser.add_argument("--out", default="Zygosity_Blocks_Graph.png", help="Output graph filename")
    parser.add_argument("--key", default="Zygosity_Sample_Key.csv", help="Output mapping key CSV")
    args = parser.parse_args()

    if args.blastn_dir and not args.blast_db:
        print("ERROR: You provided --blastn_dir but forgot to provide a --blast_db path!")
        return

    os.makedirs(args.outdir, exist_ok=True)
    if args.blastn_dir:
        os.makedirs(args.blastn_dir, exist_ok=True)
        
    pre_trim_asm_dir = os.path.abspath(args.pre_trimmed_dir)
    trim_asm_dir = os.path.abspath(args.assembly_dir)
    pre_gff_dir = os.path.abspath(args.pre_trimmed_gff_dir)
    reann_gff_dir = os.path.abspath(args.reannotated_gff_dir)
    sr_dir = os.path.abspath(args.short_read_dir)
    out_dir = os.path.abspath(args.outdir)
    
    blastn_dir_path = os.path.abspath(args.blastn_dir) if args.blastn_dir else None

    print("Loading SNP annotations and LOH Data...")
    loh_df = pd.read_csv(args.loh)
    snp_df = pd.read_excel(args.snps)

    with open(args.samples, 'r') as f:
        samples = [line.strip() for line in f if line.strip() and line.strip() != "BSH"]

    print(f"Distributing {len(samples)} samples across {args.jobs} parallel jobs...")

    assembly_data = {}
    longest_sample = None
    max_len = 0

    tasks = []
    for s in samples:
        pre_trim_fasta = os.path.join(pre_trim_asm_dir, f"{s}_trimmed_MHC.fasta")
        region = calculate_trim_region(s, pre_trim_fasta, pre_gff_dir, loh_df, snp_df)
        if not region: continue
        tasks.append((s, region, pre_trim_asm_dir, trim_asm_dir, reann_gff_dir, sr_dir, out_dir, blastn_dir_path, args.blast_db, args.window_size, args.hetz_per_window, args.min_coverage, args.MAF, args.min_mapQ, args.threads))
    
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(process_sample_wrapper, task): task[0] for task in tasks}
        
        for future in as_completed(futures):
            sample_name = futures[future]
            try:
                result = future.result()
                if result:
                    assembly_data[result['sample']] = result
                    if result['length'] > max_len:
                        max_len = result['length']
                        longest_sample = result['sample']
            except Exception as exc:
                print(f"[{sample_name}] ❌ FATAL ERROR: {exc}")

    if not longest_sample:
        print("No valid assemblies processed. Exiting.")
        return

    longest_data = assembly_data[longest_sample]
    longest_genes = longest_data['genes']
    print(f"\n Longest Assembly Selected: {longest_sample} ({longest_data['length']:,} bp)")

    target_genes = ['HLA-G', 'HLA-A', 'MUC22', 'HLA-C', 'HLA-B', 'MICA', 'C4', 'HLA-DRA', 'HLA-DRB1', 'HLA-DQ', 'HLA-DP']
    gene_markers = {}
    for tg in target_genes:
        if tg in longest_genes:
            gene_markers[tg] = longest_genes[tg][0]
        else:
            for mg in longest_genes.keys():
                if tg in mg:
                    gene_markers[tg] = longest_genes[mg][0]
                    break

    for sample, data in assembly_data.items():
        if sample == longest_sample:
            data['anchor_used'] = "LONGEST (REFERENCE)"
            continue
            
        genes = data['genes']
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

    sorted_samples = sorted(
        assembly_data.keys(),
        key=lambda s: (
            0 if assembly_data[s]['is_full'] else 1, 
            assembly_data[s]['shift'], 
            -assembly_data[s]['length']
        )
    )

    het_details_list = []
    
    for sample in sorted_samples:
        data = assembly_data[sample]
        genes = data['genes']
        het_blocks = [b for b in data['blocks'] if b['is_het']]
        
        if not het_blocks:
            het_details_list.append("None")
            continue
            
        sample_het_strings = []
        for b in het_blocks:
            b_start = b['start']
            b_end = b['end']
            
            overlapping_genes = []
            for gene_name, coords in genes.items():
                g_start, g_end = coords
                if g_start <= b_end and g_end >= b_start:
                    overlapping_genes.append(gene_name)
                    
            if overlapping_genes:
                gene_str = ", ".join(overlapping_genes)
                sample_het_strings.append(f"{b_start}-{b_end} (Genes: {gene_str})")
            else:
                sample_het_strings.append(f"{b_start}-{b_end} (Intergenic)")
                
        het_details_list.append(" | ".join(sample_het_strings))

    print("Generating Plot...")

    fig, ax = plt.subplots(figsize=(12, 8))
    y_positions = list(range(len(sorted_samples), 0, -1))
    
    min_x, max_x = float('inf'), float('-inf')

    for i, sample in enumerate(sorted_samples):
        y_pos = y_positions[i]
        data = assembly_data[sample]
        shift = data['shift']
        
        het_blocks = []
        hom_blocks = []
        
        for b in data['blocks']:
            shifted_start = b['start'] + shift
            if b['is_het']:
                het_blocks.append((shifted_start, b['width']))
            else:
                hom_blocks.append((shifted_start, b['width']))
                
            min_x = min(min_x, shifted_start)
            max_x = max(max_x, shifted_start + b['width'])

        if hom_blocks:
            ax.broken_barh(hom_blocks, (y_pos - 0.3, 0.6), facecolors='royalblue', edgecolor='none')
        if het_blocks:
            ax.broken_barh(het_blocks, (y_pos - 0.4, 0.8), facecolors='black', edgecolor='black', linewidth=0.5)

    ax.set_ylim(bottom=0)
    y_min, y_max = ax.get_ylim()

    for gene, pos in gene_markers.items():
        ax.axvline(x=pos, color='red', linestyle=':', alpha=0.6, zorder=3)
        ax.text(x=pos, y=y_max, s=gene, color='black', fontsize=10, rotation=90, 
                verticalalignment='bottom', horizontalalignment='center')

    ax.set_title(f'MHC Zygosity Blocks ({args.window_size//1000}kb Windows, {args.hetz_per_window}+ SNPs/Win = HET)', fontsize=14, pad=60)
    ax.set_xlabel('Coordinates Relative to Longest Assembly (Mb)', fontsize=12)
    ax.set_ylabel('Sample Number', fontsize=12)

    mapping_df = pd.DataFrame({
        'Y_Axis_Position': y_positions,
        'Sample_Name': sorted_samples,
        'Length_bp': [assembly_data[s]['length'] for s in sorted_samples],
        'Is_Full_Assembly': [assembly_data[s]['is_full'] for s in sorted_samples],
        'Aligned_Via': [assembly_data[s]['anchor_used'] for s in sorted_samples],
        'Het_Blocks_Details': het_details_list
    })
    mapping_df.to_csv(os.path.join(args.outdir, args.key), index=False)
    print(f"Decoder Key saved as {args.key}")

    hom_patch = mpatches.Patch(color='royalblue', label='Homozygous Block')
    het_patch = mpatches.Patch(color='black', label='Heterozygous Block')
    ax.legend(handles=[hom_patch, het_patch], loc='lower left', framealpha=0.9, edgecolor='black')

    ax.xaxis.set_major_locator(ticker.MultipleLocator(500000))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x / 1_000_000:.2f}'))

    ax.grid(axis='y', linestyle='-', alpha=0.2, zorder=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    padding = 200000
    ax.set_xlim(min_x - padding, max_x + padding)

    plot_out = os.path.join(args.outdir, args.out)
    plt.tight_layout()
    plt.savefig(plot_out, dpi=300, bbox_inches='tight')
    print(f"Graph successfully saved to {plot_out}")

if __name__ == "__main__":
    main()
