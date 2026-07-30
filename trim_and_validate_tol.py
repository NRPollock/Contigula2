import argparse
import os
import pandas as pd
import subprocess
import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==========================================
# 1. PATH CONFIGURATION
# ==========================================
PRE_TRIMMED_DIR = "/home/nick/more_shizz/pre_trimmed_assemblies"
SHORT_READ_DIR = "/home/nick/more_shizz/short_reads"
LONG_READ_DIR = "/home/nick/more_shizz/long_reads"
ANNOTATION_BASE = "/home/nick/asstest/merged"

def run_command(cmd, log_file=None):
    try:
        if log_file:
            with open(log_file, "w") as f:
                subprocess.run(cmd, shell=True, check=True, stdout=f, stderr=subprocess.STDOUT)
        else:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        err_msg = f"Error running command: {cmd}"
        if not log_file:
            err_msg += f"\n{e.stderr.decode()}"
        print(err_msg)
        raise

# ==========================================
# 2. VALIDATION LOGIC
# ==========================================
def validate_assembly(sample, assembly_fasta, threads, outdir):
    print(f"[{sample}] 🔍 Validating CDS assembly...")
    
    base_name = f"{sample}_CDS_val"
    
    sr1 = os.path.join(SHORT_READ_DIR, f"{sample}_mhc_sr_R1.fastq")
    sr2 = os.path.join(SHORT_READ_DIR, f"{sample}_mhc_sr_R2.fastq")
    lr = os.path.join(LONG_READ_DIR, f"{sample}_mhc_longreads.fastq")
    
    sr_bam = os.path.join(outdir, f"{base_name}_sr.bam")
    lr_bam = os.path.join(outdir, f"{base_name}_lr.bam")
    
    bed_out = os.path.join(outdir, f"{base_name}_errors.bed")
    csv_out = os.path.join(outdir, f"{base_name}_errors.csv")
    
    # Run mapping (Skipped if BAMs already exist AND are not empty)
    if not os.path.exists(f"{assembly_fasta}.bwt") or os.path.getsize(f"{assembly_fasta}.bwt") == 0:
        run_command(f"bwa index {assembly_fasta}")
        
    if not os.path.exists(sr_bam) or os.path.getsize(sr_bam) == 0:
        run_command(f"bwa mem -t {threads} {assembly_fasta} {sr1} {sr2} | samtools sort -@ {threads} -o {sr_bam}")
        run_command(f"samtools index {sr_bam}")
        
    if not os.path.exists(lr_bam) or os.path.getsize(lr_bam) == 0:
        run_command(f"minimap2 -ax map-ont -t {threads} {assembly_fasta} {lr} | samtools sort -@ {threads} -o {lr_bam}")
        run_command(f"samtools index {lr_bam}")

    # BED/CSV Generation
    flagged_regions = []
    window_size = 1000
    
    with pysam.AlignmentFile(lr_bam, "rb") as lr_sam, pysam.AlignmentFile(sr_bam, "rb") as sr_sam:
        for ref_name in lr_sam.references:
            ref_len = lr_sam.get_reference_length(ref_name)
            
            for start in range(0, ref_len, window_size):
                end = min(start + window_size, ref_len)
                
                # Check Long Read Coverage
                lr_cov = lr_sam.count_coverage(ref_name, start, end, quality_threshold=0)
                lr_win_bases = sum(sum(arr) for arr in lr_cov)
                if lr_win_bases == 0:
                    flagged_regions.append((ref_name, start, end, "Zero LR Coverage"))
                
                # Check Short Read Coverage
                sr_cov = sr_sam.count_coverage(ref_name, start, end, quality_threshold=0)
                sr_win_bases = sum(sum(arr) for arr in sr_cov)
                if sr_win_bases == 0:
                    flagged_regions.append((ref_name, start, end, "Zero SR Coverage"))

    if flagged_regions:
        with open(bed_out, "w") as f_bed:
            for chrom, start, end, reason in flagged_regions:
                f_bed.write(f"{chrom}\t{start}\t{end}\t{reason}\n")
                
        with open(csv_out, "w") as f_csv:
            f_csv.write("Contig,Start_Position,End_Position,Error_Reason\n")
            for chrom, start, end, reason in flagged_regions:
                f_csv.write(f"{chrom},{start},{end},{reason}\n")
                
        print(f"[{sample}] Found {len(flagged_regions)} zero-coverage gaps. Saved to .bed and .csv")
    else:
        print(f"[{sample}] CDS Assembly intact! No coverage gaps found.")

# ==========================================
# 3. HELPER FUNCTIONS
# ==========================================
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

def write_fasta(filepath, header, seq):
    with open(filepath, 'w') as f:
        f.write(f"{header}\n")
        for i in range(0, len(seq), 80):
            f.write(f"{seq[i:i+80]}\n")

def get_cds_coords_from_gff(gff_file, target_gene):
    coords = []
    if not target_gene:
        return coords
        
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) > 8 and parts[2] == 'CDS':
                attr = parts[8]
                if target_gene in attr:
                    coords.extend([int(parts[3]), int(parts[4])])
    return coords

def process_sample_wrapper(args):
    return process_sample(*args)

def process_sample(sample, loh_df, snp_df, threads, continue_run, outdir):
    print(f"\n[{sample}] Starting Processing...")
    
    cds_fasta_out = os.path.join(outdir, f"{sample}_CDS_trimmed_MHC.fasta")
    
    if continue_run and os.path.exists(cds_fasta_out):
        if os.path.getsize(cds_fasta_out) > 0:
            print(f"[{sample}] --continue flag active. Existing valid FASTA found. Skipping trim and jumping to validation.")
            validate_assembly(sample, cds_fasta_out, threads, outdir)
            return sample
        else:
            print(f"[{sample}] --continue flag active, but existing FASTA is empty (0 bytes). Overwriting...")

    sample_loh = loh_df[loh_df['Sample'] == sample]
    if sample_loh.empty or sample_loh['LOH_Start_Position'].values[0] == 'ERROR':
        print(f"[{sample}] Skipping: Not found in LOH summary.")
        return sample

    hg38_start = int(sample_loh['LOH_Start_Position'].values[0])
    hg38_end = int(sample_loh['LOH_End_Position'].values[0])
    
    gff_file = os.path.join(ANNOTATION_BASE, f"{sample}_assembly", f"{sample}_mhca", f"{sample}_MHC_GPX5_ZBTB9.gff")
    fasta_file = os.path.join(PRE_TRIMMED_DIR, f"{sample}_trimmed_MHC.fasta")
    
    if not all(os.path.exists(f) for f in [gff_file, fasta_file]):
        print(f"[{sample}] Skipping: Missing GFF or FASTA file.")
        return sample

    seqs = read_fasta(fasta_file)
    header, full_seq = list(seqs.items())[0]
    total_len = len(full_seq)

    # ==========================================
    # CDS TRIMMING
    # ==========================================
    annotated_snps = snp_df.dropna(subset=['gene/annotation']).copy()
    invalid_annotations = ['0', '0.0', 'NA', 'N/A', 'NAN', 'NONE', '']
    valid_mask = ~annotated_snps['gene/annotation'].astype(str).str.strip().str.upper().isin(invalid_annotations)
    clean_snps = annotated_snps[valid_mask]
    
    bypass_end = hg38_end > 33377424
    if bypass_end:
        print(f"[{sample}] LOH End ({hg38_end}) is past 33377424. Bypassing downstream trim.")
        end_gene = "ZBTB9" 
        end_cds_coords = []
    else:
        downstream_snps = clean_snps[clean_snps['hg38_position'] >= hg38_end].sort_values('hg38_position', ascending=True)
        end_gene = None
        end_cds_coords = []
        
        for index, row in downstream_snps.iterrows():
            candidate_gene = row['gene/annotation']
            coords = get_cds_coords_from_gff(gff_file, candidate_gene)
            if coords:
                end_gene = candidate_gene
                end_cds_coords = coords
                break
            else:
                print(f"[{sample}] End anchor '{candidate_gene}' not found in GFF. Trying next...")
                
        if not end_gene:
            print(f"[{sample}] WARNING: Exhausted all end anchor candidates. Skipping.")
            return sample

    upstream_snps = clean_snps[clean_snps['hg38_position'] <= hg38_start].sort_values('hg38_position', ascending=False)
    start_gene = None
    start_cds_coords = []
    
    for index, row in upstream_snps.iterrows():
        candidate_gene = row['gene/annotation']
        
        if 'GPX5' in candidate_gene.upper():
            print(f"[{sample}] Anchor is GPX5. Bypassing upstream trim.")
            start_gene = candidate_gene
            bypass_start = True
            break
            
        coords = get_cds_coords_from_gff(gff_file, candidate_gene)
        if coords:
            start_gene = candidate_gene
            start_cds_coords = coords
            bypass_start = False
            break
        else:
            print(f"[{sample}] Start anchor '{candidate_gene}' not found in GFF. Trying next...")
            
    if not start_gene:
        print(f"[{sample}] WARNING: Exhausted all start anchor candidates. Skipping.")
        return sample

    print(f"[{sample}] Final Resolved Anchors -> Start: {start_gene}, End: {end_gene}")

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

    print(f"[{sample}] Final Assembly Slice: {cds_trim_start} to {cds_trim_end}")

    cds_seq = full_seq[cds_trim_start - 1 : cds_trim_end]
    
    dynamic_header = f">{sample}_MHC_{start_gene}_{end_gene}"
    write_fasta(cds_fasta_out, dynamic_header, cds_seq)
    
    validate_assembly(sample, cds_fasta_out, threads, outdir)

    print(f"[{sample}] Finished successfully.")
    return sample

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim MHC assemblies by GFF CDS annotations and automatically validate (Parallel).")
    parser.add_argument("--samples", required=True, help="TXT file with one sample per line")
    parser.add_argument("--loh", required=True, help="Path to the LOH Summary CSV")
    parser.add_argument("--snps", required=True, help="Path to the MEGA_MHC_SNPs.xlsx file")
    parser.add_argument("--jobs", type=int, default=2, help="Number of samples to process concurrently")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads per sample for BWA/Minimap2")
    parser.add_argument("--continue", dest="continue_run", action="store_true", help="Skip trimming if output FASTA exists and is not empty.")
    parser.add_argument("--outdir", default="/home/nick/more_shizz/trimmed_CDS", help="Directory to save output files")
    args = parser.parse_args()

    out_dir = os.path.abspath(args.outdir)
    os.makedirs(out_dir, exist_ok=True)

    loh_df = pd.read_csv(args.loh)
    print("Loading SNP annotations...")
    snp_df = pd.read_excel(args.snps)

    with open(args.samples, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]

    print(f"Found {len(samples)} samples to process.")
    print(f"Outputs will be saved to: {out_dir}")
    
    tasks = [(sample, loh_df, snp_df, args.threads, args.continue_run, out_dir) for sample in samples]
    
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(process_sample_wrapper, task): task[0] for task in tasks}
        
        for future in as_completed(futures):
            sample_name = futures[future]
            try:
                future.result()
            except Exception as exc:
                print(f"[{sample_name}] FATAL ERROR: {exc}")
                
    print("\n Trimming and validation complete!")
