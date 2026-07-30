import argparse
import os
import subprocess
import pysam
import pandas as pd
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_cmd(cmd):
    subprocess.run(cmd, shell=True, check=True)

def process_sample(sample, ref, read_dir, outdir, region, mq, min_coverage, threads):
    """Aligns a sample and creates an individual VCF."""
    bam = os.path.join(outdir, f"{sample}.bam")
    vcf = os.path.join(outdir, f"{sample}.vcf.gz")
    
    # 1. Locate Reads
    r1_matches = glob.glob(os.path.join(read_dir, f"{sample}_filtered_R1.fastq*")) + glob.glob(os.path.join(read_dir, f"{sample}_filtered_R1.fq*"))
    r2_matches = glob.glob(os.path.join(read_dir, f"{sample}_filtered_R2.fastq*")) + glob.glob(os.path.join(read_dir, f"{sample}_filtered_R2.fq*"))
    
    if not r1_matches or not r2_matches:
        print(f"❌ Skipping {sample}: Could not find matching _filtered_ reads in {read_dir}")
        return sample, None, False
        
    r1, r2 = r1_matches[0], r2_matches[0]
    
    try:
        # 2. Align (if BAM missing)
        if not os.path.exists(bam):
            print(f"🧬 Aligning {sample}...")
            rg_string = f"@RG\\tID:{sample}\\tSM:{sample}"
            run_cmd(f"bwa mem -R '{rg_string}' -t {threads} {ref} {r1} {r2} | samtools sort -@ {threads} -o {bam}")
            run_cmd(f"samtools index {bam}")
            
        # 3. Call Individual VCF (if VCF missing)
        if not os.path.exists(vcf):
            print(f"✂️ Calling VCF for {sample}...")
            run_cmd(f"bcftools mpileup -q {mq} -r {region} -a FORMAT/AD,FORMAT/DP -Ou -f {ref} {bam} | "
                    f"bcftools call -mv -Oz -o {vcf}")
            run_cmd(f"bcftools index {vcf}")
            
        return sample, vcf, True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {sample}: {e}")
        return sample, None, False

def main():
    parser = argparse.ArgumentParser(description="Parallel Individual SNP Calling and Novelty Matrix Generation.")
    parser.add_argument("--samples", required=True, help="TXT file of ALL sample names (1 per line).")
    parser.add_argument("--previous", required=True, help="TXT file of PREVIOUS sample names (1 per line).")
    parser.add_argument("--ref", required=True, help="Path to Reference FASTA.")
    parser.add_argument("--read_dir", required=True, help="Directory containing sample reads.")
    parser.add_argument("--outdir", default="Individual_SNP_Output", help="Output directory.")
    
    parser.add_argument("--region", default="chr6:28000000-34000000", help="Target coordinate block for calling (Default: chr6:28000000-34000000)")
    parser.add_argument("--mq", type=int, default=20, help="Minimum mapping quality (Default: 20)")
    parser.add_argument("--min_coverage", type=int, default=30, help="Minimum depth of coverage (Default: 30)")
    parser.add_argument("--jobs", type=int, default=12, help="Number of concurrent samples to process (Default: 12)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads per sample for BWA/Samtools (Default: 8)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    vcf_dir = os.path.join(args.outdir, "VCFs")
    os.makedirs(vcf_dir, exist_ok=True)

    # 1. Load Samples
    with open(args.samples, 'r') as f:
        all_samples = [line.strip() for line in f if line.strip()]
    with open(args.previous, 'r') as f:
        previous_samples = [line.strip() for line in f if line.strip()]
    
    new_samples = [s for s in all_samples if s not in previous_samples]
    print(f"📊 Loaded {len(all_samples)} total samples ({len(previous_samples)} Previous, {len(new_samples)} New).")
    print(f"🎯 Target Region: {args.region}")

    if not os.path.exists(f"{args.ref}.bwt"):
        print("📇 Indexing reference genome...")
        run_cmd(f"bwa index {args.ref}")

    # 2. Process all samples in PARALLEL (Align -> VCF)
    print(f"🚀 Starting parallel processing ({args.jobs} jobs, {args.threads} threads per job)...")
    valid_vcfs = {}
    
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(process_sample, sample, args.ref, args.read_dir, vcf_dir, args.region, args.mq, args.min_coverage, args.threads) for sample in all_samples]
        
        for future in as_completed(futures):
            sample, vcf, success = future.result()
            if success and vcf:
                valid_vcfs[sample] = vcf

    if not valid_vcfs:
        print("❌ No valid VCF files were generated. Exiting.")
        return

    # 3. First Pass: Build the Master List of Polymorphic Sites
    print("🔍 Pass 1: Building master list of polymorphic positions...")
    master_positions = set()
    ref_map = {}

    for sample, vcf_file in valid_vcfs.items():
        with pysam.VariantFile(vcf_file) as vcf:
            for rec in vcf:
                # Skip complex indels
                if len(rec.ref) != 1 or any(len(a) != 1 for a in rec.alts): 
                    continue
                
                # Check if variant passes filters for this sample
                fmt = rec.samples[sample]
                dp = fmt.get('DP')
                gt = fmt.get('GT')
                
                if dp is None or dp < args.min_coverage or gt is None or None in gt:
                    continue
                if gt[0] != gt[1] or gt[0] == 0: # Skip het or hom-ref
                    continue

                pos_id = f"{rec.chrom}:{rec.pos}"
                master_positions.add(pos_id)
                ref_map[pos_id] = rec.ref

    master_positions = sorted(list(master_positions), key=lambda x: int(x.split(':')[1]))
    
    if not master_positions:
        print("⚠️ No SNPs found passing criteria across any sample. Exiting.")
        return

    print(f"✅ Found {len(master_positions)} unique polymorphic sites across the cohort.")

    # 4. Second Pass: Build the Data Matrix
    print("🧮 Pass 2: Generating Cohort Matrix...")
    matrix_data = {sample: {} for sample in all_samples}
    
    # Initialize all matrix cells to the Reference base mathematically
    for sample in all_samples:
        for pos in master_positions:
            matrix_data[sample][pos] = ref_map[pos]

    # Populate actual variants from individual VCFs
    for sample, vcf_file in valid_vcfs.items():
        with pysam.VariantFile(vcf_file) as vcf:
            for rec in vcf:
                if len(rec.ref) != 1 or any(len(a) != 1 for a in rec.alts): continue
                
                pos_id = f"{rec.chrom}:{rec.pos}"
                if pos_id not in master_positions: continue
                
                fmt = rec.samples[sample]
                dp = fmt.get('DP')
                gt = fmt.get('GT')
                
                if dp is None or dp < args.min_coverage or gt is None or None in gt:
                    matrix_data[sample][pos_id] = 'N' # Insufficient coverage
                    continue
                
                if gt[0] != gt[1]:
                    matrix_data[sample][pos_id] = 'N' # Het
                elif gt[0] != 0:
                    alt_idx = gt[0] - 1
                    matrix_data[sample][pos_id] = rec.alts[alt_idx] # Hom-Alt

    df_all = pd.DataFrame(matrix_data).T 
    df_all = df_all[master_positions] 
    df_all.to_csv(os.path.join(args.outdir, "datatable_all_positions.csv"))

    # 5. Create the 7-Row Novelty Table
    print("🔬 Isolating Novel SNPs (7-Row Table)...")
    novelty_rows = []
    novel_counts_per_sample = {s: 0 for s in new_samples}
    total_snps_per_sample = {s: 0 for s in new_samples}

    for pos in master_positions:
        ref_b = ref_map[pos]
        
        valid_prev = [s for s in previous_samples if s in df_all.index]
        prev_bases = set(df_all.loc[valid_prev, pos]) - {ref_b, 'N'} if valid_prev else set()
        prev_bases = sorted(list(prev_bases))
        
        valid_new = [s for s in new_samples if s in df_all.index]
        new_bases = set(df_all.loc[valid_new, pos]) - {ref_b, 'N'} if valid_new else set()
        
        novel_bases = sorted(list(new_bases - set(prev_bases)))
        
        col_data = [
            ref_b,
            prev_bases[0] if len(prev_bases) > 0 else "-",
            prev_bases[1] if len(prev_bases) > 1 else "-",
            prev_bases[2] if len(prev_bases) > 2 else "-",
            novel_bases[0] if len(novel_bases) > 0 else "-",
            novel_bases[1] if len(novel_bases) > 1 else "-",
            novel_bases[2] if len(novel_bases) > 2 else "-"
        ]
        novelty_rows.append(col_data)

        for s in valid_new:
            s_base = df_all.loc[s, pos]
            if s_base != 'N' and s_base != ref_b:
                total_snps_per_sample[s] += 1
                if s_base in novel_bases:
                    novel_counts_per_sample[s] += 1

    row_labels = ["Ref_Base", "Previous_Alt_1", "Previous_Alt_2", "Previous_Alt_3", "New_Alt_1", "New_Alt_2", "New_Alt_3"]
    df_7row = pd.DataFrame(novelty_rows, index=master_positions, columns=row_labels).T
    df_7row.to_csv(os.path.join(args.outdir, "summary_novelty_positions.csv"))

    # 6. Create High-Level Cohort Summary with TOTAL row
    print("📝 Generating Summary Report...")
    summary_data = []
    
    for s in new_samples:
        if s in total_snps_per_sample:
            summary_data.append({
                "Sample": s,
                "Total_Hom_SNPs_vs_Ref": total_snps_per_sample[s],
                "Novel_SNPs_vs_Previous": novel_counts_per_sample[s]
            })
    
    df_summary = pd.DataFrame(summary_data)
    
    if not df_summary.empty:
        total_row = pd.DataFrame([{
            "Sample": "TOTAL",
            "Total_Hom_SNPs_vs_Ref": df_summary["Total_Hom_SNPs_vs_Ref"].sum(),
            "Novel_SNPs_vs_Previous": df_summary["Novel_SNPs_vs_Previous"].sum()
        }])
        df_summary = pd.concat([df_summary, total_row], ignore_index=True)
        df_summary.to_csv(os.path.join(args.outdir, "per_sample_novel_counts.csv"), index=False)

    print("\n✅ Pipeline Complete!")
    print(f"Outputs saved to: {args.outdir}/")

if __name__ == "__main__":
    main()
