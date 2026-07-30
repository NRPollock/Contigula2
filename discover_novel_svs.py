import argparse
import os
import subprocess
import pysam
import pandas as pd
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_cmd(cmd):
    subprocess.run(cmd, shell=True, check=True)

def align_and_sniffle(sample, ref, read_dir, outdir, threads, region):
    """Aligns Nanopore reads, slices to region, and generates a Sniffles (.snf) file."""
    bam = os.path.join(outdir, f"{sample}.bam")
    snf = os.path.join(outdir, f"{sample}.snf")
    
    if os.path.exists(snf):
        return sample, snf, True

    # Find long reads (supports .fastq, .fastq.gz, .fq)
    lr_matches = glob.glob(os.path.join(read_dir, f"{sample}_filtered_lr.fastq")) + glob.glob(os.path.join(read_dir, f"{sample}*.fq*"))
    
    if not lr_matches:
        print(f"❌ Skipping {sample}: Could not find long reads in {read_dir}")
        return sample, None, False
        
    lr = lr_matches[0]
    
    try:
        # 1. Align and optionally Slice the BAM
        if not os.path.exists(bam):
            print(f"🧬 Aligning {sample} with minimap2...")
            
            if region and region.lower() != "all":
                full_bam = os.path.join(outdir, f"{sample}_full.bam")
                # Map to full genome
                run_cmd(f"minimap2 -ax map-ont -t {threads} --MD {ref} {lr} | samtools sort -@ {threads} -o {full_bam}")
                run_cmd(f"samtools index {full_bam}")
                
                print(f"✂️ Slicing BAM to {region}...")
                run_cmd(f"samtools view -b {full_bam} {region} > {bam}")
                run_cmd(f"samtools index {bam}")
                
                # Cleanup the large temporary BAM to save disk space
                if os.path.exists(full_bam): os.remove(full_bam)
                if os.path.exists(f"{full_bam}.bai"): os.remove(f"{full_bam}.bai")
            else:
                # Standard full genome mapping
                run_cmd(f"minimap2 -ax map-ont -t {threads} --MD {ref} {lr} | samtools sort -@ {threads} -o {bam}")
                run_cmd(f"samtools index {bam}")
            
        # 2. Run Sniffles in single-sample mode
        print(f"👃 Running Sniffles2 on {sample}...")
        run_cmd(f"sniffles --input {bam} --snf {snf} --threads {threads}")
        return sample, snf, True
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {sample}: {e}")
        return sample, None, False

def main():
    parser = argparse.ArgumentParser(description="Discover Novel Structural Variations using Nanopore & Sniffles2.")
    parser.add_argument("--samples", required=True, help="TXT file of ALL sample names (1 per line).")
    parser.add_argument("--previous", required=True, help="TXT file of PREVIOUS sample names (1 per line).")
    parser.add_argument("--ref", required=True, help="Path to Reference FASTA.")
    parser.add_argument("--read_dir", required=True, help="Directory containing Nanopore sample reads.")
    parser.add_argument("--outdir", default="Novel_SV_Output", help="Output directory.")
    
    # New Region Parameter
    parser.add_argument("--region", default="chr6:28000000-34000000", help="Target coordinate block for calling (Default: chr6:28000000-34000000. Use 'all' for whole genome).")
    
    parser.add_argument("--min_sv_len", type=int, default=50, help="Minimum length of SV to report (Default: 50bp)")
    parser.add_argument("--min_support", type=int, default=10, help="Minimum read support to call an SV (Default: 10)")
    parser.add_argument("--jobs", type=int, default=4, help="Concurrent samples to process (Default: 4)")
    parser.add_argument("--threads", type=int, default=8, help="Threads per sample for minimap2/sniffles (Default: 8)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Load Samples
    with open(args.samples, 'r') as f:
        all_samples = [line.strip() for line in f if line.strip()]
    with open(args.previous, 'r') as f:
        previous_samples = [line.strip() for line in f if line.strip()]
    
    new_samples = [s for s in all_samples if s not in previous_samples]
    print(f"📊 Loaded {len(all_samples)} samples ({len(previous_samples)} Previous, {len(new_samples)} New).")
    if args.region.lower() != "all":
        print(f"🎯 Target Region: {args.region}")

    # 2. Align and Generate SNF files in PARALLEL
    print(f"🚀 Starting parallel alignments and SNF generation...")
    snf_files = []
    
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(align_and_sniffle, sample, args.ref, args.read_dir, args.outdir, args.threads, args.region) for sample in all_samples]
        for future in as_completed(futures):
            sample, snf, success = future.result()
            if success and snf:
                snf_files.append(snf)

    if not snf_files:
        print("❌ No valid SNF files were generated. Exiting.")
        return

    # 3. Joint SV Calling (Combines all .snf files into one cohort VCF)
    joint_vcf = os.path.join(args.outdir, "cohort_joint_SVs.vcf")
    snf_list_file = os.path.join(args.outdir, "snf_list.tsv")
    
    with open(snf_list_file, 'w') as f:
        for snf in snf_files:
            f.write(f"{snf}\n")

    if not os.path.exists(joint_vcf):
        print("🌍 Performing Joint SV Calling across the entire cohort...")
        run_cmd(f"sniffles --input {snf_list_file} --vcf {joint_vcf} --threads {args.threads*args.jobs} --minsupport {args.min_support} --minsvlen {args.min_sv_len}")

    # 4. Parse Joint VCF to find NOVEL SVs
    print("🔍 Parsing Joint VCF for Novel Structural Variations...")
    
    sv_matrix_data = []
    novel_svs_detailed = []
    
    total_svs_per_sample = {s: 0 for s in new_samples}
    novel_svs_per_sample = {s: 0 for s in new_samples}

    with pysam.VariantFile(joint_vcf) as vcf:
        for rec in vcf:
            # Get SV attributes
            sv_type = rec.info.get('SVTYPE', 'UNKNOWN')
            sv_len = rec.info.get('SVLEN', 0)
            if isinstance(sv_len, tuple): sv_len = sv_len[0] # Handle tuple returns
            
            # Filter by size
            if abs(sv_len) < args.min_sv_len:
                continue

            sv_id = f"{rec.chrom}:{rec.pos}_{sv_type}_{abs(sv_len)}bp"
            
            prev_has_sv = False
            new_samples_with_sv = []
            row_dict = {"SV_ID": sv_id, "CHROM": rec.chrom, "POS": rec.pos, "TYPE": sv_type, "LENGTH": abs(sv_len)}

            for sample in all_samples:
                if sample not in rec.samples:
                    row_dict[sample] = 0
                    continue
                    
                fmt = rec.samples[sample]
                gt = fmt.get('GT')
                
                # Check if sample actually has the SV (0/1 or 1/1)
                has_variant = 1 if gt and (1 in gt) else 0
                row_dict[sample] = has_variant

                if has_variant:
                    if sample in previous_samples:
                        prev_has_sv = True
                    elif sample in new_samples:
                        new_samples_with_sv.append(sample)
                        total_svs_per_sample[sample] += 1

            sv_matrix_data.append(row_dict)

            # 5. The Core Logic: Is it strictly Novel?
            if len(new_samples_with_sv) > 0 and not prev_has_sv:
                novel_svs_detailed.append({
                    "SV_ID": sv_id,
                    "CHROM": rec.chrom,
                    "POS": rec.pos,
                    "SV_TYPE": sv_type,
                    "SV_LENGTH_BP": abs(sv_len),
                    "Found_In_Samples": ", ".join(new_samples_with_sv)
                })
                
                for s in new_samples_with_sv:
                    novel_svs_per_sample[s] += 1

    # 6. Export the Full Presence/Absence Matrix
    print("🧮 Generating Cohort SV Matrix...")
    df_all_svs = pd.DataFrame(sv_matrix_data)
    df_all_svs.to_csv(os.path.join(args.outdir, "datatable_all_SVs_presence.csv"), index=False)

    # 7. Export the Detailed Novel SV Table
    print("🔬 Isolating Novel Structural Variations...")
    df_novel = pd.DataFrame(novel_svs_detailed)
    if not df_novel.empty:
        df_novel.to_csv(os.path.join(args.outdir, "detailed_novel_SVs_list.csv"), index=False)
    else:
        print("⚠️ No strictly novel SVs found!")

    # 8. Create High-Level Summary
    print("📝 Generating Summary Report...")
    summary_data = []
    total_novel_svs = len(df_novel) if not df_novel.empty else 0
    
    for s in new_samples:
        if s in total_svs_per_sample:
            summary_data.append({
                "Sample": s,
                "Total_SVs_Found": total_svs_per_sample[s],
                "Strictly_Novel_SVs": novel_svs_per_sample[s]
            })
    
    df_summary = pd.DataFrame(summary_data)
    
    if not df_summary.empty:
        total_row = pd.DataFrame([{
            "Sample": "TOTAL",
            "Total_SVs_Found": df_summary["Total_SVs_Found"].sum(),
            "Strictly_Novel_SVs": df_summary["Strictly_Novel_SVs"].sum()
        }])
        df_summary = pd.concat([df_summary, total_row], ignore_index=True)
        df_summary.to_csv(os.path.join(args.outdir, "per_sample_novel_SV_counts.csv"), index=False)

    # Terminal Output
    print("\n✅ SV Pipeline Complete!")
    print(f"Total 'New' Samples Analyzed: {len(new_samples)}")
    print(f"Total Strictly Novel SVs Discovered: {total_novel_svs}")
    print(f"Outputs saved to: {args.outdir}/")

if __name__ == "__main__":
    main()
