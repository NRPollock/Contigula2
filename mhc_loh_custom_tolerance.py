import argparse
import os
import pysam
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

def create_target_bed(excel_path):
    """Reads the SNP list natively from Excel and extracts hg38_position."""
    df = pd.read_excel(excel_path)
    df["hg38_position"] = df["hg38_position"].astype(int)
    df = df.sort_values(by="hg38_position")
    return df["hg38_position"].tolist()

def filter_balanced_het(in_vcf, out_vcf, min_maf=0.30):
    """Parses a VCF and filters for heterozygous sites with balanced allele coverage."""
    balanced_positions = []
    
    with pysam.VariantFile(in_vcf) as vcf_in:
        header = vcf_in.header
        with pysam.VariantFile(out_vcf, 'w', header=header) as vcf_out:
            for rec in vcf_in:
                samp = rec.samples[0]
                ad = samp.get('AD')
                
                if ad and len(ad) >= 2:
                    ref_depth = ad[0]
                    alt_depth = ad[1]
                    total_depth = ref_depth + alt_depth
                    
                    if total_depth > 0:
                        minor_allele_ratio = min(ref_depth, alt_depth) / total_depth
                        if minor_allele_ratio >= min_maf:
                            vcf_out.write(rec)
                            balanced_positions.append(rec.pos)
                            
    return balanced_positions

def get_breaks(het_snps, window_size, hets_per_window):
    """
    Filters out 'tolerated' heterozygous SNPs basedthresholds.
    If the number of SNPs in a given window exceeds the tolerated amount, 
    those SNPs are classified as breaks in homozygosity.
    """
    # If 0 hets are tolerated, every single het is a break
    if hets_per_window == 0:
        return het_snps

    breaks = set()
    # To violate the tolerance, we need (tolerated + 1) SNPs in the window
    trigger_count = hets_per_window + 1

    if len(het_snps) < trigger_count:
        return []

    # Slide a window of size 'trigger_count' across the SNPs
    for i in range(len(het_snps) - (trigger_count - 1)):
        # If the distance between the first and last SNP in this group is within the window_size
        if het_snps[i + trigger_count - 1] - het_snps[i] <= window_size:
            # Add all SNPs in this violating group to the breaks list
            for j in range(trigger_count):
                breaks.add(het_snps[i + j])
            
    return sorted(list(breaks))

def calculate_loh(all_snps, het_snps):
    """Finds the longest stretch of homozygosity between breaks."""
    stretches = []
    het_indices = [all_snps.index(h) for h in het_snps if h in all_snps]
    
    if not het_indices:
        if len(all_snps) > 0:
            return all_snps[0], all_snps[-1]
        return None, None
        
    if het_indices[0] > 0:
        stretches.append((all_snps[0], all_snps[het_indices[0] - 1]))
        
    for i in range(len(het_indices) - 1):
        idx1 = het_indices[i]
        idx2 = het_indices[i+1]
        if idx2 - idx1 > 1:
            start_snp = all_snps[idx1 + 1]
            end_snp = all_snps[idx2 - 1]
            stretches.append((start_snp, end_snp))
            
    if het_indices[-1] < len(all_snps) - 1:
        stretches.append((all_snps[het_indices[-1] + 1], all_snps[-1]))
        
    if not stretches:
        return None, None
        
    stretches.sort(key=lambda x: x[1] - x[0], reverse=True)
    longest = stretches[0]
    
    return longest[0], longest[1]

def process_sample(sample_info):
    s_name, all_snps, out_dir, window_size, hets_per_window = sample_info
    
    f4 = os.path.join(out_dir, f"{s_name}_het.vcf.gz")
    f5 = os.path.join(out_dir, f"{s_name}_balanced_het.vcf")
    
    try:
        balanced_snps = []
        
        # 1. Grab the balanced SNPs from existing files
        if os.path.exists(f5):
            with pysam.VariantFile(f5) as vcf_in:
                for rec in vcf_in:
                    balanced_snps.append(rec.pos)
        elif os.path.exists(f4):
            print(f"[{s_name}] Generating balanced het VCF...")
            balanced_snps = filter_balanced_het(f4, f5, min_maf=0.30)
        else:
            print(f"❌ SKIPPING {s_name}: No VCF files found in {out_dir}.")
            return {"Sample": s_name, "LOH_Start_Position": "ERROR", "LOH_End_Position": "Missing VCF", "LOH_Length": "N/A"}
            
        # 2. Apply the Custom Tolerance Filter
        final_breaks = get_breaks(balanced_snps, window_size, hets_per_window)
        
        # 3. Calculate the LOH Tract
        start_loh, end_loh = calculate_loh(all_snps, final_breaks)
        loh_length = (end_loh - start_loh) if (start_loh and end_loh) else 0
        
        print(f"✅ Analyzed {s_name} (Ignored {len(balanced_snps) - len(final_breaks)} tolerated SNPs)")
        return {"Sample": s_name, "LOH_Start_Position": start_loh, "LOH_End_Position": end_loh, "LOH_Length": loh_length}
        
    except Exception as e:
        print(f"❌ FAILED {s_name}: {e}")
        return {"Sample": s_name, "LOH_Start_Position": "ERROR", "LOH_End_Position": "Processing Error", "LOH_Length": "N/A"}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate LOH with Custom Heterozygous Tolerance")
    parser.add_argument("--samples", required=True, help="TXT file with one sample per line")
    parser.add_argument("--snps", required=True, help="Path to MEGA_MHC_SNPs.xlsx")
    parser.add_argument("--out_dir", default="vcf_output", help="Directory where the existing VCFs are stored")
    parser.add_argument("--window_size", type=int, default=500000, help="Window size in base pairs (default: 500000)")
    parser.add_argument("--hets_per_window", type=int, default=2, help="Number of heterozygous SNPs tolerated per window (default: 2)")
    parser.add_argument("--jobs", type=int, default=4, help="Number of samples to process in parallel")
    args = parser.parse_args()

    # Load target SNPs from Excel
    all_target_positions = create_target_bed(args.snps)
    
    with open(args.samples, 'r') as f:
        sample_names = [line.strip() for line in f if line.strip()]
        
    # Pass the new arguments into the task
    tasks = [(s, all_target_positions, args.out_dir, args.window_size, args.hets_per_window) for s in sample_names]
        
    summary_file = os.path.join(args.out_dir, f"LOH_Tol_{args.hets_per_window}per{args.window_size//1000}kb_Summary.csv")
    
    with open(summary_file, "w") as f:
        f.write("Sample,LOH_Start_Position,LOH_End_Position,LOH_Length\n")
        
    print(f"Starting LOH analysis (Tolerating {args.hets_per_window} hets per {args.window_size} bp window)...")
    
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(process_sample, task): task[0] for task in tasks}
        
        for future in as_completed(futures):
            s_name = futures[future]
            try:
                res = future.result()
                with open(summary_file, "a") as f:
                    f.write(f"{res['Sample']},{res['LOH_Start_Position']},{res['LOH_End_Position']},{res['LOH_Length']}\n")
            except Exception as exc:
                with open(summary_file, "a") as f:
                    f.write(f"{s_name},ERROR,Fatal Exception,N/A\n")
                    
    print(f"\n Analysis Complete! Summary saved to {summary_file}")
