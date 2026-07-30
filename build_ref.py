#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


# ----------------------------
# Helpers
# ----------------------------
def run(cmd: List[str], *, check: bool = True) -> None:
    print("[CMD]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=check)

def ensure_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"[ERROR] Required tool not found in PATH: {name}")

def sanitize_header_for_fasta(h: str) -> str:
    return re.sub(r"\s+", "_", h.strip()) or "custom_ref"

def fasta_to_dict(path: str) -> Dict[str, str]:
    """
    Read a (multi-)FASTA into a dict: {header_first_token: sequence}
    Use the first token of the header as the contig name to match minimap2/samtools rname.
    """
    seqs: Dict[str, str] = {}
    header: Optional[str] = None
    chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    name = header.split()[0]
                    if name in seqs:
                        raise ValueError(f"Duplicate FASTA record name (first token) '{name}' in {path}")
                    seqs[name] = "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

        if header is not None:
            name = header.split()[0]
            if name in seqs:
                raise ValueError(f"Duplicate FASTA record name (first token) '{name}' in {path}")
            seqs[name] = "".join(chunks)

    if not seqs:
        raise ValueError(f"No FASTA records found in {path}")

    return seqs

def write_fasta(path: str, header: str, seq: str, wrap: int = 60) -> None:
    with open(path, "w", encoding="utf-8") as out:
        out.write(f">{header}\n")
        if wrap and wrap > 0:
            for i in range(0, len(seq), wrap):
                out.write(seq[i:i + wrap] + "\n")
        else:
            out.write(seq + "\n")


# ----------------------------
# Coverage parsing / stats
# ----------------------------
@dataclass
class CoverageRow:
    rname: str
    length: int
    numreads: int
    covbases: int
    coverage_pct: float
    meandepth: float
    meanmapq: float

@dataclass
class DepthStats:
    rname: str
    length: int
    zero_bases: int
    ltK_bases: int
    geK_bases: int
    mean_depth: float
    breadth_geK_pct: float

@dataclass
class ScoredContig:
    rname: str
    length: int
    numreads: int
    covbases: int
    zero_bases: int
    ltK_bases: int
    geK_bases: int
    coverage_pct: float
    meandepth_cov: float
    mean_depth_depthtool: float
    meanmapq: float
    breadth_geK_pct: float


def parse_samtools_coverage(tsv_text: str) -> Dict[str, CoverageRow]:
    """
    Parse `samtools coverage` output (tab-separated).
    Typical columns:
      rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
    Return dict keyed by rname.
    """
    rows: Dict[str, CoverageRow] = {}
    for line in tsv_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            raise ValueError(f"Unexpected samtools coverage line (need >=9 fields): {line}")

        rname = parts[0]
        startpos = int(parts[1])
        endpos = int(parts[2])
        length = endpos - startpos + 1
        numreads = int(parts[3])
        covbases = int(parts[4])
        coverage_pct = float(parts[5])
        meandepth = float(parts[6])
        meanmapq = float(parts[8])

        rows[rname] = CoverageRow(
            rname=rname,
            length=length,
            numreads=numreads,
            covbases=covbases,
            coverage_pct=coverage_pct,
            meandepth=meandepth,
            meanmapq=meanmapq,
        )
    if not rows:
        raise ValueError("No rows parsed from samtools coverage output.")
    return rows


def compute_depth_stats(bam: str, min_mapq: int, min_baseq: int, K: int) -> Dict[str, DepthStats]:
    """
    Compute per-contig depth-based stats using:
      samtools depth -aa -q <min_mapq> -Q <min_baseq> <bam>

    This includes ALL positions (including those with 0 depth),
    directly count:
      - zero_bases
      - bases with depth < K
      - bases with depth >= K
      - mean depth
      - breadth(depth>=K)
    """
    cmd = [
        "samtools", "depth",
        "-aa",
        "-q", str(min_mapq),
        "-Q", str(min_baseq),
        "-d", "0",
        bam
    ]
    print("[CMD]", " ".join(cmd), flush=True)

    # streaming parse
    counts_len: Dict[str, int] = {}
    counts_zero: Dict[str, int] = {}
    counts_ltK: Dict[str, int] = {}
    counts_geK: Dict[str, int] = {}
    sums_depth: Dict[str, int] = {}

    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as p:
        assert p.stdout is not None
        for line in p.stdout:
            line = line.rstrip("\n")
            if not line:
                continue
            rname, _pos, depth_s = line.split("\t")
            d = int(depth_s)

            counts_len[rname] = counts_len.get(rname, 0) + 1
            sums_depth[rname] = sums_depth.get(rname, 0) + d

            if d == 0:
                counts_zero[rname] = counts_zero.get(rname, 0) + 1
            if d < K:
                counts_ltK[rname] = counts_ltK.get(rname, 0) + 1
            else:
                counts_geK[rname] = counts_geK.get(rname, 0) + 1

        rc = p.wait()
        if rc != 0:
            raise SystemExit(f"[ERROR] samtools depth failed with exit code {rc}")

    stats: Dict[str, DepthStats] = {}
    for rname, length in counts_len.items():
        zero = counts_zero.get(rname, 0)
        ltK = counts_ltK.get(rname, 0)
        geK = counts_geK.get(rname, 0)
        mean_depth = (sums_depth.get(rname, 0) / length) if length else 0.0
        breadth_geK_pct = (geK / length * 100.0) if length else 0.0

        stats[rname] = DepthStats(
            rname=rname,
            length=length,
            zero_bases=zero,
            ltK_bases=ltK,
            geK_bases=geK,
            mean_depth=mean_depth,
            breadth_geK_pct=breadth_geK_pct,
        )

    if not stats:
        raise ValueError("No depth stats computed (empty output from samtools depth?).")

    return stats


def coverage_for_reference(
    ref_fa: str,
    reads: str,
    out_prefix: str,
    threads: int,
    outdir: str,
    min_mapq: int,
    min_baseq: int,
    K: int
) -> Tuple[str, List[ScoredContig]]:
    """
    Align reads to ref_fa -> sorted/indexed BAM -> samtools coverage (filtered) + samtools depth (filtered)
    Return bam_path and list of per-contig scored stats.
    """
    bam = os.path.join(outdir, f"{out_prefix}.sorted.bam")

    # minimap2 -> samtools sort
    cmd = [
        "bash", "-lc",
        f"minimap2 -t {threads} -ax map-ont {ref_fa} {reads} | "
        f"samtools sort -@ {max(1, threads-1)} -o {bam} -"
    ]
    run(cmd)
    run(["samtools", "index", "-@", str(max(1, threads-1)), bam])

    # samtools coverage with MAPQ/baseQ filters for meanmapq/meandepth/covbases
    cov_cmd = ["samtools", "coverage", "-H", "-q", str(min_mapq), "-Q", str(min_baseq), bam]
    print("[CMD]", " ".join(cov_cmd), flush=True)
    cov_out = subprocess.check_output(cov_cmd, text=True)
    cov = parse_samtools_coverage(cov_out)

    # depth stats (also filtered) to get breadth>=K and ltK counts robustly
    depth_stats = compute_depth_stats(bam=bam, min_mapq=min_mapq, min_baseq=min_baseq, K=K)

    # merge
    scored: List[ScoredContig] = []
    for rname, ds in depth_stats.items():
        # if a contig has absolutely no alignments passing filters, samtools coverage may omit it
        cr = cov.get(rname, CoverageRow(rname=rname, length=ds.length, numreads=0, covbases=ds.length - ds.zero_bases,
                                       coverage_pct=(100.0*(ds.length - ds.zero_bases)/ds.length if ds.length else 0.0),
                                       meandepth=0.0, meanmapq=0.0))
        scored.append(ScoredContig(
            rname=rname,
            length=ds.length,
            numreads=cr.numreads,
            covbases=cr.covbases,
            zero_bases=ds.zero_bases,
            ltK_bases=ds.ltK_bases,
            geK_bases=ds.geK_bases,
            coverage_pct=cr.coverage_pct,
            meandepth_cov=cr.meandepth,
            mean_depth_depthtool=ds.mean_depth,
            meanmapq=cr.meanmapq,
            breadth_geK_pct=ds.breadth_geK_pct,
        ))

    return bam, scored


def pick_best(scored: List[ScoredContig]) -> ScoredContig:
    """
    Selection (MAPQ/baseQ filtered):
      1) maximize breadth_geK_pct
      2) maximize meanmapq
      3) maximize mean_depth_depthtool (robust mean across all positions)
      4) minimize zero_bases
      5) then stable by name
    """
    return sorted(
        scored,
        key=lambda r: (-r.breadth_geK_pct, -r.meanmapq, -r.mean_depth_depthtool, r.zero_bases, r.rname)
    )[0]


def write_report(scored: List[ScoredContig], path: str, K: int, min_mapq: int, min_baseq: int) -> None:
    with open(path, "w", encoding="utf-8") as out:
        out.write(f"# Filters: min_mapq={min_mapq}, min_baseq={min_baseq}\n")
        out.write(f"# Depth threshold K={K} (breadth>=K used for selection)\n")
        out.write("\t".join([
            "rname",
            "length",
            "numreads",
            "covbases",
            "zero_bases",
            f"lt{K}_bases",
            f"ge{K}_bases",
            "coverage_pct",
            "breadth_geK_pct",
            "mean_depth_allpos",
            "meandepth_cov",
            "meanmapq"
        ]) + "\n")
        for r in scored:
            out.write("\t".join([
                r.rname,
                str(r.length),
                str(r.numreads),
                str(r.covbases),
                str(r.zero_bases),
                str(r.ltK_bases),
                str(r.geK_bases),
                f"{r.coverage_pct:.6f}",
                f"{r.breadth_geK_pct:.6f}",
                f"{r.mean_depth_depthtool:.6f}",
                f"{r.meandepth_cov:.6f}",
                f"{r.meanmapq:.2f}",
            ]) + "\n")


# ----------------------------
# Main
# ----------------------------
def main() -> None:
    ensure_tool("minimap2")
    ensure_tool("samtools")
    ensure_tool("bash")

    ap = argparse.ArgumentParser(
        description="Pick best-covered Class I and Class II contigs (breadth at depth>=K after MAPQ filtering) and splice with Ns."
    )
    ap.add_argument("--classI_fasta", required=True, help="Multi-FASTA of Class I candidates.")
    ap.add_argument("--classII_fasta", required=True, help="Multi-FASTA of Class II candidates.")
    ap.add_argument("--reads", required=True, help="ONT reads (FASTQ/FASTA; .gz ok).")
    ap.add_argument("--sample", required=True, help="Sample name used in output FASTA header.")
    ap.add_argument(
        "--outname",
        help="Basename for output combined FASTA (no extension). Default: <sample>.custom_ref"
    )
    ap.add_argument("-o", "--outdir", default="best_ref_out", help="Output directory (default: best_ref_out).")
    ap.add_argument("-t", "--threads", type=int, default=16, help="Threads for minimap2 (default: 16).")
    ap.add_argument("--wrap", type=int, default=0, help="FASTA wrap length (default: 0; 0 = no wrap).")
    ap.add_argument("--gapN", type=int, default=20, help="Number of Ns between Class I and Class II (default: 20).")

    # New scoring knobs
    ap.add_argument("--min_mapq", type=int, default=10, help="Min MAPQ for counting coverage/depth (default: 10).")
    ap.add_argument("--min_baseq", type=int, default=0, help="Min base quality for counting depth (default: 0).")
    ap.add_argument("--K", type=int, default=3, help="Depth threshold K for breadth(depth>=K) (default: 3).")

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Map + score Class I
    _, scoredI = coverage_for_reference(
        ref_fa=args.classI_fasta,
        reads=args.reads,
        out_prefix=f"{args.sample}.classI",
        threads=args.threads,
        outdir=args.outdir,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        K=args.K
    )
    reportI = os.path.join(args.outdir, f"{args.sample}.classI.coverage.tsv")
    write_report(scoredI, reportI, K=args.K, min_mapq=args.min_mapq, min_baseq=args.min_baseq)

    # Map + score Class II
    _, scoredII = coverage_for_reference(
        ref_fa=args.classII_fasta,
        reads=args.reads,
        out_prefix=f"{args.sample}.classII",
        threads=args.threads,
        outdir=args.outdir,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        K=args.K
    )
    reportII = os.path.join(args.outdir, f"{args.sample}.classII.coverage.tsv")
    write_report(scoredII, reportII, K=args.K, min_mapq=args.min_mapq, min_baseq=args.min_baseq)

    bestI = pick_best(scoredI)
    bestII = pick_best(scoredII)

    # Extract sequences
    seqsI = fasta_to_dict(args.classI_fasta)
    seqsII = fasta_to_dict(args.classII_fasta)

    if bestI.rname not in seqsI:
        raise SystemExit(f"[ERROR] Best Class I contig '{bestI.rname}' not found in {args.classI_fasta}")
    if bestII.rname not in seqsII:
        raise SystemExit(f"[ERROR] Best Class II contig '{bestII.rname}' not found in {args.classII_fasta}")

    combined = seqsI[bestI.rname] + ("N" * args.gapN) + seqsII[bestII.rname]

    # Output file naming
    base_name = sanitize_header_for_fasta(args.outname) if args.outname else sanitize_header_for_fasta(args.sample) + ".custom_ref"
    out_fa = os.path.join(args.outdir, f"{base_name}.fasta")

    out_header = (
        sanitize_header_for_fasta(args.sample)
        + f"|classI={bestI.rname}|classII={bestII.rname}"
        + f"|gapN={args.gapN}|K={args.K}|minMAPQ={args.min_mapq}|minBQ={args.min_baseq}"
    )
    write_fasta(out_fa, out_header, combined, wrap=args.wrap)

    # Summary
    summary_path = os.path.join(args.outdir, f"{args.sample}.chosen.txt")
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write(f"Sample\t{args.sample}\n")
        out.write(f"Filters\tmin_mapq={args.min_mapq}\tmin_baseq={args.min_baseq}\tK={args.K}\n")
        out.write(
            "Chosen_ClassI\t{r}\tbreadth_geK_pct={b:.6f}\tltK_bases={l}\tzero_bases={z}\t"
            "mean_depth_allpos={md:.6f}\tmeanmapq={mq:.2f}\n".format(
                r=bestI.rname, b=bestI.breadth_geK_pct, l=bestI.ltK_bases, z=bestI.zero_bases,
                md=bestI.mean_depth_depthtool, mq=bestI.meanmapq
            )
        )
        out.write(
            "Chosen_ClassII\t{r}\tbreadth_geK_pct={b:.6f}\tltK_bases={l}\tzero_bases={z}\t"
            "mean_depth_allpos={md:.6f}\tmeanmapq={mq:.2f}\n".format(
                r=bestII.rname, b=bestII.breadth_geK_pct, l=bestII.ltK_bases, z=bestII.zero_bases,
                md=bestII.mean_depth_depthtool, mq=bestII.meanmapq
            )
        )
        out.write(f"Output_FASTA\t{out_fa}\n")
        out.write(f"ClassI_Report\t{reportI}\n")
        out.write(f"ClassII_Report\t{reportII}\n")

    print("\n[DONE] Selection criteria (MAPQ/baseQ filtered):")
    print(f"  Primary : maximize breadth(depth >= {args.K})")
    print("  Secondary: maximize mean MAPQ")
    print("  Tertiary : maximize mean depth (all positions)")
    print("  Then     : minimize zero-coverage bases\n")

    print("[DONE] Chosen candidates:")
    print(f"  Class I : {bestI.rname}  breadth>=K={bestI.breadth_geK_pct:.3f}%  ltK={bestI.ltK_bases}  zero={bestI.zero_bases}  meanDepth={bestI.mean_depth_depthtool:.3f}  meanMAPQ={bestI.meanmapq:.2f}")
    print(f"  Class II: {bestII.rname}  breadth>=K={bestII.breadth_geK_pct:.3f}%  ltK={bestII.ltK_bases}  zero={bestII.zero_bases}  meanDepth={bestII.mean_depth_depthtool:.3f}  meanMAPQ={bestII.meanmapq:.2f}")

    print("\n[DONE] Wrote custom reference:", out_fa)
    print("[DONE] Reports:")
    print("  ", reportI)
    print("  ", reportII)
    print("  ", summary_path)


if __name__ == "__main__":
    main()

