#!/usr/bin/env bash
set -o pipefail

# =============================================================================
# End-to-end pipeline: filter_reads + build_ref/RCCX replace + assembly
# + filter out SPAdes contigs shorter than --min-contig-length
# + (optional) include Flye contigs if Flye finishes
# + ragtag patch + pilon
#
# USAGE:
#   bash full_assembly.sh SAMPLE THREADS [WKDIR] [--wkdir PATH] [--start-step N] [--min-contig-length N] [--rerun]
#
# STEPS:
#   0 = Filter reads to MHC
#   1 = RCCX mapping + coverage check
#   2 = build_ref + replace RCCX region
#   3 = Assemble (SPAdes required; Flye optional)
#   4 = Filter contigs (minlen) + combine (SPAdes required; Flye optional)
#   5 = RagTag patch
#   6 = Pilon polishing
#
# SKIP BEHAVIOR:
#   If "${wkdir}/${name}_final.fasta" exists and is non-empty, the script exits
#   immediately with code 0 unless --rerun is provided.
#   Intermediate steps are also skipped if their expected outputs already exist
#   and are non-empty, unless --rerun is provided.
# =============================================================================

# -------------------- CLI --------------------
name="${1:?sample name required}"
threads="${2:?threads required}"
wkdir="${3:-}"

START_STEP=0
MIN_CONTIG_LENGTH=3000
RERUN=0

shift 2
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-step|--start_step)
      START_STEP="${2:?missing number for --start-step}"
      shift 2
      ;;
    --wkdir)
      wkdir="${2:?missing path for --wkdir}"
      shift 2
      ;;
    --min-contig-length|--min_contig_length)
      MIN_CONTIG_LENGTH="${2:?missing number for --min-contig-length}"
      shift 2
      ;;
    --rerun)
      RERUN=1
      shift
      ;;
    --help|-h)
      echo "Usage: $0 SAMPLE THREADS [WKDIR] [--wkdir PATH] [--start-step N] [--min-contig-length N] [--rerun]"
      echo "Steps: 0,1,2,3,4,5,6"
      echo "Note: SPAdes is required. Flye is optional; pipeline continues if Flye output is missing."
      exit 0
      ;;
    *)
      if [[ -z "${wkdir}" ]]; then
        wkdir="$1"
        shift
      else
        echo "WARNING: ignoring unknown argument: $1" >&2
        shift
      fi
      ;;
  esac
done

# -------------------- ONE-STOP CONFIG --------------------
BS_DIR="/home/nick/asstest/merged"

DEFAULT_WKDIR="${BS_DIR}/${name}_assembly"
wkdir="${wkdir:-$DEFAULT_WKDIR}"
mkdir -p "${wkdir}"

LR_FASTQ="/home/nick/asstest/merged/long_read/${name}.fastq"

SR_R1_GZ=( "/home/nick/asstest/merged/short_read/${name}_R1.fastq.gz" )
SR_R2_GZ=( "/home/nick/asstest/merged/short_read/${name}_R2.fastq.gz" )

for f in "${SR_R1_GZ[@]}" "${SR_R2_GZ[@]}"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing SR FASTQ: $f" >&2
    exit 2
  fi
done

if [[ ! -s "${LR_FASTQ}" ]]; then
  echo "ERROR: missing LR FASTQ: ${LR_FASTQ}" >&2
  exit 2
fi

if (( ${#SR_R1_GZ[@]} != ${#SR_R2_GZ[@]} )); then
  echo "ERROR: Unequal number of R1 and R2 files for sample '${name}'" >&2
  echo "  R1 (${#SR_R1_GZ[@]}): ${SR_R1_GZ[*]}" >&2
  echo "  R2 (${#SR_R2_GZ[@]}): ${SR_R2_GZ[*]}" >&2
  exit 2
fi

IFS=$'\n' SR_R1_GZ=( $(printf "%s\n" "${SR_R1_GZ[@]}" | sort) )
IFS=$'\n' SR_R2_GZ=( $(printf "%s\n" "${SR_R2_GZ[@]}" | sort) )
unset IFS

CUSTOM_HG38="/home/nick/asstest/refs/chimaric/refwalker20241209/bashgenerator/hg38p14.customMHC.fa"
FILTER_BED="${BS_DIR}/filter_padded.bed"
RCCX_REFS="/home/nick/asstest/refs/RCCX/RCCX_refs.fasta"

CLASSI_FASTA="${BS_DIR}/refs/split/MHC_classI.fasta"
CLASSII_FASTA="${BS_DIR}/refs/split/MHC_classII.fasta"

PILON_JAR="/home/nick/micromamba/envs/contigula/share/pilon-1.24-0/pilon.jar"

TRUSTED_CONTIGS="/home/nick/asstest/trusted_contigs.fa"
# TRUSTED_CONTIGS=""

# -------------------- Pilon settings --------------------
MAX_PILON_ROUNDS=1
PILON_MINDEPTH=15
PILON_TRACKS=1

ENABLE_SR_MAPQ_FILTER=1
SR_MIN_MAPQ=5

ENABLE_LR_SOFTCLIP=0
DISABLE_LR_SECONDARY=1
JAVA_MEM="250G"

log="${wkdir}/${name}.log"

# Canonical pipeline file names
SR1="${wkdir}/${name}_mhc_sr_R1.fastq"
SR2="${wkdir}/${name}_mhc_sr_R2.fastq"
LR="${wkdir}/${name}_mhc_longreads.fastq"

FINAL_OUT="${wkdir}/${name}_final.fasta"
BEST_RCCX="${wkdir}/${name}_best_RCCX.fasta"
SPLIT_RCCX_MHC="${wkdir}/${name}_split_RCCX_MHC.fasta"
SPADES_OUT="${wkdir}/${name}_spades"
FLYE_OUT="${wkdir}/${name}_flye"
SPADES_CTG="${SPADES_OUT}/contigs.fasta"
FLYE_CTG="${FLYE_OUT}/assembly.fasta"
SPADES_FILT="${wkdir}/${name}_spades.contigs.min${MIN_CONTIG_LENGTH}.fasta"
FLYE_FILT="${wkdir}/${name}_flye.assembly.min${MIN_CONTIG_LENGTH}.fasta"
COMBINED="${wkdir}/${name}_combined.min${MIN_CONTIG_LENGTH}.fa"
RAGTAG_OUT="${wkdir}/${name}_ragtag_patch"
PATCHED="${RAGTAG_OUT}/ragtag.patch.fasta"
REF="${SPLIT_RCCX_MHC}"

# build_ref.py output directory depends on sample name, keep explicit here
BEST_REF_DIR="${wkdir}/${name}_best_ref"
BEST_REF_FASTA="${BEST_REF_DIR}/${name}_split_MHC.fasta"

ts() { date +"%Y-%m-%d %H:%M:%S"; }
say() { echo "[$(ts)] $*" | tee -a "${log}"; }
require_file() { [[ -s "$1" ]] || { echo "ERROR: missing file: $1" >&2; exit 2; }; }

step_outputs_exist() {
  [[ "${RERUN}" -eq 1 ]] && return 1
  local f
  for f in "$@"; do
    [[ -s "${f}" ]] || return 1
  done
  return 0
}

# -------------------- Early exit if final output already exists --------------------
if [[ "${RERUN}" -ne 1 && -s "${FINAL_OUT}" ]]; then
  echo "[$(ts)] SKIP: final output already exists and is non-empty: ${FINAL_OUT}"
  echo "[$(ts)] SKIP: use --rerun to force reprocessing" | tee -a "${log}" >/dev/null
  exit 0
fi

start_time=$(date +%s)
elapsed_time() {
  local now=$(date +%s)
  local elapsed=$((now - start_time))
  echo "Elapsed Time: $((elapsed/60))m $((elapsed%60))s"
}

say "START_STEP=${START_STEP}  wkdir=${wkdir}  sample=${name}  threads=${threads}"
say "MIN_CONTIG_LENGTH=${MIN_CONTIG_LENGTH}  (contigs < this will be dropped)"
say "RERUN=${RERUN}"
say "Flye is OPTIONAL; SPAdes is REQUIRED"

# =============================================================================
# 0) Filter reads to MHC
# =============================================================================
if (( START_STEP <= 0 )); then
  say "Step 0: Filter reads to MHC (LR + SR)  |  ${name}"
  elapsed_time

  require_file "${LR_FASTQ}"
  for f in "${SR_R1_GZ[@]}" "${SR_R2_GZ[@]}"; do
    require_file "$f"
  done
  require_file "${CUSTOM_HG38}"
  require_file "${FILTER_BED}"

  # -------------------- LR branch --------------------
  if [[ "${RERUN}" -ne 1 && -s "${LR}" ]]; then
    say "  0A-0C: LR output already present; skipping minimap2 LR filtering"
  else
    say "  0A: map long reads to custom hg38"
    minimap2 -ax map-ont -t "${threads}" \
      "${CUSTOM_HG38}" "${LR_FASTQ}" \
      | samtools sort -@ "${threads}" \
      | samtools view -Sb -@ "${threads}" -h \
      > "${wkdir}/${name}_lr_on_custom_hg38.bam"
    samtools index -@ "${threads}" "${wkdir}/${name}_lr_on_custom_hg38.bam"

    say "  0B: filter long reads to MHC bed"
    samtools view -@ "${threads}" -S -h -L "${FILTER_BED}" \
      "${wkdir}/${name}_lr_on_custom_hg38.bam" \
      > "${wkdir}/${name}_mhc_lr.sam"

    samtools view -@ "${threads}" -Sb -h "${wkdir}/${name}_mhc_lr.sam" > "${wkdir}/${name}_mhc_lr.bam"
    samtools index -@ "${threads}" "${wkdir}/${name}_mhc_lr.bam"

    say "  0C: LR bam -> fastq"
    PicardCommandLine SamToFastq -F "${LR}" -I "${wkdir}/${name}_mhc_lr.bam"

    require_file "${LR}"
  fi

  # -------------------- SR branch --------------------
  if [[ "${RERUN}" -ne 1 && -s "${SR1}" && -s "${SR2}" ]]; then
    say "  0D-0F: SR outputs already present; skipping bwa mem SR filtering"
  else
    say "  0D: map short reads to custom hg38"
    bwa mem -t "${threads}" "${CUSTOM_HG38}" "${SR_R1_GZ[@]}" "${SR_R2_GZ[@]}" \
      | samtools sort -@ "${threads}" \
      | samtools view -Sb -@ "${threads}" -h \
      > "${wkdir}/${name}_sr_on_custom_hg38.bam"
    samtools index -@ "${threads}" "${wkdir}/${name}_sr_on_custom_hg38.bam"

    say "  0E: filter short reads to MHC bed"
    samtools view -@ "${threads}" -Sb -h -L "${FILTER_BED}" \
      "${wkdir}/${name}_sr_on_custom_hg38.bam" \
      > "${wkdir}/${name}_mhc_sr.bam"
    samtools index -@ "${threads}" "${wkdir}/${name}_mhc_sr.bam"

    say "  0F: SR bam -> fastq pairs (LENIENT validation)"
    PicardCommandLine SamToFastq \
      -F  "${SR1}" \
      -F2 "${SR2}" \
      -I  "${wkdir}/${name}_mhc_sr.bam" \
      --VALIDATION_STRINGENCY LENIENT

    require_file "${SR1}"
    require_file "${SR2}"
  fi
fi

# =============================================================================
# 1) RCCX depth check + best RCCX selection
# =============================================================================
if (( START_STEP <= 1 )); then
  if step_outputs_exist "${BEST_RCCX}"; then
    say "Step 1: output already present; skipping"
  else
    say "Step 1: RCCX mapping + coverage check"
    elapsed_time

    require_file "${RCCX_REFS}"
    require_file "${LR}"

    minimap2 -t "${threads}" --secondary=no -a -x map-ont \
      "${RCCX_REFS}" "${LR}" \
      | samtools sort -@ "${threads}" \
      | samtools view -Sb -@ "${threads}" \
      > "${wkdir}/${name}_RCCX_lr.bam"

    samtools depth -a -o "${wkdir}/${name}_RCCX_depths_lr.txt" "${wkdir}/${name}_RCCX_lr.bam"

    python "${BS_DIR}/coverage_check.py" \
      "${wkdir}/${name}_RCCX_depths_lr.txt" \
      "${wkdir}/${name}_coverage_zeros.txt" \
      "${BEST_RCCX}"

    require_file "${BEST_RCCX}"
  fi
fi

# =============================================================================
# 2) Build best ref (class I/II selection) + replace RCCX region
# =============================================================================
if (( START_STEP <= 2 )); then
  if step_outputs_exist "${SPLIT_RCCX_MHC}"; then
    say "Step 2: output already present; skipping"
    require_file "${SPLIT_RCCX_MHC}"
  else
    say "Step 2: build_ref.py + replace_region.py"
    elapsed_time

    require_file "${CLASSI_FASTA}"
    require_file "${CLASSII_FASTA}"
    require_file "${LR}"
    require_file "${BEST_RCCX}"

    python "${BS_DIR}/build_ref.py" \
      --classI_fasta "${CLASSI_FASTA}" \
      --classII_fasta "${CLASSII_FASTA}" \
      --reads "${LR}" \
      --sample "${name}" \
      --outname "${name}_split_MHC" \
      -o "${BEST_REF_DIR}" \
      -t "${threads}" \
      --min_mapq 10 \
      --K 3

    require_file "${BEST_REF_FASTA}"

    python "${BS_DIR}/replace_region.py" \
      -t "${BEST_REF_FASTA}" \
      -i "${BEST_RCCX}" \
      -5 "GTCTGACACAAGCATTAGTGAGATGCTCCCCTCGAAGAATAGTCTTGTTTCTTCTAAGGACTGATTCTCACCCCGGCTTTGGCTCTCCTAATTTTAGAGG" \
      -3 "GCGGCGTCTCAGGGCAGGACAGGGAAGTCTCCCTCACTTGTCCCCTGCAACAGGGGCTGAGCCACAACCGACTGTGGATCTCGGCAGCGACAGTGAGGAG" \
      -o "${SPLIT_RCCX_MHC}"

    require_file "${SPLIT_RCCX_MHC}"
  fi
fi

# =============================================================================
# 3) Assemble: SPAdes required; Flye optional
# =============================================================================
if (( START_STEP <= 3 )); then
  if step_outputs_exist "${SPADES_CTG}"; then
    say "Step 3: SPAdes output already present; skipping"
    if [[ -s "${FLYE_CTG}" ]]; then
      say "  Flye assembly found: ${FLYE_CTG}"
    else
      say "  WARNING: Flye output missing/empty: ${FLYE_CTG}. Step 4 will continue without Flye contigs."
    fi
  else
    say "Step 3: SPAdes + (optional) Flye assembly"
    elapsed_time

    require_file "${SR1}"
    require_file "${SR2}"
    require_file "${LR}"

    mkdir -p "${SPADES_OUT}"
    mkdir -p "${FLYE_OUT}"

    spades_args=( -o "${SPADES_OUT}" -t "${threads}" --cov-cutoff auto -1 "${SR1}" -2 "${SR2}" --nanopore "${LR}" )
    if [[ -n "${TRUSTED_CONTIGS}" && -s "${TRUSTED_CONTIGS}" ]]; then
      spades_args+=( --trusted-contigs "${TRUSTED_CONTIGS}" )
    fi

    spades "${spades_args[@]}" 2>&1 | tee -a "${log}"

    require_file "${SPADES_CTG}"

    # Flye is best-effort: do not fail pipeline if it errors or produces no assembly.fasta
    say "  Running Flye (best-effort; pipeline continues if Flye fails)"
    set +e
    flye --nano-hq "${LR}" -o "${FLYE_OUT}" -t "${threads}" -i 3 2>&1 | tee -a "${log}"
    flye_rc=${PIPESTATUS[0]}
    set -e

    if [[ "${flye_rc}" -ne 0 ]]; then
      say "  WARNING: Flye exited non-zero (rc=${flye_rc}). Continuing without Flye contigs."
    fi
    if [[ ! -s "${FLYE_CTG}" ]]; then
      say "  WARNING: Flye output missing/empty: ${FLYE_CTG}. Continuing without Flye contigs."
    else
      say "  Flye assembly found: ${FLYE_CTG}"
    fi
  fi
fi

# =============================================================================
# 4) Filter contigs (SPAdes required; Flye optional) + combine
# =============================================================================
if (( START_STEP <= 4 )); then
  if step_outputs_exist "${COMBINED}"; then
    say "Step 4: combined contigs already present; skipping"
  else
    say "Step 4: Filter contigs (minlen=${MIN_CONTIG_LENGTH}) + combine"
    elapsed_time

    # SPAdes must exist
    require_file "${SPADES_CTG}"

    say "  Filtering SPAdes contigs to >= ${MIN_CONTIG_LENGTH} bp"
    python - <<'PY' "${SPADES_CTG}" "${SPADES_FILT}" "${MIN_CONTIG_LENGTH}" 2>&1 | tee -a "${log}"
import sys
inf, outf, minlen = sys.argv[1], sys.argv[2], int(sys.argv[3])

def fasta_iter(path):
    h, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq)
                h = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if h is not None:
            yield h, "".join(seq)

kept = dropped = 0
with open(outf, "w") as out:
    for h, s in fasta_iter(inf):
        s = s.replace(" ", "").upper()
        if len(s) >= minlen:
            out.write(f">{h}\n{s}\n")
            kept += 1
        else:
            dropped += 1

print(f"[minlen] {inf}: kept={kept} dropped={dropped} (minlen={minlen})", file=sys.stderr)
PY

    # Flye optional: if missing, create empty flye_filt so cat always works
    if [[ -s "${FLYE_CTG}" ]]; then
      say "  Filtering Flye contigs to >= ${MIN_CONTIG_LENGTH} bp"
      python - <<'PY' "${FLYE_CTG}" "${FLYE_FILT}" "${MIN_CONTIG_LENGTH}" 2>&1 | tee -a "${log}"
import sys
inf, outf, minlen = sys.argv[1], sys.argv[2], int(sys.argv[3])

def fasta_iter(path):
    h, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq)
                h = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if h is not None:
            yield h, "".join(seq)

kept = dropped = 0
with open(outf, "w") as out:
    for h, s in fasta_iter(inf):
        s = s.replace(" ", "").upper()
        if len(s) >= minlen:
            out.write(f">{h}\n{s}\n")
            kept += 1
        else:
            dropped += 1

print(f"[minlen] {inf}: kept={kept} dropped={dropped} (minlen={minlen})", file=sys.stderr)
PY
    else
      say "  WARNING: Flye contigs not available; creating empty ${FLYE_FILT}"
      : > "${FLYE_FILT}"
    fi

    say "  SPAdes filtered size: $(stat -c%s "${SPADES_FILT}" 2>/dev/null || echo 0) bytes"
    say "  Flye filtered size:   $(stat -c%s "${FLYE_FILT}" 2>/dev/null || echo 0) bytes"

    if [[ ! -s "${SPADES_FILT}" ]]; then
      echo "ERROR: SPAdes produced 0 contigs >= ${MIN_CONTIG_LENGTH} bp (SPAdes is required)" >&2
      exit 2
    fi
    if [[ ! -s "${FLYE_FILT}" ]]; then
      say "  WARNING: Flye produced 0 contigs >= ${MIN_CONTIG_LENGTH} bp (continuing with SPAdes only)"
    fi

    cat "${SPADES_FILT}" "${FLYE_FILT}" > "${COMBINED}"

    require_file "${COMBINED}"
    say "  Combined contigs written: ${COMBINED}"
  fi
fi

# =============================================================================
# 5) RagTag patch
# =============================================================================
if (( START_STEP <= 5 )); then
  if step_outputs_exist "${PATCHED}"; then
    say "Step 5: RagTag output already present; skipping"
  else
    say "Step 5: RagTag patch"
    elapsed_time

    say "  RagTag inputs:"
    say "    REF=${REF}"
    say "    combined=${COMBINED}"

    require_file "${REF}"
    require_file "${COMBINED}"

    rm -rf "${RAGTAG_OUT}"
    mkdir -p "${RAGTAG_OUT}"

    ragtag.py patch -f 250 -u -w -t "${threads}" -o "${RAGTAG_OUT}" "${REF}" "${COMBINED}" 2>&1 | tee -a "${log}"

    require_file "${PATCHED}"
  fi
fi

# =============================================================================
# 6) Pilon polishing (SR + ONT each round, early-stop)
# =============================================================================
if (( START_STEP <= 6 )); then
  if step_outputs_exist "${FINAL_OUT}"; then
    say "Step 6: final output already present; skipping"
  else
    say "Step 6: Pilon polishing (SR + ONT)"
    elapsed_time

    require_file "${PILON_JAR}"
    require_file "${SR1}"
    require_file "${SR2}"
    require_file "${LR}"

    pilon_dir="${wkdir}/${name}_pilon"
    mkdir -p "${pilon_dir}"

    pilon_common_args=()
    if [[ "${PILON_TRACKS}" -eq 1 ]]; then
      pilon_common_args+=( --tracks )
    fi
    if [[ "${PILON_MINDEPTH}" -gt 0 ]]; then
      pilon_common_args+=( --mindepth "${PILON_MINDEPTH}" )
    fi

    LR_MM2_OPTS=( -t "${threads}" -a -x map-ont )
    if [[ "${DISABLE_LR_SECONDARY}" -eq 1 ]]; then
      LR_MM2_OPTS+=( --secondary=no )
    fi
    if [[ "${ENABLE_LR_SOFTCLIP}" -eq 1 ]]; then
      LR_MM2_OPTS+=( -Y )
    fi

    require_file "${PATCHED}"

    prev="${PATCHED}"
    final=""

    for r in $(seq 1 "${MAX_PILON_ROUNDS}"); do
      outprefix="${pilon_dir}/${name}_pilon_r${r}"
      outfa="${outprefix}.fasta"

      say "  Pilon round ${r}: remap SR+ONT to current genome"
      bwa index "${prev}" >/dev/null 2>&1 || true

      sr_bam_raw="${pilon_dir}/${name}_sr_r${r}.bam"
      bwa mem -t "${threads}" "${prev}" "${SR1}" "${SR2}" \
        | samtools sort -@ "${threads}" -o "${sr_bam_raw}"
      samtools index "${sr_bam_raw}"

      sr_bam="${sr_bam_raw}"
      if [[ "${ENABLE_SR_MAPQ_FILTER}" -eq 1 ]]; then
        say "    SR MAPQ filter enabled: >=${SR_MIN_MAPQ}"
        sr_bam="${pilon_dir}/${name}_sr_r${r}.q${SR_MIN_MAPQ}.bam"
        samtools view -b -q "${SR_MIN_MAPQ}" "${sr_bam_raw}" > "${sr_bam}"
        samtools index "${sr_bam}"
      fi

      ont_bam="${pilon_dir}/${name}_ont_r${r}.bam"
      minimap2 "${LR_MM2_OPTS[@]}" "${prev}" "${LR}" \
        | samtools sort -@ "${threads}" | samtools view -@ "${threads}" -b -q "${SR_MIN_MAPQ}" > "${ont_bam}"
      samtools index "${ont_bam}"

      say "  Pilon round ${r}: run pilon (mindepth=${PILON_MINDEPTH}, tracks=${PILON_TRACKS})"
      java -Xmx"${JAVA_MEM}" -jar "${PILON_JAR}" \
        --genome "${prev}" \
        --frags "${sr_bam}" \
        --nanopore "${ont_bam}" \
        --output "${outprefix}" \
        --threads "${threads}" \
        --fix all \
        --changes --vcf \
        "${pilon_common_args[@]}" 2>&1 | tee -a "${log}"

      require_file "${outfa}"

      prev_md5="$(grep -v '^>' "${prev}" | tr -d '\n' | md5sum | awk '{print $1}')"
      out_md5="$(grep -v '^>' "${outfa}" | tr -d '\n' | md5sum | awk '{print $1}')"
      if [[ "${prev_md5}" == "${out_md5}" ]]; then
        say "  Early stop: round ${r} produced no sequence changes."
        final="${outfa}"
        break
      fi

      prev="${outfa}"
      final="${outfa}"
    done

    require_file "${final}"
    cp -f "${final}" "${FINAL_OUT}"
  fi
fi

say "=== COMPLETE ==="
say "Final polished assembly: ${FINAL_OUT}"
say "Log: ${log}"

if [[ -s "${BS_DIR}/${name}.fa" ]]; then
  fsa --maxram 250000000 "${BS_DIR}/${name}.fa" "${FINAL_OUT}" > "${BS_DIR}/${name}_fullrefstest.msa"
  say "Wrote MSA: ${BS_DIR}/${name}_fullrefstest.msa"
else
  say "NOTE: ${BS_DIR}/${name}.fa not found; skipping fsa alignment"
fi

end_time=$(date +%s)
elapsed=$((end_time - start_time))
say "Total time elapsed: $((elapsed/3600))h $(((elapsed%3600)/60))m $((elapsed%60))s"

mkdir -p "${BS_DIR}/final"

cp "${FINAL_OUT}" "${BS_DIR}/final/"
cd "${BS_DIR}/final/" || exit 2

python "${BS_DIR}/fastareformat.py" "${name}_final.fasta"
python "${BS_DIR}/rename_header.py" "${name}_final_reformatted.fasta"
rm "${name}_final.fasta"
rm "${name}_final_reformatted.fasta"
cd "${BS_DIR}" || exit 2
