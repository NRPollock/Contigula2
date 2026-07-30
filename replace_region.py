#!/usr/bin/env python3
"""
replace_region.py – replace the region between a 5′ and 3′ motif
                     with a sequence taken from a second FASTA file
"""

import argparse
import sys
from textwrap import wrap


# ----------------------------------------------------------------------
# FASTA helpers (no external libraries needed)
# ----------------------------------------------------------------------
def read_single_fasta(path):
    """
    Return (header, sequence) from a FASTA file that contains exactly one
    sequence record.  Raises ValueError otherwise.
    """
    header = None
    seq_chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"More than one record found in {path}")
                header = line[1:]
            else:
                seq_chunks.append(line.upper())

    if header is None:
        raise ValueError(f"No FASTA record found in {path}")

    return header, "".join(seq_chunks)


def write_fasta(header, sequence, handle, width=60):
    """Write FASTA to an already opened file handle."""
    handle.write(f">{header}\n")
    for chunk in wrap(sequence, width):
        handle.write(chunk + "\n")


# ----------------------------------------------------------------------
# Core logic
# ----------------------------------------------------------------------
def replace_region(template_seq, insert_seq, five, three):
    five = five.upper()
    three = three.upper()

    start = template_seq.find(five)
    if start == -1:
        sys.exit("ERROR: 5′ motif not found in template sequence")

    stop = template_seq.find(three, start + len(five))
    if stop == -1:
        sys.exit("ERROR: 3′ motif not found downstream of 5′ motif")

    # slice: keep everything before the 5′ motif and after the 3′ motif
    new_seq = template_seq[:start] + insert_seq + template_seq[stop + len(three):]
    return new_seq


# ----------------------------------------------------------------------
# Command-line interface
# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Replace the region between 5′ and 3′ motifs "
                    "with a sequence from another FASTA file")
    p.add_argument("-t", "--template", required=True,
                   help="template FASTA (region will be replaced here)")
    p.add_argument("-i", "--insert", required=True,
                   help="FASTA containing the replacement sequence "
                        "(should include BOTH motifs)")
    p.add_argument("-5", "--five", required=True,
                   help="literal 5′ motif sequence")
    p.add_argument("-3", "--three", required=True,
                   help="literal 3′ motif sequence")
    p.add_argument("-o", "--output", default="edited.fasta",
                   help="output FASTA (default: edited.fasta)")
    args = p.parse_args()

    # ------------------------------------------------------------------
    # Read input files
    # ------------------------------------------------------------------
    tmpl_head, tmpl_seq = read_single_fasta(args.template)
    ins_head,  ins_seq  = read_single_fasta(args.insert)

    # ------------------------------------------------------------------
    # Perform replacement
    # ------------------------------------------------------------------
    new_seq = replace_region(tmpl_seq, ins_seq, args.five, args.three)

    # ------------------------------------------------------------------
    # Save result
    # ------------------------------------------------------------------
    with open(args.output, "w") as out_fh:
        write_fasta(tmpl_head + "_edited", new_seq, out_fh)

    print(f"✓   Written: {args.output}")


if __name__ == "__main__":
    main()
