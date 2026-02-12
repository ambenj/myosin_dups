#!/usr/bin/env python3

import sys
from Bio import SeqIO
import re

############################################################
# FASTA writer with proper line wrapping
############################################################
def write_fasta(header, seq, handle, width=60):
    handle.write(f">{header}\n")
    seq_str = str(seq)
    for i in range(0, len(seq_str), width):
        handle.write(seq_str[i:i+width] + "\n")


############################################################
# Parse GFF for features by Name=
############################################################
def parse_gff_features(gff_file, gene_names):
    features = {}
    for line in open(gff_file):
        if line.startswith("#"):
            continue
        parts = line.rstrip("\t\n").split("\t")
        if len(parts) < 9:
            continue
        contig, _, ftype, start, end, _, strand, _, attrs = parts
        for g in gene_names:
            # exact Name= match
            if re.search(r"Name=" + re.escape(g) + r"(;|$)", attrs):
                features[g] = (contig, int(start), int(end), strand)
    return features


############################################################
# Extract a region and subset GFF
############################################################
def extract_region(seq_record, gff_file, region_start, region_end, region_contig, new_contig_name):
    if region_start > region_end:
        region_start, region_end = region_end, region_start

    # Extract subsequence
    subseq = seq_record.seq[region_start - 1 : region_end]

    # Subset GFF and convert coordinates
    new_gff_entries = []
    for line in open(gff_file):
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        contig, source, ftype, start, end, score, strand, phase, attrs = parts

        if contig != region_contig:
            continue

        start, end = int(start), int(end)
        if start >= region_start and end <= region_end:
            new_start = start - region_start + 1
            new_end = end - region_start + 1

            parts[0] = new_contig_name
            parts[3] = str(new_start)
            parts[4] = str(new_end)

            new_gff_entries.append("\t".join(parts))

    return subseq, new_gff_entries


############################################################
# Orient sequence + GFF entries
############################################################
def orient(subseq, gff_entries, keep_forward=True):
    if keep_forward:
        return subseq, gff_entries

    # reverse complement
    rc = subseq.reverse_complement()
    L = len(rc)

    new_entries = []
    for line in gff_entries:
        parts = line.split("\t")
        old_start = int(parts[3])
        old_end = int(parts[4])

        new_start = L - old_end + 1
        new_end = L - old_start + 1

        parts[3] = str(new_start)
        parts[4] = str(new_end)

        # flip strand
        if parts[6] == "+":
            parts[6] = "-"
        elif parts[6] == "-":
            parts[6] = "+"

        new_entries.append("\t".join(parts))

    return rc, new_entries


############################################################
# MAIN
############################################################
def main():
    if len(sys.argv) != 5:
        print("Usage: script.py <fasta> <gff1> <gff2> <assembly_name>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    gff1 = sys.argv[2]
    gff2 = sys.argv[3]
    assembly = sys.argv[4]

    # Load FASTA
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    ############################################################
    # REGION 1: nlrc5 → syt19
    ############################################################
    genes1 = parse_gff_features(gff1, ["nlrc5", "syt19"])
    if ("nlrc5" not in genes1) or ("syt19" not in genes1):
        sys.exit("ERROR: Could not find nlrc5 and/or syt19 in GFF1")

    contig_n, n_start, n_end, n_strand = genes1["nlrc5"]
    contig_s, s_start, s_end, s_strand = genes1["syt19"]

    if contig_n != contig_s:
        sys.exit("ERROR: nlrc5 and syt19 are on different contigs in GFF1")

    # Get nlrc5 3'
    nlrc5_3 = n_end if n_strand == "+" else n_start
    # Get syt19 5'
    syt19_5 = s_start if s_strand == "+" else s_end

    region1_start = nlrc5_3
    region1_end = syt19_5

    # orientation logic
    forward1 = region1_start < region1_end

    header1_short = f"{assembly}_chrXIX_MYH3C"
    header1 = f"{header1_short} [Gasterosteus aculeatus] Gasterosteus aculeatus chromosome XIX MYH3C region, genomic sequence"

    seq_record1 = seqs[contig_n]

    region1_seq, g1_subset = extract_region(seq_record1, gff1, region1_start, region1_end, contig_n, header1_short)
    region1_seq, g1_subset = orient(region1_seq, g1_subset, forward1)

    ############################################################
    # REGION 2: LOC120812805 → LOC120812204
    ############################################################
    genes2 = parse_gff_features(gff2, ["LOC120812805", "LOC120812204"])
    if ("LOC120812805" not in genes2) or ("LOC120812204" not in genes2):
        sys.exit("ERROR: Could not find LOC120812805 and/or LOC120812204 in GFF2")

    contig_a, a_start, a_end, a_strand = genes2["LOC120812805"]
    contig_b, b_start, b_end, b_strand = genes2["LOC120812204"]

    if contig_a != contig_b:
        sys.exit("ERROR: LOC120812805 and LOC120812204 are on different contigs in GFF2")

    locA_5 = a_start if a_strand == "+" else a_end
    locB_3 = b_end if b_strand == "+" else b_start

    region2_start = locA_5
    region2_end = locB_3

    forward2 = region2_start < region2_end

    header2_short = f"{assembly}_chrY_MYH3C"
    header2 = f"{header2_short} [Gasterosteus aculeatus] Gasterosteus aculeatus chromosome Y MYH3C region, genomic sequence"

    seq_record2 = seqs[contig_a]

    region2_seq, g2_subset = extract_region(seq_record2, gff2, region2_start, region2_end, contig_a, header2_short)
    region2_seq, g2_subset = orient(region2_seq, g2_subset, forward2)

    ############################################################
    # WRITE OUTPUTS
    ############################################################

    fasta_out = f"{assembly}_MYH3C_regions.fasta"
    gff_out = f"{assembly}_MYH3C_regions.gff"

    # FASTA
    with open(fasta_out, "w") as outfa:
        write_fasta(header1, region1_seq, outfa)
        write_fasta(header2, region2_seq, outfa)

    # GFF (merged)
    with open(gff_out, "w") as outg:
        for entry in g1_subset:
            outg.write(entry + "\n")
        for entry in g2_subset:
            outg.write(entry + "\n")

    print("Done.")
    print(f"FASTA written to: {fasta_out}")
    print(f"GFF written to:   {gff_out}")


if __name__ == "__main__":
    main()
