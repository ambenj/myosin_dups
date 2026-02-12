#!/usr/bin/env python3
import sys
import re

infile = sys.argv[1]
outfile = sys.argv[2]

def clean_coord(coord, is_start, seq_len):
    """Remove < or > unless coordinate is truly at sequence boundary."""
    raw = coord
    stripped = coord.lstrip("<>").strip()
    pos = int(stripped)

    # Keep boundary indicators at real boundaries
    if is_start and pos == 1:
        return raw
    if (not is_start) and pos == seq_len:
        return raw

    # Otherwise remove < or >
    return stripped


def is_ncRNA(ftype):
    """Determine if feature type is a non-coding RNA."""
    return ("RNA" in ftype and ftype != "mRNA")


# ----------------------
# Parse input
# ----------------------
with open(infile) as f:
    lines = f.readlines()

# Determine sequence length from REFERENCE line
seq_len = None
for line in lines:
    if "\tREFERENCE" in line:
        start, end, _ = line.strip().split("\t", 2)
        seq_len = int(end)
        break

if seq_len is None:
    sys.exit("ERROR: Could not determine contig length from REFERENCE line.")


# --------------------------
# Process features
# --------------------------
out = []
i = 0

last_gene_indices = []     # indices in output list for gene feature intervals
last_gene_coords = []      # (start, end) pairs

current_feature_type = None
current_feature_interval_indices = []  # for multi-interval features

while i < len(lines):
    line = lines[i]

    # Match lines like: start end type   OR   start end (interval continuation)
    m = re.match(r"^\s*([<>]?\d+)\s+([<>]?\d+)(?:\s+(\S+))?", line)
    if m:
        start, end, ftype = m.groups()

        # New feature begins when ftype is present
        if ftype:
            current_feature_type = ftype
            current_feature_interval_indices = []

            # If feature is a gene, store indices for fixing later
            if ftype == "gene":
                last_gene_indices = []
                last_gene_coords = []

        # Record interval for potential fixing
        interval_index = len(out)
        current_feature_interval_indices.append(interval_index)

        # If this is a gene interval, store it
        if current_feature_type == "gene":
            last_gene_indices.append(interval_index)
            last_gene_coords.append((start, end))

        # If it's an ncRNA interval
        if is_ncRNA(current_feature_type):
            # Fix coordinates for this interval
            new_start = clean_coord(start, True, seq_len)
            new_end   = clean_coord(end, False, seq_len)
            # Overwrite line with corrected interval
            if ftype:
                line = f"{new_start}\t{new_end}\t{current_feature_type}\n"
            else:
                line = f"{new_start}\t{new_end}\n"

            # Also update parent gene intervals
            for idx, (gstart, gend) in zip(last_gene_indices, last_gene_coords):
                fixed_start = clean_coord(gstart, True, seq_len)
                fixed_end   = clean_coord(gend, False, seq_len)
                out[idx] = f"{fixed_start}\t{fixed_end}\tgene\n"

                # Update stored coords for future ncRNA intervals
                last_gene_coords = [
                    (clean_coord(s, True, seq_len), clean_coord(e, False, seq_len))
                    for (s, e) in last_gene_coords
                ]

        out.append(line)
        i += 1
        continue

    # Non-feature line, copy as-is
    out.append(line)
    i += 1


# --------------------------
# Write output
# --------------------------
with open(outfile, "w") as f:
    f.writelines(out)

print(f"Finished. Updated ncRNA + parent gene + multi-interval features written to {outfile}")
