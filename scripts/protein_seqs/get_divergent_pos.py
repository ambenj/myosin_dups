from Bio import SeqIO
import pandas as pd
import argparse

def extract_divergent_positions(fasta_file, output_file):
    # Read aligned sequences into a dictionary
    records = list(SeqIO.parse(fasta_file, "fasta"))
    seq_dict = {record.id: str(record.seq) for record in records}

    # Ensure all sequences are the same length
    seq_lengths = {len(seq) for seq in seq_dict.values()}
    if len(seq_lengths) != 1:
        raise ValueError("All sequences must be the same length in the alignment.")

    seq_length = seq_lengths.pop()
    sequence_ids = list(seq_dict.keys())

    # Initialize list of rows for the output
    divergent_rows = []

    # Iterate over positions
    for pos in range(seq_length):
        column = [seq_dict[seq_id][pos] for seq_id in sequence_ids]
        non_gap_residues = [aa for aa in column if aa != '-']
        if len(set(non_gap_residues)) > 1:
            # Keep only truly variable (non-gap divergent) positions
            row = [pos + 1] + column  # 1-based position
            divergent_rows.append(row)

    # Create dataframe
    df = pd.DataFrame(divergent_rows, columns=["Position"] + sequence_ids)

    # Output to file
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Divergent positions saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract divergent positions from a protein alignment.")
    parser.add_argument("-i", "--input", required=True, help="Input protein alignment in FASTA format.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file with divergent positions.")

    args = parser.parse_args()
    extract_divergent_positions(args.input, args.output)
