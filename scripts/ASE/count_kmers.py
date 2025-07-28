import gzip
import csv
import argparse
from collections import defaultdict
from Bio import SeqIO

def open_fastq(filename):
    return gzip.open(filename, 'rt') if filename.endswith('.gz') else open(filename, 'r')

def load_kmers(kmer_file):
    kmers = []
    with open(kmer_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            kmers.append({
                'name': row['kmer_name'],
                'kmer': row['kmer'],
                'kmer_rc': row['kmer_rc'],
                'match_pairs': set()  # store normalized read ID prefixes
            })
    return kmers

def normalize_read_id(read_id):
    # Normalize /1 /2 or space-separated identifiers
    return read_id.split('/')[0].split()[0]

def match_kmers(fastq_file, kmers, match_field):
    with open_fastq(fastq_file) as handle:
        for record in SeqIO.parse(handle, 'fastq'):
            seq = str(record.seq)
            norm_id = normalize_read_id(record.id)

            for k in kmers:
                if k[match_field] in seq:
                    k['match_pairs'].add(norm_id)


def write_match_output(kmers, output_match_file):
    with open(output_match_file, 'w') as f:
        for k in kmers:
            for rid in sorted(k['match_pairs']):
                f.write(f"{k['name']}\t{rid}\n")

def write_kmer_count_output(kmers, output_kmer_file):
    with open(output_kmer_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['kmer_name', 'kmer', 'kmer_rc', 'match_count'])
        for k in kmers:
            writer.writerow([k['name'], k['kmer'], k['kmer_rc'], len(k['match_pairs'])])

def main():
    parser = argparse.ArgumentParser(description="Search paired FASTQ files for kmer matches.")
    parser.add_argument("kmer_file", help="Tab-delimited file with columns: kmer_name, kmer, kmer_rc")
    parser.add_argument("forward_fastq", help="Forward reads FASTQ file (can be .gz)")
    parser.add_argument("reverse_fastq", help="Reverse reads FASTQ file (can be .gz)")
    parser.add_argument("--prefix", default="", help="Prefix to add to output file names")
    args = parser.parse_args()

    # Compose output filenames with prefix
    match_reads_file = f"{args.prefix}_matching_reads.txt" if args.prefix else "matching_reads.txt"
    updated_kmer_file = f"{args.prefix}_kmers_with_counts.tsv" if args.prefix else "kmers_with_counts.tsv"


    # Load kmers
    kmers = load_kmers(args.kmer_file)

    # Search forward reads for reverse complement kmer matches
    match_kmers(args.forward_fastq, kmers, 'kmer_rc')

    # Search reverse reads for forward kmer matches
    match_kmers(args.reverse_fastq, kmers, 'kmer')

    # Write outputs
    write_match_output(kmers, match_reads_file)
    write_kmer_count_output(kmers, updated_kmer_file)

if __name__ == '__main__':
    main()
