#!/bin/bash
#SBATCH --job-name=kmer
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=24G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Check if correct number of arguments is provided
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 input_kmer_file sample_2.fastq.gz sample_1.fastq.gz sample output_file"
  exit 1
fi

# Assign input arguments to variables
input_file=$1
fileA=$2
fileB=$3
sample=$4
output_file=$5

# Check if files exist
if [ ! -f "$input_file" ]; then
  echo "Input file not found!"
  exit 1
fi

if [ ! -f "$fileA" ]; then
  echo "File A not found!"
  exit 1
fi

if [ ! -f "$fileB" ]; then
  echo "File B not found!"
  exit 1
fi

# Read the header line
header=$(head -n 1 "$input_file")

# Output the header with additional columns for the counts
echo -e "$header\tcount_kmer\tcount_kmerRC\tsample" > $output_file

# Read the input file line by line, skipping the header
tail -n +2 "$input_file" | while IFS=$'\t' read -r kmer_name kmer exon spans_intron sequence_source length min max min_original max_original allele source comments kmer_rc; do
  # Search for seq4 in fileA and count the hits
  echo -e "Searching $kmer in $fileA"
  count_kmer=$(zgrep -e "$kmer" "$fileA" | wc -l)
  echo -e "Kmer count: $count_kmer"

  # Search for seq5 in fileB and count the hits
  echo -e "Searching $kmer_rc in $fileB"
  count_kmerRC=$(zgrep -e "$kmer_rc" "$fileB" | wc -l)
  echo -e "Kmer count: $count_kmerRC"

  # Output the original columns along with the hit counts
  echo -e "$kmer_name\t$kmer\t$exon\t$spans_intron\t$sequence_source\t$length\t$min\t$max\t$min_original\t$max_original\t$allele\t$source\t$comments\t$kmer_rc\t$count_kmer\t$count_kmerRC\t$sample" >> $output_file
done