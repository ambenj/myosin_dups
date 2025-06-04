#!/bin/bash
#SBATCH --job-name=blat
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --account=kingsley

# Run with:
# sbatch blat.sh <ref.fa>

# Get user input
REF=$1
QUERY=/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/stickleback_v5_myosin_genes.fasta

# Run blat
blat $REF $QUERY ${REF}.BLAT.psl

# Sort output file
head -n 5 ${REF}.BLAT.psl > ${REF}.BLAT.sorted.psl
tail -n +6 ${REF}.BLAT.psl | sort -nr -k1 >> ${REF}.BLAT.sorted.psl
