#!/bin/bash
#SBATCH --job-name=kmer
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=24G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Run kmer counting python script
# Run with:
# sbatch count_kmers.sh <kmers.txt> <reads_1.fq.gz> <reads_2.fq.gz> <output_prefix>

module load biopython/1.81

KMERS=$1
FOR_READS=$2
REV_READS=$3
PREFIX=$4

python /labs/kingsley/ambenj/myosin_dups/scripts/ASE/count_kmers.py $KMERS $FOR_READS $REV_READS --prefix $PREFIX