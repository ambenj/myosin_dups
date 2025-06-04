#!/bin/bash
#SBATCH --job-name=get_roi
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --account=kingsley

# Run with:
# sbatch extract_reads_ROI.sh <chr:start-end> <input_bam> <output_fasta>

POS=$1    #Ex: chrXIX:2416650-2916650
BAM=$2
FASTA=$3

# Extract reads from region of interest
samtools view -u $BAM $POS | samtools fasta -0 $FASTA -