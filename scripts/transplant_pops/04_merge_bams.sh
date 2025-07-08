#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --ntasks=4
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --account=kingsley

# Run with: sbatch 04_merge_bams.sh <BAM_PREFIX>

module load samtools/1.19

BAM_PREFIX=$1
THREADS=4

mkdir -p 04_merge_bams
samtools merge -@ $THREADS 04_merge_bams/${BAM_PREFIX}.realignGA5_C4masked.sort.merged.bam 03_markdups_addRG/${BAM_PREFIX}_*.realignGA5_C4masked.sort.picard.bam
