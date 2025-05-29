#!/bin/bash
#SBATCH --job-name=bamfastq
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=kingsley

# Convert unmapped bam with pacbio reads to fastq format
# Run with:
# sbatch bam_to_fastq.sh <input.bam>

INPUT_BAM=$1
BASE=${INPUT_BAM%.bam}

samtools fastq -@ ${SLURM_CPUS_PER_TASK} $INPUT_BAM > ${BASE}.fastq