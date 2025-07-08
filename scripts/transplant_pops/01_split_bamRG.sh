#!/bin/bash
#SBATCH --job-name=split
#SBATCH --ntasks=8
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --account=kingsley

# Run with: sbatch 01_split_bamRG.sh <BAM>

module load samtools/1.19

BAM=$1
THREADS=8

mkdir -p 01_split_bams
cd 01_split_bams
samtools split -@ $THREADS $BAM

