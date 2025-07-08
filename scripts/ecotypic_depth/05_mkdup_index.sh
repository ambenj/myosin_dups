#!/bin/bash
#SBATCH --job-name=dup_index
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --account=kingsley

# Run with: sbatch 05_mkdup_index.sh <BAM>

module load samtools/1.19
module load picard/3.1.1

BAM=$1

# Get name for output file
FILE=${BAM##*/}
BASE=${FILE%.bam}
BAM_OUT="05_mkdup_index/${BASE}.mkdup.bam"
METRICS="05_mkdup_index/${BASE}.mkdup.metrics.txt"

# Make output directory if it doesn't exist
mkdir -p 05_mkdup_index

# Mark duplicates in merged bam
picard MarkDuplicates -I ${BAM} -O ${BAM_OUT} -M ${METRICS}

# Index bam
samtools index $BAM_OUT
