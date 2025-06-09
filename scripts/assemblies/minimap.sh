#!/bin/bash
#SBATCH --job-name=minimap
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --account=kingsley

# Map oxford nanopre reads to reference using minimap2
# Run with:
# sbatch minimap.sh <ref.fa> <reads.fastq> <out_dir> <base>

REF=$1
INPUT_FASTQ=$2
DIR=$3
BASE=$4

# Map reads
minimap2 -ax map-ont -L -t 7 $REF $INPUT_FASTQ > ${DIR}/${BASE}.sam 

# Convert to bam
samtools view -@ 8 -b ${DIR}/${BASE}.sam | samtools sort -@ 8 > ${DIR}/${BASE}.sorted.bam 

# Index
samtools index -@ 8 ${DIR}/${BASE}.sorted.bam