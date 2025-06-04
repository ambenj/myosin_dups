#!/bin/bash
#SBATCH --job-name=pb_align
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --account=kingsley

# Align pacbio read to reference
# Run with sbatch pbmm2_align.sh <ref.fasta|.mmi> <reads.fastq|.fofn> <out.bam> <CCS|SUBREAD>

REF=$1
FASTQ=$2
OUT_BAM=$3
PRESET=$4 #CCS/HIFI or SUBREAD

pbmm2 align $REF $FASTQ $OUT_BAM --preset $PRESET --sort --log-level INFO