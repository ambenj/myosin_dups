#!/bin/bash
#SBATCH --job-name=samCov
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Run with: sbatch 06_samtools_coverage_wg.sh <bam> <outdir>

# Load modules
module load samtools/1.19

BAM=$1
OUTDIR=$2

# Get name for output file
FILE=${BAM##*/}
BASE=${FILE%.bam}
OUT_FILE="${OUTDIR}/${BASE}_mapq3_wg_coverage.txt"

# Make out directory if it does not exist
mkdir -p $OUTDIR

printf "chr\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" > $OUT_FILE

# Get coverage for whole genome
samtools coverage -H -q 3 --ff 260 $BAM >> $OUT_FILE
