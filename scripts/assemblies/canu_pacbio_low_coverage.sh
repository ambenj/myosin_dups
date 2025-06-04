#!/bin/bash
#SBATCH --job-name=canu_low
#SBATCH --mem=64G
#SBATCH --time=2-0:00:00
#SBATCH --cpus-per-task=4
#SBATCH --account=kingsley

# Run with:
# sbatch canu_pacbio_low_coverage.sh <input.fastq> <base_name> <genome_size>

INPUT_FASTQ=$1
BASE=$2
SIZE=$3
DIR=${INPUT_FASTQ%/*}/canuLow_${BASE}

mkdir -p $DIR

# Run canu
canu -p $BASE -d $DIR gridOptions='--account=kingsley --time=2-0:00:00' genomeSize=${SIZE} correctedErrorRate=0.075 stopOnLowCoverage=2 minInputCoverage=2 -pacbio $INPUT_FASTQ