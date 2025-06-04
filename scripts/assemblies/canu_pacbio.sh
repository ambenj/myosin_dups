#!/bin/bash
#SBATCH --job-name=canu
#SBATCH --mem=64G
#SBATCH --time=2-0:00:00
#SBATCH --cpus-per-task=4
#SBATCH --account=kingsley

# Run with:
# sbatch canu.sh <input.fastq> <base> <genome_size>

INPUT_FASTQ=$1
BASE=$2
SIZE=$3
DIR=${INPUT_FASTQ%/*}/canu_${BASE}

# Run hicanu
canu -p $BASE -d $DIR gridOptions='--account=kingsley --time=4-0:00:00' genomeSize=${SIZE} -pacbio $INPUT_FASTQ