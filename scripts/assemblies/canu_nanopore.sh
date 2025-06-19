#!/bin/bash
#SBATCH --job-name=canu
#SBATCH --mem=64G
#SBATCH --time=2-0:00:00
#SBATCH --cpus-per-task=4
#SBATCH --account=kingsley

# Run with:
# sbatch canu_nanopore.sh <input.fastq> <out_dir> <base> <genome_size>

INPUT_FASTQ=$1
DIR=$2
BASE=$3
SIZE=$4

# Run hicanu
canu -p $BASE -d ${DIR}/${BASE} gridOptions='--account=kingsley --time=2-0:00:00' genomeSize=${SIZE}  cnsMemory=20 gridOptionsOVS="--mem-per-cpu=64g" -nanopore $INPUT_FASTQ