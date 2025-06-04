#!/bin/bash
#SBATCH --job-name=canu
#SBATCH --mem=64G
#SBATCH --time=4-0:00:00
#SBATCH --cpus-per-task=8
#SBATCH --account=kingsley

# Run with:
# sbatch hicanu.sh <input.fastq> <base>

INPUT_FASTQ=$1
BASE=$2
DIR=/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/$BASE/hicanu/

# Run hicanu
canu -p $BASE -d $DIR gridOptions='--account=kingsley --time=4-0:00:00' genomeSize=500m -pacbio-hifi $INPUT_FASTQ