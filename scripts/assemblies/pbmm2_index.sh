#!/bin/bash
#SBATCH --job-name=index
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --account=kingsley

# Index reference for pbmm2 alignment
# Run with pbmm2_index.sh <ref.fa> <SUBREAD|HIFI|CCS|ISOSEQ>

# Get user input
REF=$1
PRESET=$2

# Get name for output index file
BASE=${REF%.fasta}
BASE=${BASE%.fa}
INDEX=${BASE}.${PRESET}.mmi

echo $INDEX

# Index reference
pbmm2 index $REF $INDEX --preset $PRESET