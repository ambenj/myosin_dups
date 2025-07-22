#!/bin/bash
#SBATCH --job-name=blat
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=24G
#SBATCH --account=kingsley

# Run with:
# sbatch blat.sh <ref.fa> <query.fa> <prefix>

# Get user input
REF=$1
QUERY=$2
PREFIX=$3

# Run blat
blat $REF $QUERY ${PREFIX}.BLAT.psl -tileSize=8 -stepSize=1 -minMatch=1 -minScore=20 -oneOff=1

# Sort output file
head -n 5 ${PREFIX}.BLAT.psl > ${PREFIX}.BLAT.sorted.psl
tail -n +6 ${PREFIX}.BLAT.psl | sort -nr -k1 >> ${PREFIX}.BLAT.sorted.psl