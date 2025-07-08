#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --ntasks=16
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Get input files
ALIGN_FILE=$1
BASE=${ALIGN_FILE%.fasta}
OUT_PREFIX="${BASE}_iqtree_boot100"

# Run IQ-TREE
## Default is using ModelFinderPlus to find best fit model
## Using standard bootstrap (minimum recommended is 100)
iqtree -s $ALIGN_FILE -b 100 -T AUTO -ntmax 16 --prefix $OUT_PREFIX