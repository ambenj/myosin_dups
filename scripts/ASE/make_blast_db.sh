#!/bin/bash
#SBATCH --job-name=blastdb
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --account=kingsley

# Create blast db from input fasta file
# Run with:
# sbatch make_blast_db.sh <ref.fa>


module load ncbi-blast/2.15.0+ 

# Get input fasta reference file
INPUT_FA=$1

# Make blast database
makeblastdb -in $INPUT_FA -dbtype nucl