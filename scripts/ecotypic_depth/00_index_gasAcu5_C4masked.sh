#!/bin/bash
#SBATCH --job-name=index
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --account=kingsley

# Run with: sbatch 00_index_gasAcu5_C4masked.sh 

module load bwa/0.7.17

bwa index /labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly_MYH3C4dup_hardmasked_noChrY.fa

