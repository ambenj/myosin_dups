#!/bin/bash
#SBATCH --job-name=index
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --account=kingsley

# Index reference genome for pbmm2 alignment

# Run with: sbatch index_isoseq.sh <ref.fa>

REF=$1

pbmm2 index --preset ISOSEQ $REF ${REF}.isoseq.mmi 