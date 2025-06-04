#!/bin/bash
#SBATCH --job-name=liftoff
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --account=kingsley

# Liftoff gene annotations from one assembly to another

# Run with:
# sbatch liftoff_NCBI.sh <assembly.fasta>

# Get target assembly
TARGET=$1

# Source annotation files
SOURCE_GFF="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/ncbi_annotations_genes_build100_converted.MYHrevised.sorted.gff3"
SOURCE_FA="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly.fa"

# Get name for output file
BASE=${TARGET%.fasta}
BASE=${BASE%.fa}
DIR=${TARGET%/*}


# Run liftoff on NCBI annotations with revised myosin genes
liftoff -g $SOURCE_GFF -o ${BASE}.NCBI.MYHrevised.gff -u ${BASE}.NCBI.MYHrevised.gff -dir ${BASE}_liftoff.NCBI.MYHrevised_intermediate_files -copies -sc 0.98 $TARGET $SOURCE_FA
