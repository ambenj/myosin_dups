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
SOURCE_GFF_all="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/ncbi_annotations_genes_build100_converted.MYHrevised.sorted.gff3"
SOURCE_GFF_chrXIX="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/ncbi_annotations_genes_build100_converted.MYHrevised.sorted.chrXIX.gff3"
SOURCE_GFF_chrY="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/ncbi_annotations_genes_build100_converted.MYHrevised.sorted.chrY.gff3"

SOURCE_FA="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly.fa"

# Get name for output file
BASE=${TARGET%.fasta}
BASE=${BASE%.fa}
DIR=${TARGET%/*}


# Run liftoff on all NCBI annotations with revised myosin genes
liftoff -g $SOURCE_GFF_all -o ${BASE}.NCBI.MYHrevised.gff -u ${BASE}.NCBI.MYHrevised.unlifted.gff -dir ${BASE}_liftoff.NCBI.MYHrevised_intermediate_files -copies -sc 0.92 $TARGET $SOURCE_FA
sed -E 's/CDS_[0-9]+/CDS/g' - ${BASE}.NCBI.MYHrevised.gff > ${BASE}.NCBI.MYHrevised_geneious.gff


# Run liftoff on chrXIX NCBI annotations with revised myosin genes
liftoff -g $SOURCE_GFF_chrXIX -o ${BASE}.NCBI.MYHrevised_chrXIX.gff -u ${BASE}.NCBI.MYHrevised_chrXIX.unlifted.gff -dir ${BASE}_liftoff.NCBI.MYHrevised_chrXIX_intermediate_files -copies -sc 0.92 $TARGET $SOURCE_FA
sed -E 's/CDS_[0-9]+/CDS/g' - ${BASE}.NCBI.MYHrevised_chrXIX.gff > ${BASE}.NCBI.MYHrevised_chrXIX_geneious.gff

# Run liftoff on chrY NCBI annotations with revised myosin genes
liftoff -g $SOURCE_GFF_chrY -o ${BASE}.NCBI.MYHrevised_chrY.gff -u ${BASE}.NCBI.MYHrevised_chrY.unlifted.gff -dir ${BASE}_liftoff.NCBI.MYHrevised_chrY_intermediate_files -copies -sc 0.92 $TARGET $SOURCE_FA
sed -E 's/CDS_[0-9]+/CDS/g' - ${BASE}.NCBI.MYHrevised_chrY.gff > ${BASE}.NCBI.MYHrevised_chrY_geneious.gff
