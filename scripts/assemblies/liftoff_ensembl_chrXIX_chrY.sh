#!/bin/bash
#SBATCH --job-name=liftoff
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --account=kingsley

# Liftoff gene annotations from one assembly to another

# Run with:
# sbatch liftoff_ensembl_chrXIX_chrY.sh <assembly.fasta>

# Get target assembly
TARGET=$1

# Source annotation files
SOURCE_GTF1="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/stickleback_v5_ensembl_genes_chrXIX_agat_MYHrevised.gtf"
SOURCE_GTF2="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/stickleback_v5_maker_genes_chrY_agat_MYHrevised.gtf"
SOURCE_GTF3="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/stickleback_v5_maker_genes_chrY_agat_MYHrevised.gtf"
SOURCE_FA="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly.fa"

# Get name for output file
BASE=${TARGET%.fasta}
DIR=${TARGET%/*}

echo ${BASE}
echo ${DIR}


# Run liftoff on ensembl genes 
#liftoff -g $SOURCE_GTF1 -o ${BASE}.liftoffEnsemblchrXIX.gtf -u ${BASE}.unliftedEnsemblchrXIX.gtf -dir ${DIR}/liftoffEnsemblchrXIX_intermediate_files -copies -sc 0.98 $TARGET $SOURCE_FA

# Run liftoff on make ChrY genes
liftoff -g $SOURCE_GTF2 -o ${BASE}.liftoffchrY.gtf -u ${BASE}.unliftedchY.gtf -dir ${DIR}/liftoffchrY_intermediate_files -copies -sc 0.98 $TARGET $SOURCE_FA

# Run liftoff on make ChrY genes
liftoff -g $SOURCE_GTF3 -o ${BASE}.liftoffchrY.gtf -u ${BASE}.unliftedchY.gtf -dir ${DIR}/liftoffchrY_intermediate_files -copies -sc 0.98 $TARGET $SOURCE_FA