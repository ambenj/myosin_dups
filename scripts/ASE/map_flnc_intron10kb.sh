#!/bin/bash
#SBATCH --job-name=align
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --account=kingsley

# Align FLNC reads to reference genome 

# Run with sbatch map_flnc_intron10kb.sh <base_name> <ref.fa> <ref_sub_dir>

BASE=$1
REF=$2
RAW_DIR="/labs/kingsley/ambenj/myosin_dups/raw/kinnex/temperature_ase/flnc_reads"
DIR="/labs/kingsley/ambenj/myosin_dups/analysis/ase/flnc_mapping"
REF_SUBDIR=$3

module load samtools/1.20

# Make reference sub directory if it does not exist
mkdir -p ${DIR}/${REF_SUBDIR}

# Align FLNC reads with intron size limit of 10kb
pbmm2 align --preset ISOSEQ -G 10000 --sort ${RAW_DIR}/${BASE}.flnc.bam ${REF} ${DIR}/${REF_SUBDIR}/${BASE}.mapped.bam

# Filter alignments to primary alignments with MAPQ > 3
samtools view -@ 8  -F 256 -q 3 -b ${DIR}/${REF_SUBDIR}/${BASE}.mapped.bam > ${DIR}/${REF_SUBDIR}/${BASE}.mapped.primary.mapq3.bam

# Index filtered alignments
samtools index ${DIR}/${REF_SUBDIR}/${BASE}.mapped.primary.mapq3.bam 

# Downsample reads to 10% for easier visibility in IGV
samtools view -@ 8 -b -s 0.1 ${DIR}/${REF_SUBDIR}/${BASE}.mapped.primary.mapq3.bam > ${DIR}/${REF_SUBDIR}/${BASE}.mapped.primary.mapq3.10.bam 

# Index dowsampled alignments
samtools index ${DIR}/${REF_SUBDIR}/${BASE}.mapped.primary.mapq3.10.bam 