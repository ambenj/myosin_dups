#!/bin/bash
#SBATCH --job-name=gff_to_tbl
#SBATCH --ntasks=1
#SBATCH --time=30:00
#SBATCH --mem=4G
#SBATCH --account=kingsley

# Convert gff annotations from liftoff to genbank feature table format for genbank submission

# Run with:
# sbatch gff_to_tbl.sh <genome.fasta> <annotations.gff>

# Get input
FASTA=$1
GFF=$2

# Get basename
GFF_BASE="${GFF%.*}"

# Tidy gff and fix phase
module load genometools/1.6.1
gt gff3 -sort -tidy -retainids $GFF > ${GFF_BASE}.tidied.gff

# Replace lnc_RNA with ncRNA
sed -i 's/lnc_RNA/ncRNA/g' ${GFF_BASE}.tidied.gff

# Convert to tbl format for genbank submission
python2.7 /labs/kingsley/ambenj/tools/genomeannotation-GAG-997e384/gag.py --fix_start_stop --fasta $FASTA --gff ${GFF_BASE}.tidied.gff --out ${GFF_BASE}.tidied_gag_output

# Replace %2C with ,; Need to do this after tbl format conversion otherwise splits products on multiple lines
sed -i 's/%2C/,/g' ${GFF_BASE}.tidied_gag_output/genome.tbl

# Replace NCBI database label
sed -i 's/|ncbi|/|KingsleySU|/g' ${GFF_BASE}.tidied_gag_output/genome.tbl

# Remove false incomplete notation from non-coding RNAs
python3 /labs/kingsley/ambenj/myosin_dups/scripts/genbank_submission/clean_tbl_ncRNA.py ${GFF_BASE}.tidied_gag_output/genome.tbl ${GFF_BASE}.tidied_gag_output/genome.ncRNAfixed.tbl