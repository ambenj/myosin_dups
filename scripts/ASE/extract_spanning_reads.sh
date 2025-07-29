#!/bin/bash
#SBATCH --job-name=extract
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=24G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Run with:
# sbatch extract_spanning_reads.sh <sample.bam>

module load biopython/1.81

BAM=$1
OUTDIR=/labs/kingsley/ambenj/myosin_dups/analysis/ase/flnc_mapping/align_RABS-3copy_BEPA-5copy/exons3-35_counts
REGIONS=/labs/kingsley/ambenj/myosin_dups/analysis/ase/flnc_mapping/align_RABS-3copy_BEPA-5copy/RABS_BEPA_exons3-35_positions.txt

# Make output directory if it does not already exist
mkdir -p $OUTDIR

# Run script to extract and count reads spanning regions of interest
python /labs/kingsley/ambenj/myosin_dups/scripts/ASE/extract_spanning_reads.py \
  -b $BAM \
  -r $REGIONS \
  -o $OUTDIR \