#!/bin/bash
#SBATCH --job-name=samCov
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem=16G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Run with: sbatch 06_samtools_coverage_roi.sh <bam> <regions.bed> <outdir>

# Load modules
module load samtools/1.19

BAM=$1
REGIONS=$2
OUTDIR=$3

# Get name for output file
FILE=${BAM##*/}
BASE=${FILE%.bam}
OUT_FILE="${OUTDIR}/${BASE}_mapq3_roi_coverage.txt" 

# Make out directory if it does not exist
mkdir -p $OUTDIR

printf "desc\tregion\tchr\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tbam\n" > $OUT_FILE

# Get coverage for rois
while IFS=$'\t' read -r -a region_list
do
	region="${region_list[0]}:${region_list[1]}-${region_list[2]}"
	result=$(samtools coverage -H -q 3 --ff 260 -r $region $BAM)
	printf "${region_list[3]}\t${region}\t${result}\t${BAM}\n" >> $OUT_FILE
 
done < $REGIONS

