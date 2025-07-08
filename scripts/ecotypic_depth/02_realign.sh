#!/bin/bash
#SBATCH --job-name=realign
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --account=kingsley

# Run with: sbatch 02_realign.sh <BAM>

module load samtools/1.19
module load bwa/0.7.17

BAM=$1
THREADS=8
REF="/labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly_MYH3C4dup_hardmasked_noChrY.fa"

# Get name for output file
FILE=${BAM##*/}
BASE=${FILE%.bam}
BAM_OUT="02_realign/${BASE}.realignGA5_C4masked.sort.bam"
echo "Realigned bam name is $BAM_OUT"

# Make output directory if it doesn't exist
mkdir -p 02_realign

# Extract and realign reads
samtools sort -n -@ $THREADS $BAM | \
samtools fastq -@ $THREADS -t - | \
bwa mem -p -t $THREADS $REF - | \
samtools view -buh -@ $THREADS - | \
samtools sort -@ $THREADS -o $BAM_OUT -

# Index bam
# samtools index -@ $THREADS -b $BAM_OUT ${BAM_OUT}.bai
