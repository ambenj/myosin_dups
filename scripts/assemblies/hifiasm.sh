#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --mem=64G
#SBATCH --time=3-0:00:00
#SBATCH --cpus-per-task=8
#SBATCH --account=kingsley

# Perform de novo assembly of pacbio hifi reads using hifiasm, then convert output files to fasta format
# Run with:
# sbatch hifiasm.sh <input.fastq.gz> <base_name>

INPUT_FASTQ=$1
BASE=$2
DIR=/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/$BASE/hifiasm/

mkdir -p $DIR
cd $DIR

# Assemble using hifiasm
hifiasm -o $BASE -t ${SLURM_CPUS_PER_TASK} $INPUT_FASTQ

# Convert gfa to fa
awk '/^S/{print ">"$2;print $3}' ${BASE}.bp.hap1.p_ctg.gfa > ${BASE}.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${BASE}.bp.hap2.p_ctg.gfa > ${BASE}.bp.hap2.p_ctg.fa 
awk '/^S/{print ">"$2;print $3}' ${BASE}.bp.p_ctg.gfa > ${BASE}.bp.p_ctg.fa  
