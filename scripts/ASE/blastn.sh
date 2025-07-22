#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH --time=12:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --account=kingsley

# Run blast (megablast) search on a set of input sequences and output a table of blast hit
# Run with:
# sbatch blastn.sh <query.fa> <database.fa> <outfile.txt>

QUERY=$1
DB=$2 
OUTFILE=$3


module load ncbi-blast/2.15.0+
blastn -task blastn-short -query $QUERY -db $DB -num_threads ${SLURM_CPUS_PER_TASK} -out ${OUTFILE} -evalue 1 \
        -outfmt "6 qseqid sseqid evalue stitle pident length mismatch gapopen qstart qend sstart send bitscore"