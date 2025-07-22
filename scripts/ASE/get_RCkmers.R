library(Biostrings)
library(tidyverse)

# Run with:
# Rscript get_RCkmers.sh <input_table_file> <output_table_file> <output_forward.fasta> <output_reverse.fasta>
# Requires column called kmer that includes the kmer sequence to be reverse complimented

# Get input and output file names from command line
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[[1]]
output_file <- args[[2]]
output_fasta_file <- args[[3]]
output_fasta_RCfile <- args[[4]]

# Function to get reverse complement
get_reverse_complement <- function(seq) {  
	dna_seq <- DNAString(seq)
  rev_comp_seq <- reverseComplement(dna_seq)
  return(as.character(rev_comp_seq))
}

# Read data filed
df <- read_tsv(input_file)

# Add a new column with reverse complement sequences using dplyr
df <- df %>%
  dplyr::mutate(kmer_rc = map_chr(kmer, get_reverse_complement))
  
df %>%
   write_tsv(output_file)

# Write fasta file with forward kmers
sequences_set <- DNAStringSet(df$kmer)
names(sequences_set) <- df$kmer_name
writeXStringSet(sequences_set, filepath = output_fasta_file)

# Write fasta file with reverse complement kmers
sequences_set_RC <- DNAStringSet(df$kmer_rc)
names(sequences_set_RC) <- df$kmer_name
writeXStringSet(sequences_set_RC, filepath = output_fasta_RCfile)