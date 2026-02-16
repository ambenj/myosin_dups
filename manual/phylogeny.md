# Myosin phylogeny
*Methods for phylogenetic analysis of MYH3C genes*

**Relevant Figure**: Figure 2B

Run phylogenetic analysis on CDS nucleotide alignment (MAFFT):
```bash
# Load iqtree (IQ-TREE multicore version 2.3.6)
conda activate tree

# Run IQ-TREE
sbatch scripts/phylogeny/iqtree.sh /labs/kingsley/ambenj/myosin_dups/analysis/phylogeny/cds/All_CDS_w3outgroups_mafft_alignment.fasta
```