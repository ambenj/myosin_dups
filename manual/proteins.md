# Protein sequence analysis

*Methods for protein sequence analysis*

**Relevant Figures:**:
* Figure 2C
* Table S3

## Convert alignments to divergent sites
* Protein sequences were aligned with MAFFT in Geneious
* Edited alignment file so that position of frameshift is represented by a “&” and all following positions are “-”

Get divergent sequences:
```
python3 scripts/protein_seqs/get_divergent_pos.py -i analysis/protein_msa/all_FSmod_mafft_align.fasta -o analysis/protein_msa/all_protein_FSmod_mafft_aligndivergent_sites.tsv
```

## Plot divergent sites
Plot divergent sites and perform tests for domain enrichment:  
&nbsp;&nbsp;&nbsp;&nbsp;R notebook: [scripts/R/protein_alignment_figure_final.Rmd](../scripts/R/protein_alignment_figure_final.Rmd)  
&nbsp;&nbsp;&nbsp;&nbsp;Html output: [scripts/R/protein_alignment_figure_final.html](../scripts/R/protein_alignment_figure_final.html)
