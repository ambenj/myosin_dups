# NAHR Breakpoint Analysis
*Methods for determining and visualizing NAHR breakpoints*

**Relevant Figures**: Figure 5C

Extract divergent positions from alignment:
```
python3 scripts/protein_seqs/get_divergent_pos.py -i analysis/nahr_analysis/5-6-copy_dup_mafft_alignment_no_gaps.fasta -o analysis/nahr_analysis/5-6-copy_dup_nogaps_divergent_sites.txt
```

Assign divergent positions to matching copy (A or L) and plot:  
&nbsp;&nbsp;&nbsp;&nbsp;R notebook: [scripts/R/NAHR_figure_final.Rmd](../scripts/R/NAHR_figure_final.Rmd)  
&nbsp;&nbsp;&nbsp;&nbsp;Html output: [scripts/R/NAHR_figure_final.html](../scripts/R/NAHR_figure_final.html)