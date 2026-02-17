# Rapid and repeated evolution of myosin copy number in threespine stickleback
[![DOI](https://zenodo.org/badge/992221541.svg)](https://doi.org/10.5281/zenodo.18666765)

This repository contains the code used for following paper:

>Yoxsimer AM, Daugherty RR, Hare EE, Chan YF, Jones FC, Kingman GAR, Offenberg EG, Howes TR, Zhang H, Pollen AA, Brady SD, Xie KT, Chen HI, Lowe CB, Au EH, Grimwood J, Schmutz J, Myers RM, Schluter D, Heins DC, Reyes ML, Baker JA, Jónsson B, Reimchen TE, Bell MA, Kingsley DM. Rapid and repeated evolution of myosin copy number in threespine stickleback. bioRxiv [Preprint]. 2025 Dec 25:2025.12.22.696110. doi: [10.64898/2025.12.22.696110](https://doi.org/10.64898/2025.12.22.696110).

*Note: these scripts were configured for analysis on the SCG Cluster at Stanford University and most would require adjusting paths and software in order to run elsewhere.*

Detailed instructions and commands for analysis are available as listed below:

| Section | Relevant Figures | Description |
|---|---|---|
|[Genome tracks](manual/genometracks.md)|Figure 1A| Methods for lifting over annotations and visualizing genome tracks|
|[Ecotypic depth](manual/ecotypic_depth.md)|Figures 1B-D<br>Figure S1<br>Figure S2<br>Table S1| Methods for generating simulated data, mapping short read data, and analyzing read depth in global stickleback populations|
|[Transplant pool-seq](manual/poolseq.md)|Figures 1E-F<br>Figure S3|Methods for mapping short read pool-seq data and analyzing read depth in transplant stickleback populations|
|[Assemblies](manual/assemblies.md) | Figure 2A<br>Table S2 | Methods for *de novo* genome assembly from long reads and evaluation of assembly quality |
|[Myosin phylogeny](manual/phylogeny.md)|Figure 2B| Methods for phylogenetic analysis of *MYH3C* genes|
|[Protein analysis](manual/proteins.md)|Figure 2C<br>Table S3| Methods for protein sequence analysis|
|[Allele-specific expression](manual/ase.md)|Figure 3<br>Figure S6<br>Table S4|Methods for processing RNA-seq and Kinnex reads and evaluating allele-specific expression|
|[Enhancer](manual/enhancer.md)|Figure 4|Methods for processing images and analyzing GFP intensity scores from transgenic fish|
|[NAHR analysis](manual/nahr.md)|Figure 5C|Methods for determining and visualizing NAHR breakpoints|
|[QTL](manual/qtl.md)|Figure S5|Methods for performing QTL analysis|
