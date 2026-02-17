# Allele-specific expression analysis
*Methods for processing RNA-seq and Kinnex reads and evaluating allele-specific expression*

**Relevant Figures**:
* Figure 3
* Figure S6
* Table S4

## RNA-seq preprocessing and *k*-mer counting
BLAT kmers to RABS reference to check for off-targets and variants:
```bash
# BLAT kmers to new RABS reference
sbatch scripts/ASE/blat.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/RABS_Duke_GAcu_1.0/GCA_046562415.1_Duke_GAcu_1.0_genomic.fasta analysis/ase/kmer_analysis/MYH_BEPA_RABS_27mers.fasta analysis/ase/kmer_analysis/blat/MYH_BEPA_RABS_27mers_RABSblat

# BLAT reverse complement kmers
sbatch scripts/ASE/blat.sh /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/RABS_Duke_GAcu_1.0/GCA_046562415.1_Duke_GAcu_1.0_genomic.fasta analysis/ase/kmer_analysis/MYH_BEPA_RABS_27mers_RC.fasta analysis/ase/kmer_analysis/blat/MYH_BEPA_RABS_27mers_RC_RABSblat
```
BLAT kmers to BEPA reference to check for off-targets and variants:
```bash
# BLAT kmers to stickleback_v5 assembly without chrY and with 4th copy masked
sbatch scripts/ASE/blat.sh /labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly_MYH3C4dup_hardmasked_noChrY.fa analysis/ase/kmer_analysis/MYH_BEPA_RABS_27mers.fasta analysis/ase/kmer_analysis/blat/MYH_BEPA_RABS_27mers_BEPAblat

# BLAT reverse complement kmers
sbatch scripts/ASE/blat.sh /labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly_MYH3C4dup_hardmasked_noChrY.fa analysis/ase/kmer_analysis/MYH_BEPA_RABS_27mers_RC.fasta analysis/ase/kmer_analysis/blat/MYH_BEPA_RABS_27mers_RC_BEPAblat
```

Trim and QC raw RNA-seq reads:
```bash
module load snakemake
snakemake --configfile configs/myosin_temp_ASE_RNA_preprocess.yaml --snakefile scripts/ASE/RNA_preprocess.smk --profile scg --jobs 300 --restart-times 0 --rerun-incomplete
```

Run *k*-mer counting script on all samples:
```bash
for s in 4 34 54 55-58 192 193-196 2 29 137-140 138 139-142 185-188 186 187-190 18 65-68 66 67-70 227-230 228 229-232 11 25 113-116 114 115-118 6 24 35-38 36 37-40 155-158 156 157-160 3 28 56 173-176 174 175-178 194 209-212 210 211-214; do sbatch scripts/ASE/count_kmers.sh analysis/ase/kmer_analysis/MYH_BEPA_RABS_27mers_RC.txt analysis/ase/RNA_preprocess/01_fastp_trim/${s}_R1.trimmed.fastq.gz analysis/ase/RNA_preprocess/01_fastp_trim/${s}_R2.trimmed.fastq.gz analysis/ase/kmer_analysis/kmer_counts/${s}; done
```
Compile the results into one file:
```bash
head -n 1 113-116_kmers_with_counts.tsv > 27mer_counts_all.txt
for f in *_kmers_with_counts.tsv; do tail -n +2 ${f} >> 27mer_counts_all.txt; done
```

## Kinnex read mapping and mapped read counting

Create concatenated RABS reference and stickleback_v5 modified to have 5-copy BEPA assemblies region (using same sequence from depth simulations):
```bash
cat /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/RABS_Duke_GAcu_1.0/GCA_046562415.1_Duke_GAcu_1.0_genomic.fasta \
/labs/kingsley/ambenj/myosin_dups/analysis/ecotypic_depth/depth_simulations/genome_subsets/stickleback_v5_no_chrXIX_chrY.fa \
/labs/kingsley/ambenj/myosin_dups/analysis/ecotypic_depth/depth_simulations/genome_subsets/stickleback_v5-BEPA_XAB031-5copy_chrXIX.fa > \
/labs/kingsley/ambenj/myosin_dups/analysis/ecotypic_depth/depth_simulations/genome_subsets/Duke_GAcu_1.0-RABS-3copy.stickleback_v5-BEPA_XAB031-5copy.female.fa
```

Index concatenated RABSxBEPA reference:
```bash
sbatch scripts/ASE/index_isoseq.sh /labs/kingsley/ambenj/myosin_dups/analysis/ecotypic_depth/depth_simulations/genome_subsets/Duke_GAcu_1.0-RABS-3copy.stickleback_v5-BEPA_XAB031-5copy.female.fa
```
Run script to map FLNC reads to concatenated RABSxBEPA reference:
```bash
for s in 113-116_F1_23wpf_3C_14_Ad 155-158_F1_23wpf_18C_03_Ad 228_F1_23wpf_18C_15_Ab 35-38_F1_23wpf_3C_01_Ad 114_F1_23wpf_3C_14_Ab 156_F1_23wpf_18C_03_Ab 229-232_F1_23wpf_18C_15_F 36_F1_23wpf_3C_01_Ab 115-118_F1_23wpf_3C_14_F 157-160_F1_23wpf_18C_03_F 24_F1_31dpf_18C_08 37-40_F1_23wpf_3C_01_RF 11_F1_31dpf_3C_11 227-230_F1_23wpf_18C_15_Ad 29_F1_31dpf_18C_13 6_F1_31dpf_3C_06; do sbatch scripts/ASE/map_flnc_intron10kb.sh $s /labs/kingsley/ambenj/myosin_dups/analysis/ecotypic_depth/depth_simulations/genome_subsets/Duke_GAcu_1.0-RABS-3copy.stickleback_v5-BEPA_XAB031-5copy.female.fa.isoseq.mmi align_RABS-3copy_BEPA-5copy; done
```
Run script to extract and count reads that span exons 3 to 35:
```bash
for s in 113-116_F1_23wpf_3C_14_Ad 155-158_F1_23wpf_18C_03_Ad 228_F1_23wpf_18C_15_Ab 35-38_F1_23wpf_3C_01_Ad 114_F1_23wpf_3C_14_Ab 156_F1_23wpf_18C_03_Ab 229-232_F1_23wpf_18C_15_F 36_F1_23wpf_3C_01_Ab 115-118_F1_23wpf_3C_14_F 157-160_F1_23wpf_18C_03_F 24_F1_31dpf_18C_08 37-40_F1_23wpf_3C_01_RF 11_F1_31dpf_3C_11 227-230_F1_23wpf_18C_15_Ad 29_F1_31dpf_18C_13 6_F1_31dpf_3C_06; do sbatch scripts/ASE/extract_spanning_reads.sh analysis/ase/flnc_mapping/align_RABS-3copy_BEPA-5copy/${s}.mapped.primary.mapq3.bam; done
```

Compile results into one file:
```bash
head -n 1 6_F1_31dpf_3C_06.mapped.primary.mapq3.counts.txt > exon3-35_counts_all.txt
for f in *counts.txt; do tail -n +2 ${f} >> exon3-35_counts_all.txt; done
```

## *K*-mer normalization and plotting
Normalize *k*-mers and plot both *k*-mer and kinnex results:  
&nbsp;&nbsp;&nbsp;&nbsp;R notebook: [scripts/R/temperature_ase_final.Rmd](../scripts/R/temperature_ase_final.Rmd)  
&nbsp;&nbsp;&nbsp;&nbsp;Html output: [scripts/R/temperature_ase_final.html](../scripts/R/temperature_ase_final.html)