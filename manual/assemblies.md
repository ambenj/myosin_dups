# Long read assembly analysis

*Methods for *de novo* genome assembly from long reads and evaluation of assembly quality*

**Relevant Figures**: 
* Figure 2A
* Table S2

This section shows methods used to assemble, annotate, and quality check the assemblies that were used in the main figures. For assemblies that did not meet quality minimums for use in the paper, see [this section](./assemblies_unused.md).

## *De novo* assembly
Perform *de novo* assembly using Hifiasm:
```bash
conda activate myosin_assembly

# Start hifiasm assembly for BOULxBDGB
sbatch scripts/assemblies/hifiasm.sh /labs/kingsley/ambenj/myosin_dups/raw/sra/SRR16093180.fastq BOULxBDGB

# Start hifiasm assembly for BLAU22_12
sbatch scripts/assemblies/hifiasm.sh "/labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/BLAU22_12/BLAU22_12.AY/DEMUX_10091/demultiplex.bc2007--bc2007.hifi_reads.fastq.gz /labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/BLAU22_12/BLAU22_12.AY/DEMUX_10196/demultiplex.bc2007--bc2007.hifi_reads.fastq.gz" BLAU22_12

# Start hifiasm assembly for KFSY21_29
sbatch scripts/assemblies/hifiasm.sh /labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool1/Voxsimer_Pool01.HFSS/r84046_20240524_201927_1_C01/hifi_reads/m84046_240525_002538_s3.hifi_reads.bc2045.fastq.gz KFSY21_29

# Start hifiasm assembly for DNSE21_15
sbatch scripts/assemblies/hifiasm.sh /labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool1/Voxsimer_Pool01.HFSS/r84046_20240524_201927_1_C01/hifi_reads/m84046_240525_002538_s3.hifi_reads.bc2046.fastq.gz DNSE21_15

# Start hifiasm assembly for CMCB21_08
sbatch scripts/assemblies/hifiasm.sh /labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool1/Voxsimer_Pool01.HFSS/r84046_20240524_201927_1_C01/hifi_reads/m84046_240525_002538_s3.hifi_reads.bc2047.fastq.gz CMCB21_08

# Start hifiasm assembly for JADE21_11
sbatch scripts/assemblies/hifiasm.sh /labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool2/Voxsimer_Pool02.HFSS/r84046_20240605_185212_1_C01/hifi_reads/m84046_240605_230202_s1.hifi_reads.bc2023.fastq.gz JADE21_11

# Start hifiasm assembly for BEPA_XAB031
sbatch scripts/assemblies/hifiasm.sh /labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool2/Voxsimer_Pool02.HFSS/r84046_20240605_185212_1_C01/hifi_reads/m84046_240605_230202_s1.hifi_reads.bc2024.fastq.gz BEPA_XAB031

# Start hifiasm assembly for ECHO21_44
sbatch scripts/assemblies/hifiasm.sh /labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool2/Voxsimer_Pool02.HFSS/r84046_20240605_185212_1_C01/hifi_reads/m84046_240605_230202_s1.hifi_reads.bc2048.fastq.gz ECHO21_44
```

BLAT hifiasm assemblies with Ensembl *MYH3C* sequences to find contigs with *MYH3C* genes:
```bash
conda activate myosin_assembly
for f in analysis/assemblies/*/hifiasm/*.fa; do sbatch scripts/assemblies/blat.sh $f; done
```

## Annotate assemblies
Liftoff NCBI annotations and revised myosin annotations to assemblies:
```bash
conda activate myosin_annotate

for f in /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/SALR_CH213/SALR_10367_final_Contig6.fasta /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/*/hifiasm/*.hifiasm.myosin.fasta /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/fGasAcu3/GCA_964276395.1_fGasAcu3.hap1.1_genomic_chrX_chrY.fasta /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/RABS_Duke_GAcu_1.0/GCA_046562415.1_Duke_GAcu_1.0_genomic_chr19.fasta; do sbatch scripts/assemblies/liftoff_NCBI.sh $f; done
```

## Evaluation of assemblies
### Hifiasm assemblies generated from this paper:
Combine the two haplotypes into one assembly:
```bash
for s in BEPA_XAB031 BLAU22_12 BOULxBDGB CMCB21_08 DNSE21_15 ECHO21_44 JADE21_11 KFSY21_29; do cat analysis/assemblies/${s}/hifiasm/${s}.bp.hap?.p_ctg.fa > analysis/assemblies/${s}/hifiasm/${s}.bp.hap12.p_ctg.fa; done
```
Align HiFi reads to each assembly:
```bash
# Align BOULxBDGB pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BOULxBDGB/hifiasm/BOULxBDGB.bp.hap12.p_ctg.fa \
/labs/kingsley/ambenj/myosin_dups/raw/sra/SRR16093180.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BOULxBDGB/hifiasm/BOULxBDGB.bp.hap12.p_ctg.bam \
HIFI

# Align BLAU22_12 pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BLAU22_12/hifiasm/BLAU22_12.bp.hap12.p_ctg.fa \
/labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/BLAU22_12/BLAU22_12.AY/BLAU22_12_pacbio.fofn \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BLAU22_12/hifiasm/BLAU22_12.bp.hap12.p_ctg.bam \
HIFI

# Align KFSY21_29 pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/KFSY21_29/hifiasm/KFSY21_29.bp.hap12.p_ctg.fa \
/labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool1/Voxsimer_Pool01.HFSS/r84046_20240524_201927_1_C01/hifi_reads/m84046_240525_002538_s3.hifi_reads.bc2045.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/KFSY21_29/hifiasm/KFSY21_29.bp.hap12.p_ctg.bam \
HIFI

# Align DNSE21_15 pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/DNSE21_15/hifiasm/DNSE21_15.bp.hap12.p_ctg.fa \
/labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool1/Voxsimer_Pool01.HFSS/r84046_20240524_201927_1_C01/hifi_reads/m84046_240525_002538_s3.hifi_reads.bc2046.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/DNSE21_15/hifiasm/DNSE21_15.bp.hap12.p_ctg.bam \
HIFI

# Align CMCB21_08 pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/CMCB21_08/hifiasm/CMCB21_08.bp.hap12.p_ctg.fa \
/labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool1/Voxsimer_Pool01.HFSS/r84046_20240524_201927_1_C01/hifi_reads/m84046_240525_002538_s3.hifi_reads.bc2047.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/CMCB21_08/hifiasm/CMCB21_08.bp.hap12.p_ctg.bam \
HIFI

# Align JADE21_11 pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/JADE21_11/hifiasm/JADE21_11.bp.hap12.p_ctg.fa \
/labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool2/Voxsimer_Pool02.HFSS/r84046_20240605_185212_1_C01/hifi_reads/m84046_240605_230202_s1.hifi_reads.bc2023.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/JADE21_11/hifiasm/JADE21_11.bp.hap12.p_ctg.bam \
HIFI

# Align BEPA_XAB031 pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BEPA_XAB031/hifiasm/BEPA_XAB031.bp.hap12.p_ctg.fa \
/labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool2/Voxsimer_Pool02.HFSS/r84046_20240605_185212_1_C01/hifi_reads/m84046_240605_230202_s1.hifi_reads.bc2024.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/BEPA_XAB031/hifiasm/BEPA_XAB031.bp.hap12.p_ctg.bam \
HIFI

# Align ECHO21_44 pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/ECHO21_44/hifiasm/ECHO21_44.bp.hap12.p_ctg.fa \
/labs/kingsley/ambenj/myosin_dups/raw/pacbio_hifi/AK21_pool2/Voxsimer_Pool02.HFSS/r84046_20240605_185212_1_C01/hifi_reads/m84046_240605_230202_s1.hifi_reads.bc2048.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/ECHO21_44/hifiasm/ECHO21_44.bp.hap12.p_ctg.bam \
HIFI
```

### Publicly-available assemblies
Align long reads back to each genome assembly:
```bash
conda activate myosin_assembly

# Align RABS/Duke_GAcu_1.0 pacbio reads to self genome assembly
sbatch scripts/assemblies/pbmm2_align.sh  \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/RABS_Duke_GAcu_1.0/GCA_046562415.1_Duke_GAcu_1.0_genomic.fasta \
/labs/kingsley/ambenj/myosin_dups/raw/sra/Duke_GAcu_1.0_pacbio.fofn \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/RABS_Duke_GAcu_1.0/GCA_046562415.1_Duke_GAcu_1.0_genomic.bam \
SUBREAD

# Align fGasAcu3 pacbio reads to self genome assembly 
sbatch scripts/assemblies/pbmm2_align.sh \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/fGasAcu3/GCA_964276395.1_fGasAcu3.hap1.1_genomic.fasta \
/labs/kingsley/ambenj/myosin_dups/raw/sra/ERR11458793.fastq \
/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/fGasAcu3/GCA_964276395.1_fGasAcu3.hap1.1_genomic.bam \
HIFI
```

## Visualization
Plot diagram of representative copy number haplotypes:  
&nbsp;&nbsp;&nbsp;&nbsp;R notebook: [scripts/R/CNV_alleles_figure.Rmd](../scripts/R/CNV_alleles_figure.Rmd)    
&nbsp;&nbsp;&nbsp;&nbsp;Html output: [scripts/R/CNV_alleles_figure.html](../scripts/R/CNV_alleles_figure.html)  