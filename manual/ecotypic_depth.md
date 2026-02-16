# Ecotypic depth analysis

*Methods for generating simulated data, mapping short read data, and analyzing read depth in global stickleback populations*

**Relevant Figures**:
* Figures 1B-D
* Figure S1
* Figure S2
* Table S1

## Read depth from Roberts Kingman et al. (2021) global samples
We used the global whole genome sequence dataset described in Roberts Kingman et al. (2021) which was mapped to *gasAcu1-4*. Unmapped reads can also be found at NCBI SRA: PRJNA247503.

### Process and map reads
Index *stickleback v. 5* reference genome (with masked duplication region 2 and no chrY):
```
sbatch 00_index_gasAcu5_C4masked.sh
```

Split bams by read group:
```
for i in /labs/kingsley/shared_data/stickleback_data/227_genomes_gasAcu1-4/bams/*_*_*.recal.bam; do sbatch 01_split_bamRG.sh ${i}; done
```

Realign bams to *stickleback v. 5* reference genome (with masked duplication region 2 and no ChrY):
```
# Realign bams for 206 split bams
for i in 01_split_bams/*.bam; do sbatch 02_realign.sh $i; done

# Realign bams for 21 original genomes (no split needed as these to not have multiple read groups)
for i in /labs/kingsley/shared_data/stickleback_data/227_genomes_gasAcu1-4/bams/????.recal.bam /labs/kingsley/shared_data/stickleback_data/227_genomes_gasAcu1-4/bams/???.recal.bam /labs/kingsley/shared_data/stickleback_data/227_genomes_gasAcu1-4/bams/????_?.recal.bam; do sbatch 02_realign.sh $i; done
```

Assign read groups and mark duplicates:
```
# Start interactive session
sdev

# Get list of read group data
module load samtools/1.19
for i in 01_split_bams/*_*_*.recal_*.bam; do printf $i; printf "\t"; samtools view -H $i | grep '^@RG' - ; done > 206genomes_read_group_data.txt

# Prepare commands for individual jobs
python2 03_mark_dups_add_rgs.py 206genomes_read_group_data.txt > 03_mark_dups_add_rgs_206genomes.sh

# Make output folder
mkdir 03_markdups_addRG

# Run mark duplicates and assign read groups
bash 03_mark_dups_add_rgs_206genomes.sh
```

Merge separate read group files into one per sample:
```
# Merge bams (for 206 genomes)
for i in /labs/kingsley/shared_data/stickleback_data/227_genomes_gasAcu1-4/bams/*_*_*.recal.bam; do FILE=${i##*/}; BASE=${FILE%.bam}; sbatch 04_merge_bams.sh $BASE; done
```

Mark duplicates and index:
```
# Mark dups and index (for 206 genomes)
for i in 04_merge_bams/*.bam; do sbatch 05_mkdup_index.sh $i; done

# Mark dups and index (for 21 original genomes)
for i in 02_realign/????.recal.realignGA5_C4masked.sort.bam 02_realign/???.recal.realignGA5_C4masked.sort.bam 02_realign/????_?.recal.realignGA5_C4masked.sort.bam; do sbatch 05_mkdup_index.sh $i; done
```
### Read depth analysis
Calculate read depth in regions of interest (i.e. *MYH3C* region):
```
# Calculate roi read depth for each sample
for i in 05_mkdup_index/*.bam; do sbatch 06_samtools_coverage_roi.sh $i /labs/kingsley/ambenj/myosin_dups/analysis/ecotypic_depth/stickleback_v5_ncbi_myosin_annotations_noChrY.bed 06_samtools_coverage_roi; done

# Gather roi results in one file
cd 06_samtools_coverage_roi
head -n 1 ZERO_X_2011_05.recal.realignGA5_C4masked.sort.merged.mkdup_mapq3_roi_coverage.txt > roi_coverage_mapq3.txt
for f in *mapq3_roi_coverage.txt; do tail -n +2 $f >> roi_coverage_mapq3.txt; done
```
Calculate read depth for whole genome:
```
for i in 05_mkdup_index/*.bam; do sbatch 06_samtools_coverage_wg.sh $i 06_samtools_coverage_wg; done
```


## Read depth from simulated genomes

### Prepare *MYH3C* region modified reference genomes
Make freshwater modified reference genomes:
```
# Get stickleback_v5 chrY and chrXIX
module load samtools/1.19
samtools faidx /labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly.fa chrXIX > stickleback_v5_chrXIX.fa
samtools faidx /labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly.fa chrY > stickleback_v5_chrY.fa

# remove chrXIX and chrY from stickleback_v5
module load seqkit/2.3.1
seqkit grep -vrp "^chrXIX|^chrY" /labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/stickleback_v5_assembly.fa > stickleback_v5_no_chrXIX_chrY.fa
```
```
# Make modified stickleback_v5 chrXIX containing 5-copy BEPA_XAB031 sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/BEPA_XAB031_chrXIX_h2tg000056l_MYHregion_RC.fasta stickleback_v5-BEPA_XAB031-5copy_chrXIX.fa

# Make modified stickleback_v5 chrXIX containing 5-copy BLAU22_12 sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/BLAU22_12_chrXIX_h1tg000036l_MYHregion_RC.fasta stickleback_v5-BLAU22_12-6copy_chrXIX.fa

# Make modified stickleback_v5 chrXIX containing 4-copy BOULxBDGB sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/BOULxBDGB_chrXIX_h1tg000075l_MYHregion.fasta stickleback_v5-BOULxBDGB-4copy_chrXIX.fa

# Make modified stickleback_v5 chrXIX containing 5-copy CMCB21_08 sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/CMCB21_08_chrXIX_h1tg000028l_MYHregion.fasta stickleback_v5-CMCB21_08-5copy_chrXIX.fa

# Make modified stickleback_v5 chrXIX containing 5-copy DNSE21_15 sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/DNSE21_15_chrXIX_h1tg000257l_MYHregion_RC.fasta stickleback_v5-DNSE21_15-5copy_chrXIX.fa

# Make modified stickleback_v5 chrXIX containing 5-copy ECHO21_44 sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/ECHO21_44_chrXIX_h1tg000119l_MYHregion_RC.fasta stickleback_v5-ECHO21_44-5copy_chrXIX.fa

# Make modified stickleback_v5 chrXIX containing 4-copy fGasAcu3 sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/fGasAcu3_chrX_OZ193841.1_MYHregion.fasta stickleback_v5-fGasAcu3-4copy_chrXIX.fa

# Make modified stickleback_v5 chrXIX containing 5-copy JADE21_11 sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/JADE21_11_chrXIX_h2tg000028l_MYHregion_RC.fasta stickleback_v5-JADE21_11-5copy_chrXIX.fa

# Make modified stickleback_v5 chrXIX containing 5-copy KFSY21_29 sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py chrXIX 2601499 2730636 stickleback_v5_chrXIX.fa ../nlrc5-syt19_regions/KFSY21_29_chrXIX_h1tg000035l_MYHregion_RC.fasta stickleback_v5-KFSY21_29-5copy_chrXIX.fa
```
```
# Make female genome reference for stickleback_v5 with freshwater assembly MYH regions (Note: fasta has two copies of each chr so should simulate half desired depth)
for s in stickleback_v5-BEPA_XAB031-5copy stickleback_v5-DNSE21_15-5copy stickleback_v5-JADE21_11-5copy stickleback_v5-BLAU22_12-6copy stickleback_v5-ECHO21_44-5copy stickleback_v5-KFSY21_29-5copy stickleback_v5-BOULxBDGB-4copy stickleback_v5-CMCB21_08-5copy stickleback_v5-fGasAcu3-4copy; do cat stickleback_v5_no_chrXIX_chrY.fa stickleback_v5_no_chrXIX_chrY.fa ${s}_chrXIX.fa ${s}_chrXIX.fa > ../${s}.female.fa; done
```
Make marine modified references:
```
# Get Duke GAcu_1.0 (RABS) chr19
samtools faidx /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/RABS_Duke_GAcu_1.0/GCA_046562415.1_Duke_GAcu_1.0_genomic.fasta CM102094.1 > Duke_GAcu_1.0_chr19.fa

# remove chr19 from Duke GAcu_1.0 (RABS)
seqkit grep -vrp "^CM102094.1" /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/RABS_Duke_GAcu_1.0/GCA_046562415.1_Duke_GAcu_1.0_genomic.fasta > Duke_GAcu_1.0_no_chr19.fa
```
```
# Make modified Duke GAcu_1.0 (RABS) chrXIX containing SALR sequence in myosin region
python3 ../../../../scripts/depth_simulations/replace_seq.py CM102094.1 2756462 2854692 Duke_GAcu_1.0_chr19.fa ../nlrc5-syt19_regions/SALR_10367_final_Contig6_226555-bp_MYHregion_RC.fasta Duke_GAcu_1.0-SALR-3copy_chr19.fa
```
```
# Make female genome reference for Duke_GAcu_1.0 with marine assembly MYH regions (Note: fasta has two copies of each chr so should simulate half desired depth)
cat Duke_GAcu_1.0_no_chr19.fa Duke_GAcu_1.0_no_chr19.fa Duke_GAcu_1.0_chr19.fa Duke_GAcu_1.0_chr19.fa > ../Duke_GAcu_1.0-RABS-3copy.female.fa

cat Duke_GAcu_1.0_no_chr19.fa Duke_GAcu_1.0_no_chr19.fa Duke_GAcu_1.0-SALR-3copy_chr19.fa Duke_GAcu_1.0-SALR-3copy_chr19.fa > ../Duke_GAcu_1.0-SALR-3copy.female.fa
```
### Simulate short read data from modified references and align
Simulate short reads from assemblies and map to stickleback_v5 (duplication region 2 masked):
```
cd /labs/kingsley/ambenj/myosin_dups
module load snakemake
snakemake --configfile configs/206_genomes_stickleback_v5masked_noChrY_sim_and_align.yaml --snakefile scripts/depth_simulations/sim_and_align.smk --profile scg --jobs 100 --restart-times 0 --rerun-incomplete
```

### Read depth analysis
Gather roi results in one file:
```
head -n 1 Duke_GAcu_1.0-RABS-3copy.female.fa_sim4X.rep0.mapq3_roi_coverage.txt > roi_coverage_mapq3.txt
for f in *mapq3_roi_coverage.txt; do tail -n +2 $f >> roi_coverage_mapq3.txt; done
```

## Neutral region phylogenetic analysis
Source files: 
* Vcf from Roberts Kingman et al. 2021: `227_genomes.final.filtered.vcf`
* Sensitive Pacific EcoPeaks from Roberts Kingman et al. 2021: `c150.sensitive.50kb.final.peaks`
* Specific Pacific EcoPeaks from Roberts Kingman et al. 2021: `c155.sensitive.50kb.final.peaks`
* Text file with names of samples in Table S1 global comparison group: `file_names.txt`

Combine EcoPeak files and remove extraneous data:
```
cat c150.sensitive.50kb.final.peaks c155.sensitive.50kb.final.peaks | cut -f1-3 - | sort -k1,1 -k2,2n - | bedtools merge -i - > combined_ecopeaks.bed
```

Make vcf subset with only selected samples and with MAF > 0.02:
```
bcftools view 227_genomes.final.filtered.vcf -S file_names.txt -i 'MAF>=0.02' -o subset_samples.vcf 
```

Define EcoPeak exclusion regions (5 Mb exclusion from any EcoPeak):
```
bedtools slop -i combined_ecopeaks.bed -b 5000000 -g  gasAcu1-4.fa.fai | bedtools merge -i - > eco_exclusion.bed
bedtools intersect -v -header -a subset_samples.vcf -b eco_exclusion.bed > eco_exclusion.vcf
```

Randomly subsample to 100k SNPs:
``` 
bcftools view -H eco_exclusion.vcf | shuf -n 100000 - > eco.tmp.vcf
bcftools view -h eco_exclusion.vcf | cat - eco.tmp.vcf > subsample.eco.vcf
```
Convert from VCF to PHYLIP (vcf2phylip by Ortiz 2019):
```
python3 vcf2phylip.py -i subsample.eco.vcf -o eco.100k.phy
```

Build tree (iqtree by Wong et al 2025 and Minh et al 2020):
```
iqtree3 -s subsample.eco.min4.phy -m GTR+G -nt AUTO -bb 1000
```
## Normalize read depth and plot results

Normalize and plot read depths:  
&nbsp;&nbsp;&nbsp;&nbsp;R notebook: [scripts/R/ecotypic_read_depth_sticklebackv5_noChrY_final.Rmd](../scripts/R/ecotypic_read_depth_sticklebackv5_noChrY_final.Rmd)  
&nbsp;&nbsp;&nbsp;&nbsp;Html output: [scripts/R/ecotypic_read_depth_sticklebackv5_noChrY_final.html](../scripts/R/ecotypic_read_depth_sticklebackv5_noChrY_final.html)