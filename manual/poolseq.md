# Transplant Pool-seq depth and SNP analysis

*Methods for mapping short read pool-seq data and analyzing read depth in transplant stickleback populations*

**Relevant Figures**:
* Figures 1E-F
* Figure S3

## Process and map reads
We used the transplant pool-seq dataset described in Roberts Kingman et al. (2021) which was mapped to *gasAcu1-4*. Unmapped reads can also be found at NCBI SRA: PRJNA671824.

Index *stickleback v. 5* reference genome (with masked duplication region 2):
```
sbatch 00_index_gasAcu5_C4masked.sh 
```

Split bams by read group:
```
for i in /labs/kingsley/ambenj/myosin/transplant_poolseq_depth/gasAcuv5_C4masked/original_bams/*.bam; do sbatch 01_split_bamRG.sh ${i}; done
```

Realign bams to *stickleback v. 5* reference genome (with masked duplication region 2):
```
for i in 01_split_bams/*.bam; do sbatch 02_realign.sh $i; done
```

Assign read groups and mark duplicates:
```
# Start interactive session
sdev

# Get list of read group data
module load samtools/1.19
for i in 01_split_bams/*.bam; do printf $i; printf "\t"; samtools view -H $i | grep '^@RG' - ; done > transplant_poolseq_read_group_data.txt

# Prepare commands for individual jobs
python2 03_mark_dups_add_rgs.py transplant_poolseq_read_group_data.txt > 03_mark_dups_add_rgs_transplant_poolseq.sh

# Make output folder
mkdir 03_markdups_addRG

# Run mark duplicates and assign read groups
bash 03_mark_dups_add_rgs_transplant_poolseq.sh
```
Merge separate read group files into one per sample:
```
# Merge bams
for i in original_bams/*.bam; do FILE=${i##*/}; BASE=${FILE%.bam}; sbatch 04_merge_bams.sh $BASE; done
```
Mark duplicates and index:
```
for i in 04_merge_bams/*.bam; do sbatch 05_mkdup_index.sh $i; done
```

## Read depth analysis
Calculate read depth in regions of interest (i.e. *MYH3C* region):
```
# Calculate roi read depth for each sample
for i in 05_mkdup_index/*.bam; do sbatch 06_samtools_coverage_roi.sh $i /labs/kingsley/ambenj/ref_genomes/Nath2020_sticklebackv5/annotations/stickleback_v5_myosin_manual_annotations_v2.bed 06_samtools_coverage_roi; done

# Gather roi results in one file
cd 06_samtools_coverage_roi
head -n 1 CH-2009.PE_ME_merged.sort.Mkdup.realign.BQrecal.realignGA5_C4masked.sort.merged.mkdup_mapq3_roi_coverage.txt > roi_coverage_mapq3.txt
for f in *mapq3_roi_coverage.txt; do tail -n +2 $f >> roi_coverage_mapq3.txt; done
```

Calculate read depth for whole genome:
```
for i in 05_mkdup_index/*.bam; do sbatch 06_samtools_coverage_wg.sh $i 06_samtools_coverage_wg; done
```

Normalize read depths; plot read depths and allele frequencies:  
&nbsp;&nbsp;&nbsp;&nbsp;R notebook: [scripts/R/transplant_read_depth_sticklebackv5_final.Rmd](../scripts/R/transplant_read_depth_sticklebackv5_final.Rmd)  
&nbsp;&nbsp;&nbsp;&nbsp;Html output: [scripts/R/transplant_read_depth_sticklebackv5_final.html](../scripts/R/transplant_read_depth_sticklebackv5_final.html)