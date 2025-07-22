#!/bin/bash
#SBATCH --job-name=vcftools
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# module load bedtools/2.27.1
module load vcftools/0.1.16-20-gd511f46

DIR="/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/gasAcu1-4"

# Filter to all variants in Sensitive global ecopeak overlapping myosin region
#vcftools --vcf /labs/kingsley/shared_data/stickleback_data/227_genomes_gasAcu1-4/vcfs/227_genomes.final.filtered.vcf --bed ${DIR}/global_sensitive_EcoPeak_MYH.bed --recode --out ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak

# Remove variants overlapping duplication region
vcftools --vcf ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak.recode.vcf --exclude-bed ${DIR}/MYH_duplication_regions.bed --recode --out ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak.noDup

# Calculate heterozygosity
vcftools --vcf ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak.noDup.recode.vcf --out ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak.noDup --het



# Filter to all variants in Sensitive global ecopeak overlapping myosin region
# bedtools intersect -a /labs/kingsley/shared_data/stickleback_data/227_genomes_gasAcu1-4/vcfs/227_genomes.final.filtered.vcf -b ${DIR}/global_sensitive_EcoPeak_MYH.bed > ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak.vcf

# # Remove variants overlapping duplication region
# bedtools intersect -v -a ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak.vcf -b ${DIR}/MYH_duplication_regions.bed > ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak.noDup.vcf

# # Calculate heterozygosity
# vcftools --het --vcf ${DIR}/227_genomes.final.filtered.MYHSensitiveEcopeak.noDup.vcf --out /labs/kingsley/ambenj/myosin_dups/analysis/assemblies/gasAcu1-4/227_genomes.final.filtered.MYHSensitiveEcopeak.noDup
