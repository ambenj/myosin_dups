#!/bin/bash
#SBATCH --job-name=dl
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --account=default
#SBATCH --partition=interactive


#Run with:
#sbatch 00_download.sh <str_to_match>

module load rclone

rclone copy --progress --include $1 kingsley_gdrive:Data\ backup/Veeramah\ Pool\ Seq/ original_bams/
