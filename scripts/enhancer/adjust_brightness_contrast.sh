#!/bin/bash
#SBATCH --job-name=image
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --account=kingsley
#SBATCH --partition=batch

# Adjust image and contrast levels for photos

module load imagemagick/7.0.5-7

cd /labs/kingsley/ambenj/myosin_dups/analysis/enhancer_assay/images_for_pub
mkdir -p b25c35adj


for j in *.JPG
do 
  convert -brightness-contrast 25x35 "$j" b25c35adj/b25c35adj_"$j"
done