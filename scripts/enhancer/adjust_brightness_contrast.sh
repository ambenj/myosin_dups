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
mkdir -p b40c20adj
mkdir -p b20c20adj
mkdir -p b10c20adj
mkdir -p b15c20adj
mkdir -p b15c15adj
mkdir -p b15c25adj


for j in *.JPG
do 
  #convert -brightness-contrast 40x20 "$j" b40c20adj/b40c20adj_"$j"
  #convert -brightness-contrast 20x20 "$j" b20c20adj/b20c20adj_"$j"
  #convert -brightness-contrast 10x20 "$j" b10c20adj/b10c20adj_"$j"
  #convert -brightness-contrast 15x20 "$j" b15c20adj/b15c20adj_"$j"
  #convert -brightness-contrast 15x15 "$j" b15c15adj/b15c15adj_"$j"
  convert -brightness-contrast 15x25 "$j" b15c25adj/b15c25adj_"$j"
done