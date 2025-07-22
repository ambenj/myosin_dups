#!/bin/bash
#SBATCH --job-name=liftover
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --account=kingsley

# Liftover tracks from gasAcu1-4 to stickleback v5
# Run with:
# sbatch liftover.sh

OLD_DIR="/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/gasAcu1-4"
NEW_DIR="/labs/kingsley/ambenj/myosin_dups/analysis/assemblies/stickleback_v5/gasAcu1-4_liftover"
CHAIN="${OLD_DIR}/v4_to_v5.chain.txt"

module load ucsc_tools/469

# liftOver ${OLD_DIR}/global_sensitive_EcoPeaks.bed ${CHAIN} ${NEW_DIR}/global_sensitive_EcoPeaks_lifted.bed ${NEW_DIR}/global_sensitive_EcoPeaks_unlifted.bed
# liftOver ${OLD_DIR}/global_specific_EcoPeaks.bed ${CHAIN} ${NEW_DIR}/global_specific_EcoPeaks_lifted.bed ${NEW_DIR}/global_specific_EcoPeaks_unlifted.bed
# liftOver ${OLD_DIR}/CH_SC_LB_sensitive_TempoPeaks.bed ${CHAIN} ${NEW_DIR}/CH_SC_LB_sensitive_TempoPeaks_lifted.bed ${NEW_DIR}/CH_SC_LB_sensitive_TempoPeaks_unlifted.bed
# liftOver ${OLD_DIR}/CH_SC_LB_specific_TempoPeaks.bed ${CHAIN} ${NEW_DIR}/CH_SC_LB_specific_TempoPeaks_lifted.bed ${NEW_DIR}/CH_SC_LB_specific_TempoPeaks_unlifted.bed
# liftOver ${OLD_DIR}/gasAcu1-4_gapTrack.bed ${CHAIN} ${NEW_DIR}/gasAcu1-4_gapTrack_lifted.bed ${NEW_DIR}/gasAcu1-4_gapTrack_unlifted.bed
liftOver ${OLD_DIR}/CH_SC_LB_specific_TempoPeaks_MYHregion_ends.bed ${CHAIN} ${NEW_DIR}/CH_SC_LB_specific_TempoPeaks_MYHregion_ends_lifted.bed ${NEW_DIR}/CH_SC_LB_specific_TempoPeaks_MYHregion_ends_unlifted.bed