#!/bin/bash
#SBATCH --job-name=owl_pca_runA_nomac
#SBATCH --account=nn8013k
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_pca_runA_nomac_%j.out
#SBATCH --error=owl_pca_runA_nomac_%j.err

set -euo pipefail

module purge
module load StdEnv
module load PLINK/2.00a3.7-foss-2022a

########################################
# INPUT
########################################

VCF="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/autosomes.final.vcf.gz"

OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/pca_runA_qc_relatedness_no_mac"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

########################################
# STEP 0: Define sample-removal lists
#
# Run A matches the final ANGSD PCA sample logic:
# qc_flags == OK, relatedness_flag == keep, contamination_flag == keep.
#
# LD pruning is estimated from the larger non-contaminated dataset, because
# PLINK refuses LD pruning when fewer than 50 samples remain. The stricter
# Run A sample filter is then applied at the PCA step.
########################################

cat <<EOF > remove_samples_runA_contaminated_only_no_mac.txt
GRL04_merged_X14 GRL04_merged_X14
SKW02_merged_X22 SKW02_merged_X22
EOF

cat <<EOF > remove_samples_runA_qc_relatedness_no_mac.txt
FNM09_merged_X2 FNM09_merged_X2
FNM11_merged_X4 FNM11_merged_X4
FNM12_5-X5_S119_L004 FNM12_5-X5_S119_L004
FNM15_merged_X8 FNM15_merged_X8
GRL04_merged_X14 GRL04_merged_X14
GRL05_merged_X15 GRL05_merged_X15
SKW02_merged_X22 SKW02_merged_X22
SKW03_merged_X23 SKW03_merged_X23
SKW10_merged_X30 SKW10_merged_X30
SKW11_merged_X31 SKW11_merged_X31
SKW16_merged_X36 SKW16_merged_X36
SKW17_merged_X37 SKW17_merged_X37
WRG01_1-WRG1_S28_L001 WRG01_1-WRG1_S28_L001
WRG02_2-WRG2_S33_L001 WRG02_2-WRG2_S33_L001
WRG03_3-WRG3_S39_L001 WRG03_3-WRG3_S39_L001
WRG05_4-WRG5_S4_L001 WRG05_4-WRG5_S4_L001
WRG08_5-WRG8_S10_L001 WRG08_5-WRG8_S10_L001
WRG11_8-WRG11_S24_L001 WRG11_8-WRG11_S24_L001
WRG48_20-X53_S27_L002 WRG48_20-X53_S27_L002
EOF

########################################
echo "==== Step 1: Convert VCF to PLINK without MAC filtering in larger non-contaminated dataset ===="
########################################

plink \
  --vcf "${VCF}" \
  --double-id \
  --allow-extra-chr \
  --remove remove_samples_runA_contaminated_only_no_mac.txt \
  --biallelic-only strict \
  --set-missing-var-ids '@:#:$1:$2' \
  --make-bed \
  --out autosomes_clean_runA_pruning_base_no_mac

########################################
echo "==== Step 2: LD pruning in larger non-contaminated dataset ===="
########################################

plink \
  --bfile autosomes_clean_runA_pruning_base_no_mac \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out autosomes_clean_runA_no_mac_pruned

########################################
echo "==== Step 3: PCA with strict Run A sample filter ===="
########################################

plink \
  --bfile autosomes_clean_runA_pruning_base_no_mac \
  --allow-extra-chr \
  --remove remove_samples_runA_qc_relatedness_no_mac.txt \
  --extract autosomes_clean_runA_no_mac_pruned.prune.in \
  --pca 20 \
  --out autosomes_clean_runA_no_mac_pca

echo "==== DONE ===="
