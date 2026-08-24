#!/bin/bash
#SBATCH --job-name=owl_pca_runB_nomac
#SBATCH --account=nn8013k
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_pca_runB_nomac_%j.out
#SBATCH --error=owl_pca_runB_nomac_%j.err

set -euo pipefail

module purge
module load StdEnv
module load PLINK/2.00a3.7-foss-2022a

########################################
# INPUT
########################################

VCF="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/autosomes.final.vcf.gz"

OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/pca_runB_contaminated_only_no_mac"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

########################################
# STEP 0: Remove contaminated / strong QC outlier samples only
########################################

cat <<EOF > remove_samples_runB_contaminated_only_no_mac.txt
GRL04_merged_X14 GRL04_merged_X14
SKW02_merged_X22 SKW02_merged_X22
EOF

########################################
echo "==== Step 1: Convert VCF to PLINK without MAC filtering ===="
########################################

plink \
  --vcf "${VCF}" \
  --double-id \
  --allow-extra-chr \
  --remove remove_samples_runB_contaminated_only_no_mac.txt \
  --biallelic-only strict \
  --set-missing-var-ids '@:#:$1:$2' \
  --make-bed \
  --out autosomes_clean_runB_no_mac

########################################
echo "==== Step 2: LD pruning ===="
########################################

plink \
  --bfile autosomes_clean_runB_no_mac \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out autosomes_clean_runB_no_mac_pruned

########################################
echo "==== Step 3: PCA ===="
########################################

plink \
  --bfile autosomes_clean_runB_no_mac \
  --allow-extra-chr \
  --extract autosomes_clean_runB_no_mac_pruned.prune.in \
  --pca 20 \
  --out autosomes_clean_runB_no_mac_pca

echo "==== DONE ===="
