#!/bin/bash
#SBATCH --job-name=pixy_ALL_noMAF
#SBATCH --account=nn8013k
#SBATCH --time=1-00:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --output=pixy_ALL_noMAF_%j.out
#SBATCH --error=pixy_ALL_noMAF_%j.err

set -euo pipefail

source /cluster/projects/nn9244k/olekto/miniconda3/etc/profile.d/conda.sh
conda activate pixy2
module load BCFtools/1.22-GCC-14.2.0

which pixy
pixy --version
which tabix

WORKDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL_noMAF"
OLD_PIXY_DIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_new"

VCF="${WORKDIR}/SnowyOwl_autosomes_allSites_noContam_noQUAL_DP5_noMAF.vcf.gz"
POPFILE="${WORKDIR}/population_all.txt"
OUTDIR="${WORKDIR}/01_pixy_all"

mkdir -p "${OUTDIR}"

if [[ ! -s "${POPFILE}" ]]; then
  cp "${OLD_PIXY_DIR}/population_all.txt" "${POPFILE}"
fi

pixy \
  --stats pi tajima_d \
  --vcf "${VCF}" \
  --populations "${POPFILE}" \
  --window_size 50000 \
  --n_cores 8 \
  --output_folder "${OUTDIR}" \
  --output_prefix SnowyOwl_ALL_noQUAL_DP5_noMAF

echo "DONE"
