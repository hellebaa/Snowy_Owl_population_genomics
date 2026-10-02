#!/bin/bash
#SBATCH --job-name=pixy_regions_noQUAL
#SBATCH --account=nn8013k
#SBATCH --time=2-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=pixy_regions_noQUAL_%j.out
#SBATCH --error=pixy_regions_noQUAL_%j.err

set -euo pipefail

source /cluster/projects/nn9244k/olekto/miniconda3/etc/profile.d/conda.sh
conda activate pixy2
module load BCFtools/1.22-GCC-14.2.0

which pixy
pixy --version
which tabix

WORKDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL"
OLD_PIXY_DIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_new"

VCF="${WORKDIR}/SnowyOwl_autosomes_allSites_noContam_noQUAL.vcf.gz"
POPFILE="${WORKDIR}/population_regions.txt"
OUTDIR="${WORKDIR}/02_pixy_regions"

mkdir -p "${OUTDIR}"

if [[ ! -s "${POPFILE}" ]]; then
  cp "${OLD_PIXY_DIR}/population_regions.txt" "${POPFILE}"
fi

pixy \
  --stats pi tajima_d fst dxy \
  --vcf "${VCF}" \
  --populations "${POPFILE}" \
  --window_size 50000 \
  --n_cores 8 \
  --output_folder "${OUTDIR}" \
  --output_prefix SnowyOwl_regions_noQUAL

echo "DONE"
