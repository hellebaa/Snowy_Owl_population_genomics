#!/bin/bash
#SBATCH --job-name=owl_inv_DP5
#SBATCH --account=nn8013k
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-32
#SBATCH --output=owl_inv_DP5_%A_%a.out
#SBATCH --error=owl_inv_DP5_%A_%a.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

# Regenerate invariant sites for pixy without a site-level QUAL filter.
# Genotypes with depth < 5 are set to missing. RGQ is not used after the
# SUPER_10 sensitivity test showed it strongly reduced the denominator and
# inflated pi relative to noGTfilter and DP-only runs.

INDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4"
OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4/02_non_variant/07_invariant_noQUAL_DP5"
RENAME="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4/01_snps/sample_rename.txt"

CHR="SUPER_${SLURM_ARRAY_TASK_ID}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "Processing ${CHR}"

bcftools reheader \
  -s "${RENAME}" \
  "${INDIR}/${CHR}.vcf.gz" \
| bcftools view \
  -s ^GRL04_merged_X14,SKW02_merged_X22 \
  -i 'ALT="."' \
  -Ou \
| bcftools filter \
  -S . \
  -e 'FMT/DP<5' \
  -Oz \
  --threads "${THREADS}" \
  -o "${OUTDIR}/${CHR}.invariant.noQUAL.DP5.vcf.gz"

tabix -f -p vcf "${OUTDIR}/${CHR}.invariant.noQUAL.DP5.vcf.gz"

N=$(bcftools index -n "${OUTDIR}/${CHR}.invariant.noQUAL.DP5.vcf.gz")
echo "${CHR}: ${N} invariant records retained before pixy missing-data accounting"
echo -e "${CHR}\t${N}" > "${OUTDIR}/${CHR}.invariant.noQUAL.DP5.count.tsv"

echo "DONE ${CHR}"
