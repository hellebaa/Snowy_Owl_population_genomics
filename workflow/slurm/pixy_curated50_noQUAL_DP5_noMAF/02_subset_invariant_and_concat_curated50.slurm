#!/bin/bash
#SBATCH --job-name=owl_concat_cur50
#SBATCH --account=nn8013k
#SBATCH --time=1-00:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_concat_cur50_%j.out
#SBATCH --error=owl_concat_cur50_%j.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_curated50_noQUAL_DP5_noMAF"
INV_SOURCE="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL/SnowyOwl_autosomes_invariant_noQUAL_DP5.vcf.gz"
INV_CURATED="${OUTDIR}/SnowyOwl_autosomes_invariant_curated50_noQUAL_DP5.vcf.gz"
VARIANTS="${OUTDIR}/SnowyOwl_autosomes_snps_curated50_missing20_noMAF.recalc.vcf.gz"
OUTVCF="${OUTDIR}/SnowyOwl_autosomes_allSites_curated50_noQUAL_DP5_noMAF.vcf.gz"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

test -s curated50_samples.txt
test -s "${INV_SOURCE}"
test -s "${INV_SOURCE}.tbi" || test -s "${INV_SOURCE}.csi"
test -s "${VARIANTS}"
test -s "${VARIANTS}.tbi" || test -s "${VARIANTS}.csi"

echo "Subsetting invariant VCF to curated 50 samples"
bcftools view \
  --threads "${THREADS}" \
  -S curated50_samples.txt \
  -Oz \
  -o "${INV_CURATED}" \
  "${INV_SOURCE}"

tabix -f -p vcf "${INV_CURATED}"

cat > vcf_list_curated50_noQUAL_DP5_noMAF.txt << EOF
${INV_CURATED}
${VARIANTS}
EOF

echo "Concatenating curated invariant and no-MAF variant sites"
bcftools concat \
  -a \
  -f vcf_list_curated50_noQUAL_DP5_noMAF.txt \
  -Oz \
  --threads "${THREADS}" \
  -o "${OUTVCF}"

tabix -f -p vcf "${OUTVCF}"

echo "Counts"
echo -n "Curated invariant records: "
bcftools index -n "${INV_CURATED}"
echo -n "Curated no-MAF variant records: "
bcftools index -n "${VARIANTS}"
echo -n "Curated all-sites records: "
bcftools index -n "${OUTVCF}"

echo "DONE"
