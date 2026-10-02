#!/bin/bash
#SBATCH --job-name=owl_pixy_concat_noQUAL
#SBATCH --account=nn8013k
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_pixy_concat_noQUAL_%j.out
#SBATCH --error=owl_pixy_concat_noQUAL_%j.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

INV_DIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4/02_non_variant/06_invariant_noQUAL"
VARIANTS="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/SnowyOwl_autosomes_finalFilters_noContam.recalc.vcf.gz"
OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "Checking invariant chromosome files"
for i in $(seq 1 32)
do
  test -s "${INV_DIR}/SUPER_${i}.invariant.noQUAL.vcf.gz"
  test -s "${INV_DIR}/SUPER_${i}.invariant.noQUAL.vcf.gz.tbi"
done

echo "Creating autosomal invariant VCF list"
for i in $(seq 1 32)
do
  echo "${INV_DIR}/SUPER_${i}.invariant.noQUAL.vcf.gz"
done > invariant_autosomes_noQUAL.list

echo "Concatenating autosomal invariant VCF"
bcftools concat \
  -f invariant_autosomes_noQUAL.list \
  -Oz \
  --threads "${THREADS}" \
  -o SnowyOwl_autosomes_invariant_noQUAL.vcf.gz

tabix -f -p vcf SnowyOwl_autosomes_invariant_noQUAL.vcf.gz

echo "Preparing all-sites VCF list"
cat > vcf_list_noQUAL.txt << EOF
${OUTDIR}/SnowyOwl_autosomes_invariant_noQUAL.vcf.gz
${VARIANTS}
EOF

echo "Concatenating invariant and variant sites"
bcftools concat \
  -a \
  -f vcf_list_noQUAL.txt \
  -Oz \
  --threads "${THREADS}" \
  -o SnowyOwl_autosomes_allSites_noContam_noQUAL.vcf.gz

tabix -f -p vcf SnowyOwl_autosomes_allSites_noContam_noQUAL.vcf.gz

echo "Counts"
echo -n "Invariant records: "
bcftools index -n SnowyOwl_autosomes_invariant_noQUAL.vcf.gz
echo -n "Variant records: "
bcftools index -n "${VARIANTS}"
echo -n "All-sites records: "
bcftools index -n SnowyOwl_autosomes_allSites_noContam_noQUAL.vcf.gz

echo "DONE"
