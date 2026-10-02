#!/bin/bash
#SBATCH --job-name=owl_snp_noMAF
#SBATCH --account=nn8013k
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_snp_noMAF_%j.out
#SBATCH --error=owl_snp_noMAF_%j.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

# Genotype-filtered autosomal SNP VCF before the final MAF >= 0.02 site filter.
# This is the preferred SNP input for pixy diversity/Tajima's D because rare
# segregating sites must be retained.
SOURCE_VCF="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/autosomes.step5.GTfiltered.vcf.gz"
OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL_noMAF"
OUTVCF="${OUTDIR}/SnowyOwl_autosomes_snps_noContam_missing20_noMAF.recalc.vcf.gz"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

test -s "${SOURCE_VCF}"
test -s "${SOURCE_VCF}.tbi" || test -s "${SOURCE_VCF}.csi"

echo "Source SNP VCF:"
echo "${SOURCE_VCF}"
echo -n "Source records: "
bcftools index -n "${SOURCE_VCF}"

echo "Creating SNP VCF without MAF filtering"
bcftools view \
  --threads "${THREADS}" \
  -s ^GRL04_merged_X14,SKW02_merged_X22 \
  -m2 -M2 -v snps \
  -Ou \
  "${SOURCE_VCF}" \
| bcftools +fill-tags \
  -Ou \
  -- -t AC,AN,MAF,F_MISSING \
| bcftools view \
  --threads "${THREADS}" \
  -i 'F_MISSING<=0.20 && AC>0 && AC<AN' \
  -Ou \
| bcftools +fill-tags \
  -Ou \
  -- -t AC,AN,MAF,F_MISSING \
| bcftools view \
  --threads "${THREADS}" \
  -Oz \
  -o "${OUTVCF}"

tabix -f -p vcf "${OUTVCF}"

echo -n "No-MAF SNP records retained: "
bcftools index -n "${OUTVCF}"

echo "DONE"
