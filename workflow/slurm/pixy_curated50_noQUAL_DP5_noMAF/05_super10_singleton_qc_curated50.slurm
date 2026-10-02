#!/bin/bash
#SBATCH --job-name=owl_cur50_singQC
#SBATCH --account=nn8013k
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=owl_cur50_singQC_%j.out
#SBATCH --error=owl_cur50_singQC_%j.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

WORKDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_curated50_noQUAL_DP5_noMAF"
SNP_VCF="${WORKDIR}/SnowyOwl_autosomes_snps_curated50_missing20_noMAF.recalc.vcf.gz"

cd "${WORKDIR}"

bcftools view -r SUPER_10 \
  -i 'INFO/AC=1 || INFO/AN-INFO/AC=1' \
  -Ou "${SNP_VCF}" |
bcftools query -i 'GT="het"' \
  -f '[%CHROM\t%POS\t%SAMPLE\t%DP\t%GQ\t%AD\n]' \
  > singleton_quality_curated50.tsv

bcftools view -r SUPER_10 \
  -i 'INFO/AC/INFO/AN>=0.1 && INFO/AC/INFO/AN<=0.9' \
  -Ou "${SNP_VCF}" |
bcftools query -i 'GT="het"' \
  -f '[%CHROM\t%POS\t%SAMPLE\t%DP\t%GQ\t%AD\n]' \
  > common_het_quality_curated50.tsv

echo "Line counts"
wc -l singleton_quality_curated50.tsv common_het_quality_curated50.tsv

echo "DONE"
