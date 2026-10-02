#!/bin/bash
#SBATCH --job-name=owl_concat_noMAF
#SBATCH --account=nn8013k
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_concat_noMAF_%j.out
#SBATCH --error=owl_concat_noMAF_%j.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

INV_VCF="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL/SnowyOwl_autosomes_invariant_noQUAL_DP5.vcf.gz"
VARIANTS="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL_noMAF/SnowyOwl_autosomes_snps_noContam_missing20_noMAF.recalc.vcf.gz"
OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL_noMAF"
OUTVCF="${OUTDIR}/SnowyOwl_autosomes_allSites_noContam_noQUAL_DP5_noMAF.vcf.gz"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

test -s "${INV_VCF}"
test -s "${INV_VCF}.tbi" || test -s "${INV_VCF}.csi"
test -s "${VARIANTS}"
test -s "${VARIANTS}.tbi" || test -s "${VARIANTS}.csi"

cat > vcf_list_noQUAL_DP5_noMAF.txt << EOF
${INV_VCF}
${VARIANTS}
EOF

echo "Concatenating invariant and no-MAF variant sites"
bcftools concat \
  -a \
  -f vcf_list_noQUAL_DP5_noMAF.txt \
  -Oz \
  --threads "${THREADS}" \
  -o "${OUTVCF}"

tabix -f -p vcf "${OUTVCF}"

echo "Counts"
echo -n "Invariant records: "
bcftools index -n "${INV_VCF}"
echo -n "No-MAF variant records: "
bcftools index -n "${VARIANTS}"
echo -n "All-sites records: "
bcftools index -n "${OUTVCF}"

echo "DONE"
