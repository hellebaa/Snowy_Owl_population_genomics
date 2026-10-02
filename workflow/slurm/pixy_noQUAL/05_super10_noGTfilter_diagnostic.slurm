#!/bin/bash
#SBATCH --job-name=pixy_s10_noGTfilter
#SBATCH --account=nn8013k
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=pixy_s10_noGTfilter_%j.out
#SBATCH --error=pixy_s10_noGTfilter_%j.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

source /cluster/projects/nn9244k/olekto/miniconda3/etc/profile.d/conda.sh
conda activate pixy2
module load BCFtools/1.22-GCC-14.2.0

which bcftools
which tabix
which pixy
pixy --version

CHR="SUPER_10"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

INDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4"
RENAME="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4/01_snps/sample_rename.txt"
VARIANTS="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/SnowyOwl_autosomes_finalFilters_noContam.recalc.vcf.gz"
WORKDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL"
OLD_PIXY_DIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_new"
TESTDIR="${WORKDIR}/SUPER_10_diagnostic_noGTfilter"

mkdir -p "${TESTDIR}/pixy_all"
cd "${TESTDIR}"

echo "================================="
echo "1. Extract ${CHR} invariant sites with no QUAL filter and no genotype-level DP/RGQ filter"
echo "================================="

bcftools reheader \
  -s "${RENAME}" \
  "${INDIR}/${CHR}.vcf.gz" \
| bcftools view \
  -s ^GRL04_merged_X14,SKW02_merged_X22 \
  -i 'ALT="."' \
  -Oz \
  --threads "${THREADS}" \
  -o "${CHR}.invariant.noQUAL.noGTfilter.vcf.gz"

tabix -f -p vcf "${CHR}.invariant.noQUAL.noGTfilter.vcf.gz"

echo -n "Invariant records without genotype filter: "
bcftools index -n "${CHR}.invariant.noQUAL.noGTfilter.vcf.gz"

echo "================================="
echo "2. Extract ${CHR} variant sites from final SNP VCF"
echo "================================="

bcftools view \
  -r "${CHR}" \
  "${VARIANTS}" \
  -Oz \
  --threads "${THREADS}" \
  -o "${CHR}.variants.finalFilters.noContam.vcf.gz"

tabix -f -p vcf "${CHR}.variants.finalFilters.noContam.vcf.gz"

echo -n "Variant records: "
bcftools index -n "${CHR}.variants.finalFilters.noContam.vcf.gz"

echo "================================="
echo "3. Concatenate ${CHR} invariant and variant sites"
echo "================================="

cat > "${CHR}.vcf_list.noGTfilter.txt" << EOF
${TESTDIR}/${CHR}.invariant.noQUAL.noGTfilter.vcf.gz
${TESTDIR}/${CHR}.variants.finalFilters.noContam.vcf.gz
EOF

bcftools concat \
  -a \
  -f "${CHR}.vcf_list.noGTfilter.txt" \
  -Oz \
  --threads "${THREADS}" \
  -o "${CHR}.allSites.noQUAL.noGTfilter.vcf.gz"

tabix -f -p vcf "${CHR}.allSites.noQUAL.noGTfilter.vcf.gz"

echo -n "All-sites records: "
bcftools index -n "${CHR}.allSites.noQUAL.noGTfilter.vcf.gz"

echo "================================="
echo "4. Run pixy ALL on ${CHR} diagnostic all-sites VCF"
echo "================================="

POPFILE="${TESTDIR}/population_all.txt"
if [[ ! -s "${POPFILE}" ]]; then
  cp "${OLD_PIXY_DIR}/population_all.txt" "${POPFILE}"
fi

pixy \
  --stats pi tajima_d \
  --vcf "${CHR}.allSites.noQUAL.noGTfilter.vcf.gz" \
  --populations "${POPFILE}" \
  --window_size 50000 \
  --n_cores 8 \
  --output_folder "${TESTDIR}/pixy_all" \
  --output_prefix "SnowyOwl_ALL_${CHR}_noQUAL_noGTfilter"

echo "================================="
echo "5. Summarize pixy no_sites for diagnostic run"
echo "================================="

awk '
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="no_sites") s=i;
    if ($i=="count_missing") m=i;
    if ($i=="count_comparisons") c=i;
  }
}
NR>1 {
  windows++;
  sites += $s;
  missing += $m;
  comps += $c;
  if ($s==0) zero++;
}
END {
  print "windows", windows;
  print "zero_site_windows", zero;
  print "callable_sites_sum", sites;
  print "count_missing_sum", missing;
  print "count_comparisons_sum", comps;
}
' "${TESTDIR}/pixy_all/SnowyOwl_ALL_${CHR}_noQUAL_noGTfilter_pi.txt"

echo "DONE"
