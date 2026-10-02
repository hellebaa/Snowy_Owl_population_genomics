#!/bin/bash
#SBATCH --job-name=pixy_s10_filter_cmp
#SBATCH --account=nn8013k
#SBATCH --time=1-00:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --output=pixy_s10_filter_cmp_%j.out
#SBATCH --error=pixy_s10_filter_cmp_%j.err

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
TESTDIR="${WORKDIR}/SUPER_10_filter_comparison"
POPFILE="${TESTDIR}/population_all.txt"

mkdir -p "${TESTDIR}"
cd "${TESTDIR}"

if [[ ! -s "${POPFILE}" ]]; then
  cp "${OLD_PIXY_DIR}/population_all.txt" "${POPFILE}"
fi

echo "================================="
echo "Extract ${CHR} variant sites once"
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

run_mode () {
  local mode="$1"
  local filter_expr="$2"

  echo "================================="
  echo "Mode: ${mode}"
  echo "Filter expression: ${filter_expr:-none}"
  echo "================================="

  if [[ -z "${filter_expr}" ]]; then
    bcftools reheader \
      -s "${RENAME}" \
      "${INDIR}/${CHR}.vcf.gz" \
    | bcftools view \
      -s ^GRL04_merged_X14,SKW02_merged_X22 \
      -i 'ALT="."' \
      -Oz \
      --threads "${THREADS}" \
      -o "${CHR}.invariant.${mode}.vcf.gz"
  else
    bcftools reheader \
      -s "${RENAME}" \
      "${INDIR}/${CHR}.vcf.gz" \
    | bcftools view \
      -s ^GRL04_merged_X14,SKW02_merged_X22 \
      -i 'ALT="."' \
      -Ou \
    | bcftools filter \
      -S . \
      -e "${filter_expr}" \
      -Oz \
      --threads "${THREADS}" \
      -o "${CHR}.invariant.${mode}.vcf.gz"
  fi

  tabix -f -p vcf "${CHR}.invariant.${mode}.vcf.gz"

  echo -n "Invariant records (${mode}): "
  bcftools index -n "${CHR}.invariant.${mode}.vcf.gz"

  cat > "${CHR}.vcf_list.${mode}.txt" << EOF
${TESTDIR}/${CHR}.invariant.${mode}.vcf.gz
${TESTDIR}/${CHR}.variants.finalFilters.noContam.vcf.gz
EOF

  bcftools concat \
    -a \
    -f "${CHR}.vcf_list.${mode}.txt" \
    -Oz \
    --threads "${THREADS}" \
    -o "${CHR}.allSites.${mode}.vcf.gz"

  tabix -f -p vcf "${CHR}.allSites.${mode}.vcf.gz"

  echo -n "All-sites records (${mode}): "
  bcftools index -n "${CHR}.allSites.${mode}.vcf.gz"

  mkdir -p "${TESTDIR}/pixy_${mode}"

  pixy \
    --stats pi tajima_d \
    --vcf "${CHR}.allSites.${mode}.vcf.gz" \
    --populations "${POPFILE}" \
    --window_size 50000 \
    --n_cores 8 \
    --output_folder "${TESTDIR}/pixy_${mode}" \
    --output_prefix "SnowyOwl_ALL_${CHR}_${mode}"

  awk -v mode="${mode}" '
  NR==1 {
    for (i=1; i<=NF; i++) {
      if ($i=="avg_pi") pi=i;
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
    if ($pi!="NA") {
      pi_weighted_num += $pi * $s;
      pi_weighted_den += $s;
    }
  }
  END {
    mean_pi = (pi_weighted_den > 0 ? pi_weighted_num / pi_weighted_den : "NA");
    print mode, windows, zero, sites, missing, comps, mean_pi;
  }
  ' "${TESTDIR}/pixy_${mode}/SnowyOwl_ALL_${CHR}_${mode}_pi.txt" \
    >> "${TESTDIR}/SUPER_10_filter_comparison_summary.tsv"
}

echo -e "mode\twindows\tzero_site_windows\tcallable_sites_sum\tcount_missing_sum\tcount_comparisons_sum\tweighted_mean_pi" \
  > "${TESTDIR}/SUPER_10_filter_comparison_summary.tsv"

run_mode "noGTfilter" ""
run_mode "DP5" 'FMT/DP<5'
run_mode "DP5_RGQ20" 'FMT/DP<5 || FMT/RGQ<20'

echo "================================="
echo "Summary"
echo "================================="
column -t "${TESTDIR}/SUPER_10_filter_comparison_summary.tsv" || cat "${TESTDIR}/SUPER_10_filter_comparison_summary.tsv"

echo "DONE"
