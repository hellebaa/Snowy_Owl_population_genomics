#!/bin/bash
#SBATCH --job-name=pixy_s10_sens
#SBATCH --account=nn8013k
#SBATCH --time=1-00:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --output=pixy_s10_sens_%j.out
#SBATCH --error=pixy_s10_sens_%j.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

source /cluster/projects/nn9244k/olekto/miniconda3/etc/profile.d/conda.sh
conda activate pixy2

which bcftools
which pixy
pixy --version
which tabix

CHR="SUPER_10"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

SOURCE_SNP_VCF="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/autosomes.step5.GTfiltered.vcf.gz"
SOURCE_INV_VCF="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL/SnowyOwl_autosomes_invariant_noQUAL_DP5.vcf.gz"
OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_SUPER10_sensitivity"

mkdir -p "${OUTDIR}"/{sample_lists,vcfs,pixy,summary}
cd "${OUTDIR}"

test -s "${SOURCE_SNP_VCF}"
test -s "${SOURCE_SNP_VCF}.tbi" || test -s "${SOURCE_SNP_VCF}.csi"
test -s "${SOURCE_INV_VCF}"
test -s "${SOURCE_INV_VCF}.tbi" || test -s "${SOURCE_INV_VCF}.csi"

echo "Preparing sample lists"

# Historical pixy sample set: all samples except the two contamination-remove samples.
bcftools query -l "${SOURCE_SNP_VCF}" \
  | grep -vx 'GRL04_merged_X14' \
  | grep -vx 'SKW02_merged_X22' \
  > sample_lists/all66.samples.txt

# Diagnostic set: remove only WRG01 from the historical pixy sample set.
grep -vx 'WRG01_1-WRG1_S28_L001' sample_lists/all66.samples.txt \
  > sample_lists/noWRG01_65.samples.txt

# Curated set: QC OK, contamination keep, and one representative from the FNM related group.
cat > sample_lists/curated50.samples.txt << 'EOF'
FNM01_9-X42_S17_L002
FNM02_10-X43_S20_L002
FNM03_11-X44_S23_L002
FNM04_12-X45_S26_L002
FNM05_13-X46_S29_L002
FNM06_14-X47_S32_L002
FNM07_15-X48_S35_L002
FNM08_merged_X1
FNM10_merged_X3
FNM12_5-X5_S119_L004
FNM13_merged_X6
FNM16_16-X49_S37_L002
FNM17_17-X50_S18_L002
GRL01_merged_X11
GRL02_merged_X12
GRL03_merged_X13
GRL06_merged_X16
GRL07_merged_X17
GRL08_merged_X18
GRL09_merged_X19
GRL10_merged_X20
NY01_DAG-BS-1_S1
NY02_DAG-BS-2_S2
NY03_DAG-BS-3_S3
NY04_DAG-BS-4_S4
NY05_DAG-BS-5_S5
NY06_DAG-BS-6_S6
NY07_DAG-BS-7_S7
NY08_DAG-BS-8_S8
NY09_DAG-BS-9_S9
NY10_DAG-BS-10_S10
NY11_DAG-BS-11_S11
NY12_DAG-BS-12_S12
NY13_DAG-BS-13_S13
SKW01_merged_X21
SKW04_merged_X24
SKW05_merged_X25
SKW06_merged_X26
SKW07_merged_X27
SKW08_merged_X28
SKW09_merged_X29
SKW12_merged_X32
SKW13_merged_X33
SKW14_merged_X34
SKW15_merged_X35
WRG06_21-X54_S30_L002
WRG09_6-WRG9_S15_L001
WRG10_7-WRG10_S19_L001
WRG32_18-X51_S21_L002
WRG46_19-X52_S24_L002
EOF

for sample_set in all66 noWRG01_65 curated50
do
  echo "${sample_set} sample count: $(wc -l < sample_lists/${sample_set}.samples.txt)"
done

cat > summary/run_manifest.tsv << 'EOF'
sample_set	variant_filter	description
all66	noMAF	Historical pixy sample set; no minor-allele-frequency filter.
all66	maf005	Historical pixy sample set; MAF >= 0.005.
all66	maf01	Historical pixy sample set; MAF >= 0.01.
all66	maf02	Historical pixy sample set; MAF >= 0.02.
noWRG01_65	noMAF	Historical pixy sample set excluding only WRG01; no minor-allele-frequency filter.
noWRG01_65	maf005	Historical pixy sample set excluding only WRG01; MAF >= 0.005.
noWRG01_65	maf01	Historical pixy sample set excluding only WRG01; MAF >= 0.01.
noWRG01_65	maf02	Historical pixy sample set excluding only WRG01; MAF >= 0.02.
curated50	noMAF	Curated QC-passing sample set; no minor-allele-frequency filter.
curated50	maf005	Curated QC-passing sample set; MAF >= 0.005.
curated50	maf01	Curated QC-passing sample set; MAF >= 0.01.
curated50	maf02	Curated QC-passing sample set; MAF >= 0.02.
EOF

filter_expr() {
  case "$1" in
    noMAF)
      echo 'F_MISSING<=0.20 && AC>0 && AC<AN'
      ;;
    maf005)
      echo 'F_MISSING<=0.20 && MAF>=0.005 && AC>0 && AC<AN'
      ;;
    maf01)
      echo 'F_MISSING<=0.20 && MAF>=0.01 && AC>0 && AC<AN'
      ;;
    maf02)
      echo 'F_MISSING<=0.20 && MAF>=0.02 && AC>0 && AC<AN'
      ;;
    *)
      echo "Unknown filter: $1" >&2
      exit 1
      ;;
  esac
}

echo -e "sample_set\tvariant_filter\tn_samples\tn_variant_records\tn_invariant_records\tn_all_sites_records\tpi_windows\tpi_callable_sites\tpi_weighted_mean\ttajima_windows\ttajima_callable_sites\ttajima_weighted_mean" \
  > summary/SUPER10_pixy_sensitivity_summary.tsv

for sample_set in all66 noWRG01_65 curated50
do
  sample_file="sample_lists/${sample_set}.samples.txt"
  n_samples="$(wc -l < "${sample_file}")"

  inv_vcf="vcfs/${CHR}.${sample_set}.invariant.noQUAL.DP5.vcf.gz"
  echo "Subsetting invariant ${CHR} for ${sample_set}"
  bcftools view \
    --threads "${THREADS}" \
    -r "${CHR}" \
    -S "${sample_file}" \
    -Oz \
    -o "${inv_vcf}" \
    "${SOURCE_INV_VCF}"
  tabix -f -p vcf "${inv_vcf}"
  n_inv="$(bcftools index -n "${inv_vcf}")"

  popfile="sample_lists/${sample_set}.population_all.txt"
  awk '{print $1 "\tALL"}' "${sample_file}" > "${popfile}"

  for variant_filter in noMAF maf005 maf01 maf02
  do
    expr="$(filter_expr "${variant_filter}")"
    prefix="${CHR}.${sample_set}.${variant_filter}"
    snp_vcf="vcfs/${prefix}.snps.vcf.gz"
    allsites_vcf="vcfs/${prefix}.allSites.vcf.gz"
    vcf_list="vcfs/${prefix}.vcf_list.txt"
    pixy_outdir="pixy/${sample_set}/${variant_filter}"

    mkdir -p "${pixy_outdir}"

    echo "Creating SNP VCF: ${prefix}"
    bcftools view \
      --threads "${THREADS}" \
      -r "${CHR}" \
      -S "${sample_file}" \
      -m2 -M2 -v snps \
      -Ou \
      "${SOURCE_SNP_VCF}" \
    | bcftools +fill-tags \
      -Ou \
      -- -t AC,AN,MAF,F_MISSING \
    | bcftools view \
      --threads "${THREADS}" \
      -i "${expr}" \
      -Ou \
    | bcftools +fill-tags \
      -Ou \
      -- -t AC,AN,MAF,F_MISSING \
    | bcftools view \
      --threads "${THREADS}" \
      -Oz \
      -o "${snp_vcf}"

    tabix -f -p vcf "${snp_vcf}"
    n_var="$(bcftools index -n "${snp_vcf}")"

    cat > "${vcf_list}" << EOF
${inv_vcf}
${snp_vcf}
EOF

    echo "Concatenating all-sites VCF: ${prefix}"
    bcftools concat \
      -a \
      -f "${vcf_list}" \
      -Oz \
      --threads "${THREADS}" \
      -o "${allsites_vcf}"
    tabix -f -p vcf "${allsites_vcf}"
    n_all="$(bcftools index -n "${allsites_vcf}")"

    echo "Running pixy: ${prefix}"
    pixy \
      --stats pi tajima_d \
      --vcf "${allsites_vcf}" \
      --populations "${popfile}" \
      --window_size 50000 \
      --n_cores "${THREADS}" \
      --output_folder "${pixy_outdir}" \
      --output_prefix "${prefix}"

    pi_file="${pixy_outdir}/${prefix}_pi.txt"
    tajima_file="${pixy_outdir}/${prefix}_tajima_d.txt"

    pi_summary="$(awk '
      BEGIN { FS = OFS = "\t" }
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if ($i == "avg_pi") pi = i
          if ($i == "no_sites") s = i
        }
      }
      NR > 1 && $pi != "NA" && $s > 0 {
        n++
        sites += $s
        w += $pi * $s
      }
      END {
        if (sites > 0) {
          print n, sites, w / sites
        } else {
          print 0, 0, "NA"
        }
      }
    ' "${pi_file}")"

    tajima_summary="$(awk '
      BEGIN { FS = OFS = "\t" }
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if ($i == "tajima_d") td = i
          if ($i == "no_sites") s = i
        }
      }
      NR > 1 && $td != "NA" && $s > 0 {
        n++
        sites += $s
        w += $td * $s
      }
      END {
        if (sites > 0) {
          print n, sites, w / sites
        } else {
          print 0, 0, "NA"
        }
      }
    ' "${tajima_file}")"

    echo -e "${sample_set}\t${variant_filter}\t${n_samples}\t${n_var}\t${n_inv}\t${n_all}\t${pi_summary}\t${tajima_summary}" \
      >> summary/SUPER10_pixy_sensitivity_summary.tsv
  done
done

echo "Summary:"
cat summary/SUPER10_pixy_sensitivity_summary.tsv

echo "DONE"
