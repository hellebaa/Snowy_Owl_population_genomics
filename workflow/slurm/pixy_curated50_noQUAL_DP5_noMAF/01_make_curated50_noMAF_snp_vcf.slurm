#!/bin/bash
#SBATCH --job-name=owl_snp_cur50
#SBATCH --account=nn8013k
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_snp_cur50_%j.out
#SBATCH --error=owl_snp_cur50_%j.err

set -euo pipefail

module purge
module load BCFtools/1.22-GCC-14.2.0

SOURCE_VCF="/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/autosomes.step5.GTfiltered.vcf.gz"
OUTDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_curated50_noQUAL_DP5_noMAF"
OUTVCF="${OUTDIR}/SnowyOwl_autosomes_snps_curated50_missing20_noMAF.recalc.vcf.gz"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

cat > curated50_samples.txt << 'EOF'
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

echo "Curated sample count:"
wc -l curated50_samples.txt

test -s "${SOURCE_VCF}"
test -s "${SOURCE_VCF}.tbi" || test -s "${SOURCE_VCF}.csi"

echo "Source SNP VCF:"
echo "${SOURCE_VCF}"
echo -n "Source records: "
bcftools index -n "${SOURCE_VCF}"

echo "Creating curated no-MAF SNP VCF"
bcftools view \
  --threads "${THREADS}" \
  -S curated50_samples.txt \
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

echo -n "Curated no-MAF SNP records retained: "
bcftools index -n "${OUTVCF}"

echo "DONE"
