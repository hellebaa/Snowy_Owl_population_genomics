#!/bin/bash
#SBATCH --job-name=pixy_regions_cur50
#SBATCH --account=nn8013k
#SBATCH --time=2-00:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=pixy_regions_cur50_%j.out
#SBATCH --error=pixy_regions_cur50_%j.err

set -euo pipefail

source /cluster/projects/nn9244k/olekto/miniconda3/etc/profile.d/conda.sh
conda activate pixy2
module load BCFtools/1.22-GCC-14.2.0

which pixy
pixy --version
which tabix

WORKDIR="/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_curated50_noQUAL_DP5_noMAF"
VCF="${WORKDIR}/SnowyOwl_autosomes_allSites_curated50_noQUAL_DP5_noMAF.vcf.gz"
POPFILE="${WORKDIR}/population_regions_curated50.txt"
OUTDIR="${WORKDIR}/02_pixy_regions"

mkdir -p "${OUTDIR}"
cd "${WORKDIR}"

cat > "${POPFILE}" << 'EOF'
FNM01_9-X42_S17_L002	FNM
FNM02_10-X43_S20_L002	FNM
FNM03_11-X44_S23_L002	FNM
FNM04_12-X45_S26_L002	FNM
FNM05_13-X46_S29_L002	FNM
FNM06_14-X47_S32_L002	FNM
FNM07_15-X48_S35_L002	FNM
FNM08_merged_X1	FNM
FNM10_merged_X3	FNM
FNM12_5-X5_S119_L004	FNM
FNM13_merged_X6	FNM
FNM16_16-X49_S37_L002	FNM
FNM17_17-X50_S18_L002	FNM
GRL01_merged_X11	GRL
GRL02_merged_X12	GRL
GRL03_merged_X13	GRL
GRL06_merged_X16	GRL
GRL07_merged_X17	GRL
GRL08_merged_X18	GRL
GRL09_merged_X19	GRL
GRL10_merged_X20	GRL
NY01_DAG-BS-1_S1	NYS
NY02_DAG-BS-2_S2	NYS
NY03_DAG-BS-3_S3	NYS
NY04_DAG-BS-4_S4	NYS
NY05_DAG-BS-5_S5	NYS
NY06_DAG-BS-6_S6	NYS
NY07_DAG-BS-7_S7	NYS
NY08_DAG-BS-8_S8	NYS
NY09_DAG-BS-9_S9	NYS
NY10_DAG-BS-10_S10	NYS
NY11_DAG-BS-11_S11	NYS
NY12_DAG-BS-12_S12	NYS
NY13_DAG-BS-13_S13	NYS
SKW01_merged_X21	SKW
SKW04_merged_X24	SKW
SKW05_merged_X25	SKW
SKW06_merged_X26	SKW
SKW07_merged_X27	SKW
SKW08_merged_X28	SKW
SKW09_merged_X29	SKW
SKW12_merged_X32	SKW
SKW13_merged_X33	SKW
SKW14_merged_X34	SKW
SKW15_merged_X35	SKW
WRG06_21-X54_S30_L002	WRG
WRG09_6-WRG9_S15_L001	WRG
WRG10_7-WRG10_S19_L001	WRG
WRG32_18-X51_S21_L002	WRG
WRG46_19-X52_S24_L002	WRG
EOF

pixy \
  --stats pi tajima_d fst dxy \
  --vcf "${VCF}" \
  --populations "${POPFILE}" \
  --window_size 50000 \
  --n_cores 8 \
  --output_folder "${OUTDIR}" \
  --output_prefix SnowyOwl_regions_curated50_noQUAL_DP5_noMAF

echo "DONE"
