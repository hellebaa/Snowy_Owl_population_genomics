#!/bin/bash
#SBATCH --job-name=owl_pca_gatk_full50
#SBATCH --account=nn8013k
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_pca_gatk_full50_%j.out
#SBATCH --error=owl_pca_gatk_full50_%j.err
set -euo pipefail
# Submit from the directory containing this script and its companion files.
INPUT_DIR="${PCA_INPUT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
KEEP="${INPUT_DIR}/curated50_samples.txt"
[[ -s "$KEEP" ]] || { echo "Missing $KEEP" >&2; exit 1; }
[[ $(sort -u "$KEEP" | wc -l) -eq 50 ]] || { echo "Expected 50 unique sample IDs" >&2; exit 1; }
[[ $(wc -l < "$KEEP") -eq 50 ]] || { echo "Expected exactly 50 rows" >&2; exit 1; }
export LC_ALL=C
module purge
module load StdEnv

module load PLINK/2.00a3.7-foss-2022a
# Preserve the existing Run A MAC3 marker set, including its original
# MAC filtering and LD pruning in the 66-sample non-contaminated set.
BASE="${PCA_GATK_BASE:-/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/pca_runA_qc_relatedness_mac3}"
OUTDIR="${PCA_GATK_OUT:-/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/pca_curated50_mac3}"
BFILE="${BASE}/autosomes_clean_runA_pruning_base_mac3"
PRUNE="${BASE}/autosomes_clean_runA_mac3_pruned.prune.in"
for file in "$BFILE.bed" "$BFILE.bim" "$BFILE.fam" "$PRUNE"; do
    [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
mkdir -p "$OUTDIR"
cd "$OUTDIR"
cp "$KEEP" curated50_samples.txt
awk '{print $1, $1}' "$KEEP" > curated50.keep
# Check every requested individual exists, before calling PLINK.
comm -23 <(sort "$KEEP") <(awk '{print $2}' "$BFILE.fam" | sort) > missing_samples.txt
[[ ! -s missing_samples.txt ]] || { cat missing_samples.txt >&2; exit 1; }
plink --bfile "$BFILE" --allow-extra-chr --keep curated50.keep \
  --extract "$PRUNE" --pca 50 --threads "${SLURM_CPUS_PER_TASK:-4}" \
  --out autosomes_curated50_mac3_full_pca
awk '{print $2}' autosomes_curated50_mac3_full_pca.eigenvec | sort > actual_samples.txt
sort "$KEEP" > expected_samples.txt
diff -u expected_samples.txt actual_samples.txt
[[ $(wc -l < actual_samples.txt) -eq 50 ]]
echo "Completed GATK PCA: 50 individuals; existing Run A MAC3/LD-pruned markers."
