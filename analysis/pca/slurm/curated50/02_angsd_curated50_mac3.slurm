#!/bin/bash
#SBATCH --job-name=owl_pca_angsd_50
#SBATCH --account=nn8013k
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=owl_pca_angsd_50_%j.out
#SBATCH --error=owl_pca_angsd_50_%j.err
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

# Use the cluster's default R module, or set PCA_R_MODULE to a known version.
module load "${PCA_R_MODULE:-R}"
HAPLO="${PCA_HAPLO:-${INPUT_DIR}/snowyowl_autosomes_thinned.haplo}"
SAMPLES="${PCA_ANGSD_SAMPLES:-${INPUT_DIR}/sample_names_pca.txt}"
OUTDIR="${PCA_ANGSD_OUT:-/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/pca_angsd_curated50_mac3}"
for file in "$HAPLO" "$SAMPLES" "$INPUT_DIR/angsd_curated50.R"; do
    [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done
mkdir -p "$OUTDIR"
Rscript "$INPUT_DIR/angsd_curated50.R" "$HAPLO" "$SAMPLES" "$KEEP" "$OUTDIR"
echo "Completed ANGSD PCA: 50 individuals; MAC >= 3 recalculated after subsetting."
