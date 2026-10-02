# Current curated50 PCA update

The final autosomal PCA sample set now matches pixy: 50 individuals, retaining FNM12. Run `Rscript scripts/02e_pca_curated50.R` from `new/` to import downloaded GATK results, regenerate ANGSD locally and update Figure S2, Tables S2/S4 and sample inclusion. GATK: 564,409 SNPs; ANGSD: 39,231 sites, PC1 2.203%, PC2 2.187%. Z/W diagnostics are unchanged.

The full-spectrum export has now been downloaded and verified: all 50 eigenvalues are available, and the leading coordinates agree with the earlier export. GATK PC1–PC4 explain 2.145%, 2.102%, 2.095% and 2.089% of total variance. Figure S2 and Table S4 now use the full-spectrum denominator. The figure shows two GATK PC1 outliers, GRL01 and WRG09, from different regions; their cause is unresolved.

The following records describe historical runs, not the current final sample set.

# PCA

This analysis summarizes population structure using autosomal PCA, ANGSD
pseudohaploid PCA, and KING kinship outputs.

Expected local inputs are stored outside Git under `data/external/pca/`:

- `gatk/historical/autosomes_pca.eigenvec`: historical all-sample PLINK/GATK-style eigenvectors, with FID, IID, and PC columns.
- `gatk/historical/autosomes_pca.eigenval`: optional eigenvalues for percent-variance axis labels.
- `gatk/historical/autosomal_rm_ind.eigenvec` and `gatk/historical/autosomal_rm_ind.eigenval`: historical filtered PCA output after removal of PCA/QC outliers.
- `gatk/runA_qc_relatedness_mac3/`, `gatk/runA_qc_relatedness_no_mac/`, `gatk/runB_contaminated_only_mac3/`, and `gatk/runB_contaminated_only_no_mac/`: folders for the regenerated GATK/PLINK PCA runs.
- `king/autosomes_clean_pruned.king`: KING kinship output from the LD-pruned autosomal SNP dataset.
- `king/autosomes_clean_pruned.king.id`: sample IDs for the KING output.
- `angsd/sample_names_pca.txt` and thinned ANGSD haplo inputs can be kept here for provenance or future ANGSD/PCAngsd reruns.

During migration from the historical project folder, `scripts/02_pca.R` falls back to the corresponding files in `../old/pca/` if the `data/external/pca/` copies are not present. Generated figures and tables are written to `results/figures/pca/` and `results/tables/pca/`; these outputs are intentionally ignored by Git.

Primary generated outputs:

- `SnowyOwl_PCA_autosomes_filtered_PC1_PC2.png/pdf`
- `SnowyOwl_PCA_autosomes_filtered_PC1_PC3.png/pdf`
- `SnowyOwl_PCA_autosomes_filtered_PC3_PC4.png/pdf`
- `SnowyOwl_PCA_autosomes_filtered_scree.png/pdf`
- `SnowyOwl_PCA_autosomes_all_samples_PC1_PC2.png/pdf`
- `SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_PC1_PC2.png/pdf`
- `SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_PC1_PC2.png/pdf`
- `SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_mac3_PC1_PC2.png/pdf`
- `SnowyOwl_PCA_ANGSD_pseudohaploid_Z_PC1_PC2.png/pdf`
- `SnowyOwl_PCA_ANGSD_pseudohaploid_W_PC1_PC2.png/pdf`
- `SnowyOwl_KING_kinship_heatmap.png/pdf`

## GATK/PLINK PCA provenance status

The historical filtered GATK/PLINK PCA input in
`data/external/pca/gatk/historical/` is `autosomal_rm_ind.eigenvec` with 63
samples. Mapping the old sequence IDs in that file through
`../old/sample_rename.txt` shows that the five samples absent from the filtered
PCA are:

- `GRL04_merged_X14`
- `SKW02_merged_X22`
- `FNM11_merged_X4`
- `FNM15_merged_X8`
- `WRG01_1-WRG1_S28_L001`

The historical all-sample GATK/PLINK PCA input is `autosomes_pca.eigenvec` with
68 samples. The exact upstream SLURM script that produced
`gatk/historical/autosomal_rm_ind.eigenvec` should still be confirmed from the
HPC project history. A candidate `pca_remove_FNM` script removes
`FNM12_5-X5_S119_L004` instead of `FNM15_merged_X8` and
`WRG01_1-WRG1_S28_L001`, and therefore does not exactly match the historical
filtered PCA file used during migration.

A separate contaminated-sample-only PCA script was found in
`08_pop_analysis/01_pca/02_pca_helle/pca_snowy_autosomes_removed_ind.slurm`.
This script removes only:

- `GRL04_merged_X14`
- `SKW02_merged_X22`

It converts the filtered autosomal VCF to PLINK format, performs LD pruning
with `--indep-pairwise 50 10 0.2`, and runs PCA on the pruned variants:

```bash
plink \
  --vcf /cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/autosomes.final.vcf.gz \
  --double-id \
  --allow-extra-chr \
  --remove remove_samples.txt \
  --make-bed \
  --out autosomes_clean

plink \
  --bfile autosomes_clean \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out autosomes_clean_pruned

plink \
  --bfile autosomes_clean \
  --allow-extra-chr \
  --extract autosomes_clean_pruned.prune.in \
  --pca 20 \
  --out autosomes_clean_pca
```

This contaminated-sample-only workflow likely corresponds to the
`autosomes_clean_pca.*` outputs in the HPC folder, not to the 63-sample
`autosomal_rm_ind.eigenvec` file currently used for the filtered manuscript PCA
plot.

ANGSD haplo inputs are recoded site-by-site as major allele = 0, minor allele = 1, and missing = `NA`; invariant and multi-allelic sites are excluded. Missing genotypes are mean-imputed after centering, which is equivalent to setting missing values to zero contribution in the centered genotype matrix.

For manuscript Figure S2, the current preferred ANGSD PCA candidate is `SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_mac3_PC1_PC2.png/pdf`. This uses the curated metadata to retain filtered-PCA samples with `qc_flags == OK`, `relatedness_flag == keep`, and `contamination_flag == keep`, and applies a minor-allele-count threshold of 3. The all-sample and MAC >= 1 outputs are kept as diagnostics because they are dominated by individual high-missingness/private-variant axes.

Pairwise kinship estimates were generated on the HPC cluster with PLINK
2.00a3.7. First, the LD-pruned autosomal binary PLINK dataset was generated
from `autosomes_clean` using the retained variants listed in
`autosomes_clean_pruned.prune.in`. The square KING kinship matrix was then
calculated from the pruned dataset:

```bash
module load PLINK/2.00a3.7-foss-2022a

plink2 \
  --bfile autosomes_clean \
  --allow-extra-chr \
  --extract autosomes_clean_pruned.prune.in \
  --make-bed \
  --out autosomes_clean_pruned

plink2 \
  --bfile autosomes_clean_pruned \
  --allow-extra-chr \
  --make-king square \
  --out autosomes_clean_pruned
```

The resulting KING kinship estimates were used for relatedness screening and
for the relatedness heatmap generated by `scripts/02_pca.R`.

The PCA colors use the manuscript population palette defined in `scripts/00_setup.R`.
