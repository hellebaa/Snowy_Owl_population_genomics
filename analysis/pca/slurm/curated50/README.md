# Autosomal PCA with the pixy 50 individuals

Both jobs retain exactly the 50 IDs in `curated50_samples.txt`, copied from the curated50 pixy workflow. FNM12 is included. Existing PCA outputs and manuscript files are not overwritten.

## Upload from your Mac

Upload the scripts, sample lists and R helper together:

```bash
rsync -ravz /Users/hellebaalsrud/Documents/Projects/Snowy_owl_population/new/analysis/pca/slurm/curated50/ hellb@login-1.saga.sigma2.no:/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/curated50_scripts/
rsync -ravz /Users/hellebaalsrud/Documents/Projects/Snowy_owl_population/old/pca/PCA_INPUT/snowyowl_autosomes_thinned.haplo hellb@login-1.saga.sigma2.no:/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/curated50_scripts/
```

The matrix is the existing 55 MB thinned ANGSD input used locally. Uploading this known input avoids guessing its current cluster location; `sample_names_pca.txt` supplies its original column order.

## Submit on Saga

```bash
cd /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/curated50_scripts
sbatch 01_gatk_curated50_mac3.slurm
sbatch 02_angsd_curated50_mac3.slurm
```

The jobs are independent. Account `nn8013k` and the PLINK module follow the existing PCA SLURMs; storage remains under `nn9244k`. ANGSD uses base R only and loads the default `R` module. If necessary, set `PCA_R_MODULE` to an available version before submission. Cluster modules and input availability have not been checked remotely.

## Methods and outputs

GATK reuses the Run A binary genotype input and existing MAC3/LD-pruned marker list. Its MAC threshold was originally applied in the 66-individual non-contaminated set, not recalculated among 50. This deliberately preserves the old marker selection to isolate the sample change. Results go to `pca_curated50_mac3` alongside the old Run A directory. The eigenvectors are checked against all 50 requested IDs.

ANGSD recalculates biallelic polymorphism and MAC >=3 after subsetting the existing thinned matrix to 50 individuals. The recoding, centering, mean-imputation and covariance PCA match `new/scripts/02_pca.R`. Results go to `pca_angsd_curated50_mac3`, including score/eigenvalue CSVs, actual sample IDs, retained sites and R session details. No genotype calling or new random allele sampling is performed.

After both jobs complete, download their outputs:

```bash
rsync -ravz hellb@login-1.saga.sigma2.no:/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/pca_curated50_mac3/ /Users/hellebaalsrud/Documents/Projects/Snowy_owl_population/new/data/external/pca/gatk/curated50_mac3/
rsync -ravz hellb@login-1.saga.sigma2.no:/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/pca_angsd_curated50_mac3/ /Users/hellebaalsrud/Documents/Projects/Snowy_owl_population/new/data/external/pca/angsd/curated50_mac3/
```

The plotting scripts still point to the previous 49-individual outputs. Update them and the manuscript after checking the rerun results. Z/W diagnostic PCAs are unchanged.

## Export the full GATK spectrum for variance percentages

The first rerun exported only 20 eigenvalues. These support the plotted coordinates, but their sum is not total variance. `03_gatk_curated50_full_spectrum.slurm` uses the same samples and markers and requests all 50 components, with a separate `autosomes_curated50_mac3_full_pca` output prefix. Upload this script alongside the existing companion files, then submit from that directory:

```bash
sbatch 03_gatk_curated50_full_spectrum.slurm
```

Download the `pca_curated50_mac3` directory again using the command above. Running `Rscript scripts/02e_pca_curated50.R` locally from `new/` will detect the full spectrum and populate the GATK percentages in the supporting figure and Table S4.

### Full-spectrum upload from your Mac

Upload the full-spectrum SLURM and its required sample list together:

```bash
rsync -ravz /Users/hellebaalsrud/Documents/Projects/Snowy_owl_population/new/analysis/pca/slurm/curated50/03_gatk_curated50_full_spectrum.slurm /Users/hellebaalsrud/Documents/Projects/Snowy_owl_population/new/analysis/pca/slurm/curated50/curated50_samples.txt hellb@login-1.saga.sigma2.no:/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/curated50_scripts/
```

### Full-spectrum submission on Saga

```bash
cd /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/curated50_scripts
sbatch 03_gatk_curated50_full_spectrum.slurm
```

### Full-spectrum download to your Mac

After successful completion, download the output directory, including the new full-spectrum eigenvalues and eigenvectors:

```bash
rsync -ravz hellb@login-1.saga.sigma2.no:/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/01_pca/02_pca_helle/pca_curated50_mac3/ /Users/hellebaalsrud/Documents/Projects/Snowy_owl_population/new/data/external/pca/gatk/curated50_mac3/
```

Then regenerate the supporting outputs locally:

```bash
cd /Users/hellebaalsrud/Documents/Projects/Snowy_owl_population/new
Rscript scripts/02e_pca_curated50.R
```
