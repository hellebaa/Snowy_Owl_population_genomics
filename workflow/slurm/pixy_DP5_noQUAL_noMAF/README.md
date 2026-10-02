# pixy all-sites workflow: corrected invariant sites and no MAF filter on SNPs

This workflow is the preferred version for final pixy diversity analyses.

It keeps the corrected invariant-site treatment from `pixy_DP5_noQUAL`:

* invariant records are retained as `ALT=.`
* no site-level `QUAL` filter is applied to invariant records
* invariant genotypes with `FMT/DP < 5` are set to missing

It differs from the previous pixy input by regenerating the SNP VCF without a
minor-allele-frequency filter. This is important for Tajima's D because the
statistic depends on the site-frequency spectrum, especially rare variants.

## SNP input

Expected source VCF on Saga:

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/autosomes.step5.GTfiltered.vcf.gz
```

This is the genotype-filtered autosomal SNP VCF before the final
`MAF >= 0.02` site filter. The script removes contaminated samples,
recalculates `AC`, `AN`, `MAF`, and `F_MISSING`, then retains biallelic SNPs
with:

* `F_MISSING <= 0.20`
* `AC > 0`
* `AC < AN`

No MAF filter is applied.

## Run order

```bash
sbatch 01_make_noMAF_snp_vcf.slurm
sbatch 02_concat_allSites_noQUAL_DP5_noMAF.slurm
sbatch 03_pixy_ALL_noQUAL_DP5_noMAF.slurm
sbatch 04_pixy_regions_noQUAL_DP5_noMAF.slurm
```

The concatenation step reuses the corrected invariant VCF from:

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL/SnowyOwl_autosomes_invariant_noQUAL_DP5.vcf.gz
```

## Expected outputs

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL_noMAF/SnowyOwl_autosomes_snps_noContam_missing20_noMAF.recalc.vcf.gz
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL_noMAF/SnowyOwl_autosomes_allSites_noContam_noQUAL_DP5_noMAF.vcf.gz
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL_noMAF/01_pixy_all/
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL_noMAF/02_pixy_regions/
```
