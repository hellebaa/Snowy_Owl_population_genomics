# Corrected pixy all-sites workflow without invariant-site QUAL filtering

These SLURM scripts regenerate the pixy all-sites VCF after removing the
site-level `QUAL >= 30` filter from invariant sites. This corrects the low
callable-site denominator observed in the previous pixy input VCF.

The old files are not overwritten. Corrected files are written on Saga under:

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4/02_non_variant/06_invariant_noQUAL
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL
```

## Order to run

From any folder on Saga, submit:

```bash
sbatch 01_extract_invariant_noQUAL_array.slurm
```

Wait until all 32 array tasks finish successfully, then submit:

```bash
sbatch 02_concat_allSites_noQUAL.slurm
```

After the corrected all-sites VCF has been created, run:

```bash
sbatch 03_pixy_ALL_noQUAL.slurm
sbatch 04_pixy_regions_noQUAL.slurm
```

The pixy scripts copy `population_all.txt` and `population_regions.txt` from
the previous pixy folder if they are not already present in the corrected
workflow folder.

## Important filtering choice

Invariant sites are retained as `ALT="."` without a site-level `QUAL` filter.
Low-confidence genotypes are set to missing using:

```text
FMT/DP < 5 or FMT/RGQ < 20
```

This follows the pixy recommendation to avoid applying variant-style `QUAL`
filters to invariant sites. Variant sites remain filtered separately using the
existing GATK SNP-filtered VCF:

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/SnowyOwl_autosomes_finalFilters_noContam.recalc.vcf.gz
```

## Quick checks after each step

After step 1:

```bash
grep -i "error\|done\|invariant records" owl_inv_noQUAL_*.out owl_inv_noQUAL_*.err
cat /cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4/02_non_variant/06_invariant_noQUAL/*.count.tsv
```

After step 2:

```bash
bcftools index -n /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL/SnowyOwl_autosomes_invariant_noQUAL.vcf.gz
bcftools index -n /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL/SnowyOwl_autosomes_allSites_noContam_noQUAL.vcf.gz
```

After pixy:

```bash
ls -lh /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL/01_pixy_all
ls -lh /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL/02_pixy_regions
```

## SUPER_10 genotype-filter diagnostic

If the corrected all-sites VCF is genome-scale but pixy still reports a much
smaller `no_sites` denominator, run:

```bash
sbatch 05_super10_noGTfilter_diagnostic.slurm
```

This diagnostic extracts `SUPER_10` invariant sites with no site-level `QUAL`
filter and no genotype-level `DP/RGQ` filter, concatenates those sites with the
final filtered `SUPER_10` variant sites, and runs pixy on that chromosome only.
Compare its `callable_sites_sum` against the full noQUAL workflow. If `no_sites`
jumps strongly, the `DP/RGQ` genotype filter is responsible for much of the
remaining denominator reduction.

To explicitly separate the effects of depth and reference genotype quality,
run:

```bash
sbatch 06_super10_filter_comparison.slurm
```

This compares three `SUPER_10` invariant-site treatments:

```text
noGTfilter  = no genotype-level filtering
DP5         = set genotypes with DP < 5 to missing
DP5_RGQ20   = set genotypes with DP < 5 or RGQ < 20 to missing
```

The output summary is written to:

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_noQUAL/SUPER_10_filter_comparison/SUPER_10_filter_comparison_summary.tsv
```
