# pixy workflow: invariant sites without QUAL, DP-only genotype filter

This workflow is the corrected final candidate for pixy analyses after the
`SUPER_10` sensitivity test:

```text
noGTfilter   pi = 0.000216663
DP5          pi = 0.000227859
DP5_RGQ20    pi = 0.00112283
```

The depth-only filter had little effect on `SUPER_10` diversity estimates,
whereas adding `RGQ < 20` greatly reduced the callable denominator and inflated
pi. This workflow therefore:

* retains invariant sites as `ALT="."`
* does not apply a site-level invariant-site `QUAL` filter
* sets invariant genotypes with `FMT/DP < 5` to missing
* does not filter invariant genotypes on `RGQ`
* uses the existing final filtered autosomal SNP VCF for variant sites

Cluster output folders:

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/01_gatk4-4/02_non_variant/07_invariant_noQUAL_DP5
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL
```

## Run order

```bash
sbatch 01_extract_invariant_noQUAL_DP5_array.slurm
```

After all 32 array jobs finish:

```bash
sbatch 02_concat_allSites_noQUAL_DP5.slurm
```

Then:

```bash
sbatch 03_pixy_ALL_noQUAL_DP5.slurm
sbatch 04_pixy_regions_noQUAL_DP5.slurm
```

## Quick checks

After concatenation:

```bash
bcftools index -n /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL/SnowyOwl_autosomes_invariant_noQUAL_DP5.vcf.gz
bcftools index -n /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL/SnowyOwl_autosomes_allSites_noContam_noQUAL_DP5.vcf.gz
```

After pixy ALL:

```bash
awk '
NR==1 {
  for (i=1; i<=NF; i++) {
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
}
END {
  print "windows", windows;
  print "zero_site_windows", zero;
  print "callable_sites_sum", sites;
  print "count_missing_sum", missing;
  print "count_comparisons_sum", comps;
}
' /cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_DP5_noQUAL/01_pixy_all/SnowyOwl_ALL_noQUAL_DP5_pi.txt
```
