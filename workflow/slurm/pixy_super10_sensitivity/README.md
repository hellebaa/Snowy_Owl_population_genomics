# SUPER_10 pixy sensitivity checks

This workflow runs a focused `SUPER_10` pixy sensitivity analysis to document
how nucleotide diversity and Tajima's D respond to sample filtering and rare
variant filtering.

It compares three sample sets:

* `all66`: the previous pixy sample set, excluding only `GRL04` and `SKW02`
* `noWRG01_65`: the previous pixy sample set excluding only `WRG01`
* `curated50`: high-quality curated sample set, retaining `FNM12` as the
  representative from the related Fennoscandia group

For each sample set it tests four SNP filters after subsetting and recalculating
`AC`, `AN`, `MAF`, and `F_MISSING`:

* `noMAF`: `F_MISSING <= 0.20`, no minor-allele-frequency filter
* `maf005`: `F_MISSING <= 0.20` and `MAF >= 0.005`
* `maf01`: `F_MISSING <= 0.20` and `MAF >= 0.01`
* `maf02`: `F_MISSING <= 0.20` and `MAF >= 0.02`

The invariant VCF is subset to each sample set before concatenation.

## Run

```bash
sbatch 01_super10_pixy_sample_filter_maf_sensitivity.slurm
```

Main output:

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/08_pop_analysis/03_pixy_SUPER10_sensitivity/summary/SUPER10_pixy_sensitivity_summary.tsv
```

This is intended as a diagnostic/audit analysis, not the final genome-wide
dataset.
