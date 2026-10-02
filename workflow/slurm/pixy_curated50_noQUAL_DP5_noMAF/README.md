# Curated 50-sample pixy workflow

This workflow reruns pixy after removing samples with high missingness,
contamination-review/remove flags, and close relatives except for the retained
representative `FNM12_5-X5_S119_L004`.

Retained samples:

* `qc_flags == OK`
* `contamination_flag == keep`
* `relatedness_flag == keep`, except `FNM12` retained as the representative
  individual from the related Fennoscandia group

Final sample count: 50.

Population counts:

* FNM: 13
* GRL: 8
* NYS: 13
* SKW: 11
* WRG: 5

The SNP VCF is regenerated without MAF filtering because Tajima's D requires
rare variants. It keeps biallelic polymorphic SNPs after subsetting to the
curated 50 samples and applying `F_MISSING <= 0.20`.

The invariant VCF is subset to the same curated 50 samples before concatenation
so sample columns match exactly.

## Run order

```bash
sbatch 01_make_curated50_noMAF_snp_vcf.slurm
sbatch 02_subset_invariant_and_concat_curated50.slurm
sbatch 03_pixy_ALL_curated50_noQUAL_DP5_noMAF.slurm
sbatch 04_pixy_regions_curated50_noQUAL_DP5_noMAF.slurm
sbatch 05_super10_singleton_qc_curated50.slurm
```

After step 2 completes, steps 3 and 4 can be run independently. Step 5 can run
after step 1 and is a quick QC check for singleton heterozygote support on
`SUPER_10`.

## Sample-size note

WRG has only five retained individuals, so FST and dXY estimates involving WRG
will be noisier than comparisons among larger regional samples. This is not a
fatal problem for genome-wide summaries, but window-based FST scans should keep
minimum SNP/site thresholds and should be interpreted cautiously.
