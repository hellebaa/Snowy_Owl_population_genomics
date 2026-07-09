# GONE Demographic Analysis

## Input dataset

GONE analyses were run on genotype-filtered SNPs from syntenic autosomal
microchromosomes with inferred recombination rates. No linkage pruning was
performed, following recommendations for GONE.

The primary analysis uses 24 high-depth individuals. This dataset was selected
after sensitivity analyses showed that the full 64-individual polymorphic-site
dataset produced implausibly high effective population size estimates. The
inflation was consistent with sensitivity to low-depth individuals and very rare
variants. In the 24-individual high-depth subset, the lowest observable allele
frequency is one allele copy in 48 chromosomes (0.0208), and GONE trajectories
were stable across the tested variant-filtering settings.

Sensitivity analyses were retained for transparency:

- all 64 unrelated/non-contaminated individuals, polymorphic sites
- all 64 individuals, MAF >= 0.005
- all 64 individuals, MAF >= 0.01
- all 64 individuals, MAF >= 0.02
- 24 high-depth individuals, polymorphic sites
- 24 high-depth individuals, MAF >= 0.005
- 24 high-depth individuals, MAF >= 0.01
- 24 high-depth individuals, MAF >= 0.02

The all-64 polymorphic and MAF >= 0.005 runs are treated as diagnostic
sensitivity analyses only, because they inferred biologically implausible
effective population sizes exceeding one million individuals over parts of the
recent past. The all-64 MAF >= 0.01 and MAF >= 0.02 runs were more consistent
with the high-depth analyses.

## Samples

Samples removed before the 64-individual sensitivity analyses:
- GRL04_merged_X14
- SKW02_merged_X22
- FNM11_merged_X4
- FNM15_merged_X8

Retained related sample:
- FNM12_5-X5_S119_L004

The primary 24-individual analysis further restricted the dataset to individuals
with high sequencing depth. The exact sample list should be stored with the HPC
GONE input files or project sample metadata.

## Recombination rate estimate

Inference of recent effective population size using GONE requires an estimate of
the recombination rate. Because no recombination map is currently available for
Snowy Owl (*Bubo scandiacus*), recombination rates were approximated using
chromosome-specific estimates from chicken (*Gallus gallus*).

To minimize bias associated with differences in chromosome structure between
species, analyses were restricted to a subset of Snowy Owl autosomes satisfying
four criteria:

- One-to-one chromosomal homology with chicken based on whole-genome synteny
  analyses.
- No evidence of chromosome fusion or fission relative to the homologous chicken
  chromosome.
- Comparable chromosome morphology between species, with both chromosomes
  classified as acrocentric.
- Availability of published chromosome-specific recombination rate estimates for
  the homologous chicken chromosome.

Macrochromosomes were excluded because recombination rates differ substantially
between macrochromosomes and microchromosomes in birds. Chromosomes with evidence
of fusion, uncertain homology, or lacking published recombination estimates were
also excluded.

Application of these criteria resulted in retention of 17 Snowy Owl
microchromosomes:

- SUPER_11
- SUPER_13
- SUPER_14
- SUPER_15
- SUPER_17
- SUPER_18
- SUPER_19
- SUPER_20
- SUPER_21
- SUPER_22
- SUPER_23
- SUPER_24
- SUPER_25
- SUPER_26
- SUPER_27
- SUPER_28
- SUPER_29

For each retained chromosome, the chromosome-specific recombination rate reported
for the homologous chicken chromosome was assigned to the Snowy Owl chromosome. A
length-weighted mean recombination rate was then calculated across all retained
chromosomes:

```text
6.63 cM/Mb
```

This value was used in all GONE analyses.

## Interpretation

GONE estimates recent effective population size from linkage disequilibrium.
For this study, inference is interpreted most strongly over the most recent
approximately 200 generations, corresponding to roughly 1,600 years assuming an
eight-year generation time. Older trajectories are shown as exploratory
sensitivity analyses and should be interpreted cautiously because LD-based
signals become weaker, smoother, and more sensitive to recombination-rate
assumptions, marker filtering, sample composition, and genotyping error farther
back in time.

The main plotting script is:

```text
scripts/09_gone_demography.R
```

Key outputs:

- `results/figures/gone/SnowyOwl_GONE_Ne_timeseries.png`
- `results/figures/gone/SnowyOwl_GONE_sensitivity_Ne_by_sample_set.png`
- `results/figures/gone/SnowyOwl_GONE_sensitivity_Ne_timeseries_log10.png`
- `results/figures/gone/SnowyOwl_GONE_sensitivity_Ne_by_sample_set_longterm.png`
- `results/tables/gone/SnowyOwl_GONE_sensitivity_key_generations.csv`
- `results/tables/gone/SnowyOwl_GONE_sensitivity_input_diagnostics_summary.csv`
