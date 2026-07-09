# Runs of Homozygosity (ROH) Analysis

## Overview

Runs of homozygosity (ROH) were used to characterize individual inbreeding levels and recent demographic history in Snowy Owls.

ROHs are long homozygous genomic segments inherited identically by descent and provide information about both recent and historical inbreeding. Short ROHs generally reflect older demographic processes, whereas long ROHs indicate more recent shared ancestry.

Analyses were performed using PLINK v1.9.

---

## Input Dataset

Input VCF:

```text
/cluster/projects/nn9244k/for_emily/SnowyOwl/05_genotyped/03_filtered_autosomes/SnowyOwl_autosomes_finalFilters_noContam.recalc.vcf.gz
```

This dataset contains:

- Autosomal SNPs only
- Genotype filtering:
  - DP ≥ 5
  - GQ ≥ 20
- Site filtering:
  - MAF ≥ 0.02
  - F_MISSING ≤ 0.20
- Contaminated samples removed:
  - GRL04_merged_X14
  - SKW02_merged_X22

Final SNP count:

- 1,324,912 SNPs

All remaining individuals were retained for ROH analyses.

Unlike GONE analyses, no additional samples were removed based on relatedness because ROH is an individual-level metric and related individuals do not bias ROH estimation.

---

## PLINK Conversion

The filtered VCF was converted to binary PLINK format:

```bash
plink \
  --vcf SnowyOwl_autosomes_finalFilters_noContam.recalc.vcf.gz \
  --make-bed
```

Chromosome names were retained using:

```bash
--allow-extra-chr
--chr-set 70 no-xy
```

---

## ROH Calling Parameters

ROHs were identified using PLINK's sliding-window approach with the following settings:

```bash
--homozyg
--homozyg-window-snp 50
--homozyg-window-het 3
--homozyg-window-missing 10
--homozyg-window-threshold 0.05

--homozyg-snp 50
--homozyg-kb 300
--homozyg-density 50
--homozyg-gap 1000
--homozyg-het 3
```

These settings define an ROH as:

- Minimum 50 SNPs
- Minimum length 300 kb
- Maximum SNP spacing 50 kb
- Maximum gap between consecutive SNPs 1 Mb
- Up to 3 heterozygous sites allowed per ROH
- Up to 10 missing genotypes allowed per sliding window

Allowing a small number of heterozygous sites helps account for genotyping errors and sequencing artefacts.

---

## Output Files

Main ROH file:

```text
SnowyOwl_autosomes_noContam_ROH.hom
```

Per-individual summary:

```text
SnowyOwl_autosomes_noContam_ROH.hom.indiv
```

Length-bin summary:

```text
SnowyOwl_autosomes_noContam_ROH.ROH_length_bins.tsv
```

---

## ROH Length Categories

ROHs were classified into three length classes:

| Category | Length |
|-----------|----------|
| Short ROH | 0.3–1 Mb |
| Intermediate ROH | 1–5 Mb |
| Long ROH | >5 Mb |

### Interpretation

- Short ROHs primarily reflect older demographic history.
- Intermediate ROHs reflect more recent shared ancestry.
- Long ROHs indicate recent inbreeding and mating among close relatives.

---

## Downstream Analyses

For each individual, total ROH length is calculated as:

```text
Sum of all ROH segment lengths
```

Genomic inbreeding (FROH) is estimated as:

```text
FROH = Total ROH length / Total autosomal genome length
```

Population-level summaries can be calculated by sampling region and visualized using boxplots and ROH length distributions.

---

## Relationship to Other Analyses

### pixy

pixy requires both variant and invariant sites and therefore used a separate all-sites VCF.

### GONE

GONE used a different SNP dataset:

- No MAF filtering
- Related individuals removed
- Subset of syntenic microchromosomes only

### ROH

ROH analyses use the final filtered autosomal SNP dataset because ROH detection relies on high-quality polymorphic markers rather than invariant sites.

---

## Software

- PLINK v1.9

## References

Purcell S. et al. (2007) PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. *American Journal of Human Genetics* 81:559–575.

Chang C.C. et al. (2015) Second-generation PLINK: Rising to the challenge of larger and richer datasets. *GigaScience* 4:7.