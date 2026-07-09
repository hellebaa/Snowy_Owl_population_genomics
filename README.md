# Snowy Owl Population Genomics

Reproducible downstream analyses for snowy owl population genomics.

This repository is intended to publish analysis code, metadata documentation,
and manuscript-ready tables/figures. Large genomic data files are not tracked
in Git.

## Project Status

The repository is being rebuilt from cleaned HPC analyses. Historical and
exploratory files are kept outside this repository in the local `old/`
directory.

## Repository Structure

```text
analysis/                 Notes and method records for major analyses
config/                   Version-controlled path template
data/metadata/            Canonical sample metadata and data dictionary
docs/                     Data inventories and provenance notes
manuscript/               Manuscript-facing notes, methods, and figure plans
results/figures/          Regenerated manuscript figures, ignored by Git
results/tables/           Regenerated manuscript tables, ignored by Git
scripts/                  R scripts for downstream analyses
workflow/                 HPC commands, SLURM scripts, and provenance records
```

## Data Policy

Tracked:

- Analysis scripts
- Metadata templates and non-sensitive curated metadata
- Documentation of input files, filtering rules, and software versions
- Small derived tables required for manuscript transparency, when appropriate

Not tracked:

- VCF/BCF files
- BAM/CRAM files
- PLINK binary files
- Large pixy/GONE/window-level outputs
- Local machine or HPC paths

Large inputs should be placed locally in `data/external/` or
`data/processed/`. These folders are ignored by Git except for `.gitkeep`
placeholders.

## Rebuilding the Analyses

1. Open `Snowy_Owl_population_genomics.Rproj` in RStudio.
2. Download required HPC outputs into `data/external/`.
3. Record filenames, source paths, checksums, and generation steps in
   `docs/input_data_inventory.md`.
4. Fill `data/metadata/samples_curated.csv` using
   `data/metadata/samples_dictionary.csv` as the column guide.
5. Run downstream R scripts from `scripts/` in numerical order.

## Current Key Inputs

The downstream analyses expect these primary inputs once downloaded from HPC:

- `SnowyOwl_autosomes_finalFilters_noContam.recalc.vcf.gz`
- `SnowyOwl_autosomes_allSites_noContam.vcf.gz`
- pixy output tables for pi, Tajima's D, FST, and dXY
- PCA outputs
- ROH/heterozygosity outputs
- GONE outputs
- curated sample metadata

## Existing Analysis Notes

- `analysis/SNP_pipeline/README.md`
- `analysis/pixy/README.md`
- `analysis/pca/README.md`
