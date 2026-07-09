# Data

This repository publishes code and metadata templates, not large genomic data files.

## Tracked in Git

- `metadata/samples_curated.csv`: canonical sample metadata table.
- `metadata/samples_dictionary.csv`: definitions for every metadata column.
- Small derived tables may be committed when they are needed to reproduce manuscript figures and do not contain sensitive information.

## Not Tracked in Git

Large files such as VCF/BCF, BAM/CRAM, PLINK files, and pixy window outputs should be stored in `data/external/` or `data/processed/` locally, but are ignored by Git. Record their filenames, checksums, software versions, and generation steps in `docs/input_data_inventory.md`.
