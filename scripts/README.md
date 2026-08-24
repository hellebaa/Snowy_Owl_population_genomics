# Scripts

Downstream analysis scripts should be numbered in execution order.

Current scaffold:

- `00_setup.R`: shared package checks, paths, and plotting theme.
- `01_metadata_qc.R`: rebuilds the curated sample metadata table from historical sample metadata, sample renaming, sex/QC metrics, pixy population files, ROH output, and PCA IDs.
- `02_pca.R`: plots all-sample and filtered GATK autosomal PCA axes, ANGSD pseudohaploid autosomal/Z/W PCA axes, plus a KING kinship heatmap from PCA/KING outputs stored under `data/external/pca/`, with fallback to `../old/pca/` and `../old/re/` during project migration.
- `02b_gatk_pca_regenerated.R`: plots regenerated GATK/PLINK autosomal PCA sensitivity runs from `data/external/pca/gatk/`, including strict QC/relatedness-filtered Run A and contaminated-sample-only Run B, each with and without MAC >= 3 filtering.
- `02c_pca_supporting_figure.R`: assembles the selected four-panel supporting PCA figure from GATK Run A MAC >= 3, ANGSD autosomal MAC >= 3, ANGSD Z, and ANGSD W PCA outputs.
- `02d_angsd_w_female_pca.R`: recomputes ANGSD W chromosome PCA using female samples only, avoiding male samples that should not carry W-linked sequence.
- `03_pixy_all_genomewide.R`: genome-wide pi and Tajima's D plots and summaries for the `ALL` pixy analysis. Also compares macrochromosomes (`SUPER_1`, `SUPER_2`, `SUPER_3`, `SUPER_4`, `SUPER_5`, `SUPER_6`, `SUPER_7`, `SUPER_8`, `SUPER_9`, `SUPER_10`, `SUPER_11`) and microchromosomes based on chicken synteny.
- `04_pixy_centromere_regions.R`: compares pi and Tajima's D in annotated putative centromere intervals, 1 Mb pericentromeric flanks, and the remaining genomic background.
- `05_pixy_regions_differentiation.R`: explores regional pairwise FST and dXY, including pairwise summaries, top FST windows, recurrent outlier windows, and Manhattan plots with pair-specific high-FST outliers highlighted.
- `06_fst_outlier_annotation.R`: annotates recurrent FST outlier windows with overlapping genes, repeat classes, and distance to putative centromeres.
- `07_pixy_all_genomic_features.R`: relates ALL-sample pi and Tajima's D to gene density, repeat content, repeat classes, centromere distance, and macro/microchromosome class.
- `08_roh_summary.R`: summarizes individual ROH burden, FROH, ROH length classes, and regional ROH patterns.
- `09_gone_demography.R`: plots GONE effective population size through time using the syntenic microchromosome recombination-rate analysis, with the 24 high-depth individuals as the primary run and all64/high-depth MAF sensitivity runs as diagnostics.
- `10_roh_landscape.R`: plots ROH landscapes, including long-ROH tracks for high/low FROH individuals, ROH burden by length class, and genome-wide ROH density in 500 kb windows.
- `11_supporting_tables.R`: prepares compact manuscript supporting-table drafts from curated metadata and analysis outputs, with large machine-readable tables left as repository/archive outputs.
- `12_ena_fastq_overview.py`: gathers historical sequencing IDs, submission-table records, lane statistics, and old alignment QC into draft FASTQ/ENA metadata tables under `data/metadata/ena/`.

Planned scripts:

- `13_figures.R`
