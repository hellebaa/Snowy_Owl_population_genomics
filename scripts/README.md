# Scripts

Downstream analysis scripts should be numbered in execution order.

Current scaffold:

- `00_setup.R`: shared package checks, paths, and plotting theme.
- `03_pixy_all_genomewide.R`: genome-wide pi and Tajima's D plots and summaries for the `ALL` pixy analysis. Also compares macrochromosomes (`SUPER_1`, `SUPER_2`, `SUPER_3`, `SUPER_4`, `SUPER_5`, `SUPER_6`, `SUPER_7`, `SUPER_8`, `SUPER_9`, `SUPER_10`, `SUPER_11`) and microchromosomes based on chicken synteny.
- `04_pixy_centromere_regions.R`: compares pi and Tajima's D in annotated putative centromere intervals, 1 Mb pericentromeric flanks, and the remaining genomic background.
- `05_pixy_regions_differentiation.R`: explores regional pairwise FST and dXY, including pairwise summaries, top FST windows, recurrent outlier windows, and Manhattan plots with pair-specific high-FST outliers highlighted.
- `06_fst_outlier_annotation.R`: annotates recurrent FST outlier windows with overlapping genes, repeat classes, and distance to putative centromeres.
- `07_pixy_all_genomic_features.R`: relates ALL-sample pi and Tajima's D to gene density, repeat content, repeat classes, centromere distance, and macro/microchromosome class.
- `08_roh_summary.R`: summarizes individual ROH burden, FROH, ROH length classes, and regional ROH patterns.
- `09_gone_demography.R`: plots GONE effective population size through time using the syntenic microchromosome recombination-rate analysis, with the 24 high-depth individuals as the primary run and all64/high-depth MAF sensitivity runs as diagnostics.
- `10_roh_landscape.R`: plots ROH landscapes, including long-ROH tracks for high/low FROH individuals, ROH burden by length class, and genome-wide ROH density in 500 kb windows.

Planned scripts:

- `01_metadata_qc.R`
- `02_pca.R`
- `11_figures.R`
