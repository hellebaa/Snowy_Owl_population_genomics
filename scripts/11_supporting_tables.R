# Prepare compact manuscript supporting tables from reproducible analysis outputs.

source("scripts/00_setup.R")

supporting_dir <- file.path("manuscript", "supporting_tables")
dir.create(supporting_dir, recursive = TRUE, showWarnings = FALSE)

write_supporting_csv <- function(x, filename) {
  readr::write_csv(x, file.path(supporting_dir, filename), na = "")
}

round_numeric <- function(x, digits = 6) {
  dplyr::mutate(x, dplyr::across(where(is.numeric), ~ round(.x, digits)))
}

metadata <- readr::read_csv(paths$metadata, show_col_types = FALSE)

# S1: Sample metadata overview.
table_s1 <- metadata |>
  dplyr::transmute(
    sample_id,
    sample_short_id,
    individual_id,
    country,
    location,
    population_id = population_final,
    collection_year,
    sex_field = dplyr::recode(sex_field, Male = "M", Female = "F", Unknown = "NA"),
    age,
    tissue = dplyr::recode(tissue, Muscle = "MU", Blood = "BL", Blodfjær = "BF", Feather = "FE", .default = tissue)
  )
write_supporting_csv(table_s1, "Table_S1_sample_metadata_and_analysis_inclusion.csv")

# S2: Sequencing, genotype QC, and analysis inclusion summary.
table_s2 <- metadata |>
  dplyr::transmute(
    sample_short_id,
    mean_depth,
    missingness,
    heterozygosity,
    sex_field = dplyr::recode(sex_field, Male = "M", Female = "F", Unknown = "NA"),
    sex_genomic = dplyr::recode(sex_genomic, Male = "M", Female = "F", Unknown = "NA"),
    sex_final = dplyr::recode(sex_final, Male = "M", Female = "F", Unknown = "NA"),
    qc_flags,
    relatedness_flag,
    contamination_flag,
    include_pca_filtered = dplyr::if_else(include_pca_filtered, "T", "F"),
    include_pixy = dplyr::if_else(include_pixy, "T", "F"),
    include_roh = dplyr::if_else(include_roh, "T", "F"),
    include_gone = dplyr::if_else(include_gone, "T", "F"),
    notes
  ) |>
  round_numeric(6)
write_supporting_csv(table_s2, "Table_S2_sequencing_and_genotype_QC.csv")

# S3: Pairwise KING kinship estimates above the reporting threshold.
king_pairs <- tibble::tibble()
king_file <- paths$pca_king_matrix
king_id_file <- paths$pca_king_ids

if (file.exists(king_file) && file.exists(king_id_file)) {
  king_ids <- readr::read_table(
    king_id_file,
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    dplyr::rename(sample_id = IID)

  king <- as.matrix(readr::read_table(
    king_file,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_double())
  ))
  rownames(king) <- king_ids$sample_id
  colnames(king) <- king_ids$sample_id
  idx <- which(upper.tri(king) & king > 0.03, arr.ind = TRUE)

  if (nrow(idx) > 0) {
    king_pairs <- tibble::tibble(
      sample_id_1 = rownames(king)[idx[, "row"]],
      sample_id_2 = colnames(king)[idx[, "col"]],
      king_kinship = king[idx]
    ) |>
      dplyr::left_join(
        metadata |> dplyr::select(sample_id, sample_short_id_1 = sample_short_id, population_1 = population_final, relatedness_flag_1 = relatedness_flag),
        by = c("sample_id_1" = "sample_id")
      ) |>
      dplyr::left_join(
        metadata |> dplyr::select(sample_id, sample_short_id_2 = sample_short_id, population_2 = population_final, relatedness_flag_2 = relatedness_flag),
        by = c("sample_id_2" = "sample_id")
      ) |>
      dplyr::transmute(
        sample_id_1, sample_id_2, sample_short_id_1, sample_short_id_2,
        population_1, population_2, king_kinship, relatedness_flag_1,
        relatedness_flag_2,
        notes = "Pairwise KING kinship > 0.03 in autosomal pruned SNP dataset."
      )
  }
}

table_s3 <- king_pairs |>
  dplyr::arrange(dplyr::desc(king_kinship), sample_id_1, sample_id_2) |>
  round_numeric(6)
write_supporting_csv(table_s3, "Table_S3_sample_exclusions_and_KING_relatedness.csv")

# S4: PCA dataset summaries.
count_rows <- function(file) if (file.exists(file)) nrow(readr::read_csv(file, show_col_types = FALSE)) else NA_integer_
first_pve <- function(file, pc_name) {
  if (!file.exists(file)) return(NA_real_)
  x <- readr::read_csv(file, show_col_types = FALSE)
  x$pve[match(pc_name, x$pc)]
}

gatk_eigen <- if (file.exists(paths$pca_filtered_eigenval)) scan(paths$pca_filtered_eigenval, quiet = TRUE) else numeric()
gatk_pve <- gatk_eigen / sum(gatk_eigen) * 100

table_s4 <- tibble::tribble(
  ~dataset, ~input_type, ~sample_filter, ~site_filter, ~n_samples, ~n_sites_or_snps, ~pc1_variance_percent, ~pc2_variance_percent, ~recommended_use,
  "GATK autosomes, filtered", "diploid SNP genotypes", "Historical filtered PCA set; high-missingness and close-related samples removed", "Filtered autosomal SNPs, LD-pruned before PCA", count_rows(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_autosomes_filtered_scores.csv")), NA_real_, gatk_pve[1], gatk_pve[2], "Figure S1 candidate",
  "ANGSD pseudohaploid autosomes", "pseudohaploid alleles", "QC-passing samples; relatedness_flag == keep; contamination_flag == keep", "Biallelic sites, MAC >= 3", count_rows(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_mac3_scores.csv")), 38800, first_pve(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_mac3_eigenvalues.csv"), "PC1"), first_pve(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_mac3_eigenvalues.csv"), "PC2"), "Figure S2 candidate",
  "ANGSD pseudohaploid Z", "pseudohaploid alleles", "All available ANGSD samples", "Biallelic sites, MAC >= 1", count_rows(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_Z_scores.csv")), NA_real_, first_pve(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_Z_eigenvalues.csv"), "PC1"), first_pve(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_Z_eigenvalues.csv"), "PC2"), "Diagnostic only",
  "ANGSD pseudohaploid W", "pseudohaploid alleles", "All available ANGSD samples", "Biallelic sites, MAC >= 1", count_rows(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_W_scores.csv")), NA_real_, first_pve(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_W_eigenvalues.csv"), "PC1"), first_pve(file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_W_eigenvalues.csv"), "PC2"), "Diagnostic only"
) |>
  round_numeric(3)
write_supporting_csv(table_s4, "Table_S4_PCA_dataset_summaries.csv")

# S5: Pairwise differentiation and divergence among sampling regions.
fst_summary <- readr::read_csv(
  file.path(paths$results_tables, "pixy", "SnowyOwl_regions_pairwise_fst_summary.csv"),
  show_col_types = FALSE
) |>
  dplyr::mutate(
    pop1_label = region_display(pop1),
    pop2_label = region_display(pop2),
    pair_label = paste(pop1_label, pop2_label, sep = "_vs_")
  )

dxy_summary <- readr::read_tsv(paths$pixy_regions_dxy, show_col_types = FALSE) |>
  dplyr::filter(!is.na(avg_dxy), no_sites > 0) |>
  dplyr::mutate(
    pop1_label = region_display(pop1),
    pop2_label = region_display(pop2),
    pair_label = paste(pop1_label, pop2_label, sep = "_vs_")
  ) |>
  dplyr::group_by(pair_label) |>
  dplyr::summarise(
    mean_dxy_weighted_by_sites = stats::weighted.mean(avg_dxy, no_sites, na.rm = TRUE),
    median_dxy = stats::median(avg_dxy, na.rm = TRUE),
    .groups = "drop"
  )

table_s5 <- fst_summary |>
  dplyr::left_join(dxy_summary, by = "pair_label") |>
  dplyr::transmute(
    comparison = pair_label,
    population_1 = pop1_label,
    population_2 = pop2_label,
    n_windows,
    mean_fst_weighted_by_snps,
    median_fst,
    q95_fst,
    q99_fst,
    max_fst,
    mean_dxy_weighted_by_sites,
    median_dxy
  ) |>
  round_numeric(6)
write_supporting_csv(table_s5, "Table_S5_pairwise_FST_and_dXY.csv")

# S10: Chromosomes used for GONE demographic inference.
gone_chromosomes <- c(
  "SUPER_11", "SUPER_13", "SUPER_14", "SUPER_15", "SUPER_17", "SUPER_18",
  "SUPER_19", "SUPER_20", "SUPER_21", "SUPER_22", "SUPER_23", "SUPER_24",
  "SUPER_25", "SUPER_26", "SUPER_27", "SUPER_28", "SUPER_29"
)

table_s10 <- tibble::tibble(
  snow_owl_chromosome = gone_chromosomes,
  snow_owl_chromosome_number = stringr::str_remove(gone_chromosomes, "SUPER_"),
  homologous_chicken_chromosome = NA_character_,
  chicken_recombination_rate_cm_per_mb = NA_real_,
  selection_criteria = "One-to-one synteny with chicken; no fusion/fission relative to chicken homolog; comparable acrocentric morphology; available chicken chromosome-specific recombination-rate estimate.",
  notes = "Manual curation required: add homologous chicken chromosome and chromosome-specific recombination rate from source table/literature."
)
write_supporting_csv(table_s10, "Table_S10_GONE_chromosomes_used_for_demographic_inference.csv")

# S11: GONE input datasets and sensitivity runs.
table_s11 <- readr::read_csv(
  file.path(paths$results_tables, "gone", "SnowyOwl_GONE_sensitivity_summary.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(
    run_id, sample_set, filter, n_generations, generation_time_years,
    present_year, min_ne_full_output, max_ne_full_output,
    most_recent_ne_full_output, oldest_ne_full_output,
    max_years_before_present_full_output, most_recent_ne_plotted
  ) |>
  round_numeric(3)
write_supporting_csv(table_s11, "Table_S11_GONE_dataset_and_sensitivity_summary.csv")

# S8: ROH summary by individual.
table_s8 <- readr::read_csv(
  file.path(paths$results_tables, "roh", "SnowyOwl_ROH_individual_summary.csv"),
  show_col_types = FALSE
) |>
  dplyr::left_join(
    metadata |> dplyr::select(sample_id, sample_short_id),
    by = c("IID" = "sample_id")
  ) |>
  dplyr::transmute(
    sample_short_id,
    region = region_label,
    roh_count = ROH_count,
    roh_total_mb = roh_total_mb,
    froh_total,
    froh_0.3_1Mb,
    froh_1_5Mb,
    froh_gt5Mb
  ) |>
  round_numeric(6)
write_supporting_csv(table_s8, "Table_S8_ROH_individual_summary.csv")

# S9: ROH summary by sampling region.
table_s9 <- readr::read_csv(
  file.path(paths$results_tables, "roh", "SnowyOwl_ROH_region_summary.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(
    region_label, n_individuals, mean_froh, median_froh, max_froh,
    mean_roh_count, median_roh_count, mean_roh_total_mb, median_roh_total_mb
  ) |>
  dplyr::rename(region = region_label) |>
  round_numeric(6)
write_supporting_csv(table_s9, "Table_S9_ROH_region_summary.csv")

# S6: Genome-wide diversity summaries.
pi_macro_micro <- readr::read_csv(
  file.path(paths$results_tables, "pixy", "SnowyOwl_ALL_pi_macro_micro_summary.csv"),
  show_col_types = FALSE
) |>
  dplyr::mutate(metric = "pi")
td_macro_micro <- readr::read_csv(
  file.path(paths$results_tables, "pixy", "SnowyOwl_ALL_tajima_d_macro_micro_summary.csv"),
  show_col_types = FALSE
) |>
  dplyr::mutate(metric = "Tajima_D")

table_s6 <- dplyr::bind_rows(pi_macro_micro, td_macro_micro) |>
  dplyr::select(metric, chromosome_class, n_windows, callable_sites, mean_unweighted, mean_weighted_by_sites, median) |>
  round_numeric(6)
write_supporting_csv(table_s6, "Table_S6_genome_wide_diversity_macro_micro.csv")

# S7: Centromere/pericentromere diversity summary.
table_s7 <- readr::read_csv(
  file.path(paths$results_tables, "pixy", "SnowyOwl_ALL_pixy_centromere_region_summary_min_sites.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(
    metric, centromere_region, min_no_sites, n_windows,
    n_windows_with_value, callable_sites, mean_unweighted,
    mean_weighted_by_sites, median, min, max
  ) |>
  round_numeric(6)
write_supporting_csv(table_s7, "Table_S7_centromere_pericentromere_diversity.csv")

table_index <- tibble::tribble(
  ~table_id, ~filename, ~title, ~recommended_destination, ~manual_status,
  "Table S1", "Table_S1_sample_metadata_and_analysis_inclusion.csv", "Sample metadata and analysis inclusion overview.", "Supporting information", "Needs manual check: verify sample identifiers, locations, collection years, observed sex, age, tissue, and public metadata fields.",
  "Table S2", "Table_S2_sequencing_and_genotype_QC.csv", "Sequencing and genotype QC summary for all Snowy Owl samples.", "Supporting information", "Needs manual check: confirm QC thresholds, sex-call wording, and exclusion language.",
  "Table S3", "Table_S3_sample_exclusions_and_KING_relatedness.csv", "Pairwise KING kinship estimates above the reporting threshold.", "Supporting information", "Confirm that 0.03 is the preferred reporting threshold; full all-pair KING results are available as online data.",
  "Table S4", "Table_S4_PCA_dataset_summaries.csv", "Summary of GATK and ANGSD PCA datasets and filters.", "Supporting information", "Needs manual check: confirm GATK SNP count/LD-pruning details if exact SNP numbers should be reported.",
  "Table S5", "Table_S5_pairwise_FST_and_dXY.csv", "Pairwise mean FST and dXY among sampling regions.", "Supporting information", "Likely ready after checking population labels and rounding.",
  "Table S6", "Table_S6_genome_wide_diversity_macro_micro.csv", "Genome-wide diversity summaries by macro- and microchromosomes.", "Supporting information", "Manual check: confirm macro/microchromosome definition in final text.",
  "Table S7", "Table_S7_centromere_pericentromere_diversity.csv", "Diversity in centromeric, pericentromeric, and background regions.", "Supporting information", "Manual check: emphasize that centromere annotations are available for only a subset of chromosomes.",
  "Table S8", "Table_S8_ROH_individual_summary.csv", "Individual runs of homozygosity summary.", "Supporting information", "Manual check: confirm autosomal genome size and ROH length classes.",
  "Table S9", "Table_S9_ROH_region_summary.csv", "Regional runs of homozygosity summary.", "Supporting information", "Manual check: decide whether regional summaries or all-sample row should be emphasized.",
  "Table S10", "Table_S10_GONE_chromosomes_used_for_demographic_inference.csv", "Chromosomes used for GONE demographic inference.", "Supporting information", "Manual curation required: add homologous chicken chromosome and chromosome-specific recombination rate; confirm whether 6.63 cM/Mb is arithmetic or length-weighted mean.",
  "Table S11", "Table_S11_GONE_dataset_and_sensitivity_summary.csv", "GONE primary and sensitivity run summary.", "Supporting information", "Needs manual check: confirm primary run and whether implausible all64 runs should remain in the formal table."
)
write_supporting_csv(table_index, "supporting_table_index.csv")

message("Wrote supporting table drafts to ", supporting_dir)
