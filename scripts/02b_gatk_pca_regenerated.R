# Plot regenerated GATK/PLINK autosomal PCA runs.

source("scripts/00_setup.R")

gatk_dir <- file.path(paths$pca_dir, "gatk")
fig_dir <- file.path(paths$results_figures, "pca")
tab_dir <- file.path(paths$results_tables, "pca")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- readr::read_csv(paths$metadata, show_col_types = FALSE) |>
  dplyr::mutate(population_final = factor(population_final, levels = names(pop_colors)))

gatk_runs <- tibble::tribble(
  ~run_id, ~run_label, ~run_dir, ~eigenvec, ~eigenval, ~pca_log, ~prune_file,
  "runA_mac3", "Run A: QC/relatedness filter, MAC >= 3",
  "runA_qc_relatedness_mac3", "autosomes_clean_runA_mac3_pca.eigenvec",
  "autosomes_clean_runA_mac3_pca.eigenval", "autosomes_clean_runA_mac3_pca.log",
  "autosomes_clean_runA_mac3_pruned.prune.in",
  "runA_no_mac", "Run A: QC/relatedness filter, no MAC filter",
  "runA_qc_relatedness_no_mac", "autosomes_clean_runA_no_mac_pca.eigenvec",
  "autosomes_clean_runA_no_mac_pca.eigenval", "autosomes_clean_runA_no_mac_pca.log",
  "autosomes_clean_runA_no_mac_pruned.prune.in",
  "runB_mac3", "Run B: contaminated samples removed, MAC >= 3",
  "runB_contaminated_only_mac3", "autosomes_clean_runB_mac3_pca.eigenvec",
  "autosomes_clean_runB_mac3_pca.eigenval", "autosomes_clean_runB_mac3_pca.log",
  "autosomes_clean_runB_mac3_pruned.prune.in",
  "runB_no_mac", "Run B: contaminated samples removed, no MAC filter",
  "runB_contaminated_only_no_mac", "autosomes_clean_runB_no_mac_pca.eigenvec",
  "autosomes_clean_runB_no_mac_pca.eigenval", "autosomes_clean_runB_no_mac_pca.log",
  "autosomes_clean_runB_no_mac_pruned.prune.in"
)

read_plink_pca <- function(run_row) {
  eigenvec_file <- file.path(gatk_dir, run_row$run_dir, run_row$eigenvec)
  eigenval_file <- file.path(gatk_dir, run_row$run_dir, run_row$eigenval)

  if (!file.exists(eigenvec_file) || !file.exists(eigenval_file)) {
    stop("Missing PCA files for ", run_row$run_id, call. = FALSE)
  }

  pca_raw <- readr::read_table(
    eigenvec_file,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_double(), X1 = readr::col_character(), X2 = readr::col_character())
  )
  pc_cols <- paste0("PC", seq_len(ncol(pca_raw) - 2))
  names(pca_raw) <- c("fid", "sample_id", pc_cols)

  scores <- pca_raw |>
    dplyr::left_join(metadata, by = "sample_id") |>
    dplyr::mutate(
      population_final = factor(population_final, levels = names(pop_colors)),
      pca_label = stringr::str_extract(sample_id, "^[^_]+"),
      run_id = run_row$run_id,
      run_label = run_row$run_label
    )

  eigen <- readr::read_table(
    eigenval_file,
    col_names = "eigenvalue",
    col_types = readr::cols(eigenvalue = readr::col_double())
  ) |>
    dplyr::mutate(
      pc = paste0("PC", dplyr::row_number()),
      pve = eigenvalue / sum(eigenvalue) * 100,
      run_id = run_row$run_id,
      run_label = run_row$run_label
    )

  list(scores = scores, eigen = eigen, pc_cols = pc_cols)
}

count_lines <- function(file) {
  if (!file.exists(file)) {
    return(NA_integer_)
  }
  length(readr::read_lines(file, progress = FALSE))
}

count_unique_lines <- function(file) {
  if (!file.exists(file)) {
    return(NA_integer_)
  }
  length(unique(readr::read_lines(file, progress = FALSE)))
}

log_value <- function(file, pattern, fallback = NA_integer_) {
  if (!file.exists(file)) {
    return(fallback)
  }
  lines <- readr::read_lines(file, progress = FALSE)
  hit <- stringr::str_match(lines, pattern)
  hit <- hit[!is.na(hit[, 2]), , drop = FALSE]
  if (nrow(hit) == 0) {
    return(fallback)
  }
  as.integer(hit[1, 2])
}

axis_label <- function(pc, eigen) {
  pve <- eigen$pve[match(pc, eigen$pc)]
  paste0(pc, " (", round(pve, 1), "%)")
}

plot_single_run <- function(scores, eigen, x_pc, y_pc, filename_stub) {
  p <- ggplot(scores, aes(x = .data[[x_pc]], y = .data[[y_pc]], color = population_final)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey80") +
    geom_point(size = 2.8, alpha = 0.9) +
    scale_color_manual(values = pop_colors, drop = FALSE, na.value = "grey45") +
    labs(
      x = axis_label(x_pc, eigen),
      y = axis_label(y_pc, eigen),
      color = "Region"
    ) +
    coord_equal() +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    )

  ggsave(file.path(fig_dir, paste0(filename_stub, ".png")), p, width = 5.8, height = 4.6, dpi = 350)
  ggsave(file.path(fig_dir, paste0(filename_stub, ".pdf")), p, width = 5.8, height = 4.6)
  invisible(p)
}

plot_faceted <- function(scores_all, eigen_all, x_pc, y_pc, filename_stub) {
  axis_df <- eigen_all |>
    dplyr::filter(pc %in% c(x_pc, y_pc)) |>
    dplyr::select(run_id, pc, pve) |>
    tidyr::pivot_wider(names_from = pc, values_from = pve)

  plot_data <- scores_all |>
    dplyr::left_join(axis_df, by = "run_id") |>
    dplyr::mutate(
      facet_label = paste0(
        run_label, "\n",
        x_pc, " ", round(.data[[x_pc.y]], 1), "%; ",
        y_pc, " ", round(.data[[y_pc.y]], 1), "%"
      )
    )
}

# Build combined score/eigenvalue tables.
pca_objects <- lapply(seq_len(nrow(gatk_runs)), function(i) read_plink_pca(gatk_runs[i, ]))
names(pca_objects) <- gatk_runs$run_id

scores_all <- dplyr::bind_rows(lapply(pca_objects, `[[`, "scores"))
eigen_all <- dplyr::bind_rows(lapply(pca_objects, `[[`, "eigen"))

readr::write_csv(
  scores_all |>
    dplyr::select(run_id, run_label, sample_id, population_final, dplyr::starts_with("PC")),
  file.path(tab_dir, "SnowyOwl_PCA_GATK_regenerated_scores.csv")
)
readr::write_csv(
  eigen_all,
  file.path(tab_dir, "SnowyOwl_PCA_GATK_regenerated_eigenvalues.csv")
)

for (run_id in names(pca_objects)) {
  obj <- pca_objects[[run_id]]
  plot_single_run(
    obj$scores, obj$eigen, "PC1", "PC2",
    paste0("SnowyOwl_PCA_GATK_", run_id, "_PC1_PC2")
  )
  plot_single_run(
    obj$scores, obj$eigen, "PC1", "PC3",
    paste0("SnowyOwl_PCA_GATK_", run_id, "_PC1_PC3")
  )
  plot_single_run(
    obj$scores, obj$eigen, "PC3", "PC4",
    paste0("SnowyOwl_PCA_GATK_", run_id, "_PC3_PC4")
  )
}

axis_summary <- eigen_all |>
  dplyr::filter(pc %in% c("PC1", "PC2")) |>
  dplyr::select(run_id, pc, pve) |>
  tidyr::pivot_wider(names_from = pc, values_from = pve)

facet_scores <- scores_all |>
  dplyr::left_join(axis_summary, by = "run_id") |>
  dplyr::mutate(
    facet_label = paste0(run_label, "\nPC1 ", round(PC1.y, 1), "%; PC2 ", round(PC2.y, 1), "%")
  )

p_faceted <- ggplot(facet_scores, aes(x = PC1.x, y = PC2.x, color = population_final)) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey82") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey82") +
  geom_point(size = 2.1, alpha = 0.9) +
  facet_wrap(~ facet_label, scales = "free", ncol = 2) +
  scale_color_manual(values = pop_colors, drop = FALSE, na.value = "grey45") +
  labs(x = "PC1", y = "PC2", color = "Region") +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.text = element_text(size = 8.5)
  )

ggsave(file.path(fig_dir, "SnowyOwl_PCA_GATK_regenerated_PC1_PC2_faceted.png"), p_faceted, width = 8.2, height = 7.0, dpi = 350)
ggsave(file.path(fig_dir, "SnowyOwl_PCA_GATK_regenerated_PC1_PC2_faceted.pdf"), p_faceted, width = 8.2, height = 7.0)

summary_table <- gatk_runs |>
  dplyr::rowwise() |>
  dplyr::mutate(
    pca_log_file = file.path(gatk_dir, run_dir, pca_log),
    prune_file_path = file.path(gatk_dir, run_dir, prune_file),
    retained_samples = count_lines(file.path(gatk_dir, run_dir, eigenvec)),
    retained_pcs = count_lines(file.path(gatk_dir, run_dir, eigenval)),
    pruned_variants = count_lines(prune_file_path),
    unique_pruned_variant_ids = count_unique_lines(prune_file_path),
    pca_variants = log_value(pca_log_file, "([0-9]+) variants and [0-9]+ people pass filters and QC"),
    pca_samples = log_value(pca_log_file, "[0-9]+ variants and ([0-9]+) people pass filters and QC"),
    pc1_pve = eigen_all$pve[eigen_all$run_id == run_id & eigen_all$pc == "PC1"][1],
    pc2_pve = eigen_all$pve[eigen_all$run_id == run_id & eigen_all$pc == "PC2"][1]
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    run_id, run_label, retained_samples, retained_pcs, pca_variants, pca_samples,
    pruned_variants, unique_pruned_variant_ids, pc1_pve, pc2_pve
  )

readr::write_csv(summary_table, file.path(tab_dir, "SnowyOwl_PCA_GATK_regenerated_run_summary.csv"))

message("Wrote regenerated GATK PCA plots to ", fig_dir)
message("Wrote regenerated GATK PCA tables to ", tab_dir)
