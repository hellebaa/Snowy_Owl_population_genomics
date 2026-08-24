# Plot PCA and KING relatedness summaries from external PCA outputs.

source("scripts/00_setup.R")

required_pca_packages <- c("data.table", "pheatmap")
missing_pca_packages <- required_pca_packages[
  !vapply(required_pca_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pca_packages) > 0) {
  stop(
    "Install missing packages before running PCA plotting: ",
    paste(missing_pca_packages, collapse = ", "),
    call. = FALSE
  )
}

find_input <- function(primary, fallback, required = TRUE) {
  if (!is.null(primary) && file.exists(primary)) {
    primary
  } else if (!is.null(fallback) && file.exists(fallback)) {
    fallback
  } else if (required) {
    stop("Could not find input: ", primary, " or fallback: ", fallback, call. = FALSE)
  } else {
    NA_character_
  }
}

old_dir <- file.path("..", "old")
old_pca_dir <- file.path(old_dir, "pca")

fig_dir <- file.path(paths$results_figures, "pca")
tab_dir <- file.path(paths$results_tables, "pca")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- readr::read_csv(paths$metadata, show_col_types = FALSE) |>
  dplyr::mutate(population_final = factor(population_final, levels = names(pop_colors)))

rename_map <- readr::read_table(
  file.path(old_dir, "sample_rename.txt"),
  col_names = c("old_id", "sample_id"),
  col_types = readr::cols(.default = readr::col_character())
)

read_pca <- function(eigenvec_file, eigenval_file = NA_character_, use_rename_map = FALSE) {
  pca_raw <- readr::read_table(
    eigenvec_file,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_double(), X1 = readr::col_character(), X2 = readr::col_character())
  )

  pc_cols <- paste0("PC", seq_len(ncol(pca_raw) - 2))
  names(pca_raw) <- c("fid", "input_id", pc_cols)

  pca <- pca_raw |>
    dplyr::mutate(old_id = input_id)

  if (use_rename_map) {
    pca <- pca |>
      dplyr::left_join(rename_map, by = "old_id") |>
      dplyr::mutate(sample_id = dplyr::coalesce(sample_id, old_id))
  } else {
    pca <- pca |>
      dplyr::mutate(sample_id = input_id)
  }

  pca <- pca |>
    dplyr::left_join(metadata, by = "sample_id") |>
    dplyr::mutate(
      population_final = dplyr::coalesce(as.character(population_final), region_display(stringr::str_sub(sample_id, 1, 3))),
      population_final = factor(population_final, levels = names(pop_colors)),
      pca_label = stringr::str_extract(sample_id, "^[^_]+")
    )

  if (!is.na(eigenval_file) && file.exists(eigenval_file)) {
    eigenval <- readr::read_table(
      eigenval_file,
      col_names = "eigenvalue",
      col_types = readr::cols(eigenvalue = readr::col_double())
    ) |>
      dplyr::mutate(pc = paste0("PC", dplyr::row_number()), pve = eigenvalue / sum(eigenvalue) * 100)
  } else {
    eigenval <- tibble::tibble(pc = pc_cols, pve = NA_real_)
  }

  list(scores = pca, eigenval = eigenval, pc_cols = pc_cols)
}

axis_label <- function(pc, eigenval) {
  pve <- eigenval$pve[match(pc, eigenval$pc)]
  if (length(pve) == 1 && !is.na(pve)) {
    paste0(pc, " (", round(pve, 1), "%)")
  } else {
    pc
  }
}

plot_pca_pair <- function(pca_obj, x_pc, y_pc, filename_stub) {
  plot_data <- pca_obj$scores |>
    dplyr::filter(!is.na(.data[[x_pc]]), !is.na(.data[[y_pc]]))

  p <- ggplot(plot_data, aes(x = .data[[x_pc]], y = .data[[y_pc]], color = population_final)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey80") +
    geom_point(size = 2.8, alpha = 0.9) +
    scale_color_manual(values = pop_colors, drop = FALSE, na.value = "grey45") +
    labs(
      x = axis_label(x_pc, pca_obj$eigenval),
      y = axis_label(y_pc, pca_obj$eigenval),
      color = "Region"
    ) +
    coord_equal() +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    )

  ggsave(file.path(fig_dir, paste0(filename_stub, ".png")), p, width = 5.8, height = 4.6, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename_stub, ".pdf")), p, width = 5.8, height = 4.6)
}

plot_scree <- function(pca_obj, filename_stub) {
  if (!any(!is.na(pca_obj$eigenval$pve))) {
    return(invisible(NULL))
  }

  scree <- pca_obj$eigenval |>
    dplyr::slice_head(n = 10) |>
    dplyr::mutate(pc = factor(pc, levels = pc))

  p_scree <- ggplot(scree, aes(x = pc, y = pve)) +
    geom_col(fill = "grey45", width = 0.7) +
    labs(x = NULL, y = "Variance explained (%)")

  ggsave(file.path(fig_dir, paste0(filename_stub, ".png")), p_scree, width = 5, height = 3.5, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename_stub, ".pdf")), p_scree, width = 5, height = 3.5)
}

all_pca <- read_pca(
  find_input(paths$pca_eigenvec, file.path(old_pca_dir, "autosomes_pca.eigenvec")),
  find_input(paths$pca_eigenval, NULL, required = FALSE),
  use_rename_map = FALSE
)

plot_pca_pair(all_pca, "PC1", "PC2", "SnowyOwl_PCA_autosomes_all_samples_PC1_PC2")
plot_pca_pair(all_pca, "PC1", "PC3", "SnowyOwl_PCA_autosomes_all_samples_PC1_PC3")
plot_pca_pair(all_pca, "PC3", "PC4", "SnowyOwl_PCA_autosomes_all_samples_PC3_PC4")
plot_scree(all_pca, "SnowyOwl_PCA_autosomes_all_samples_scree")

readr::write_csv(
  all_pca$scores |> dplyr::select(sample_id, population_final, dplyr::all_of(all_pca$pc_cols)),
  file.path(tab_dir, "SnowyOwl_PCA_autosomes_all_samples_scores.csv")
)

filtered_pca_file <- find_input(
  paths$pca_filtered_eigenvec,
  file.path(old_dir, "re", "autosomal_rm_ind.eigenvec"),
  required = FALSE
)

if (!is.na(filtered_pca_file)) {
  filtered_pca <- read_pca(
    filtered_pca_file,
    find_input(paths$pca_filtered_eigenval, file.path(old_dir, "re", "autosomal_rm_ind.eigenval"), required = FALSE),
    use_rename_map = TRUE
  )

  plot_pca_pair(filtered_pca, "PC1", "PC2", "SnowyOwl_PCA_autosomes_filtered_PC1_PC2")
  plot_pca_pair(filtered_pca, "PC1", "PC3", "SnowyOwl_PCA_autosomes_filtered_PC1_PC3")
  plot_pca_pair(filtered_pca, "PC3", "PC4", "SnowyOwl_PCA_autosomes_filtered_PC3_PC4")
  plot_scree(filtered_pca, "SnowyOwl_PCA_autosomes_filtered_scree")

  readr::write_csv(
    filtered_pca$scores |> dplyr::select(sample_id, population_final, dplyr::all_of(filtered_pca$pc_cols)),
    file.path(tab_dir, "SnowyOwl_PCA_autosomes_filtered_scores.csv")
  )
}

king_file <- find_input(paths$pca_king_matrix, file.path(old_pca_dir, "autosomes_clean_pruned.king"))
king_id_file <- find_input(paths$pca_king_ids, file.path(old_pca_dir, "autosomes_clean_pruned.king.id"))

king_ids <- readr::read_table(
  king_id_file,
  col_types = readr::cols(.default = readr::col_character())
) |>
  dplyr::rename(fid = `#FID`, sample_id = IID)

king <- as.matrix(readr::read_table(
  king_file,
  col_names = FALSE,
  col_types = readr::cols(.default = readr::col_double())
))

stopifnot(nrow(king) == nrow(king_ids), ncol(king) == nrow(king_ids))

sample_order <- king_ids |>
  dplyr::left_join(metadata, by = "sample_id") |>
  dplyr::mutate(
    population_final = factor(population_final, levels = names(pop_colors)),
    short_id = stringr::str_extract(sample_id, "^[^_]+")
  ) |>
  dplyr::arrange(population_final, short_id) |>
  dplyr::pull(sample_id)

rownames(king) <- king_ids$sample_id
colnames(king) <- king_ids$sample_id
king <- king[sample_order, sample_order]

king_labels <- stringr::str_extract(rownames(king), "^[^_]+")
rownames(king) <- king_labels
colnames(king) <- king_labels

ann <- metadata |>
  dplyr::filter(sample_id %in% sample_order) |>
  dplyr::mutate(short_id = stringr::str_extract(sample_id, "^[^_]+")) |>
  dplyr::arrange(factor(population_final, levels = names(pop_colors)), short_id) |>
  dplyr::select(short_id, Region = population_final) |>
  tibble::column_to_rownames("short_id")

png(file.path(fig_dir, "SnowyOwl_KING_kinship_heatmap.png"), width = 2200, height = 2000, res = 300)
pheatmap::pheatmap(
  king,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann,
  annotation_col = ann,
  annotation_colors = list(Region = pop_colors),
  color = colorRampPalette(c("white", "lightgoldenrod2", "firebrick3"))(100),
  border_color = NA,
  fontsize = 5,
  main = "KING kinship"
)
dev.off()

pdf(file.path(fig_dir, "SnowyOwl_KING_kinship_heatmap.pdf"), width = 7.2, height = 6.8)
pheatmap::pheatmap(
  king,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann,
  annotation_col = ann,
  annotation_colors = list(Region = pop_colors),
  color = colorRampPalette(c("white", "lightgoldenrod2", "firebrick3"))(100),
  border_color = NA,
  fontsize = 5,
  main = "KING kinship"
)
dev.off()

message("Wrote PCA and KING figures to ", fig_dir)

angsd_sample_file <- find_input(
  paths$angsd_pca_sample_names,
  file.path(old_pca_dir, "PCA_INPUT", "sample_names_pca.txt"),
  required = FALSE
)

read_angsd_haplo_pca <- function(
    haplo_file,
    sample_file,
    dataset_label,
    output_stub,
    keep_sample_ids = NULL,
    min_mac = 1,
    n_sites_report = TRUE) {
  if (is.na(haplo_file) || is.na(sample_file)) {
    return(invisible(NULL))
  }

  samples <- readr::read_lines(sample_file)

  haplo <- data.table::fread(
    haplo_file,
    header = FALSE,
    data.table = FALSE,
    showProgress = FALSE
  )

  if (ncol(haplo) != length(samples) + 3) {
    stop(
      "ANGSD haplo sample count mismatch for ", dataset_label, ": ",
      ncol(haplo) - 3, " genotype columns but ", length(samples), " sample names.",
      call. = FALSE
    )
  }

  geno_char <- as.matrix(haplo[, -(1:3), drop = FALSE])
  storage.mode(geno_char) <- "character"

  if (!is.null(keep_sample_ids)) {
    keep_samples <- samples %in% keep_sample_ids
    if (sum(keep_samples) < 3) {
      stop("Fewer than three ANGSD samples retained for ", dataset_label, call. = FALSE)
    }
    geno_char <- geno_char[, keep_samples, drop = FALSE]
    samples <- samples[keep_samples]
  }

  recode_site <- function(alleles) {
    alleles[alleles == "N" | alleles == "0" | alleles == "." | alleles == ""] <- NA_character_
    counts <- sort(table(alleles, useNA = "no"), decreasing = TRUE)
    if (length(counts) < 2) {
      return(rep(NA_real_, length(alleles)))
    }
    if (length(counts) > 2) {
      return(rep(NA_real_, length(alleles)))
    }

    major <- names(counts)[1]
    minor <- names(counts)[2]
    recoded <- rep(NA_real_, length(alleles))
    recoded[alleles == major] <- 0
    recoded[alleles == minor] <- 1
    recoded
  }

  geno <- t(apply(geno_char, 1, recode_site))
  allele_counts <- rowSums(geno == 1, na.rm = TRUE)
  called_counts <- rowSums(!is.na(geno))
  minor_counts <- pmin(allele_counts, called_counts - allele_counts)
  keep_sites <- called_counts >= 3 & minor_counts >= min_mac
  geno <- geno[keep_sites, , drop = FALSE]

  site_means <- rowMeans(geno, na.rm = TRUE)
  centered <- sweep(geno, 1, site_means, "-")
  centered[is.na(centered)] <- 0

  sample_cov <- crossprod(centered) / (nrow(centered) - 1)
  eig <- eigen(sample_cov, symmetric = TRUE)

  pve <- eig$values / sum(eig$values) * 100
  scores <- as.data.frame(eig$vectors[, seq_len(min(10, ncol(eig$vectors))), drop = FALSE])
  names(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores <- scores |>
    dplyr::mutate(sample_id = samples, .before = 1) |>
    dplyr::left_join(metadata, by = "sample_id") |>
    dplyr::mutate(
      population_final = dplyr::coalesce(as.character(population_final), region_display(stringr::str_sub(sample_id, 1, 3))),
      population_final = factor(population_final, levels = names(pop_colors))
    )

  eigen_out <- tibble::tibble(
    pc = paste0("PC", seq_along(eig$values)),
    eigenvalue = eig$values,
    pve = pve,
    dataset = dataset_label,
    retained_samples = length(samples),
    min_mac = min_mac,
    retained_biallelic_sites = nrow(centered)
  )

  readr::write_csv(
    scores,
    file.path(tab_dir, paste0(output_stub, "_scores.csv"))
  )
  readr::write_csv(
    eigen_out,
    file.path(tab_dir, paste0(output_stub, "_eigenvalues.csv"))
  )

  p <- ggplot(scores, aes(x = PC1, y = PC2, color = population_final)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey80") +
    geom_point(size = 2.8, alpha = 0.9) +
    scale_color_manual(values = pop_colors, drop = FALSE, na.value = "grey45") +
    labs(
      x = paste0("PC1 (", round(pve[1], 1), "%)"),
      y = paste0("PC2 (", round(pve[2], 1), "%)"),
      color = "Region"
    ) +
    coord_equal() +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    )

  ggsave(file.path(fig_dir, paste0(output_stub, "_PC1_PC2.png")), p, width = 5.8, height = 4.6, dpi = 300)
  ggsave(file.path(fig_dir, paste0(output_stub, "_PC1_PC2.pdf")), p, width = 5.8, height = 4.6)

  if (n_sites_report) {
    message(dataset_label, ": retained ", length(samples), " samples and ", nrow(centered), " biallelic polymorphic sites with MAC >= ", min_mac, ".")
  }

  invisible(scores)
}

angsd_autosomes_file <- find_input(
  paths$angsd_autosomes_haplo_thinned,
  file.path(old_pca_dir, "PCA_INPUT", "snowyowl_autosomes_thinned.haplo"),
  required = FALSE
)

read_angsd_haplo_pca(
  angsd_autosomes_file,
  angsd_sample_file,
  dataset_label = "ANGSD pseudohaploid autosomes",
  output_stub = "SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes",
  min_mac = 1
)

angsd_qc_sample_ids <- metadata |>
  dplyr::filter(
    include_pca_filtered,
    qc_flags == "OK",
    relatedness_flag == "keep",
    contamination_flag == "keep"
  ) |>
  dplyr::pull(sample_id)

read_angsd_haplo_pca(
  angsd_autosomes_file,
  angsd_sample_file,
  dataset_label = "ANGSD pseudohaploid autosomes, QC-filtered samples",
  output_stub = "SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered",
  keep_sample_ids = angsd_qc_sample_ids,
  min_mac = 1
)

read_angsd_haplo_pca(
  angsd_autosomes_file,
  angsd_sample_file,
  dataset_label = paste0("ANGSD pseudohaploid autosomes, QC-filtered samples, MAC >= ", paths$angsd_pca_min_mac),
  output_stub = paste0("SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_mac", paths$angsd_pca_min_mac),
  keep_sample_ids = angsd_qc_sample_ids,
  min_mac = paths$angsd_pca_min_mac
)

run_angsd_sex_chromosome_pca <- isTRUE(paths$run_angsd_sex_chromosome_pca) ||
  tolower(Sys.getenv("RUN_ANGSD_SEX_CHR_PCA", unset = "false")) %in% c("true", "1", "yes")

if (run_angsd_sex_chromosome_pca) {
  angsd_z_file <- find_input(paths$angsd_z_haplo, file.path(old_pca_dir, "PCA_INPUT", "snowyowl_Z.haplo"), required = FALSE)
  angsd_w_file <- find_input(paths$angsd_w_haplo, file.path(old_pca_dir, "PCA_INPUT", "snowyowl_W.haplo"), required = FALSE)

  read_angsd_haplo_pca(
    angsd_z_file,
    angsd_sample_file,
    dataset_label = "ANGSD pseudohaploid Z chromosome",
    output_stub = "SnowyOwl_PCA_ANGSD_pseudohaploid_Z",
    n_sites_report = FALSE
  )

  read_angsd_haplo_pca(
    angsd_w_file,
    angsd_sample_file,
    dataset_label = "ANGSD pseudohaploid W chromosome",
    output_stub = "SnowyOwl_PCA_ANGSD_pseudohaploid_W",
    n_sites_report = FALSE
  )
}
