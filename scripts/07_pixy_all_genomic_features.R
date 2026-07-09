# Relate ALL-sample pixy diversity metrics to genomic features.

source("scripts/00_setup.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

pi_file <- here(paths$pixy_all_pi %||% file.path(paths$pixy_dir, "SnowyOwl_ALL_pi.txt"))
tajima_file <- here(paths$pixy_all_tajima_d %||% file.path(paths$pixy_dir, "SnowyOwl_ALL_tajima_d.txt"))
gff_file <- here(paths$genome_gff %||% "data/external/genome/bBubSca1.1.hap1.gff")
repeatmasker_file <- here(paths$repeatmasker_out %||% "data/external/genome/bBubSca1.1.hap1.fasta.out")
centromere_file <- here(paths$centromere_bed %||% "data/metadata/centromeres/putative_centromeres.bed")

out_tab_dir <- here(paths$results_tables, "pixy")
out_fig_dir <- here(paths$results_figures, "pixy")
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

extract_gff_attr <- function(attributes, key) {
  match <- stringr::str_match(attributes, paste0("(^|;)", key, "=([^;]+)"))[, 3]
  dplyr::na_if(match, "")
}

read_genes <- function(file) {
  gff_lines <- readr::read_lines(file) |>
    stringr::str_subset("^[^#]") |>
    stringr::str_subset("\\t")

  gff_fields <- stringr::str_split_fixed(gff_lines, "\\t", 9)

  tibble::as_tibble(gff_fields, .name_repair = "minimal") |>
    rlang::set_names(c("chromosome", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) |>
    dplyr::filter(.data$type == "gene", stringr::str_detect(.data$chromosome, "^SUPER_")) |>
    dplyr::mutate(
      start = as.numeric(.data$start),
      end = as.numeric(.data$end),
      gene_id = extract_gff_attr(.data$attributes, "ID")
    ) |>
    dplyr::select(chromosome, start, end, gene_id)
}

read_repeats <- function(file) {
  readr::read_table(
    file,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ) |>
    dplyr::transmute(
      chromosome = .data$X5,
      start = as.numeric(.data$X6),
      end = as.numeric(.data$X7),
      repeat_name = .data$X10,
      repeat_class = .data$X11
    ) |>
    dplyr::filter(stringr::str_detect(.data$chromosome, "^SUPER_"), !is.na(.data$start), !is.na(.data$end))
}

read_centromeres <- function(file) {
  first_data_line <- readr::read_lines(file, n_max = 20) |>
    stringr::str_subset("^[^#]") |>
    utils::head(1)

  has_header <- length(first_data_line) == 1 &&
    stringr::str_detect(first_data_line, stringr::regex("^(Chr|chromosome)\\t", ignore_case = TRUE))

  bed <- readr::read_tsv(
    file,
    col_names = has_header,
    comment = "#",
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )

  if (has_header) {
    bed <- bed |>
      dplyr::rename(
        chromosome = dplyr::any_of(c("Chr", "chr", "chrom", "chromosome")),
        start_1based = dplyr::any_of(c("Start", "start", "start_1based")),
        end = dplyr::any_of(c("End", "end")),
        morphology = dplyr::any_of(c("Morphology_new", "Morphology", "morphology"))
      ) |>
      dplyr::mutate(start = as.numeric(.data$start_1based))
  } else {
    names(bed)[seq_len(min(ncol(bed), 5))] <- c("chromosome", "start0", "end", "centromere_id", "morphology")[seq_len(min(ncol(bed), 5))]
    bed <- bed |>
      dplyr::mutate(start = as.numeric(.data$start0) + 1)
  }

  if (!"morphology" %in% names(bed)) {
    bed$morphology <- NA_character_
  }

  bed |>
    dplyr::mutate(end = as.numeric(.data$end)) |>
    dplyr::filter(
      stringr::str_detect(.data$chromosome, "^SUPER_"),
      !is.na(.data$start),
      !is.na(.data$end),
      .data$end >= .data$start
    ) |>
    dplyr::select(chromosome, start, end, morphology)
}

read_pixy_all <- function(pi_file, tajima_file) {
  pi <- readr::read_tsv(pi_file, na = c("NA", "nan", "NaN", "."), show_col_types = FALSE) |>
    dplyr::select(pop, chromosome, window_pos_1, window_pos_2, avg_pi, pi_no_sites = no_sites) |>
    dplyr::filter(stringr::str_detect(.data$chromosome, "^SUPER_"))

  tajima <- readr::read_tsv(tajima_file, na = c("NA", "nan", "NaN", "."), show_col_types = FALSE) |>
    dplyr::select(pop, chromosome, window_pos_1, window_pos_2, tajima_d, tajima_no_sites = no_sites) |>
    dplyr::filter(stringr::str_detect(.data$chromosome, "^SUPER_"))

  pi |>
    dplyr::full_join(tajima, by = c("pop", "chromosome", "window_pos_1", "window_pos_2")) |>
    dplyr::mutate(
      chrom_number = suppressWarnings(
        as.integer(stringr::str_extract(.data$chromosome, "(?<=SUPER_)\\d+"))
      ),
      chromosome_class = chromosome_class(.data$chrom_number),
      window_id = paste(
        .data$chromosome,
        sprintf("%.0f", as.numeric(.data$window_pos_1)),
        sprintf("%.0f", as.numeric(.data$window_pos_2)),
        sep = ":"
      ),
      window_length_bp = .data$window_pos_2 - .data$window_pos_1 + 1,
      window_midpoint = (.data$window_pos_1 + .data$window_pos_2) / 2
    )
}

overlap_rows <- function(windows, features) {
  if (nrow(features) == 0 || nrow(windows) == 0) {
    return(tibble::tibble())
  }

  window_dt <- data.table::as.data.table(
    windows |>
      dplyr::transmute(
        chromosome = as.character(.data$chromosome),
        start = as.numeric(.data$window_pos_1),
        end = as.numeric(.data$window_pos_2),
        window_id = .data$window_id,
        window_start = as.numeric(.data$window_pos_1),
        window_end = as.numeric(.data$window_pos_2)
      )
  )

  feature_dt <- data.table::as.data.table(
    features |>
      dplyr::mutate(
        chromosome = as.character(.data$chromosome),
        start = as.numeric(.data$start),
        end = as.numeric(.data$end)
      )
  )

  data.table::setkey(window_dt, chromosome, start, end)

  data.table::foverlaps(
    feature_dt,
    window_dt,
    by.x = c("chromosome", "start", "end"),
    by.y = c("chromosome", "start", "end"),
    nomatch = 0
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      start = .data$i.start,
      end = .data$i.end
    ) |>
    dplyr::select(-`i.start`, -`i.end`)
}

summarise_gene_features <- function(windows, genes) {
  gene_hits <- overlap_rows(windows, genes)

  if (nrow(gene_hits) == 0) {
    return(tibble::tibble(window_id = windows$window_id, n_genes = 0))
  }

  gene_hits |>
    dplyr::mutate(
      overlap_bp = pmax(
        0,
        pmin(.data$end, .data$window_end) - pmax(.data$start, .data$window_start) + 1
      )
    ) |>
    dplyr::group_by(.data$window_id) |>
    dplyr::summarise(
      n_genes = dplyr::n_distinct(.data$gene_id),
      gene_bp = sum(.data$overlap_bp, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::right_join(tibble::tibble(window_id = windows$window_id), by = "window_id") |>
    dplyr::mutate(
      n_genes = tidyr::replace_na(.data$n_genes, 0L),
      gene_bp = tidyr::replace_na(.data$gene_bp, 0)
    )
}

summarise_repeat_features <- function(windows, repeats) {
  repeat_hits <- overlap_rows(windows, repeats)

  if (nrow(repeat_hits) == 0) {
    return(tibble::tibble(window_id = windows$window_id, n_repeat_hits = 0))
  }

  repeat_hits |>
    dplyr::mutate(
      overlap_bp = pmax(
        0,
        pmin(.data$end, .data$window_end) - pmax(.data$start, .data$window_start) + 1
      )
    ) |>
    dplyr::group_by(.data$window_id) |>
    dplyr::summarise(
      n_repeat_hits = dplyr::n(),
      repeat_bp = sum(.data$overlap_bp, na.rm = TRUE),
      top_repeat_class = names(sort(tapply(.data$overlap_bp, .data$repeat_class, sum), decreasing = TRUE))[1],
      satellite_bp = sum(.data$overlap_bp[.data$repeat_class == "Satellite"], na.rm = TRUE),
      simple_repeat_bp = sum(.data$overlap_bp[.data$repeat_class == "Simple_repeat"], na.rm = TRUE),
      cr1_bp = sum(.data$overlap_bp[stringr::str_detect(.data$repeat_class, "LINE/CR1")], na.rm = TRUE),
      ltr_bp = sum(.data$overlap_bp[stringr::str_detect(.data$repeat_class, "^LTR")], na.rm = TRUE),
      unknown_repeat_bp = sum(.data$overlap_bp[.data$repeat_class == "Unknown"], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::right_join(tibble::tibble(window_id = windows$window_id), by = "window_id") |>
    dplyr::mutate(
      dplyr::across(
        c(n_repeat_hits, repeat_bp, satellite_bp, simple_repeat_bp, cr1_bp, ltr_bp, unknown_repeat_bp),
        ~ tidyr::replace_na(.x, 0)
      ),
      top_repeat_class = tidyr::replace_na(.data$top_repeat_class, "")
    )
}

add_centromere_distance <- function(windows, centromeres) {
  split_centromeres <- split(centromeres, centromeres$chromosome)

  distance_df <- dplyr::bind_rows(lapply(seq_len(nrow(windows)), function(i) {
    chrom <- as.character(windows$chromosome[i])
    cen <- split_centromeres[[chrom]]
    if (is.null(cen)) {
      return(tibble::tibble(
        window_id = windows$window_id[i],
        nearest_centromere_distance_bp = NA_real_,
        nearest_centromere_morphology = NA_character_
      ))
    }

    midpoint <- windows$window_midpoint[i]
    distances <- dplyr::case_when(
      midpoint < cen$start ~ cen$start - midpoint,
      midpoint > cen$end ~ midpoint - cen$end,
      TRUE ~ 0
    )
    nearest <- which.min(distances)

    tibble::tibble(
      window_id = windows$window_id[i],
      nearest_centromere_distance_bp = distances[nearest],
      nearest_centromere_morphology = cen$morphology[nearest]
    )
  }))

  windows |>
    dplyr::left_join(distance_df, by = "window_id")
}

spearman_summary <- function(df, metric, features) {
  dplyr::bind_rows(lapply(features, function(feature) {
    feature_df <- df |>
      dplyr::filter(!is.na(.data[[metric]]), !is.na(.data[[feature]]))

    if (nrow(feature_df) < 3 || stats::sd(feature_df[[feature]], na.rm = TRUE) == 0) {
      return(tibble::tibble(metric = metric, feature = feature, n = nrow(feature_df), rho = NA_real_, p_value = NA_real_))
    }

    test <- suppressWarnings(stats::cor.test(feature_df[[metric]], feature_df[[feature]], method = "spearman", exact = FALSE))
    tibble::tibble(
      metric = metric,
      feature = feature,
      n = nrow(feature_df),
      rho = unname(test$estimate),
      p_value = test$p.value
    )
  }))
}

plot_feature_relationship <- function(df, metric, feature, y_label, x_label, output_stub) {
  p <- df |>
    dplyr::filter(!is.na(.data[[metric]]), !is.na(.data[[feature]])) |>
    ggplot(aes(x = .data[[feature]], y = .data[[metric]], color = .data$chromosome_class)) +
    geom_point(size = 0.35, alpha = 0.25, na.rm = TRUE) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.6, na.rm = TRUE) +
    scale_color_manual(values = c("Macrochromosome" = "lightgoldenrod3", "Microchromosome" = "lightcyan3")) +
    labs(x = x_label, y = y_label, color = NULL) +
    theme(legend.position = "bottom")

  ggsave(file.path(out_fig_dir, paste0(output_stub, ".png")), p, width = 6, height = 4.5, dpi = 350)
  ggsave(file.path(out_fig_dir, paste0(output_stub, ".pdf")), p, width = 6, height = 4.5)
}

pixy_windows <- read_pixy_all(pi_file, tajima_file)
genes <- read_genes(gff_file)
repeats <- read_repeats(repeatmasker_file)
centromeres <- read_centromeres(centromere_file)

window_features <- pixy_windows |>
  dplyr::left_join(summarise_gene_features(pixy_windows, genes), by = "window_id") |>
  dplyr::left_join(summarise_repeat_features(pixy_windows, repeats), by = "window_id") |>
  add_centromere_distance(centromeres) |>
  dplyr::mutate(
    gene_fraction = pmin(.data$gene_bp / .data$window_length_bp, 1),
    repeat_fraction = pmin(.data$repeat_bp / .data$window_length_bp, 1),
    satellite_fraction = pmin(.data$satellite_bp / .data$window_length_bp, 1),
    simple_repeat_fraction = pmin(.data$simple_repeat_bp / .data$window_length_bp, 1),
    cr1_fraction = pmin(.data$cr1_bp / .data$window_length_bp, 1),
    ltr_fraction = pmin(.data$ltr_bp / .data$window_length_bp, 1),
    unknown_repeat_fraction = pmin(.data$unknown_repeat_bp / .data$window_length_bp, 1),
    distance_to_centromere_mb = .data$nearest_centromere_distance_bp / 1e6
  )

feature_columns <- c(
  "gene_fraction",
  "n_genes",
  "repeat_fraction",
  "satellite_fraction",
  "simple_repeat_fraction",
  "cr1_fraction",
  "ltr_fraction",
  "unknown_repeat_fraction",
  "distance_to_centromere_mb",
  "pi_no_sites"
)

correlations_all <- dplyr::bind_rows(
  spearman_summary(window_features, "avg_pi", feature_columns),
  spearman_summary(window_features, "tajima_d", feature_columns)
) |>
  dplyr::mutate(chromosome_class = "All")

correlations_by_class <- window_features |>
  dplyr::group_split(.data$chromosome_class) |>
  lapply(function(class_df) {
    dplyr::bind_rows(
      spearman_summary(class_df, "avg_pi", feature_columns),
      spearman_summary(class_df, "tajima_d", feature_columns)
    ) |>
      dplyr::mutate(chromosome_class = unique(class_df$chromosome_class))
  }) |>
  dplyr::bind_rows()

correlations <- dplyr::bind_rows(correlations_all, correlations_by_class) |>
  dplyr::arrange(.data$metric, .data$chromosome_class, dplyr::desc(abs(.data$rho)))

feature_class_summary <- window_features |>
  dplyr::group_by(.data$chromosome_class) |>
  dplyr::summarise(
    n_windows = dplyr::n(),
    mean_pi = mean(.data$avg_pi, na.rm = TRUE),
    mean_tajima_d = mean(.data$tajima_d, na.rm = TRUE),
    mean_gene_fraction = mean(.data$gene_fraction, na.rm = TRUE),
    mean_repeat_fraction = mean(.data$repeat_fraction, na.rm = TRUE),
    mean_satellite_fraction = mean(.data$satellite_fraction, na.rm = TRUE),
    mean_cr1_fraction = mean(.data$cr1_fraction, na.rm = TRUE),
    median_distance_to_centromere_mb = median(.data$distance_to_centromere_mb, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  window_features,
  file.path(out_tab_dir, "SnowyOwl_ALL_pixy_window_genomic_features.csv")
)

readr::write_csv(
  correlations,
  file.path(out_tab_dir, "SnowyOwl_ALL_pixy_feature_spearman_correlations.csv")
)

readr::write_csv(
  feature_class_summary,
  file.path(out_tab_dir, "SnowyOwl_ALL_pixy_feature_macro_micro_summary.csv")
)

plot_feature_relationship(
  window_features,
  metric = "avg_pi",
  feature = "repeat_fraction",
  y_label = expression(pi),
  x_label = "Repeat fraction",
  output_stub = "SnowyOwl_ALL_pi_vs_repeat_fraction"
)

plot_feature_relationship(
  window_features,
  metric = "tajima_d",
  feature = "repeat_fraction",
  y_label = "Tajima's D",
  x_label = "Repeat fraction",
  output_stub = "SnowyOwl_ALL_tajima_d_vs_repeat_fraction"
)

plot_feature_relationship(
  window_features,
  metric = "avg_pi",
  feature = "gene_fraction",
  y_label = expression(pi),
  x_label = "Gene fraction",
  output_stub = "SnowyOwl_ALL_pi_vs_gene_fraction"
)

plot_feature_relationship(
  window_features,
  metric = "avg_pi",
  feature = "distance_to_centromere_mb",
  y_label = expression(pi),
  x_label = "Distance to nearest annotated centromere (Mb)",
  output_stub = "SnowyOwl_ALL_pi_vs_centromere_distance"
)

message("Wrote ALL pixy genomic feature tables to: ", out_tab_dir)
message("Wrote ALL pixy genomic feature figures to: ", out_fig_dir)
