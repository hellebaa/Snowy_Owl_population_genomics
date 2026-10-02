# Explore regional differentiation from pixy FST and dXY outputs.

source("scripts/00_setup.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

fst_file <- here(paths$pixy_regions_fst %||% file.path(paths$pixy_dir, "SnowyOwl_regions_fst.txt"))
dxy_file <- here(paths$pixy_regions_dxy %||% file.path(paths$pixy_dir, "SnowyOwl_regions_dxy.txt"))
min_fst_snps <- paths$min_fst_snps %||% 10
fst_outlier_quantile <- paths$fst_outlier_quantile %||% 0.999

out_fig_dir <- here(paths$results_figures, "pixy")
out_tab_dir <- here(paths$results_tables, "pixy")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)

read_pixy_tsv <- function(file) {
  if (!file.exists(file)) {
    stop("Missing expected pixy file: ", file, call. = FALSE)
  }

  readr::read_tsv(file, na = c("NA", "nan", "NaN", "."), show_col_types = FALSE) |>
    dplyr::filter(stringr::str_detect(.data$chromosome, "^SUPER")) |>
    dplyr::mutate(
      chromosome = as.character(.data$chromosome),
      chrom_number = suppressWarnings(
        as.integer(stringr::str_extract(.data$chromosome, "(?<=SUPER_)\\d+"))
      ),
      chromosome_class = chromosome_class(.data$chrom_number),
      chromosome_landscape_group = chromosome_landscape_group(.data$chrom_number),
      window_midpoint = (.data$window_pos_1 + .data$window_pos_2) / 2,
      pair = paste(.data$pop1, .data$pop2, sep = "_vs_")
    )
}

chromosome_order <- function(chromosomes) {
  tibble::tibble(chromosome = unique(chromosomes)) |>
    dplyr::mutate(
      chrom_number = suppressWarnings(
        as.integer(stringr::str_extract(.data$chromosome, "(?<=SUPER_)\\d+"))
      )
    ) |>
    dplyr::arrange(is.na(.data$chrom_number), .data$chrom_number, .data$chromosome) |>
    dplyr::pull(.data$chromosome)
}

add_genome_position <- function(df) {
  chrom_levels <- chromosome_order(df$chromosome)

  chrom_lengths <- df |>
    dplyr::mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) |>
    dplyr::group_by(.data$chromosome) |>
    dplyr::summarise(chrom_length = max(.data$window_pos_2, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$chromosome) |>
    dplyr::mutate(chrom_start = dplyr::lag(cumsum(.data$chrom_length), default = 0))

  df |>
    dplyr::mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) |>
    dplyr::left_join(chrom_lengths, by = "chromosome") |>
    dplyr::mutate(genome_midpoint = .data$chrom_start + .data$window_midpoint)
}

plot_pairwise_fst_heatmap <- function(pairwise_summary) {
  ordered_pops <- c("FNM", "GRL", "NYS", "SKW", "WRG")
  ordered_labels <- region_display(ordered_pops)

  pairwise_symmetric <- pairwise_summary |>
    dplyr::select(pop1, pop2, mean_fst_weighted_by_snps) |>
    dplyr::bind_rows(
      tibble::tibble(
        pop1 = pairwise_summary$pop2,
        pop2 = pairwise_summary$pop1,
        mean_fst_weighted_by_snps = pairwise_summary$mean_fst_weighted_by_snps
      )
    ) |>
    dplyr::bind_rows(
      tibble::tibble(
        pop1 = ordered_pops,
        pop2 = ordered_pops,
        mean_fst_weighted_by_snps = 0
      )
    )

  heatmap_df <- tidyr::expand_grid(pop1 = ordered_pops, pop2 = ordered_pops) |>
    dplyr::left_join(
      pairwise_symmetric,
      by = c("pop1", "pop2")
    ) |>
    dplyr::mutate(
      pop1_index = match(.data$pop1, ordered_pops),
      pop2_index = match(.data$pop2, ordered_pops),
      display_fst = pmax(.data$mean_fst_weighted_by_snps, 0, na.rm = TRUE),
      pop1 = factor(.data$pop1, levels = ordered_pops),
      pop2 = factor(.data$pop2, levels = ordered_pops)
    ) |>
    dplyr::filter(.data$pop2_index >= .data$pop1_index)

  heatmap_max <- max(heatmap_df$display_fst, na.rm = TRUE)

  p <- heatmap_df |>
    ggplot(aes(x = .data$pop1, y = .data$pop2, fill = .data$display_fst)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.4f", .data$display_fst)), size = 3) +
    scale_fill_gradient(
      low = "#ADD7E4",
      high = "#6D878E",
      limits = c(0, heatmap_max),
      name = "Mean FST"
    ) +
    scale_x_discrete(labels = ordered_labels, drop = FALSE) +
    scale_y_discrete(labels = ordered_labels, drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Pairwise regional FST") +
    coord_equal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(out_fig_dir, "SnowyOwl_regions_pairwise_fst_heatmap.png"), p, width = 5.5, height = 4.8, dpi = 350)
  ggsave(file.path(out_fig_dir, "SnowyOwl_regions_pairwise_fst_heatmap.pdf"), p, width = 5.5, height = 4.8)

  p
}

plot_genomewide_max_fst <- function(max_fst_windows) {
  axis_df <- max_fst_windows |>
    dplyr::distinct(.data$chromosome, .data$chrom_start, .data$chrom_length) |>
    dplyr::mutate(
      axis_midpoint = .data$chrom_start + (.data$chrom_length / 2),
      axis_label = stringr::str_remove(as.character(.data$chromosome), "^SUPER_")
    )

  p <- ggplot(
    max_fst_windows,
    aes(
      x = .data$genome_midpoint,
      y = .data$max_fst,
      color = .data$chromosome_landscape_group
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_point(size = 0.45, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(values = chromosome_landscape_colors, guide = "none") +
    scale_x_continuous(
      breaks = axis_df$axis_midpoint,
      labels = axis_df$axis_label,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    labs(
      x = "Chromosome",
      y = "Maximum pairwise FST",
      title = "Genome-wide maximum regional FST"
    ) +
    theme(
      axis.text.x = element_text(size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_fig_dir, "SnowyOwl_regions_max_fst_genomewide.png"), p, width = 11, height = 4.8, dpi = 350)
  ggsave(file.path(out_fig_dir, "SnowyOwl_regions_max_fst_genomewide.pdf"), p, width = 11, height = 4.8)

  p
}

plot_pairwise_fst_manhattan <- function(fst_windows) {
  axis_df <- fst_windows |>
    dplyr::distinct(.data$chromosome, .data$chrom_start, .data$chrom_length) |>
    dplyr::mutate(
      axis_midpoint = .data$chrom_start + (.data$chrom_length / 2),
      axis_label = stringr::str_remove(as.character(.data$chromosome), "^SUPER_")
    )

  p <- ggplot(
    fst_windows,
    aes(
      x = .data$genome_midpoint,
      y = .data$avg_wc_fst,
      color = dplyr::case_when(
        .data$is_pair_outlier ~ "Outlier",
        TRUE ~ .data$chromosome_landscape_group
      )
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_hline(
      aes(yintercept = .data$pair_outlier_threshold),
      linewidth = 0.25,
      linetype = "dashed",
      color = "#7a3b2e"
    ) +
    geom_point(size = 0.22, alpha = 0.75, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "Macrochromosome odd" = chromosome_landscape_colors[["Macrochromosome odd"]],
        "Macrochromosome even" = chromosome_landscape_colors[["Macrochromosome even"]],
        "Microchromosome odd" = chromosome_landscape_colors[["Microchromosome odd"]],
        "Microchromosome even" = chromosome_landscape_colors[["Microchromosome even"]],
        "Outlier" = chromosome_landscape_colors[["Outlier"]]
      ),
      breaks = c("Macrochromosome odd", "Macrochromosome even", "Microchromosome odd", "Microchromosome even", "Outlier"),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = axis_df$axis_midpoint,
      labels = axis_df$axis_label,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    facet_wrap(~pair, ncol = 2) +
    labs(
      x = "Chromosome",
      y = "Window FST",
      title = "Pairwise regional FST outlier scan",
      subtitle = paste0(
        "Outliers: pair-specific top ",
        signif((1 - fst_outlier_quantile) * 100, 3),
        "% of windows; SNPs >= ",
        min_fst_snps
      )
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 6),
      strip.text = element_text(size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(out_fig_dir, "SnowyOwl_regions_pairwise_fst_manhattan_faceted.png"),
    p,
    width = 12,
    height = 10,
    dpi = 350
  )
  ggsave(
    file.path(out_fig_dir, "SnowyOwl_regions_pairwise_fst_manhattan_faceted.pdf"),
    p,
    width = 12,
    height = 10
  )

  p
}

plot_single_pair_fst_manhattan <- function(fst_windows) {
  axis_df <- fst_windows |>
    dplyr::distinct(.data$chromosome, .data$chrom_start, .data$chrom_length) |>
    dplyr::mutate(
      axis_midpoint = .data$chrom_start + (.data$chrom_length / 2),
      axis_label = stringr::str_remove(as.character(.data$chromosome), "^SUPER_")
    )

  pair_names <- sort(unique(fst_windows$pair))

  for (pair_name in pair_names) {
    pair_df <- fst_windows |>
      dplyr::filter(.data$pair == pair_name)

    pair_threshold <- unique(pair_df$pair_outlier_threshold)

    p <- ggplot(
      pair_df,
      aes(
        x = .data$genome_midpoint,
        y = .data$avg_wc_fst,
        color = dplyr::case_when(
          .data$is_pair_outlier ~ "Outlier",
          TRUE ~ .data$chromosome_landscape_group
        )
      )
    ) +
      geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
      geom_hline(
        yintercept = pair_threshold[1],
        linewidth = 0.35,
        linetype = "dashed",
        color = "#7a3b2e"
      ) +
      geom_point(size = 0.45, alpha = 0.8, na.rm = TRUE) +
      scale_color_manual(
        values = c(
          "Macrochromosome odd" = chromosome_landscape_colors[["Macrochromosome odd"]],
          "Macrochromosome even" = chromosome_landscape_colors[["Macrochromosome even"]],
          "Microchromosome odd" = chromosome_landscape_colors[["Microchromosome odd"]],
          "Microchromosome even" = chromosome_landscape_colors[["Microchromosome even"]],
          "Outlier" = chromosome_landscape_colors[["Outlier"]]
        ),
        breaks = c("Macrochromosome odd", "Macrochromosome even", "Microchromosome odd", "Microchromosome even", "Outlier"),
        name = NULL
      ) +
      scale_x_continuous(
        breaks = axis_df$axis_midpoint,
        labels = axis_df$axis_label,
        expand = expansion(mult = c(0.005, 0.005))
      ) +
      labs(
        x = "Chromosome",
        y = "Window FST",
        title = paste("Regional FST:", pair_name),
        subtitle = paste0(
          "Outlier threshold = ",
          signif(pair_threshold[1], 4),
          "; SNPs >= ",
          min_fst_snps
        )
      ) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 7),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )

    safe_pair_name <- stringr::str_replace_all(pair_name, "[^A-Za-z0-9]+", "_")

    ggsave(
      file.path(out_fig_dir, paste0("SnowyOwl_regions_fst_manhattan_", safe_pair_name, ".png")),
      p,
      width = 11,
      height = 4.8,
      dpi = 350
    )
    ggsave(
      file.path(out_fig_dir, paste0("SnowyOwl_regions_fst_manhattan_", safe_pair_name, ".pdf")),
      p,
      width = 11,
      height = 4.8
    )
  }
}

fst <- read_pixy_tsv(fst_file) |>
  add_genome_position()

dxy <- read_pixy_tsv(dxy_file)

fst_filtered <- fst |>
  dplyr::filter(!is.na(.data$avg_wc_fst), .data$no_snps >= min_fst_snps)

pairwise_fst_summary <- fst_filtered |>
  dplyr::group_by(.data$pop1, .data$pop2, .data$pair) |>
  dplyr::summarise(
    n_windows = dplyr::n(),
    n_positive_windows = sum(.data$avg_wc_fst > 0, na.rm = TRUE),
    mean_fst_unweighted = mean(.data$avg_wc_fst, na.rm = TRUE),
    mean_fst_weighted_by_snps = weighted.mean(.data$avg_wc_fst, .data$no_snps, na.rm = TRUE),
    median_fst = median(.data$avg_wc_fst, na.rm = TRUE),
    q95_fst = as.numeric(stats::quantile(.data$avg_wc_fst, 0.95, na.rm = TRUE)),
    q99_fst = as.numeric(stats::quantile(.data$avg_wc_fst, 0.99, na.rm = TRUE)),
    q999_fst = as.numeric(stats::quantile(.data$avg_wc_fst, 0.999, na.rm = TRUE)),
    max_fst = max(.data$avg_wc_fst, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(.data$mean_fst_weighted_by_snps))

global_fst_summary <- fst_filtered |>
  dplyr::summarise(
    n_pairs = dplyr::n_distinct(.data$pair),
    n_pair_windows = dplyr::n(),
    n_genomic_windows = dplyr::n_distinct(paste(.data$chromosome, .data$window_pos_1, .data$window_pos_2)),
    min_fst_snps = min_fst_snps,
    mean_fst_unweighted = mean(.data$avg_wc_fst, na.rm = TRUE),
    mean_fst_weighted_by_snps = weighted.mean(.data$avg_wc_fst, .data$no_snps, na.rm = TRUE),
    median_fst = median(.data$avg_wc_fst, na.rm = TRUE),
    min_fst = min(.data$avg_wc_fst, na.rm = TRUE),
    max_fst = max(.data$avg_wc_fst, na.rm = TRUE),
    q95_fst = as.numeric(stats::quantile(.data$avg_wc_fst, 0.95, na.rm = TRUE)),
    q99_fst = as.numeric(stats::quantile(.data$avg_wc_fst, 0.99, na.rm = TRUE)),
    q999_fst = as.numeric(stats::quantile(.data$avg_wc_fst, 0.999, na.rm = TRUE))
  )

fst_with_thresholds <- fst_filtered |>
  dplyr::group_by(.data$pair) |>
  dplyr::mutate(
    pair_outlier_threshold = as.numeric(
      stats::quantile(.data$avg_wc_fst, fst_outlier_quantile, na.rm = TRUE)
    ),
    is_pair_outlier = .data$avg_wc_fst >= .data$pair_outlier_threshold
  ) |>
  dplyr::ungroup()

top_pairwise_fst_windows <- fst_with_thresholds |>
  dplyr::left_join(
    dxy |>
      dplyr::select(
        pop1,
        pop2,
        chromosome,
        window_pos_1,
        window_pos_2,
        avg_dxy,
        dxy_no_sites = no_sites
      ),
    by = c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2")
  ) |>
  dplyr::filter(.data$avg_wc_fst > 0) |>
  dplyr::arrange(dplyr::desc(.data$avg_wc_fst)) |>
  dplyr::select(
    pop1,
    pop2,
    pair,
    chromosome,
    window_pos_1,
    window_pos_2,
    avg_wc_fst,
    no_snps,
    avg_dxy,
    dxy_no_sites,
    pair_outlier_threshold,
    is_pair_outlier
  ) |>
  utils::head(200)

recurrent_fst_outlier_windows <- fst_with_thresholds |>
  dplyr::filter(.data$is_pair_outlier, .data$avg_wc_fst > 0) |>
  dplyr::group_by(.data$chromosome, .data$window_pos_1, .data$window_pos_2) |>
  dplyr::summarise(
    n_outlier_pairs = dplyr::n_distinct(.data$pair),
    outlier_pairs = paste(sort(unique(.data$pair)), collapse = ";"),
    max_fst = max(.data$avg_wc_fst, na.rm = TRUE),
    mean_fst_across_outlier_pairs = mean(.data$avg_wc_fst, na.rm = TRUE),
    min_no_snps = min(.data$no_snps, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(.data$n_outlier_pairs), dplyr::desc(.data$max_fst))

max_fst_windows <- fst_filtered |>
  dplyr::group_by(
    .data$chromosome,
    .data$chrom_number,
    .data$window_pos_1,
    .data$window_pos_2,
    .data$window_midpoint,
    .data$chrom_start,
    .data$chrom_length,
    .data$genome_midpoint,
    .data$chromosome_class,
    .data$chromosome_landscape_group
  ) |>
  dplyr::summarise(
    max_fst = max(.data$avg_wc_fst, na.rm = TRUE),
    pair_at_max_fst = .data$pair[which.max(.data$avg_wc_fst)][1],
    no_snps_at_max_fst = .data$no_snps[which.max(.data$avg_wc_fst)][1],
    mean_fst_across_pairs = mean(.data$avg_wc_fst, na.rm = TRUE),
    n_pairs = dplyr::n_distinct(.data$pair),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$chrom_number, .data$window_pos_1)

readr::write_csv(global_fst_summary, file.path(out_tab_dir, "SnowyOwl_regions_fst_global_summary.csv"))
readr::write_csv(pairwise_fst_summary, file.path(out_tab_dir, "SnowyOwl_regions_pairwise_fst_summary.csv"))
readr::write_csv(top_pairwise_fst_windows, file.path(out_tab_dir, "SnowyOwl_regions_top200_pairwise_fst_windows.csv"))
readr::write_csv(recurrent_fst_outlier_windows, file.path(out_tab_dir, "SnowyOwl_regions_recurrent_fst_outlier_windows.csv"))
readr::write_csv(max_fst_windows, file.path(out_tab_dir, "SnowyOwl_regions_max_fst_by_window.csv"))

fst_heatmap <- plot_pairwise_fst_heatmap(pairwise_fst_summary)
max_fst_plot <- plot_genomewide_max_fst(max_fst_windows)
pairwise_fst_manhattan <- plot_pairwise_fst_manhattan(fst_with_thresholds)
plot_single_pair_fst_manhattan(fst_with_thresholds)

message("Wrote regional FST summaries to: ", out_tab_dir)
message("Wrote regional FST figures to: ", out_fig_dir)
