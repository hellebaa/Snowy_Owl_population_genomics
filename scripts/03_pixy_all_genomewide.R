# Plot genome-wide nucleotide diversity and Tajima's D for all snowy owl samples.

source("scripts/00_setup.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

pixy_dir <- paths$pixy_dir %||% paths$data_external
pi_file <- here(pixy_dir, "SnowyOwl_ALL_pi.txt")
tajima_file <- here(pixy_dir, "SnowyOwl_ALL_tajima_d.txt")

out_fig_dir <- here(paths$results_figures, "pixy")
out_tab_dir <- here(paths$results_tables, "pixy")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)

read_pixy_table <- function(file) {
  if (!file.exists(file)) {
    stop("Missing expected pixy file: ", file, call. = FALSE)
  }

  readr::read_tsv(
    file,
    na = c("NA", "nan", "NaN", "."),
    show_col_types = FALSE
  ) |>
    dplyr::filter(stringr::str_detect(.data$chromosome, "^SUPER")) |>
    dplyr::mutate(
      chromosome = as.character(.data$chromosome),
      chrom_number = suppressWarnings(
        as.integer(stringr::str_extract(.data$chromosome, "(?<=SUPER_)\\d+"))
      ),
      chromosome_class = chromosome_class(.data$chrom_number),
      chromosome_landscape_group = chromosome_landscape_group(.data$chrom_number),
      window_midpoint = (.data$window_pos_1 + .data$window_pos_2) / 2
    )
}

chromosome_order <- function(chromosomes) {
  chrom_df <- tibble::tibble(chromosome = unique(chromosomes)) |>
    dplyr::mutate(
      chrom_number = suppressWarnings(
        as.integer(stringr::str_extract(.data$chromosome, "(?<=SUPER_)\\d+"))
      )
    ) |>
    dplyr::arrange(is.na(.data$chrom_number), .data$chrom_number, .data$chromosome)

  chrom_df$chromosome
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

plot_genomewide_windows <- function(df, value_col, y_label, title, output_stub) {
  axis_df <- df |>
    dplyr::distinct(.data$chromosome, .data$chrom_start, .data$chrom_length) |>
    dplyr::mutate(
      axis_midpoint = .data$chrom_start + (.data$chrom_length / 2),
      chrom_number = suppressWarnings(
        as.integer(stringr::str_extract(as.character(.data$chromosome), "(?<=SUPER_)\\d+"))
      ),
      axis_label = as.character(.data$chrom_number)
    )

  p <- ggplot(
    df,
    aes(
      x = .data$genome_midpoint,
      y = .data[[value_col]],
      color = .data$chromosome_landscape_group
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey70") +
    geom_point(size = 0.35, alpha = 0.75, na.rm = TRUE) +
    scale_color_manual(values = chromosome_landscape_colors, guide = "none") +
    scale_x_continuous(
      breaks = axis_df$axis_midpoint,
      labels = axis_df$axis_label,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    labs(x = "Chromosome", y = y_label, title = title) +
    theme(
      axis.text.x = element_text(size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(out_fig_dir, paste0(output_stub, ".png")),
    plot = p,
    width = 11,
    height = 4.8,
    dpi = 350
  )

  ggsave(
    filename = file.path(out_fig_dir, paste0(output_stub, ".pdf")),
    plot = p,
    width = 11,
    height = 4.8
  )

  p
}

plot_chromosome_class <- function(df, value_col, y_label, title, output_stub) {
  plot_df <- df |>
    dplyr::filter(!is.na(.data[[value_col]]), !is.na(.data$chromosome_class))

  class_summary <- plot_df |>
    dplyr::group_by(.data$chromosome_class) |>
    dplyr::summarise(
      n_windows = dplyr::n(),
      callable_sites = sum(.data$no_sites, na.rm = TRUE),
      mean_unweighted = mean(.data[[value_col]], na.rm = TRUE),
      mean_weighted_by_sites = weighted.mean(.data[[value_col]], .data$no_sites, na.rm = TRUE),
      median = median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )

  wilcox_p <- stats::wilcox.test(
    stats::as.formula(paste(value_col, "~ chromosome_class")),
    data = plot_df,
    exact = FALSE
  )$p.value

  test_summary <- tibble::tibble(
    metric = value_col,
    test = "Wilcoxon rank-sum test",
    grouping = "Macrochromosome vs Microchromosome",
    p_value = wilcox_p
  )

  readr::write_csv(
    class_summary,
    file.path(out_tab_dir, paste0(output_stub, "_summary.csv"))
  )
  readr::write_csv(
    test_summary,
    file.path(out_tab_dir, paste0(output_stub, "_wilcoxon_test.csv"))
  )

  mean_caption <- class_summary |>
    dplyr::arrange(.data$chromosome_class) |>
    dplyr::mutate(label = paste0(.data$chromosome_class, " = ", signif(.data$mean_weighted_by_sites, 4))) |>
    dplyr::pull(.data$label) |>
    paste(collapse = "; ")

  mean_caption <- paste("Weighted means:", mean_caption) |>
    stringr::str_wrap(width = 70)

  p <- ggplot(
    plot_df,
    aes(x = .data$chromosome_class, y = .data[[value_col]], fill = .data$chromosome_class)
  ) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, na.rm = TRUE) +
    geom_jitter(width = 0.18, size = 0.25, alpha = 0.18, na.rm = TRUE) +
    scale_fill_manual(
      values = c("Macrochromosome" = "lightgoldenrod3", "Microchromosome" = "lightcyan3"),
      guide = "none"
    ) +
    labs(
      x = NULL,
      y = y_label,
      title = title,
      subtitle = paste0("Wilcoxon p = ", format_p_value(wilcox_p)),
      caption = mean_caption
    ) +
    theme(
      axis.text.x = element_text(size = 10),
      plot.caption = element_text(hjust = 0),
      panel.grid.minor = element_blank()
    )

  ggsave(
    filename = file.path(out_fig_dir, paste0(output_stub, ".png")),
    plot = p,
    width = 5.2,
    height = 4,
    dpi = 350
  )

  ggsave(
    filename = file.path(out_fig_dir, paste0(output_stub, ".pdf")),
    plot = p,
    width = 5.2,
    height = 4
  )

  p
}

pi_windows <- read_pixy_table(pi_file) |>
  add_genome_position()

tajima_windows <- read_pixy_table(tajima_file) |>
  add_genome_position()

pi_summary <- pi_windows |>
  dplyr::summarise(
    population = dplyr::first(.data$pop),
    n_windows = dplyr::n(),
    n_windows_with_pi = sum(!is.na(.data$avg_pi)),
    callable_sites = sum(.data$no_sites, na.rm = TRUE),
    mean_pi_unweighted = mean(.data$avg_pi, na.rm = TRUE),
    mean_pi_weighted_by_sites = weighted.mean(.data$avg_pi, .data$no_sites, na.rm = TRUE),
    median_pi = median(.data$avg_pi, na.rm = TRUE),
    min_pi = min(.data$avg_pi, na.rm = TRUE),
    max_pi = max(.data$avg_pi, na.rm = TRUE)
  )

tajima_summary <- tajima_windows |>
  dplyr::summarise(
    population = dplyr::first(.data$pop),
    n_windows = dplyr::n(),
    n_windows_with_tajima_d = sum(!is.na(.data$tajima_d)),
    callable_sites = sum(.data$no_sites, na.rm = TRUE),
    mean_tajima_d_unweighted = mean(.data$tajima_d, na.rm = TRUE),
    mean_tajima_d_weighted_by_sites = weighted.mean(.data$tajima_d, .data$no_sites, na.rm = TRUE),
    median_tajima_d = median(.data$tajima_d, na.rm = TRUE),
    min_tajima_d = min(.data$tajima_d, na.rm = TRUE),
    max_tajima_d = max(.data$tajima_d, na.rm = TRUE)
  )

pi_chromosome_summary <- pi_windows |>
  dplyr::group_by(.data$chromosome, .data$chrom_number, .data$chromosome_class) |>
  dplyr::summarise(
    n_pi_windows = dplyr::n(),
    pi_callable_sites = sum(.data$no_sites, na.rm = TRUE),
    mean_pi = mean(.data$avg_pi, na.rm = TRUE),
    weighted_mean_pi = weighted.mean(.data$avg_pi, .data$no_sites, na.rm = TRUE),
    median_pi = median(.data$avg_pi, na.rm = TRUE),
    .groups = "drop"
  )

tajima_chromosome_summary <- tajima_windows |>
  dplyr::group_by(.data$chromosome, .data$chrom_number, .data$chromosome_class) |>
  dplyr::summarise(
    n_tajima_windows = sum(!is.na(.data$tajima_d)),
    tajima_callable_sites = sum(.data$no_sites, na.rm = TRUE),
    mean_tajima_d = mean(.data$tajima_d, na.rm = TRUE),
    weighted_mean_tajima_d = weighted.mean(.data$tajima_d, .data$no_sites, na.rm = TRUE),
    median_tajima_d = median(.data$tajima_d, na.rm = TRUE),
    .groups = "drop"
  )

chromosome_summary <- pi_chromosome_summary |>
  dplyr::full_join(
    tajima_chromosome_summary,
    by = c("chromosome", "chrom_number", "chromosome_class")
  ) |>
  dplyr::arrange(.data$chrom_number)

pi_chromosome_class_window_summary <- pi_windows |>
  dplyr::group_by(.data$chromosome_class) |>
  dplyr::summarise(
    n_windows = dplyr::n(),
    n_windows_with_pi = sum(!is.na(.data$avg_pi)),
    callable_sites = sum(.data$no_sites, na.rm = TRUE),
    mean_pi_unweighted = mean(.data$avg_pi, na.rm = TRUE),
    mean_pi_weighted_by_sites = weighted.mean(.data$avg_pi, .data$no_sites, na.rm = TRUE),
    median_pi = median(.data$avg_pi, na.rm = TRUE),
    min_pi = min(.data$avg_pi, na.rm = TRUE),
    max_pi = max(.data$avg_pi, na.rm = TRUE),
    .groups = "drop"
  )

tajima_chromosome_class_window_summary <- tajima_windows |>
  dplyr::group_by(.data$chromosome_class) |>
  dplyr::summarise(
    n_windows = dplyr::n(),
    n_windows_with_tajima_d = sum(!is.na(.data$tajima_d)),
    callable_sites = sum(.data$no_sites, na.rm = TRUE),
    mean_tajima_d_unweighted = mean(.data$tajima_d, na.rm = TRUE),
    mean_tajima_d_weighted_by_sites = weighted.mean(.data$tajima_d, .data$no_sites, na.rm = TRUE),
    median_tajima_d = median(.data$tajima_d, na.rm = TRUE),
    min_tajima_d = min(.data$tajima_d, na.rm = TRUE),
    max_tajima_d = max(.data$tajima_d, na.rm = TRUE),
    .groups = "drop"
  )

chromosome_class_chromosome_summary <- chromosome_summary |>
  dplyr::group_by(.data$chromosome_class) |>
  dplyr::summarise(
    n_chromosomes = dplyr::n(),
    mean_of_chromosome_mean_pi = mean(.data$mean_pi, na.rm = TRUE),
    median_of_chromosome_mean_pi = median(.data$mean_pi, na.rm = TRUE),
    mean_of_chromosome_mean_tajima_d = mean(.data$mean_tajima_d, na.rm = TRUE),
    median_of_chromosome_mean_tajima_d = median(.data$mean_tajima_d, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(pi_summary, file.path(out_tab_dir, "SnowyOwl_ALL_pi_summary.csv"))
readr::write_csv(tajima_summary, file.path(out_tab_dir, "SnowyOwl_ALL_tajima_d_summary.csv"))
readr::write_csv(chromosome_summary, file.path(out_tab_dir, "SnowyOwl_ALL_pixy_chromosome_summary.csv"))
readr::write_csv(
  pi_chromosome_class_window_summary,
  file.path(out_tab_dir, "SnowyOwl_ALL_pi_macro_micro_window_summary.csv")
)
readr::write_csv(
  tajima_chromosome_class_window_summary,
  file.path(out_tab_dir, "SnowyOwl_ALL_tajima_d_macro_micro_window_summary.csv")
)
readr::write_csv(
  chromosome_class_chromosome_summary,
  file.path(out_tab_dir, "SnowyOwl_ALL_macro_micro_chromosome_summary.csv")
)

pi_windows |>
  dplyr::filter(!is.na(.data$avg_pi)) |>
  dplyr::arrange(dplyr::desc(.data$avg_pi)) |>
  dplyr::select(
    pop,
    chromosome,
    window_pos_1,
    window_pos_2,
    avg_pi,
    no_sites,
    count_missing
  ) |>
  utils::head(50) |>
  readr::write_csv(file.path(out_tab_dir, "SnowyOwl_ALL_top50_high_pi_windows.csv"))

tajima_windows |>
  dplyr::filter(!is.na(.data$tajima_d)) |>
  dplyr::arrange(dplyr::desc(.data$tajima_d)) |>
  dplyr::select(
    pop,
    chromosome,
    window_pos_1,
    window_pos_2,
    tajima_d,
    no_sites,
    raw_pi,
    raw_watterson_theta
  ) |>
  utils::head(50) |>
  readr::write_csv(file.path(out_tab_dir, "SnowyOwl_ALL_top50_high_tajima_d_windows.csv"))

tajima_windows |>
  dplyr::filter(!is.na(.data$tajima_d)) |>
  dplyr::arrange(.data$tajima_d) |>
  dplyr::select(
    pop,
    chromosome,
    window_pos_1,
    window_pos_2,
    tajima_d,
    no_sites,
    raw_pi,
    raw_watterson_theta
  ) |>
  utils::head(50) |>
  readr::write_csv(file.path(out_tab_dir, "SnowyOwl_ALL_top50_low_tajima_d_windows.csv"))

pi_plot <- plot_genomewide_windows(
  df = pi_windows,
  value_col = "avg_pi",
  y_label = expression(pi),
  title = "Genome-wide nucleotide diversity",
  output_stub = "SnowyOwl_ALL_pi_genomewide"
)

tajima_plot <- plot_genomewide_windows(
  df = tajima_windows,
  value_col = "tajima_d",
  y_label = "Tajima's D",
  title = "Genome-wide Tajima's D",
  output_stub = "SnowyOwl_ALL_tajima_d_genomewide"
)

pi_class_plot <- plot_chromosome_class(
  df = pi_windows,
  value_col = "avg_pi",
  y_label = expression(pi),
  title = "Nucleotide diversity by chromosome class",
  output_stub = "SnowyOwl_ALL_pi_macro_micro"
)

tajima_class_plot <- plot_chromosome_class(
  df = tajima_windows,
  value_col = "tajima_d",
  y_label = "Tajima's D",
  title = "Tajima's D by chromosome class",
  output_stub = "SnowyOwl_ALL_tajima_d_macro_micro"
)

message("Wrote pixy ALL summaries to: ", out_tab_dir)
message("Wrote pixy ALL figures to: ", out_fig_dir)
