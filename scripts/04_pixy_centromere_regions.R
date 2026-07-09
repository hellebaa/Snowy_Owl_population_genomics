# Compare pixy diversity summaries in putative centromeric and pericentromeric regions.

source("scripts/00_setup.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

centromere_bed <- here(paths$centromere_bed %||% "data/metadata/centromeres/putative_centromeres.bed")
pericentromere_flank_bp <- paths$pericentromere_flank_bp %||% 1000000
min_centromere_region_sites <- paths$min_centromere_region_sites %||% 1000

pi_file <- here(paths$pixy_all_pi %||% file.path(paths$pixy_dir, "SnowyOwl_ALL_pi.txt"))
tajima_file <- here(paths$pixy_all_tajima_d %||% file.path(paths$pixy_dir, "SnowyOwl_ALL_tajima_d.txt"))

out_fig_dir <- here(paths$results_figures, "pixy")
out_tab_dir <- here(paths$results_tables, "pixy")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)

read_centromeres <- function(file) {
  if (!file.exists(file)) {
    stop("Missing centromere BED file: ", file, call. = FALSE)
  }

  first_data_line <- readr::read_lines(file, n_max = 20) |>
    stringr::str_subset("^[^#]") |>
    utils::head(1)

  has_header <- length(first_data_line) == 1 &&
    stringr::str_detect(first_data_line, regex("^(Chr|chromosome)\\t", ignore_case = TRUE))

  bed <- readr::read_tsv(
    file,
    col_names = has_header,
    comment = "#",
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )

  if (nrow(bed) == 0) {
    stop(
      "Centromere BED is empty. Add intervals to: ",
      file,
      call. = FALSE
    )
  }

  if (has_header) {
    bed <- bed |>
      dplyr::rename(
        chromosome = dplyr::any_of(c("Chr", "chr", "chrom", "chromosome")),
        start_1based = dplyr::any_of(c("Start", "start", "start_1based")),
        end = dplyr::any_of(c("End", "end")),
        chromosome_size = dplyr::any_of(c("Chr_size", "chromosome_size", "chrom_size")),
        morphology = dplyr::any_of(c("Morphology_new", "Morphology", "morphology"))
      ) |>
      dplyr::mutate(start0 = as.numeric(.data$start_1based) - 1)
  } else {
    if (ncol(bed) < 3) {
      stop("Centromere BED must contain at least chromosome, start0, and end columns.", call. = FALSE)
    }

    names(bed)[seq_len(min(ncol(bed), 7))] <- c(
      "chromosome",
      "start0",
      "end",
      "centromere_id",
      "morphology",
      "evidence",
      "notes"
    )[seq_len(min(ncol(bed), 7))]
  }

  for (optional_col in c("centromere_id", "morphology", "evidence", "notes")) {
    if (!optional_col %in% names(bed)) {
      bed[[optional_col]] <- NA_character_
    }
  }

  bed |>
    dplyr::mutate(
      start0 = as.numeric(.data$start0),
      end = as.numeric(.data$end),
      centromere_id = dplyr::if_else(
        is.na(.data$centromere_id) | .data$centromere_id == "",
        paste0(.data$chromosome, "_centromere"),
        .data$centromere_id
      )
    ) |>
    dplyr::filter(!is.na(.data$start0), !is.na(.data$end), .data$end > .data$start0)
}

read_pixy_windows <- function(file, value_col) {
  if (!file.exists(file)) {
    stop("Missing expected pixy file: ", file, call. = FALSE)
  }

  readr::read_tsv(file, na = c("NA", "nan", "NaN", "."), show_col_types = FALSE) |>
    dplyr::filter(stringr::str_detect(.data$chromosome, "^SUPER")) |>
    dplyr::mutate(
      chromosome = as.character(.data$chromosome),
      window_start0 = .data$window_pos_1 - 1,
      window_end = .data$window_pos_2,
      metric = value_col,
      value = .data[[value_col]]
    )
}

overlaps_any <- function(query, subject) {
  vapply(
    seq_len(nrow(query)),
    function(i) {
      any(
        subject$chromosome == query$chromosome[i] &
          subject$start0 < query$window_end[i] &
          subject$end > query$window_start0[i]
      )
    },
    logical(1)
  )
}

annotate_centromere_regions <- function(windows, centromeres, flank_bp) {
  pericentromeres <- centromeres |>
    dplyr::mutate(
      start0 = pmax(0, .data$start0 - flank_bp),
      end = .data$end + flank_bp
    )

  in_centromere <- overlaps_any(windows, centromeres)
  in_pericentromere_extended <- overlaps_any(windows, pericentromeres)

  windows |>
    dplyr::mutate(
      centromere_region = dplyr::case_when(
        in_centromere ~ "Centromere",
        in_pericentromere_extended ~ "Pericentromere",
        TRUE ~ "Other"
      ),
      centromere_region = factor(
        .data$centromere_region,
        levels = c("Centromere", "Pericentromere", "Other")
      )
    )
}

summarise_regions <- function(windows) {
  windows |>
    dplyr::group_by(.data$metric, .data$centromere_region) |>
    dplyr::summarise(
      n_windows = dplyr::n(),
      n_windows_with_value = sum(!is.na(.data$value)),
      callable_sites = sum(.data$no_sites, na.rm = TRUE),
      mean_unweighted = mean(.data$value, na.rm = TRUE),
      mean_weighted_by_sites = weighted.mean(.data$value, .data$no_sites, na.rm = TRUE),
      median = median(.data$value, na.rm = TRUE),
      min = min(.data$value, na.rm = TRUE),
      max = max(.data$value, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_regions <- function(windows, metric_name, y_label, output_stub) {
  plot_df <- windows |>
    dplyr::filter(.data$metric == metric_name, !is.na(.data$value), !is.na(.data$centromere_region))

  kruskal_p <- stats::kruskal.test(value ~ centromere_region, data = plot_df)$p.value

  test_summary <- tibble::tibble(
    metric = metric_name,
    test = "Kruskal-Wallis rank-sum test",
    grouping = "Centromere vs Pericentromere vs Other",
    p_value = kruskal_p
  )

  readr::write_csv(
    test_summary,
    file.path(out_tab_dir, paste0(output_stub, "_kruskal_test.csv"))
  )

  p <- plot_df |>
    ggplot(aes(x = .data$centromere_region, y = .data$value, fill = .data$centromere_region)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, na.rm = TRUE) +
    geom_jitter(width = 0.18, size = 0.25, alpha = 0.16, na.rm = TRUE) +
    scale_fill_manual(
      values = c("Centromere" = "#7a3b2e", "Pericentromere" = "#b35c44", "Other" = "#2f4f4f"),
      guide = "none"
    ) +
    labs(
      x = NULL,
      y = y_label,
      title = paste(y_label, "by putative centromere region"),
      subtitle = paste0("Kruskal-Wallis p = ", format_p_value(kruskal_p))
    ) +
    theme(axis.text.x = element_text(size = 9))

  ggsave(file.path(out_fig_dir, paste0(output_stub, ".png")), p, width = 5.5, height = 4, dpi = 350)
  ggsave(file.path(out_fig_dir, paste0(output_stub, ".pdf")), p, width = 5.5, height = 4)

  p
}

centromeres <- read_centromeres(centromere_bed)

pi_windows <- read_pixy_windows(pi_file, "avg_pi")
tajima_windows <- read_pixy_windows(tajima_file, "tajima_d")

centromere_windows <- dplyr::bind_rows(pi_windows, tajima_windows) |>
  annotate_centromere_regions(centromeres, pericentromere_flank_bp)

region_summary <- summarise_regions(centromere_windows)

region_summary_min_sites <- centromere_windows |>
  dplyr::filter(.data$no_sites >= min_centromere_region_sites) |>
  summarise_regions() |>
  dplyr::mutate(min_no_sites = min_centromere_region_sites, .before = "n_windows")

region_callability_summary <- centromere_windows |>
  dplyr::group_by(.data$metric, .data$centromere_region) |>
  dplyr::summarise(
    n_windows = dplyr::n(),
    n_windows_no_callable_sites = sum(.data$no_sites == 0, na.rm = TRUE),
    n_windows_lt_min_sites = sum(.data$no_sites < min_centromere_region_sites, na.rm = TRUE),
    median_no_sites = median(.data$no_sites, na.rm = TRUE),
    mean_no_sites = mean(.data$no_sites, na.rm = TRUE),
    min_no_sites = min(.data$no_sites, na.rm = TRUE),
    max_no_sites = max(.data$no_sites, na.rm = TRUE),
    .groups = "drop"
  )

chromosome_region_summary <- centromere_windows |>
  dplyr::group_by(.data$metric, .data$chromosome, .data$centromere_region) |>
  dplyr::summarise(
    n_windows = dplyr::n(),
    callable_sites = sum(.data$no_sites, na.rm = TRUE),
    mean_unweighted = mean(.data$value, na.rm = TRUE),
    mean_weighted_by_sites = weighted.mean(.data$value, .data$no_sites, na.rm = TRUE),
    median = median(.data$value, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  region_summary,
  file.path(out_tab_dir, "SnowyOwl_ALL_pixy_centromere_region_summary.csv")
)

readr::write_csv(
  region_summary_min_sites,
  file.path(out_tab_dir, "SnowyOwl_ALL_pixy_centromere_region_summary_min_sites.csv")
)

readr::write_csv(
  region_callability_summary,
  file.path(out_tab_dir, "SnowyOwl_ALL_pixy_centromere_region_callability_summary.csv")
)

readr::write_csv(
  chromosome_region_summary,
  file.path(out_tab_dir, "SnowyOwl_ALL_pixy_centromere_region_by_chromosome_summary.csv")
)

pi_region_plot <- plot_regions(
  centromere_windows,
  metric_name = "avg_pi",
  y_label = expression(pi),
  output_stub = "SnowyOwl_ALL_pi_centromere_regions"
)

tajima_region_plot <- plot_regions(
  centromere_windows,
  metric_name = "tajima_d",
  y_label = "Tajima's D",
  output_stub = "SnowyOwl_ALL_tajima_d_centromere_regions"
)

message("Wrote centromere-region pixy summaries to: ", out_tab_dir)
message("Wrote centromere-region pixy figures to: ", out_fig_dir)
