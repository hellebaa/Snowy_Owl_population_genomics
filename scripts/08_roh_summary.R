# Summarize and plot runs of homozygosity for snowy owl individuals.

source("scripts/00_setup.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

roh_dir <- paths$roh_dir %||% "data/external/roh"
autosomal_genome_size_bp <- paths$autosomal_genome_size_bp %||% 1173697940

hom_file <- here(roh_dir, "SnowyOwl_autosomes_noContam_ROH.hom")
indiv_file <- here(roh_dir, "SnowyOwl_autosomes_noContam_ROH.hom.indiv")
bins_file <- here(roh_dir, "SnowyOwl_autosomes_noContam_ROH.ROH_length_bins.tsv")

out_fig_dir <- here(paths$results_figures, "roh")
out_tab_dir <- here(paths$results_tables, "roh")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)

infer_region <- function(iid) {
  dplyr::case_when(
    stringr::str_detect(iid, "^FNM") ~ "FNM",
    stringr::str_detect(iid, "^GRL") ~ "GRL",
    stringr::str_detect(iid, "^NY") ~ "NYS",
    stringr::str_detect(iid, "^SKW") ~ "SKW",
    stringr::str_detect(iid, "^WRG") ~ "WRG",
    TRUE ~ "Other"
  )
}

read_plink_table <- function(file) {
  if (!file.exists(file)) {
    stop("Missing expected ROH file: ", file, call. = FALSE)
  }

  read.table(file, header = TRUE, stringsAsFactors = FALSE)
}

roh_segments <- read_plink_table(hom_file) |>
  dplyr::mutate(
    region = infer_region(.data$IID),
    region_label = region_display(.data$region),
    roh_length_mb = .data$KB / 1000,
    length_class = dplyr::case_when(
      .data$roh_length_mb < 1 ~ "0.3-1 Mb",
      .data$roh_length_mb <= 5 ~ "1-5 Mb",
      TRUE ~ ">5 Mb"
    ),
    length_class = factor(.data$length_class, levels = c("0.3-1 Mb", "1-5 Mb", ">5 Mb"))
  )

roh_indiv <- read_plink_table(indiv_file) |>
  dplyr::mutate(
    region = infer_region(.data$IID),
    region_label = region_display(.data$region),
    roh_total_mb = .data$KB / 1000,
    froh = (.data$KB * 1000) / autosomal_genome_size_bp
  )

roh_bins <- readr::read_tsv(bins_file, show_col_types = FALSE) |>
  dplyr::mutate(
    region = infer_region(.data$IID),
    region_label = region_display(.data$region),
    froh_total = (.data$ROH_total_kb * 1000) / autosomal_genome_size_bp,
    froh_0.3_1Mb = (.data$`ROH_0.3_1Mb_kb` * 1000) / autosomal_genome_size_bp,
    froh_1_5Mb = (.data$`ROH_1_5Mb_kb` * 1000) / autosomal_genome_size_bp,
    froh_gt5Mb = (.data$ROH_gt5Mb_kb * 1000) / autosomal_genome_size_bp
  )

individual_summary <- roh_indiv |>
  dplyr::left_join(
    roh_bins |>
      dplyr::select(
        IID,
        ROH_count,
        ROH_total_kb,
        `ROH_0.3_1Mb_kb`,
        ROH_1_5Mb_kb,
        ROH_gt5Mb_kb,
        froh_total,
        froh_0.3_1Mb,
        froh_1_5Mb,
        froh_gt5Mb
      ),
    by = "IID"
  )

region_summary_by_region <- individual_summary |>
  dplyr::group_by(.data$region) |>
  dplyr::summarise(
    n_individuals = dplyr::n(),
    region_label = dplyr::first(.data$region_label),
    mean_froh = mean(.data$froh, na.rm = TRUE),
    median_froh = median(.data$froh, na.rm = TRUE),
    max_froh = max(.data$froh, na.rm = TRUE),
    mean_roh_count = mean(.data$NSEG, na.rm = TRUE),
    median_roh_count = median(.data$NSEG, na.rm = TRUE),
    mean_roh_total_mb = mean(.data$roh_total_mb, na.rm = TRUE),
    median_roh_total_mb = median(.data$roh_total_mb, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$region)

all_summary <- individual_summary |>
  dplyr::summarise(
    region = "ALL",
    n_individuals = dplyr::n(),
    region_label = "All",
    mean_froh = mean(.data$froh, na.rm = TRUE),
    median_froh = median(.data$froh, na.rm = TRUE),
    max_froh = max(.data$froh, na.rm = TRUE),
    mean_roh_count = mean(.data$NSEG, na.rm = TRUE),
    median_roh_count = median(.data$NSEG, na.rm = TRUE),
    mean_roh_total_mb = mean(.data$roh_total_mb, na.rm = TRUE),
    median_roh_total_mb = median(.data$roh_total_mb, na.rm = TRUE),
    .groups = "drop"
  )

region_summary <- dplyr::bind_rows(all_summary, region_summary_by_region)

segment_summary <- roh_segments |>
  dplyr::summarise(
    n_segments = dplyr::n(),
    mean_segment_mb = mean(.data$roh_length_mb, na.rm = TRUE),
    median_segment_mb = median(.data$roh_length_mb, na.rm = TRUE),
    max_segment_mb = max(.data$roh_length_mb, na.rm = TRUE),
    individual_with_max_segment = .data$IID[which.max(.data$roh_length_mb)][1],
    chromosome_with_max_segment = .data$CHR[which.max(.data$roh_length_mb)][1],
    .groups = "drop"
  )

length_class_summary <- roh_segments |>
  dplyr::count(.data$length_class, name = "n_segments") |>
  dplyr::mutate(fraction_segments = .data$n_segments / sum(.data$n_segments))

readr::write_csv(individual_summary, file.path(out_tab_dir, "SnowyOwl_ROH_individual_summary.csv"))
readr::write_csv(region_summary, file.path(out_tab_dir, "SnowyOwl_ROH_region_summary.csv"))
readr::write_csv(segment_summary, file.path(out_tab_dir, "SnowyOwl_ROH_segment_summary.csv"))
readr::write_csv(length_class_summary, file.path(out_tab_dir, "SnowyOwl_ROH_length_class_summary.csv"))

region_order <- region_summary_by_region |>
  dplyr::arrange(.data$median_froh) |>
  dplyr::pull(.data$region_label)

froh_plot_data <- individual_summary |>
  dplyr::bind_rows(
    individual_summary |>
      dplyr::mutate(region = "ALL", region_label = "All")
  ) |>
  dplyr::mutate(region_label = factor(.data$region_label, levels = c("All", region_order)))

froh_plot_colors <- c("All" = "ivory3", pop_colors)

froh_plot <- froh_plot_data |>
  ggplot(aes(x = .data$region_label, y = .data$froh, fill = .data$region_label)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, na.rm = TRUE) +
  geom_jitter(width = 0.15, size = 1.7, alpha = 0.8, na.rm = TRUE) +
  scale_fill_manual(values = froh_plot_colors, guide = "none") +
  labs(x = "Sampling region", y = expression(F[ROH]), title = "Individual runs of homozygosity") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

roh_length_plot <- roh_segments |>
  ggplot(aes(x = .data$roh_length_mb, fill = .data$length_class)) +
  geom_histogram(binwidth = 0.25, color = "white", linewidth = 0.2, boundary = 0) +
  scale_fill_manual(values = roh_length_colors, name = NULL) +
  labs(x = "ROH segment length (Mb)", y = "Number of segments", title = "ROH length distribution") +
  theme(legend.position = "bottom")

roh_bins_long <- roh_bins |>
  dplyr::select("IID", "region", "region_label", "froh_0.3_1Mb", "froh_1_5Mb", "froh_gt5Mb") |>
  tidyr::pivot_longer(
    cols = starts_with("froh_"),
    names_to = "length_class",
    values_to = "froh"
  ) |>
  dplyr::mutate(
    length_class = dplyr::recode(
      .data$length_class,
      froh_0.3_1Mb = "0.3-1 Mb",
      froh_1_5Mb = "1-5 Mb",
      froh_gt5Mb = ">5 Mb"
    ),
    length_class = factor(.data$length_class, levels = c("0.3-1 Mb", "1-5 Mb", ">5 Mb"))
  )

froh_stacked_plot <- roh_bins_long |>
  dplyr::mutate(
    region_label = factor(.data$region_label, levels = names(pop_colors)),
    IID = factor(.data$IID, levels = roh_bins |> dplyr::arrange(.data$region_label, .data$froh_total) |> dplyr::pull(.data$IID))
  ) |>
  ggplot(aes(x = .data$IID, y = .data$froh, fill = .data$length_class)) +
  geom_col(width = 0.85) +
  facet_grid(~region_label, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = roh_length_colors, name = NULL) +
  labs(x = NULL, y = expression(F[ROH]), title = "ROH burden by length class") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = grid::unit(0.12, "lines")
  )

ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_FROH_by_region.png"), froh_plot, width = 5.5, height = 4.2, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_FROH_by_region.pdf"), froh_plot, width = 5.5, height = 4.2)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_length_distribution.png"), roh_length_plot, width = 5.5, height = 4.2, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_length_distribution.pdf"), roh_length_plot, width = 5.5, height = 4.2)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_FROH_length_classes_by_individual.png"), froh_stacked_plot, width = 9, height = 4.2, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_FROH_length_classes_by_individual.pdf"), froh_stacked_plot, width = 9, height = 4.2)

message("Wrote ROH summaries to: ", out_tab_dir)
message("Wrote ROH figures to: ", out_fig_dir)
