# Plot recent effective population size inferred by GONE and compare sensitivity runs.

source("scripts/00_setup.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

gone_dir <- paths$gone_dir %||% "data/external/gone"
generation_time_years <- paths$gone_generation_time_years %||% 8
present_year <- paths$gone_present_year %||% 2021
plot_min_generation <- paths$gone_plot_min_generation %||% 5
plot_max_years_before_present <- paths$gone_plot_max_years_before_present %||% 2000
out_fig_dir <- here(paths$results_figures, "gone")
out_tab_dir <- here(paths$results_tables, "gone")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)

parse_run_id <- function(file) {
  basename(file) |>
    stringr::str_remove("^Output_Ne_SnowyOwl_GONE_") |>
    stringr::str_remove("\\.GONE_numeric$")
}

parse_output_run_id <- function(file) {
  basename(file) |>
    stringr::str_remove("^OUTPUT_SnowyOwl_GONE_") |>
    stringr::str_remove("\\.GONE_numeric$")
}

run_sample_set <- function(run_id) {
  dplyr::case_when(
    stringr::str_detect(run_id, "^all64") ~ "All 64",
    stringr::str_detect(run_id, "^highdepth24") ~ "High-depth 24",
    TRUE ~ "Other"
  )
}

run_filter_label <- function(run_id) {
  dplyr::case_when(
    stringr::str_detect(run_id, "maf005$") ~ "MAF >= 0.005",
    stringr::str_detect(run_id, "maf01$") ~ "MAF >= 0.01",
    stringr::str_detect(run_id, "maf02$") ~ "MAF >= 0.02",
    stringr::str_detect(run_id, "poly$") ~ "Polymorphic sites",
    TRUE ~ run_id
  )
}

run_filter_order <- function(filter_label) {
  match(filter_label, c("Polymorphic sites", "MAF >= 0.005", "MAF >= 0.01", "MAF >= 0.02"))
}

read_gone_ne <- function(file) {
  if (!file.exists(file)) {
    stop("Missing expected GONE file: ", file, call. = FALSE)
  }

  run_id <- parse_run_id(file)

  readr::read_tsv(
    file,
    skip = 1,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_double())
  ) |>
    dplyr::rename(generation = "Generation", ne = "Geometric_mean") |>
    dplyr::mutate(
      run_id = run_id,
      sample_set = run_sample_set(.data$run_id),
      filter = run_filter_label(.data$run_id),
      filter_order = run_filter_order(.data$filter),
      years_before_present = .data$generation * generation_time_years,
      calendar_year = present_year - .data$years_before_present
    )
}

extract_first_number <- function(lines, pattern) {
  line <- lines[stringr::str_detect(lines, pattern)][1]

  if (is.na(line)) {
    return(NA_real_)
  }

  as.numeric(stringr::str_match(line, paste0(pattern, "\\s*=?\\s*([-+0-9.eE]+)"))[, 2])
}

read_gone_input_diagnostics <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }

  lines <- readr::read_lines(file)
  chromosome_starts <- which(stringr::str_detect(lines, "^CHROMOSOME\\s+"))

  if (length(chromosome_starts) == 0) {
    return(NULL)
  }

  run_id <- parse_output_run_id(file)
  chromosome_ends <- c(chromosome_starts[-1] - 1, length(lines))

  chromosome_summaries <- Map(function(start, end) {
    block <- lines[start:end]
    chromosome <- stringr::str_match(block[1], "^CHROMOSOME\\s+(\\S+)")[, 2]
    n_records <- extract_first_number(block, "NSNP")
    n_snps_used <- extract_first_number(block, "NSNP_calculations")
    n_monomorphic <- extract_first_number(block, "NSNP_monomorphic")

    tibble::tibble(
      run_id = run_id,
      sample_set = run_sample_set(run_id),
      filter = run_filter_label(run_id),
      filter_order = run_filter_order(.data$filter),
      gone_chromosome = chromosome,
      n_individuals = extract_first_number(block, "NIND\\(real sample\\)"),
      n_records = n_records,
      n_snps_used = n_snps_used,
      n_monomorphic = n_monomorphic,
      n_zeroes = extract_first_number(block, "NSNP_zeroes"),
      n_more_than_two_alleles = extract_first_number(block, "NSNP_\\+2alleles"),
      n_individuals_corrected = extract_first_number(block, "NIND_corrected"),
      maf_threshold = extract_first_number(block, "freq_MAF"),
      monomorphic_fraction = n_monomorphic / n_records,
      used_fraction = n_snps_used / n_records
    )
  }, chromosome_starts, chromosome_ends)

  dplyr::bind_rows(chromosome_summaries)
}

calendar_breaks <- function(df, interval = 50) {
  x_range <- range(df$calendar_year, na.rm = TRUE)
  seq(ceiling(x_range[1] / interval) * interval, floor(x_range[2] / interval) * interval, by = interval)
}

calendar_labels <- function(years, label_interval = 100) {
  dplyr::if_else(
    years %% label_interval == 0,
    dplyr::if_else(years < 0, paste0(abs(years), " BC"), as.character(years)),
    ""
  )
}

add_recent_climate_periods <- function(plot, label_y = Inf, label_vjust = 1.35, label_size = 3) {
  plot +
    annotate("rect", xmin = 400, xmax = 765, ymin = -Inf, ymax = Inf, fill = "#8FA8A0", alpha = 0.14) +
    annotate("rect", xmin = 950, xmax = 1250, ymin = -Inf, ymax = Inf, fill = "#D08A4B", alpha = 0.12) +
    annotate("rect", xmin = 1450, xmax = 1850, ymin = -Inf, ymax = Inf, fill = "#9AA8A1", alpha = 0.16) +
    annotate("text", x = 582.5, y = label_y, label = "Dark Ages Cold Period", vjust = label_vjust, size = label_size, color = "#3E5B55") +
    annotate("text", x = 1100, y = label_y, label = "Medieval Warm Period", vjust = label_vjust, size = label_size, color = "#7A4A20") +
    annotate("text", x = 1650, y = label_y, label = "Little Ice Age", vjust = label_vjust, size = label_size, color = "#4C5B57")
}

ne_files <- list.files(here(gone_dir), pattern = "^Output_Ne_SnowyOwl_GONE_.*\\.GONE_numeric$", full.names = TRUE)
ne_files <- ne_files[!stringr::str_detect(basename(ne_files), "step6_noContam_noRel")]

if (length(ne_files) == 0) {
  stop("No sensitivity GONE Ne files found in ", here(gone_dir), call. = FALSE)
}

gone_ne <- dplyr::bind_rows(lapply(ne_files, read_gone_ne)) |>
  dplyr::arrange(.data$sample_set, .data$filter_order, .data$generation)

gone_plot_data <- gone_ne |>
  dplyr::filter(
    .data$generation >= plot_min_generation,
    .data$years_before_present <= plot_max_years_before_present
  )

gone_summary <- gone_ne |>
  dplyr::group_by(.data$run_id, .data$sample_set, .data$filter, .data$filter_order) |>
  dplyr::summarise(
    n_generations = dplyr::n(),
    generation_time_years = generation_time_years,
    present_year = present_year,
    min_ne_full_output = min(.data$ne, na.rm = TRUE),
    max_ne_full_output = max(.data$ne, na.rm = TRUE),
    most_recent_ne_full_output = .data$ne[which.min(.data$generation)][1],
    oldest_ne_full_output = .data$ne[which.max(.data$generation)][1],
    max_years_before_present_full_output = max(.data$years_before_present, na.rm = TRUE),
    plot_min_generation = plot_min_generation,
    plot_max_years_before_present = plot_max_years_before_present,
    most_recent_ne_plotted = .data$ne[which(.data$generation >= plot_min_generation)[1]],
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$sample_set, .data$filter_order)

key_generations <- c(1, 5, 10, 25, 50, 100, 150, 200)

gone_key_generation_table <- gone_ne |>
  dplyr::filter(.data$generation %in% key_generations) |>
  dplyr::select(sample_set, filter, filter_order, run_id, generation, ne) |>
  dplyr::arrange(.data$sample_set, .data$filter_order, .data$generation)

log_files <- list.files(here(gone_dir), pattern = "^OUTPUT_SnowyOwl_GONE_.*\\.GONE_numeric$", full.names = TRUE)
log_files <- log_files[!stringr::str_detect(basename(log_files), "step6_noContam_noRel")]
gone_input_diagnostics <- dplyr::bind_rows(lapply(log_files, read_gone_input_diagnostics))

readr::write_csv(gone_ne, file.path(out_tab_dir, "SnowyOwl_GONE_sensitivity_Ne_timeseries.csv"))
readr::write_csv(gone_summary, file.path(out_tab_dir, "SnowyOwl_GONE_sensitivity_summary.csv"))
readr::write_csv(gone_key_generation_table, file.path(out_tab_dir, "SnowyOwl_GONE_sensitivity_key_generations.csv"))

old_gone_file <- here("..", "old", "re", "Output_Ne_CHICKEN-HOM_autosomal_rmind_PASS")
if (file.exists(old_gone_file)) {
  old_gone_ne <- readr::read_tsv(
    old_gone_file,
    skip = 1,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_double())
  ) |>
    dplyr::rename(generation = "Generation", ne = "Geometric_mean") |>
    dplyr::mutate(
      run_id = "old_CHICKEN-HOM_autosomal_rmind_PASS",
      sample_set = "Old run",
      filter = "Old run",
      filter_order = 99,
      years_before_present = .data$generation * generation_time_years,
      calendar_year = present_year - .data$years_before_present
    )

  old_new_comparison <- dplyr::bind_rows(gone_ne, old_gone_ne) |>
    dplyr::filter(.data$generation %in% key_generations) |>
    dplyr::select(sample_set, filter, run_id, generation, ne) |>
    dplyr::arrange(.data$sample_set, .data$filter, .data$generation)

  readr::write_csv(old_new_comparison, file.path(out_tab_dir, "SnowyOwl_GONE_old_sensitivity_comparison.csv"))
}

if (nrow(gone_input_diagnostics) > 0) {
  gone_input_summary <- gone_input_diagnostics |>
    dplyr::group_by(.data$run_id, .data$sample_set, .data$filter, .data$filter_order) |>
    dplyr::summarise(
      n_gone_chromosomes = dplyr::n(),
      mean_n_individuals = mean(.data$n_individuals, na.rm = TRUE),
      total_records = sum(.data$n_records, na.rm = TRUE),
      total_snps_used = sum(.data$n_snps_used, na.rm = TRUE),
      total_monomorphic = sum(.data$n_monomorphic, na.rm = TRUE),
      mean_maf_threshold = mean(.data$maf_threshold, na.rm = TRUE),
      mean_monomorphic_fraction = mean(.data$monomorphic_fraction, na.rm = TRUE),
      mean_used_fraction = mean(.data$used_fraction, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$sample_set, .data$filter_order)

  readr::write_csv(gone_input_diagnostics, file.path(out_tab_dir, "SnowyOwl_GONE_sensitivity_input_diagnostics_by_chromosome.csv"))
  readr::write_csv(gone_input_summary, file.path(out_tab_dir, "SnowyOwl_GONE_sensitivity_input_diagnostics_summary.csv"))
}

gone_sensitivity_colors <- c(
  "Polymorphic sites" = "#4B4B4B",
  "MAF >= 0.005" = "#8C6D31",
  "MAF >= 0.01" = "#C17D11",
  "MAF >= 0.02" = "#B8892D"
)

gone_log_y_limits <- range(gone_plot_data$ne[is.finite(gone_plot_data$ne) & gone_plot_data$ne > 0], na.rm = TRUE)
gone_plot <- gone_plot_data |>
  ggplot(aes(x = .data$calendar_year, y = .data$ne, color = .data$filter, linetype = .data$sample_set)) |>
  add_recent_climate_periods(label_size = 2.8) +
  geom_line(linewidth = 0.45, alpha = 0.9) +
  scale_x_continuous(
    name = "Calendar year",
    breaks = calendar_breaks(gone_plot_data, interval = 50),
    labels = calendar_labels,
    minor_breaks = NULL
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = gone_sensitivity_colors, name = "Variant set") +
  scale_linetype_manual(values = c("All 64" = "solid", "High-depth 24" = "dashed"), name = "Sample set") +
  labs(
    y = expression(N[e]),
    title = "GONE sensitivity to sample depth and variant filtering",
    subtitle = paste0("Generation time = ", generation_time_years, " years; recombination rate = 6.63 cM/Mb")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gone_facet_plot <- gone_plot_data |>
  ggplot(aes(x = .data$calendar_year, y = .data$ne, color = .data$filter)) |>
  add_recent_climate_periods(label_size = 2.8) +
  geom_line(linewidth = 0.45, alpha = 0.9) +
  facet_wrap(~sample_set, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    name = "Calendar year",
    breaks = calendar_breaks(gone_plot_data, interval = 50),
    labels = calendar_labels,
    minor_breaks = NULL
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = gone_sensitivity_colors, name = "Variant set") +
  labs(
    y = expression(N[e]),
    title = "GONE sensitivity runs by sample set",
    subtitle = paste0("Linear scale; all runs shown as inferred")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  )

gone_log_plot <- gone_plot_data |>
  ggplot(aes(x = .data$calendar_year, y = .data$ne, color = .data$filter, linetype = .data$sample_set)) +
  annotate("rect", xmin = 400, xmax = 765, ymin = gone_log_y_limits[1], ymax = gone_log_y_limits[2], fill = "#8FA8A0", alpha = 0.14) +
  annotate("rect", xmin = 950, xmax = 1250, ymin = gone_log_y_limits[1], ymax = gone_log_y_limits[2], fill = "#D08A4B", alpha = 0.12) +
  annotate("rect", xmin = 1450, xmax = 1850, ymin = gone_log_y_limits[1], ymax = gone_log_y_limits[2], fill = "#9AA8A1", alpha = 0.16) +
  annotate("text", x = 582.5, y = gone_log_y_limits[2], label = "Dark Ages Cold Period", vjust = 1.35, size = 2.8, color = "#3E5B55") +
  annotate("text", x = 1100, y = gone_log_y_limits[2], label = "Medieval Warm Period", vjust = 1.35, size = 2.8, color = "#7A4A20") +
  annotate("text", x = 1650, y = gone_log_y_limits[2], label = "Little Ice Age", vjust = 1.35, size = 2.8, color = "#4C5B57") +
  geom_line(linewidth = 0.45, alpha = 0.9) +
  scale_x_continuous(
    name = "Calendar year",
    breaks = calendar_breaks(gone_plot_data, interval = 50),
    labels = calendar_labels,
    minor_breaks = NULL
  ) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_color_manual(values = gone_sensitivity_colors, name = "Variant set") +
  scale_linetype_manual(values = c("All 64" = "solid", "High-depth 24" = "dashed"), name = "Sample set") +
  labs(
    y = expression(N[e]),
    title = "GONE sensitivity to sample depth and variant filtering",
    subtitle = "Log scale shows both inflated and plausible runs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  )

gone_long_plot_data <- gone_ne |>
  dplyr::filter(.data$generation >= plot_min_generation)

gone_long_plot <- gone_long_plot_data |>
  ggplot(aes(x = .data$calendar_year, y = .data$ne, color = .data$filter)) +
  geom_line(linewidth = 0.4, alpha = 0.9) +
  facet_wrap(~sample_set, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    name = "Calendar year",
    breaks = calendar_breaks(gone_long_plot_data, interval = 500),
    labels = function(x) calendar_labels(x, label_interval = 500),
    minor_breaks = NULL
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = gone_sensitivity_colors, name = "Variant set") +
  labs(
    y = expression(N[e]),
    title = "GONE sensitivity runs across full inferred time span",
    subtitle = paste0("Generation time = ", generation_time_years, " years; free y-axis by sample set")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_sensitivity_Ne_timeseries.png"), gone_plot, width = 8, height = 5, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_sensitivity_Ne_timeseries.pdf"), gone_plot, width = 8, height = 5)
ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_sensitivity_Ne_by_sample_set.png"), gone_facet_plot, width = 7, height = 6.5, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_sensitivity_Ne_by_sample_set.pdf"), gone_facet_plot, width = 7, height = 6.5)
ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_sensitivity_Ne_timeseries_log10.png"), gone_log_plot, width = 8, height = 5, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_sensitivity_Ne_timeseries_log10.pdf"), gone_log_plot, width = 8, height = 5)
ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_sensitivity_Ne_by_sample_set_longterm.png"), gone_long_plot, width = 7, height = 6.5, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_sensitivity_Ne_by_sample_set_longterm.pdf"), gone_long_plot, width = 7, height = 6.5)

primary_run <- gone_ne |>
  dplyr::filter(.data$sample_set == "High-depth 24", .data$filter == "Polymorphic sites")

if (nrow(primary_run) == 0) {
  primary_run <- gone_ne |>
    dplyr::filter(.data$sample_set == "High-depth 24") |>
    dplyr::arrange(.data$filter_order) |>
    dplyr::group_by(.data$run_id) |>
    dplyr::ungroup()
}

primary_plot_data <- primary_run |>
  dplyr::filter(
    .data$generation >= plot_min_generation,
    .data$years_before_present <= plot_max_years_before_present
  )

primary_summary <- primary_run |>
  dplyr::summarise(
    run_id = dplyr::first(.data$run_id),
    sample_set = dplyr::first(.data$sample_set),
    filter = dplyr::first(.data$filter),
    n_generations = dplyr::n(),
    generation_time_years = generation_time_years,
    present_year = present_year,
    min_ne_full_output = min(.data$ne, na.rm = TRUE),
    max_ne_full_output = max(.data$ne, na.rm = TRUE),
    most_recent_ne_full_output = .data$ne[which.min(.data$generation)][1],
    oldest_ne_full_output = .data$ne[which.max(.data$generation)][1],
    max_years_before_present_full_output = max(.data$years_before_present, na.rm = TRUE),
    plot_min_generation = plot_min_generation,
    plot_max_years_before_present = plot_max_years_before_present,
    most_recent_ne_plotted = primary_plot_data$ne[which.min(primary_plot_data$generation)][1],
    oldest_ne_plotted = primary_plot_data$ne[which.max(primary_plot_data$generation)][1],
    .groups = "drop"
  )

readr::write_csv(primary_run, file.path(out_tab_dir, "SnowyOwl_GONE_Ne_timeseries.csv"))
readr::write_csv(primary_summary, file.path(out_tab_dir, "SnowyOwl_GONE_summary.csv"))

primary_plot <- primary_plot_data |>
  ggplot(aes(x = .data$calendar_year, y = .data$ne)) |>
  add_recent_climate_periods(label_size = 2.8) +
  geom_point(color = "black", size = 0.35) +
  geom_line(color = "black", linewidth = 0.35, linetype = "dotted") +
  scale_x_continuous(
    name = "Calendar year",
    breaks = calendar_breaks(primary_plot_data, interval = 50),
    labels = calendar_labels,
    minor_breaks = NULL
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    y = expression(N[e]),
    title = "Recent effective population size inferred by GONE",
    subtitle = paste0("High-depth individuals; generation time = ", generation_time_years, " years; recombination rate = 6.63 cM/Mb")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank()
  )

ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_Ne_timeseries.png"), primary_plot, width = 7, height = 4.5, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_GONE_Ne_timeseries.pdf"), primary_plot, width = 7, height = 4.5)

message("Wrote GONE sensitivity summaries to: ", out_tab_dir)
message("Wrote GONE sensitivity figures to: ", out_fig_dir)
