# Plot ROH landscapes inspired by genome-wide ROH density figures.

source("scripts/00_setup.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

roh_dir <- paths$roh_dir %||% "data/external/roh"
autosomal_genome_size_bp <- paths$autosomal_genome_size_bp %||% 1173697940
roh_density_window_bp <- paths$roh_density_window_bp %||% 500000
roh_track_min_length_mb <- 1

hom_file <- here(roh_dir, "SnowyOwl_autosomes_noContam_ROH.hom")
indiv_file <- here(roh_dir, "SnowyOwl_autosomes_noContam_ROH.hom.indiv")
bins_file <- here(roh_dir, "SnowyOwl_autosomes_noContam_ROH.ROH_length_bins.tsv")
centromere_bed <- here(paths$centromere_bed %||% "data/metadata/centromeres/putative_centromeres.bed")
repeatmasker_file <- here(paths$repeatmasker_out %||% "data/external/genome/bBubSca1.1.hap1.fasta.out")

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

read_autosome_lengths <- function(file) {
  if (!file.exists(file)) {
    stop("Missing chromosome-size metadata: ", file, call. = FALSE)
  }

  readr::read_tsv(file, show_col_types = FALSE, na = c("NA", "")) |>
    dplyr::rename(
      chromosome = dplyr::any_of(c("Chr", "chromosome", "chr")),
      chromosome_size = dplyr::any_of(c("Chr_size", "chromosome_size", "chrom_size"))
    ) |>
    dplyr::mutate(
      chromosome = as.character(.data$chromosome),
      chrom_number = suppressWarnings(as.integer(stringr::str_extract(.data$chromosome, "(?<=SUPER_)\\d+"))),
      chromosome_size = as.numeric(.data$chromosome_size)
    ) |>
    dplyr::filter(
      stringr::str_detect(.data$chromosome, "^SUPER_"),
      !is.na(.data$chrom_number),
      !is.na(.data$chromosome_size)
    ) |>
    dplyr::arrange(.data$chrom_number) |>
    dplyr::select(chromosome, chrom_number, chromosome_size)
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

read_repeats <- function(file) {
  if (!file.exists(file)) {
    stop("Missing RepeatMasker file: ", file, call. = FALSE)
  }

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

overlap_rows <- function(windows, features) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install data.table before annotating ROH density windows with genomic features.", call. = FALSE)
  }

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

    midpoint <- (windows$window_pos_1[i] + windows$window_pos_2[i]) / 2
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

spearman_summary <- function(df, feature) {
  feature_df <- df |>
    dplyr::filter(!is.na(.data$roh_density), !is.na(.data[[feature]]))

  if (nrow(feature_df) < 3 || stats::sd(feature_df[[feature]], na.rm = TRUE) == 0) {
    return(tibble::tibble(test = "spearman", feature = feature, n = nrow(feature_df), estimate = NA_real_, p_value = NA_real_))
  }

  rho <- suppressWarnings(stats::cor(feature_df$roh_density, feature_df[[feature]], method = "spearman", use = "complete.obs"))
  p_value <- tryCatch(
    suppressWarnings(stats::cor.test(feature_df$roh_density, feature_df[[feature]], method = "spearman", exact = FALSE)[["p.value"]]),
    error = function(e) NA_real_
  )
  tibble::tibble(
    test = "spearman",
    feature = feature,
    n = nrow(feature_df),
    estimate = unname(rho),
    p_value = p_value
  )
}

wilcoxon_top_summary <- function(df, feature) {
  feature_df <- df |>
    dplyr::filter(!is.na(.data[[feature]]), !is.na(.data$high_roh_density_window))

  if (dplyr::n_distinct(feature_df$high_roh_density_window) < 2) {
    return(tibble::tibble(test = "wilcoxon_top50_vs_background", feature = feature, n = nrow(feature_df), estimate = NA_real_, p_value = NA_real_))
  }

  p_value <- tryCatch(
    suppressWarnings(stats::wilcox.test(feature_df[[feature]] ~ feature_df$high_roh_density_window, exact = FALSE)[["p.value"]]),
    error = function(e) NA_real_
  )
  tibble::tibble(
    test = "wilcoxon_top50_vs_background",
    feature = feature,
    n = nrow(feature_df),
    estimate = median(feature_df[[feature]][feature_df$high_roh_density_window], na.rm = TRUE) -
      median(feature_df[[feature]][!feature_df$high_roh_density_window], na.rm = TRUE),
    p_value = p_value
  )
}

add_genome_positions <- function(df, chrom_sizes) {
  chrom_offsets <- chrom_sizes |>
    dplyr::arrange(.data$chrom_number) |>
    dplyr::mutate(chrom_start = dplyr::lag(cumsum(.data$chromosome_size), default = 0))

  df |>
    dplyr::left_join(chrom_offsets, by = "chromosome") |>
    dplyr::mutate(
      genome_start = .data$chrom_start + .data$POS1,
      genome_end = .data$chrom_start + .data$POS2,
      genome_midpoint = (.data$genome_start + .data$genome_end) / 2
    )
}

make_genome_windows <- function(chrom_sizes, window_bp) {
  dplyr::bind_rows(lapply(seq_len(nrow(chrom_sizes)), function(i) {
    starts <- seq(1, chrom_sizes$chromosome_size[i], by = window_bp)
    tibble::tibble(
      chromosome = chrom_sizes$chromosome[i],
      chrom_number = chrom_sizes$chrom_number[i],
      window_pos_1 = starts,
      window_pos_2 = pmin(starts + window_bp - 1, chrom_sizes$chromosome_size[i])
    )
  }))
}

count_roh_overlap_by_window <- function(windows, segments, n_individuals) {
  overlap_counts <- vapply(seq_len(nrow(windows)), function(i) {
    overlapping_iids <- segments$IID[
      segments$chromosome == windows$chromosome[i] &
        segments$POS1 <= windows$window_pos_2[i] &
        segments$POS2 >= windows$window_pos_1[i]
    ]

    dplyr::n_distinct(overlapping_iids)
  }, integer(1))

  windows |>
    dplyr::mutate(
      n_individuals_with_roh = overlap_counts,
      roh_density = .data$n_individuals_with_roh / n_individuals
    )
}

chrom_sizes <- read_autosome_lengths(centromere_bed)
centromeres <- read_centromeres(centromere_bed)
repeats <- read_repeats(repeatmasker_file)

roh_segments <- read_plink_table(hom_file) |>
  dplyr::rename(chromosome = CHR) |>
  dplyr::mutate(
    region = infer_region(.data$IID),
    region_label = region_display(.data$region),
    roh_length_mb = .data$KB / 1000,
    length_class = dplyr::case_when(
      .data$roh_length_mb < 1 ~ "0.3-1 Mb",
      .data$roh_length_mb <= 5 ~ "1-5 Mb",
      TRUE ~ ">5 Mb"
    ),
    length_class = factor(.data$length_class, levels = names(roh_length_colors))
  ) |>
  add_genome_positions(chrom_sizes)

roh_indiv <- read_plink_table(indiv_file) |>
  dplyr::mutate(
    region = infer_region(.data$IID),
    region_label = region_display(.data$region),
    froh = (.data$KB * 1000) / autosomal_genome_size_bp
  )

roh_bins <- readr::read_tsv(bins_file, show_col_types = FALSE) |>
  dplyr::mutate(
    region = infer_region(.data$IID),
    region_label = region_display(.data$region),
    froh_0.3_1Mb = (.data$`ROH_0.3_1Mb_kb` * 1000) / autosomal_genome_size_bp,
    froh_1_5Mb = (.data$`ROH_1_5Mb_kb` * 1000) / autosomal_genome_size_bp,
    froh_gt5Mb = (.data$ROH_gt5Mb_kb * 1000) / autosomal_genome_size_bp
  )

axis_df <- chrom_sizes |>
  dplyr::mutate(
    chrom_start = dplyr::lag(cumsum(.data$chromosome_size), default = 0),
    axis_midpoint = .data$chrom_start + (.data$chromosome_size / 2),
    axis_label = as.character(.data$chrom_number)
  )

top_iids <- roh_indiv |>
  dplyr::arrange(dplyr::desc(.data$froh)) |>
  dplyr::slice_head(n = 7) |>
  dplyr::mutate(froh_rank_group = "Highest FROH")

bottom_iids <- roh_indiv |>
  dplyr::arrange(.data$froh) |>
  dplyr::slice_head(n = 7) |>
  dplyr::mutate(froh_rank_group = "Lowest FROH")

selected_iids <- dplyr::bind_rows(top_iids, bottom_iids) |>
  dplyr::mutate(
    individual_label = paste0(.data$region_label, " ", seq_len(dplyr::n())),
    plot_order = dplyr::row_number()
  )

long_roh_tracks <- roh_segments |>
  dplyr::filter(.data$roh_length_mb > roh_track_min_length_mb) |>
  dplyr::inner_join(
    selected_iids |>
      dplyr::select(IID, individual_label, plot_order, froh_rank_group, froh),
    by = "IID"
  ) |>
  dplyr::mutate(
    ymin = .data$plot_order - 0.36,
    ymax = .data$plot_order + 0.36
  )

selected_track_rows <- selected_iids |>
  dplyr::mutate(
    ymin = .data$plot_order - 0.36,
    ymax = .data$plot_order + 0.36
  )

long_roh_track_summary <- selected_iids |>
  dplyr::left_join(
    long_roh_tracks |>
      dplyr::group_by(.data$IID) |>
      dplyr::summarise(
        n_roh_gt1mb = dplyr::n(),
        total_roh_gt1mb_mb = sum(.data$roh_length_mb, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "IID"
  ) |>
  dplyr::mutate(
    n_roh_gt1mb = dplyr::coalesce(.data$n_roh_gt1mb, 0L),
    total_roh_gt1mb_mb = dplyr::coalesce(.data$total_roh_gt1mb_mb, 0)
  )

readr::write_csv(long_roh_track_summary, file.path(out_tab_dir, "SnowyOwl_ROH_top_bottom7_gt1Mb_summary.csv"))

roh_track_plot <- ggplot() +
  geom_rect(
    data = selected_track_rows,
    aes(xmin = 0, xmax = max(axis_df$chrom_start + axis_df$chromosome_size), ymin = .data$ymin, ymax = .data$ymax),
    fill = "grey95",
    color = NA
  ) +
  geom_rect(
    data = long_roh_tracks,
    aes(xmin = .data$genome_start, xmax = .data$genome_end, ymin = .data$ymin, ymax = .data$ymax, fill = .data$region_label),
    color = NA,
    alpha = 0.95
  ) +
  scale_fill_manual(values = pop_colors, name = "Region") +
  scale_x_continuous(
    breaks = axis_df$axis_midpoint,
    labels = axis_df$axis_label,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_reverse(
    breaks = selected_iids$plot_order,
    labels = selected_iids$individual_label,
    expand = expansion(add = 0.4)
  ) +
  labs(
    x = "Chromosome",
    y = NULL,
    title = "ROH longer than 1 Mb in individuals with highest and lowest FROH",
    subtitle = "Top seven individuals are shown above the bottom seven individuals"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 7),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

roh_burden_by_class <- roh_bins |>
  dplyr::select(IID, region_label, "froh_0.3_1Mb", "froh_1_5Mb", "froh_gt5Mb") |>
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
    length_class = factor(.data$length_class, levels = names(roh_length_colors))
  )

readr::write_csv(roh_burden_by_class, file.path(out_tab_dir, "SnowyOwl_ROH_individual_length_class_fraction.csv"))

roh_length_class_plot <- roh_burden_by_class |>
  ggplot(aes(x = .data$length_class, y = .data$froh, fill = .data$length_class)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, na.rm = TRUE) +
  geom_jitter(aes(color = .data$region_label), width = 0.18, size = 1.35, alpha = 0.8, na.rm = TRUE) +
  scale_fill_manual(values = roh_length_colors, guide = "none") +
  scale_color_manual(values = pop_colors, name = "Region") +
  labs(
    x = "ROH length class",
    y = expression(F[ROH]),
    title = "Distribution of ROH burden among length classes"
  ) +
  theme(legend.position = "bottom")

genome_windows <- make_genome_windows(chrom_sizes, roh_density_window_bp)

roh_density_windows <- count_roh_overlap_by_window(
  windows = genome_windows,
  segments = roh_segments,
  n_individuals = nrow(roh_indiv)
) |>
  dplyr::mutate(
    window_id = paste(.data$chromosome, .data$window_pos_1, .data$window_pos_2, sep = ":"),
    window_length_bp = .data$window_pos_2 - .data$window_pos_1 + 1
  ) |>
  dplyr::left_join(
    axis_df |>
      dplyr::select(chromosome, chrom_start),
    by = "chromosome"
  ) |>
  dplyr::mutate(
    genome_midpoint = .data$chrom_start + ((.data$window_pos_1 + .data$window_pos_2) / 2),
    chromosome_landscape_group = chromosome_landscape_group(.data$chrom_number)
  )

roh_density_features <- roh_density_windows |>
  dplyr::left_join(summarise_repeat_features(roh_density_windows, repeats), by = "window_id") |>
  add_centromere_distance(centromeres) |>
  dplyr::mutate(
    repeat_fraction = pmin(.data$repeat_bp / .data$window_length_bp, 1),
    satellite_fraction = pmin(.data$satellite_bp / .data$window_length_bp, 1),
    simple_repeat_fraction = pmin(.data$simple_repeat_bp / .data$window_length_bp, 1),
    cr1_fraction = pmin(.data$cr1_bp / .data$window_length_bp, 1),
    ltr_fraction = pmin(.data$ltr_bp / .data$window_length_bp, 1),
    unknown_repeat_fraction = pmin(.data$unknown_repeat_bp / .data$window_length_bp, 1),
    distance_to_centromere_mb = .data$nearest_centromere_distance_bp / 1e6,
    high_roh_density_window = dplyr::min_rank(dplyr::desc(.data$roh_density)) <= 50
  )

feature_columns <- c(
  "repeat_fraction",
  "satellite_fraction",
  "simple_repeat_fraction",
  "cr1_fraction",
  "ltr_fraction",
  "unknown_repeat_fraction",
  "distance_to_centromere_mb"
)

roh_density_feature_associations <- dplyr::bind_rows(
  dplyr::bind_rows(lapply(feature_columns, function(feature) spearman_summary(roh_density_features, feature))),
  dplyr::bind_rows(lapply(feature_columns, function(feature) wilcoxon_top_summary(roh_density_features, feature)))
) |>
  dplyr::mutate(p_value_label = format_p_value(.data$p_value))

roh_density_feature_summary <- roh_density_features |>
  dplyr::mutate(density_group = dplyr::if_else(.data$high_roh_density_window, "Top 50 ROH-density windows", "Remaining windows")) |>
  dplyr::group_by(.data$density_group) |>
  dplyr::summarise(
    n_windows = dplyr::n(),
    mean_roh_density = mean(.data$roh_density, na.rm = TRUE),
    median_roh_density = median(.data$roh_density, na.rm = TRUE),
    mean_repeat_fraction = mean(.data$repeat_fraction, na.rm = TRUE),
    median_repeat_fraction = median(.data$repeat_fraction, na.rm = TRUE),
    mean_satellite_fraction = mean(.data$satellite_fraction, na.rm = TRUE),
    median_satellite_fraction = median(.data$satellite_fraction, na.rm = TRUE),
    median_distance_to_centromere_mb = median(.data$distance_to_centromere_mb, na.rm = TRUE),
    n_centromere_overlapping_windows = sum(.data$distance_to_centromere_mb == 0, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(roh_density_features, file.path(out_tab_dir, "SnowyOwl_ROH_density_500kb_windows.csv"))
readr::write_csv(roh_density_features, file.path(out_tab_dir, "SnowyOwl_ROH_density_500kb_windows_genomic_features.csv"))
readr::write_csv(roh_density_feature_associations, file.path(out_tab_dir, "SnowyOwl_ROH_density_feature_associations.csv"))
readr::write_csv(roh_density_feature_summary, file.path(out_tab_dir, "SnowyOwl_ROH_density_feature_summary.csv"))

roh_density_features |>
  dplyr::arrange(dplyr::desc(.data$roh_density)) |>
  dplyr::slice_head(n = 50) |>
  readr::write_csv(file.path(out_tab_dir, "SnowyOwl_ROH_density_top50_500kb_windows.csv"))

roh_density_plot <- roh_density_features |>
  ggplot(aes(x = .data$genome_midpoint, y = .data$roh_density, color = .data$roh_density)) +
  geom_point(size = 0.45, alpha = 0.9, na.rm = TRUE) +
  scale_color_gradient(low = "#edf8e9", high = "#006d2c", name = "ROH density") +
  scale_x_continuous(
    breaks = axis_df$axis_midpoint,
    labels = axis_df$axis_label,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = "Chromosome",
    y = "Proportion of individuals",
    title = paste0("Genome-wide ROH density in ", roh_density_window_bp / 1000, " kb windows")
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 7),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

roh_density_repeat_plot <- roh_density_features |>
  ggplot(aes(x = .data$repeat_fraction, y = .data$roh_density)) +
  geom_point(aes(color = .data$high_roh_density_window), size = 0.8, alpha = 0.65, na.rm = TRUE) +
  geom_smooth(method = "loess", formula = y ~ x, color = "black", linewidth = 0.5, se = TRUE, na.rm = TRUE) +
  scale_color_manual(
    values = c("FALSE" = "grey55", "TRUE" = "#006d2c"),
    labels = c("FALSE" = "Other windows", "TRUE" = "Top 50 ROH-density windows"),
    name = NULL
  ) +
  labs(
    x = "Repeat fraction",
    y = "ROH density",
    title = "ROH density versus repeat content"
  ) +
  theme(legend.position = "bottom")

roh_density_centromere_plot <- roh_density_features |>
  dplyr::filter(!is.na(.data$distance_to_centromere_mb)) |>
  ggplot(aes(x = .data$distance_to_centromere_mb, y = .data$roh_density)) +
  geom_point(aes(color = .data$high_roh_density_window), size = 0.8, alpha = 0.65, na.rm = TRUE) +
  geom_smooth(method = "loess", formula = y ~ x, color = "black", linewidth = 0.5, se = TRUE, na.rm = TRUE) +
  scale_x_sqrt() +
  scale_color_manual(
    values = c("FALSE" = "grey55", "TRUE" = "#006d2c"),
    labels = c("FALSE" = "Other windows", "TRUE" = "Top 50 ROH-density windows"),
    name = NULL
  ) +
  labs(
    x = "Distance to nearest annotated centromere (Mb, square-root scale)",
    y = "ROH density",
    title = "ROH density versus distance to putative centromeres"
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_gt1Mb_top_bottom7_tracks.png"), roh_track_plot, width = 11, height = 5.2, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_gt1Mb_top_bottom7_tracks.pdf"), roh_track_plot, width = 11, height = 5.2)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_length_class_fraction.png"), roh_length_class_plot, width = 5.8, height = 4.4, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_length_class_fraction.pdf"), roh_length_class_plot, width = 5.8, height = 4.4)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_density_500kb_windows.png"), roh_density_plot, width = 11, height = 4.8, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_density_500kb_windows.pdf"), roh_density_plot, width = 11, height = 4.8)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_density_vs_repeat_fraction.png"), roh_density_repeat_plot, width = 5.8, height = 4.4, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_density_vs_repeat_fraction.pdf"), roh_density_repeat_plot, width = 5.8, height = 4.4)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_density_vs_centromere_distance.png"), roh_density_centromere_plot, width = 5.8, height = 4.4, dpi = 350)
ggsave(file.path(out_fig_dir, "SnowyOwl_ROH_density_vs_centromere_distance.pdf"), roh_density_centromere_plot, width = 5.8, height = 4.4)

message("Wrote ROH landscape summaries to: ", out_tab_dir)
message("Wrote ROH landscape figures to: ", out_fig_dir)
