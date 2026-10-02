# Characterize FST outlier windows by genes, repeats, and distance to centromeres.

source("scripts/00_setup.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

out_tab_dir <- here(paths$results_tables, "pixy")
out_fig_dir <- here(paths$results_figures, "pixy")
dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

fst_window_file <- file.path(out_tab_dir, "SnowyOwl_regions_max_fst_by_window.csv")
recurrent_outlier_file <- file.path(out_tab_dir, "SnowyOwl_regions_recurrent_fst_outlier_windows.csv")
top_pairwise_file <- file.path(out_tab_dir, "SnowyOwl_regions_top200_pairwise_fst_windows.csv")
gff_file <- here(paths$genome_gff %||% "data/external/genome/bBubSca1.1.hap1.gff")
repeatmasker_file <- here(paths$repeatmasker_out %||% "data/external/genome/hap1_curated.fasta.inter.fa.out")
centromere_file <- here(paths$centromere_bed %||% "data/metadata/centromeres/putative_centromeres.bed")
pericentromere_flank_bp <- paths$pericentromere_flank_bp %||% 1000000

read_required_csv <- function(file) {
  if (!file.exists(file)) {
    stop("Missing required file: ", file, call. = FALSE)
  }
  readr::read_csv(file, show_col_types = FALSE)
}

extract_gff_attr <- function(attributes, key) {
  match <- stringr::str_match(attributes, paste0("(^|;)", key, "=([^;]+)"))[, 3]
  dplyr::na_if(match, "")
}

read_genes <- function(file) {
  if (!file.exists(file)) {
    warning("Missing GFF file; gene overlap tables will be empty: ", file, call. = FALSE)
    return(tibble::tibble())
  }

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
      gene_id = extract_gff_attr(.data$attributes, "ID"),
      gene_name = extract_gff_attr(.data$attributes, "Name"),
      product = extract_gff_attr(.data$attributes, "product")
    ) |>
    dplyr::select(chromosome, start, end, gene_id, gene_name, product)
}

normalize_repeat_chromosome <- function(chromosome) {
  dplyr::case_when(
    stringr::str_detect(chromosome, "^RL_\\d+$") ~ stringr::str_replace(chromosome, "^RL_", "SUPER_"),
    chromosome == "Z" ~ "SUPER_Z",
    chromosome == "W" ~ "SUPER_W",
    TRUE ~ chromosome
  )
}

read_repeats <- function(file) {
  if (!file.exists(file)) {
    warning("Missing RepeatMasker file; repeat overlap tables will be empty: ", file, call. = FALSE)
    return(tibble::tibble())
  }

  readr::read_table(
    file,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ) |>
    dplyr::transmute(
      chromosome = normalize_repeat_chromosome(.data$X5),
      start = as.numeric(.data$X6),
      end = as.numeric(.data$X7),
      repeat_name = .data$X10,
      repeat_class = .data$X11
    ) |>
    dplyr::filter(stringr::str_detect(.data$chromosome, "^SUPER_"), !is.na(.data$start), !is.na(.data$end))
}

read_centromeres <- function(file) {
  if (!file.exists(file)) {
    warning("Missing centromere file; distance-to-centromere tables will be empty: ", file, call. = FALSE)
    return(tibble::tibble())
  }

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

  if (nrow(bed) == 0) {
    return(tibble::tibble())
  }

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

overlap_rows <- function(windows, features) {
  if (nrow(features) == 0 || nrow(windows) == 0) {
    return(tibble::tibble())
  }

  split_features <- split(features, features$chromosome)

  dplyr::bind_rows(lapply(seq_len(nrow(windows)), function(i) {
    chrom <- as.character(windows$chromosome[i])
    feat <- split_features[[chrom]]
    if (is.null(feat)) {
      return(tibble::tibble())
    }

    hits <- feat |>
      dplyr::filter(.data$start <= windows$window_pos_2[i], .data$end >= windows$window_pos_1[i])

    if (nrow(hits) == 0) {
      return(tibble::tibble())
    }

    hits |>
      dplyr::mutate(
        window_id = windows$window_id[i],
        window_start = windows$window_pos_1[i],
        window_end = windows$window_pos_2[i]
      )
  }))
}

summarise_repeat_overlap <- function(windows, repeats) {
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
      repeat_classes = paste(sort(unique(.data$repeat_class)), collapse = ";"),
      top_repeat_class = names(sort(tapply(.data$overlap_bp, .data$repeat_class, sum), decreasing = TRUE))[1],
      .groups = "drop"
    ) |>
    dplyr::right_join(tibble::tibble(window_id = windows$window_id), by = "window_id") |>
    dplyr::mutate(
      n_repeat_hits = tidyr::replace_na(.data$n_repeat_hits, 0L),
      repeat_bp = tidyr::replace_na(.data$repeat_bp, 0),
      repeat_classes = tidyr::replace_na(.data$repeat_classes, ""),
      top_repeat_class = tidyr::replace_na(.data$top_repeat_class, "")
    )
}

add_centromere_distance <- function(windows, centromeres) {
  if (nrow(centromeres) == 0) {
    return(windows)
  }

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

overlaps_any <- function(query, subject) {
  vapply(
    seq_len(nrow(query)),
    function(i) {
      any(
        subject$chromosome == query$chromosome[i] &
          subject$start <= query$window_pos_2[i] &
          subject$end >= query$window_pos_1[i]
      )
    },
    logical(1)
  )
}

annotate_centromere_regions <- function(windows, centromeres, flank_bp) {
  if (nrow(centromeres) == 0) {
    return(
      windows |>
        dplyr::mutate(
          centromere_region = factor("Other", levels = c("Centromere", "Pericentromere", "Other"))
        )
    )
  }

  pericentromeres <- centromeres |>
    dplyr::mutate(
      start = pmax(1, .data$start - flank_bp),
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

make_window_id <- function(df) {
  df |>
    dplyr::mutate(
      window_id = paste(
        .data$chromosome,
        sprintf("%.0f", as.numeric(.data$window_pos_1)),
        sprintf("%.0f", as.numeric(.data$window_pos_2)),
        sep = ":"
      )
    )
}

fst_windows <- read_required_csv(fst_window_file) |>
  make_window_id()

recurrent_outliers <- read_required_csv(recurrent_outlier_file) |>
  make_window_id()

top_pairwise <- read_required_csv(top_pairwise_file) |>
  make_window_id()

genes <- read_genes(gff_file)
repeats <- read_repeats(repeatmasker_file)
centromeres <- read_centromeres(centromere_file)

recurrent_windows <- recurrent_outliers |>
  dplyr::left_join(
    fst_windows |>
      dplyr::select(
        window_id,
        chromosome,
        window_pos_1,
        window_pos_2,
        pair_at_max_fst,
        no_snps_at_max_fst,
        mean_fst_across_pairs,
        n_pairs
      ),
    by = c("window_id", "chromosome", "window_pos_1", "window_pos_2")
  )

gene_hits <- overlap_rows(recurrent_windows, genes)
repeat_summary <- summarise_repeat_overlap(recurrent_windows, repeats)

recurrent_annotated <- recurrent_windows |>
  dplyr::left_join(
    gene_hits |>
      dplyr::group_by(.data$window_id) |>
      dplyr::summarise(
        n_genes = dplyr::n_distinct(.data$gene_id),
        gene_ids = paste(sort(unique(stats::na.omit(.data$gene_id))), collapse = ";"),
        gene_names = paste(sort(unique(stats::na.omit(.data$gene_name))), collapse = ";"),
        products = paste(sort(unique(stats::na.omit(.data$product))), collapse = " | "),
        .groups = "drop"
      ),
    by = "window_id"
  ) |>
  dplyr::left_join(repeat_summary, by = "window_id") |>
  add_centromere_distance(centromeres) |>
  dplyr::mutate(
    n_genes = tidyr::replace_na(.data$n_genes, 0L),
    distance_to_centromere_mb = .data$nearest_centromere_distance_bp / 1e6,
    repeat_fraction = pmin(.data$repeat_bp / (.data$window_pos_2 - .data$window_pos_1 + 1), 1)
  ) |>
  dplyr::arrange(dplyr::desc(.data$n_outlier_pairs), dplyr::desc(.data$max_fst))

all_windows_centromere_distance <- fst_windows |>
  add_centromere_distance(centromeres) |>
  dplyr::mutate(distance_to_centromere_mb = .data$nearest_centromere_distance_bp / 1e6)

all_windows_centromere_regions <- all_windows_centromere_distance |>
  annotate_centromere_regions(centromeres, pericentromere_flank_bp) |>
  dplyr::mutate(mean_fst_display = pmax(.data$mean_fst_across_pairs, 0, na.rm = TRUE))

centromere_region_summary <- all_windows_centromere_regions |>
  dplyr::group_by(.data$centromere_region) |>
  dplyr::summarise(
    metric = "mean_fst_across_pairs_negative_set_to_zero",
    n_windows = dplyr::n(),
    mean = mean(.data$mean_fst_display, na.rm = TRUE),
    median = median(.data$mean_fst_display, na.rm = TRUE),
    min = min(.data$mean_fst_display, na.rm = TRUE),
    max = max(.data$mean_fst_display, na.rm = TRUE),
    q25 = as.numeric(stats::quantile(.data$mean_fst_display, 0.25, na.rm = TRUE)),
    q75 = as.numeric(stats::quantile(.data$mean_fst_display, 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::select(.data$metric, .data$centromere_region, dplyr::everything())

centromere_region_test <- tibble::tibble(
  metric = "mean_fst_across_pairs_negative_set_to_zero",
  test = "Kruskal-Wallis rank-sum test",
  grouping = "Centromere vs Pericentromere vs Other",
  statistic = unname(stats::kruskal.test(mean_fst_display ~ centromere_region, data = all_windows_centromere_regions)$statistic),
  df = unname(stats::kruskal.test(mean_fst_display ~ centromere_region, data = all_windows_centromere_regions)$parameter),
  p_value = stats::kruskal.test(mean_fst_display ~ centromere_region, data = all_windows_centromere_regions)$p.value
)

centromere_distance_summary <- all_windows_centromere_distance |>
  dplyr::filter(!is.na(.data$distance_to_centromere_mb)) |>
  dplyr::mutate(
    distance_bin_mb = cut(
      .data$distance_to_centromere_mb,
      breaks = c(-Inf, 0, 1, 5, 10, 25, 50, Inf),
      labels = c("0", "0-1", "1-5", "5-10", "10-25", "25-50", ">50")
    )
  ) |>
  dplyr::group_by(.data$distance_bin_mb) |>
  dplyr::summarise(
    n_windows = dplyr::n(),
    mean_window_mean_fst = mean(.data$mean_fst_across_pairs, na.rm = TRUE),
    median_window_mean_fst = median(.data$mean_fst_across_pairs, na.rm = TRUE),
    q95_window_mean_fst = as.numeric(stats::quantile(.data$mean_fst_across_pairs, 0.95, na.rm = TRUE)),
    mean_max_fst = mean(.data$max_fst, na.rm = TRUE),
    median_max_fst = median(.data$max_fst, na.rm = TRUE),
    q95_max_fst = as.numeric(stats::quantile(.data$max_fst, 0.95, na.rm = TRUE)),
    max_fst = max(.data$max_fst, na.rm = TRUE),
    .groups = "drop"
  )

recurrent_outlier_centromere_distance_summary <- recurrent_annotated |>
  dplyr::filter(!is.na(.data$distance_to_centromere_mb)) |>
  dplyr::mutate(
    distance_bin_mb = cut(
      .data$distance_to_centromere_mb,
      breaks = c(-Inf, 0, 1, 5, 10, 25, 50, Inf),
      labels = c("0", "0-1", "1-5", "5-10", "10-25", "25-50", ">50")
    )
  ) |>
  dplyr::group_by(.data$distance_bin_mb) |>
  dplyr::summarise(
    n_recurrent_outlier_windows = dplyr::n(),
    median_max_fst = median(.data$max_fst, na.rm = TRUE),
    max_fst = max(.data$max_fst, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::right_join(
    centromere_distance_summary |>
      dplyr::select(distance_bin_mb, n_total_windows = n_windows),
    by = "distance_bin_mb"
  ) |>
  dplyr::mutate(
    n_recurrent_outlier_windows = tidyr::replace_na(.data$n_recurrent_outlier_windows, 0L),
    recurrent_outlier_fraction = .data$n_recurrent_outlier_windows / .data$n_total_windows
  )

outlier_feature_summary <- recurrent_annotated |>
  dplyr::summarise(
    n_recurrent_outlier_windows = dplyr::n(),
    n_windows_with_gene_overlap = sum(.data$n_genes > 0, na.rm = TRUE),
    n_windows_with_repeat_overlap = sum(.data$n_repeat_hits > 0, na.rm = TRUE),
    median_repeat_fraction = median(.data$repeat_fraction, na.rm = TRUE),
    median_distance_to_centromere_mb = median(.data$distance_to_centromere_mb, na.rm = TRUE)
  )

repeat_class_summary <- recurrent_annotated |>
  dplyr::filter(.data$top_repeat_class != "") |>
  dplyr::count(.data$top_repeat_class, sort = TRUE, name = "n_outlier_windows")

readr::write_csv(
  recurrent_annotated,
  file.path(out_tab_dir, "SnowyOwl_regions_recurrent_fst_outlier_windows_annotated.csv")
)

readr::write_csv(
  gene_hits,
  file.path(out_tab_dir, "SnowyOwl_regions_recurrent_fst_outlier_gene_overlaps.csv")
)

readr::write_csv(
  outlier_feature_summary,
  file.path(out_tab_dir, "SnowyOwl_regions_recurrent_fst_outlier_feature_summary.csv")
)

readr::write_csv(
  repeat_class_summary,
  file.path(out_tab_dir, "SnowyOwl_regions_recurrent_fst_outlier_repeat_class_summary.csv")
)

readr::write_csv(
  centromere_distance_summary,
  file.path(out_tab_dir, "SnowyOwl_regions_fst_distance_to_centromere_summary.csv")
)

readr::write_csv(
  centromere_region_summary,
  file.path(out_tab_dir, "SnowyOwl_regions_mean_fst_centromere_regions_summary.csv")
)

readr::write_csv(
  centromere_region_test,
  file.path(out_tab_dir, "SnowyOwl_regions_mean_fst_centromere_regions_kruskal_test.csv")
)

readr::write_csv(
  recurrent_outlier_centromere_distance_summary,
  file.path(out_tab_dir, "SnowyOwl_regions_recurrent_fst_outlier_distance_to_centromere_summary.csv")
)

distance_plot <- all_windows_centromere_distance |>
  dplyr::filter(!is.na(.data$distance_to_centromere_mb)) |>
  ggplot(aes(x = .data$distance_to_centromere_mb, y = .data$mean_fst_across_pairs)) +
  geom_point(size = 0.35, alpha = 0.25, color = "#2f4f4f", na.rm = TRUE) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, color = "#b35c44", linewidth = 0.6) +
  scale_x_continuous(trans = "sqrt") +
  labs(
    x = "Distance to nearest annotated centromere (Mb, square-root scale)",
    y = "Mean pairwise FST across regional comparisons",
    title = "Mean FST versus distance to putative centromeres"
  )

ggsave(
  file.path(out_fig_dir, "SnowyOwl_regions_mean_fst_distance_to_centromere.png"),
  distance_plot,
  width = 6,
  height = 4,
  dpi = 350
)

centromere_region_plot <- all_windows_centromere_regions |>
  dplyr::filter(!is.na(.data$centromere_region), !is.na(.data$mean_fst_display)) |>
  ggplot(aes(x = .data$centromere_region, y = .data$mean_fst_display, fill = .data$centromere_region)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, na.rm = TRUE) +
  geom_jitter(width = 0.18, size = 0.25, alpha = 0.16, na.rm = TRUE) +
  scale_fill_manual(
    values = c("Centromere" = "#7a3b2e", "Pericentromere" = "#b35c44", "Other" = "#2f4f4f"),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = "Mean pairwise FST",
    title = "Mean FST by putative centromere region",
    subtitle = paste0(
      "Negative estimates set to zero; Kruskal-Wallis p = ",
      format_p_value(centromere_region_test$p_value)
    )
  ) +
  theme(axis.text.x = element_text(size = 9))

ggsave(
  file.path(out_fig_dir, "SnowyOwl_regions_mean_fst_centromere_regions_boxplot.png"),
  centromere_region_plot,
  width = 5.5,
  height = 4,
  dpi = 350
)

ggsave(
  file.path(out_fig_dir, "SnowyOwl_regions_mean_fst_distance_to_centromere.pdf"),
  distance_plot,
  width = 6,
  height = 4
)

message("Wrote FST outlier annotation tables to: ", out_tab_dir)
message("Wrote FST-centromere distance figure to: ", out_fig_dir)
