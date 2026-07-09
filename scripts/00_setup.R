# Project setup shared by downstream analysis scripts.

required_packages <- c(
  "here",
  "readr",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "yaml"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running analyses: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(yaml)
})

paths_file <- if (file.exists(here("config", "local_paths.yml"))) {
  here("config", "local_paths.yml")
} else {
  here("config", "paths.yml")
}

paths <- yaml::read_yaml(paths_file)

pop_colors <- c(
  "SK" = "#F38400",
  "GR" = "#BE0032",
  "FN" = "#875692",
  "WR" = "#F3C300",
  "NY" = "#A1CAF1"
)

roh_length_colors <- c(
  "0.3-1 Mb" = "#C8C1B8",
  "1-5 Mb" = "#9E968C",
  ">5 Mb" = "#6F675F"
)

macrochromosomes <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)

chromosome_landscape_colors <- c(
  "Macrochromosome odd" = "lightgoldenrod4",
  "Macrochromosome even" = "lightgoldenrod2",
  "Microchromosome odd" = "lightcyan4",
  "Microchromosome even" = "lightcyan2",
  "Outlier" = "black"
)

region_display <- function(region) {
  dplyr::case_when(
    region == "FNM" ~ "FN",
    region == "GRL" ~ "GR",
    region == "NYS" ~ "NY",
    region == "SKW" ~ "SK",
    region == "WRG" ~ "WR",
    TRUE ~ region
  )
}

chromosome_class <- function(chrom_number) {
  dplyr::if_else(chrom_number %in% macrochromosomes, "Macrochromosome", "Microchromosome")
}

chromosome_landscape_group <- function(chrom_number) {
  dplyr::case_when(
    chrom_number %in% macrochromosomes & chrom_number %% 2 == 1 ~ "Macrochromosome odd",
    chrom_number %in% macrochromosomes & chrom_number %% 2 == 0 ~ "Macrochromosome even",
    !(chrom_number %in% macrochromosomes) & chrom_number %% 2 == 1 ~ "Microchromosome odd",
    !(chrom_number %in% macrochromosomes) & chrom_number %% 2 == 0 ~ "Microchromosome even",
    TRUE ~ "Microchromosome odd"
  )
}

format_p_value <- function(p_value) {
  dplyr::case_when(
    is.na(p_value) ~ "NA",
    p_value < 0.001 ~ format.pval(p_value, digits = 2, eps = 0.001),
    TRUE ~ signif(p_value, 3) |> as.character()
  )
}

theme_set(
  theme_classic(base_size = 11) +
    theme(
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
)
