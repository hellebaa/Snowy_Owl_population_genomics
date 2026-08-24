# Build a transparent sample metadata table from historical metadata and QC files.

source("scripts/00_setup.R")

required_metadata_packages <- c("readxl")
missing_metadata_packages <- required_metadata_packages[
  !vapply(required_metadata_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_metadata_packages) > 0) {
  stop(
    "Install missing packages before running metadata QC: ",
    paste(missing_metadata_packages, collapse = ", "),
    call. = FALSE
  )
}

old_dir <- normalizePath(file.path("..", "old"), mustWork = TRUE)

metadata_xlsx <- file.path(old_dir, "Metadata_curated_owls.xlsx")
qc_file <- file.path(old_dir, "QC_emilyversion_with_sex.txt")
rename_file <- file.path(old_dir, "sample_rename.txt")
location_file <- file.path(old_dir, "re", "Snowy_sample_locations.txt")

stopifnot(file.exists(metadata_xlsx), file.exists(qc_file), file.exists(rename_file))

extract_year <- function(x) {
  year <- stringr::str_extract(as.character(x), "(19|20)[0-9]{2}")
  suppressWarnings(as.integer(year))
}

normalize_sex <- function(x) {
  dplyr::case_when(
    stringr::str_to_upper(as.character(x)) %in% c("F", "FEMALE") ~ "Female",
    stringr::str_to_upper(as.character(x)) %in% c("M", "MALE") ~ "Male",
    TRUE ~ "Unknown"
  )
}

location_lookup <- tibble::tribble(
  ~region, ~population_final, ~latitude, ~longitude,
  "FNM", "FN", 69.74, 24.20,
  "GRL", "GR", 72.60, -23.20,
  "NYS", "NY", 42.81, -78.35,
  "SKW", "SK", 54.19, -106.04,
  "WRG", "WR", 71.37, 179.47
)

rename_map <- readr::read_table(
  rename_file,
  col_names = c("old_id", "sample_id"),
  col_types = readr::cols(.default = readr::col_character())
)

snowy <- readxl::read_excel(metadata_xlsx, sheet = "snowy_owl", .name_repair = "unique_quiet") |>
  dplyr::rename(
    sample_short_id = Sample_ID,
    original_individual_id = Individual,
    old_vcf_id = VCF_ID,
    old_seq_id = `Seq ID`,
    sex_field_raw = Sex_observed,
    age = Age,
    tissue = Tissue,
    first_round_depth = `Coverage 1st round`,
    first_round_heterozygosity = het_prop
  ) |>
  dplyr::mutate(old_id_for_join = dplyr::coalesce(old_vcf_id, old_seq_id)) |>
  dplyr::left_join(rename_map, by = c("old_id_for_join" = "old_id"))

qc <- readr::read_tsv(
  qc_file,
  col_types = readr::cols(.default = readr::col_guess())
) |>
  dplyr::rename(
    sample_id = Sample,
    mean_depth = AutoDepth,
    missingness = AutoMissing,
    heterozygosity = AutoHet,
    qc_flags = Flags,
    sex_genomic_raw = Sex
  )

pixy_all <- readr::read_tsv(
  paths$pixy_metadata_dir |> file.path("population_all.txt"),
  col_names = c("sample_id", "pixy_group"),
  col_types = readr::cols(.default = readr::col_character())
)

pixy_regions <- readr::read_tsv(
  paths$pixy_metadata_dir |> file.path("population_regions.txt"),
  col_names = c("sample_id", "pixy_region"),
  col_types = readr::cols(.default = readr::col_character())
)

roh_indiv <- read.table(
  paths$roh_dir |> file.path("SnowyOwl_autosomes_noContam_ROH.hom.indiv"),
  header = TRUE
) |>
  tibble::as_tibble() |>
  dplyr::transmute(sample_id = IID, roh_segments = NSEG, roh_total_kb = KB)

pca_ids <- if (file.exists(paths$pca_eigenvec)) {
  readr::read_table(
    paths$pca_eigenvec,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    dplyr::transmute(sample_id = X2)
} else if (file.exists(file.path(old_dir, "pca", "autosomes_pca.eigenvec"))) {
  readr::read_table(
    file.path(old_dir, "pca", "autosomes_pca.eigenvec"),
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    dplyr::transmute(sample_id = X2)
} else {
  tibble::tibble(sample_id = character())
}

pca_filtered_ids <- if (file.exists(paths$pca_filtered_eigenvec)) {
  readr::read_table(
    paths$pca_filtered_eigenvec,
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    dplyr::transmute(old_id = X2) |>
    dplyr::left_join(rename_map, by = "old_id") |>
    dplyr::transmute(sample_id = dplyr::coalesce(sample_id, old_id))
} else if (file.exists(file.path(old_dir, "re", "autosomal_rm_ind.eigenvec"))) {
  readr::read_table(
    file.path(old_dir, "re", "autosomal_rm_ind.eigenvec"),
    col_names = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  ) |>
    dplyr::transmute(old_id = X2) |>
    dplyr::left_join(rename_map, by = "old_id") |>
    dplyr::transmute(sample_id = dplyr::coalesce(sample_id, old_id))
} else {
  tibble::tibble(sample_id = character())
}

sample_metadata <- snowy |>
  dplyr::left_join(qc, by = "sample_id") |>
  dplyr::mutate(
    region = dplyr::case_when(
      stringr::str_starts(sample_id, "FNM") ~ "FNM",
      stringr::str_starts(sample_id, "GRL") ~ "GRL",
      stringr::str_starts(sample_id, "NY") ~ "NYS",
      stringr::str_starts(sample_id, "SKW") ~ "SKW",
      stringr::str_starts(sample_id, "WRG") ~ "WRG",
      TRUE ~ NA_character_
    ),
    collection_year = extract_year(Date),
    sex_field = normalize_sex(sex_field_raw),
    sex_genomic = normalize_sex(sex_genomic_raw),
    sex_final = dplyr::if_else(sex_genomic != "Unknown", sex_genomic, sex_field),
    relatedness_flag = dplyr::case_when(
      sample_id %in% c("GRL04_merged_X14", "SKW02_merged_X22", "FNM11_merged_X4", "FNM15_merged_X8") ~ "remove",
      sample_id == "FNM12_5-X5_S119_L004" ~ "review",
      TRUE ~ "keep"
    ),
    contamination_flag = dplyr::case_when(
      sample_id %in% c("GRL04_merged_X14", "SKW02_merged_X22") ~ "remove",
      qc_flags == "HIGH_MISSING" ~ "review",
      TRUE ~ "keep"
    ),
    include_pca = sample_id %in% pca_ids$sample_id,
    include_pca_filtered = sample_id %in% pca_filtered_ids$sample_id,
    include_pixy = sample_id %in% pixy_all$sample_id,
    include_roh = sample_id %in% roh_indiv$sample_id,
    include_gone = include_pixy & mean_depth >= 25,
    library_id = old_id_for_join,
    individual_id = dplyr::coalesce(dplyr::na_if(as.character(original_individual_id), "NA"), sample_short_id),
    population_source_rule = "Assigned from stable sample-id prefix; NY0/NY1 combined as NYS.",
    notes = dplyr::case_when(
      sample_id %in% c("GRL04_merged_X14", "SKW02_merged_X22") ~
        "Excluded from no-contaminant/no-close-relative analyses.",
      sample_id %in% c("FNM11_merged_X4", "FNM15_merged_X8") ~
        "Close relatives flagged from historical relatedness notes.",
      sample_id == "FNM12_5-X5_S119_L004" ~
        "Relatedness review flag from historical relatedness notes.",
      qc_flags == "HIGH_MISSING" ~ "High missingness flag in genomic QC.",
      TRUE ~ ""
    )
  ) |>
  dplyr::left_join(location_lookup, by = "region") |>
  dplyr::select(
    sample_id,
    individual_id,
    library_id,
    sample_short_id,
    country = Country,
    location = Location,
    region,
    population_final,
    population_source_rule,
    latitude,
    longitude,
    collection_year,
    sex_field,
    sex_genomic,
    sex_final,
    age,
    tissue,
    mean_depth,
    missingness,
    heterozygosity,
    qc_flags,
    relatedness_flag,
    contamination_flag,
    include_pca,
    include_pca_filtered,
    include_pixy,
    include_roh,
    include_gone,
    notes
  ) |>
  dplyr::arrange(population_final, sample_short_id)

readr::write_csv(sample_metadata, paths$metadata, na = "")

summary_table <- sample_metadata |>
  dplyr::count(population_final, include_pixy, include_gone, name = "n") |>
  dplyr::arrange(population_final, include_pixy, include_gone)

dir.create(paths$results_tables, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(summary_table, file.path(paths$results_tables, "sample_metadata_summary.csv"))

message("Wrote ", paths$metadata, " with ", nrow(sample_metadata), " samples.")
message("GONE high-depth samples: ", sum(sample_metadata$include_gone, na.rm = TRUE))
