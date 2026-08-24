# Recompute ANGSD W chromosome PCA using female samples only.

source("scripts/00_setup.R")

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table before running ANGSD W female-only PCA.", call. = FALSE)
}

fig_dir <- file.path(paths$results_figures, "pca")
tab_dir <- file.path(paths$results_tables, "pca")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

old_w_haplo <- file.path("..", "old", "pca", "PCA_INPUT", "snowyowl_W.haplo")
w_haplo <- if (file.exists(paths$angsd_w_haplo)) {
  paths$angsd_w_haplo
} else if (file.exists(old_w_haplo)) {
  old_w_haplo
} else {
  stop("Could not find W haplo input at ", paths$angsd_w_haplo, " or ", old_w_haplo, call. = FALSE)
}

sample_file <- paths$angsd_pca_sample_names
if (!file.exists(sample_file)) {
  stop("Could not find ANGSD sample-name file: ", sample_file, call. = FALSE)
}

metadata <- readr::read_csv(paths$metadata, show_col_types = FALSE) |>
  dplyr::mutate(population_final = factor(population_final, levels = names(pop_colors)))

samples <- readr::read_lines(sample_file)
female_samples <- metadata |>
  dplyr::filter(sex_final == "Female") |>
  dplyr::pull(sample_id)

keep_samples <- samples %in% female_samples
if (sum(keep_samples) < 3) {
  stop("Fewer than three female samples found in ANGSD sample-name file.", call. = FALSE)
}

haplo <- data.table::fread(
  w_haplo,
  header = FALSE,
  data.table = FALSE,
  showProgress = FALSE
)

if (ncol(haplo) != length(samples) + 3) {
  stop(
    "ANGSD W haplo sample count mismatch: ",
    ncol(haplo) - 3, " genotype columns but ", length(samples), " sample names.",
    call. = FALSE
  )
}

samples <- samples[keep_samples]
geno_char <- as.matrix(haplo[, -(1:3), drop = FALSE])[, keep_samples, drop = FALSE]
storage.mode(geno_char) <- "character"

recode_site <- function(alleles) {
  alleles[alleles == "N" | alleles == "0" | alleles == "." | alleles == ""] <- NA_character_
  counts <- sort(table(alleles, useNA = "no"), decreasing = TRUE)
  if (length(counts) < 2 || length(counts) > 2) {
    return(rep(NA_real_, length(alleles)))
  }

  major <- names(counts)[1]
  minor <- names(counts)[2]
  recoded <- rep(NA_real_, length(alleles))
  recoded[alleles == major] <- 0
  recoded[alleles == minor] <- 1
  recoded
}

geno <- t(apply(geno_char, 1, recode_site))
allele_counts <- rowSums(geno == 1, na.rm = TRUE)
called_counts <- rowSums(!is.na(geno))
minor_counts <- pmin(allele_counts, called_counts - allele_counts)
keep_sites <- called_counts >= 3 & minor_counts >= 1
geno <- geno[keep_sites, , drop = FALSE]

site_means <- rowMeans(geno, na.rm = TRUE)
centered <- sweep(geno, 1, site_means, "-")
centered[is.na(centered)] <- 0

sample_cov <- crossprod(centered) / (nrow(centered) - 1)
eig <- eigen(sample_cov, symmetric = TRUE)
pve <- eig$values / sum(eig$values) * 100

scores <- as.data.frame(eig$vectors[, seq_len(min(10, ncol(eig$vectors))), drop = FALSE])
names(scores) <- paste0("PC", seq_len(ncol(scores)))
scores <- scores |>
  dplyr::mutate(sample_id = samples, .before = 1) |>
  dplyr::left_join(metadata, by = "sample_id") |>
  dplyr::mutate(population_final = factor(population_final, levels = names(pop_colors)))

eigen_out <- tibble::tibble(
  pc = paste0("PC", seq_along(eig$values)),
  eigenvalue = eig$values,
  pve = pve,
  dataset = "ANGSD pseudohaploid W chromosome, females only",
  retained_samples = length(samples),
  min_mac = 1,
  retained_biallelic_sites = nrow(centered)
)

readr::write_csv(
  scores,
  file.path(tab_dir, "SnowyOwl_PCA_ANGSD_pseudohaploid_W_females_only_scores.csv")
)
readr::write_csv(
  eigen_out,
  file.path(tab_dir, "SnowyOwl_PCA_ANGSD_pseudohaploid_W_females_only_eigenvalues.csv")
)

p <- ggplot(scores, aes(x = PC1, y = PC2, color = population_final)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey80") +
  geom_point(size = 2.8, alpha = 0.9) +
  scale_color_manual(values = pop_colors, drop = FALSE, na.value = "grey45") +
  labs(
    x = paste0("PC1 (", round(pve[1], 1), "%)"),
    y = paste0("PC2 (", round(pve[2], 1), "%)"),
    color = "Region"
  ) +
  coord_equal() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave(file.path(fig_dir, "SnowyOwl_PCA_ANGSD_pseudohaploid_W_females_only_PC1_PC2.png"), p, width = 5.8, height = 4.6, dpi = 350)
ggsave(file.path(fig_dir, "SnowyOwl_PCA_ANGSD_pseudohaploid_W_females_only_PC1_PC2.pdf"), p, width = 5.8, height = 4.6)

message(
  "ANGSD W female-only PCA: retained ", length(samples),
  " female samples and ", nrow(centered), " biallelic polymorphic sites."
)
