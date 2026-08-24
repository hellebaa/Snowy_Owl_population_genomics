# Assemble supporting PCA figure from final GATK and ANGSD PCA outputs.

source("scripts/00_setup.R")

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install patchwork before assembling the PCA supporting figure.", call. = FALSE)
}

fig_dir <- file.path(paths$results_figures, "pca")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_eigen_pve <- function(file) {
  readr::read_csv(file, show_col_types = FALSE) |>
    dplyr::select(pc, pve)
}

axis_label <- function(pc, eigen) {
  pve <- eigen$pve[match(pc, eigen$pc)]
  paste0(pc, " (", round(pve, 1), "%)")
}

make_panel <- function(scores, eigen, title) {
  scores <- scores |>
    dplyr::mutate(population_final = factor(population_final, levels = names(pop_colors)))

  ggplot(scores, aes(x = PC1, y = PC2, color = population_final)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey82") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "grey82") +
    geom_point(size = 2.4, alpha = 0.9) +
    scale_color_manual(values = pop_colors, drop = FALSE, na.value = "grey45") +
    labs(
      title = title,
      x = axis_label("PC1", eigen),
      y = axis_label("PC2", eigen),
      color = "Region"
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 8, 5.5, 8)
    )
}

gatk_scores <- readr::read_csv(
  file.path(paths$results_tables, "pca", "SnowyOwl_PCA_GATK_regenerated_scores.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(run_id == "runA_mac3")

gatk_eigen <- readr::read_csv(
  file.path(paths$results_tables, "pca", "SnowyOwl_PCA_GATK_regenerated_eigenvalues.csv"),
  show_col_types = FALSE
) |>
  dplyr::filter(run_id == "runA_mac3") |>
  dplyr::select(pc, pve)

angsd_auto_scores <- readr::read_csv(
  file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_mac3_scores.csv"),
  show_col_types = FALSE
)
angsd_auto_eigen <- read_eigen_pve(
  file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_qc_filtered_mac3_eigenvalues.csv")
)

angsd_z_scores <- readr::read_csv(
  file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_Z_scores.csv"),
  show_col_types = FALSE
)
angsd_z_eigen <- read_eigen_pve(
  file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_Z_eigenvalues.csv")
)

angsd_w_scores <- readr::read_csv(
  file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_W_females_only_scores.csv"),
  show_col_types = FALSE
)
angsd_w_eigen <- read_eigen_pve(
  file.path(paths$results_tables, "pca", "SnowyOwl_PCA_ANGSD_pseudohaploid_W_females_only_eigenvalues.csv")
)

p_gatk <- make_panel(gatk_scores, gatk_eigen, "GATK autosomes")
p_angsd_auto <- make_panel(angsd_auto_scores, angsd_auto_eigen, "ANGSD autosomes")
p_z <- make_panel(angsd_z_scores, angsd_z_eigen, "ANGSD Z chromosome")
p_w <- make_panel(angsd_w_scores, angsd_w_eigen, "ANGSD W chromosome, females only")

combined <- (p_gatk + p_angsd_auto) / (p_z + p_w) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 13)
  )

ggsave(
  file.path(fig_dir, "SnowyOwl_PCA_supporting_GATK_ANGSD_Z_W.png"),
  combined,
  width = 8.8,
  height = 7.2,
  dpi = 350
)
ggsave(
  file.path(fig_dir, "SnowyOwl_PCA_supporting_GATK_ANGSD_Z_W.pdf"),
  combined,
  width = 8.8,
  height = 7.2
)

message("Wrote PCA supporting figure to ", fig_dir)
