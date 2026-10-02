# Run from new/. Import GATK results, recalculate ANGSD, and refresh PCA deliverables.
source("scripts/00_setup.R")
keep_file <- "data/external/pixy_curated50_noQUAL_DP5_noMAF/curated50_samples.txt"
keep <- readLines(keep_file)
stopifnot(length(keep) == 50, !anyDuplicated(keep))
base <- "data/external/pca/gatk/curated50_mac3/autosomes_curated50_mac3_pca"
g <- read.table(paste0(base, ".eigenvec"), stringsAsFactors = FALSE)
stopifnot(nrow(g) == 50, !anyDuplicated(g[[2]]), setequal(g[[2]], keep))
names(g) <- c("fid", "sample_id", paste0("PC", seq_len(ncol(g)-2)))
log <- readLines(paste0(base, ".log"))
hit <- regmatches(log, regexec("([0-9]+) variants and 50 people pass filters and QC", log))
hit <- Filter(function(x) length(x)>1, hit)
stopifnot(length(hit)==1)
n_snps <- as.integer(hit[[1]][2])
meta <- readr::read_csv(paths$metadata, show_col_types=FALSE)
stopifnot(all(keep %in% meta$sample_id))
meta$include_pca_filtered <- meta$sample_id %in% keep
readr::write_csv(meta, paths$metadata)
g <- dplyr::left_join(g, meta, by="sample_id")
ev <- scan(paste0(base, ".eigenval"), quiet=TRUE)
# PLINK only exported 20 components. Do not normalize a partial spectrum.
# Use all 50 eigenvalues if the optional full-spectrum rerun is downloaded.
full <- sub("_pca$", "_full_pca.eigenval", base)
if (file.exists(full)) {
  all_ev <- scan(full, quiet=TRUE)
  stopifnot(length(all_ev)==50, isTRUE(all.equal(ev, head(all_ev,length(ev)), tolerance=1e-4)))
  denominator <- sum(all_ev)
} else denominator <- if(length(ev)==50) sum(ev) else NA_real_
ge <- data.frame(pc=paste0("PC",seq_along(ev)),eigenvalue=ev,pve=100*ev/denominator,
                 variance_denominator=if(is.na(denominator)) "Unavailable: only 20 components exported" else "Full eigenvalue spectrum")
td <- file.path(paths$results_tables,"pca");dir.create(td,recursive=TRUE,showWarnings=FALSE)
readr::write_csv(g,file.path(td,"SnowyOwl_PCA_GATK_curated50_scores.csv"))
readr::write_csv(ge,file.path(td,"SnowyOwl_PCA_GATK_curated50_eigenvalues.csv"))
# Recompute locally using the exact pre-existing ANGSD thinned matrix.
status <- system2(file.path(R.home("bin"),"Rscript"),c("analysis/pca/slurm/curated50/angsd_curated50.R",
 "data/external/pca/angsd/snowyowl_autosomes_thinned.haplo",
 "data/external/pca/angsd/sample_names_pca.txt",keep_file,td))
stopifnot(status==0)
astub <- file.path(td,"SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_curated50_mac3")
a <- readr::read_csv(paste0(astub,"_scores.csv"),show_col_types=FALSE)
stopifnot(nrow(a)==50,setequal(a$sample_id,keep))
a <- dplyr::left_join(a,meta,by="sample_id")
readr::write_csv(a,paste0(astub,"_scores.csv"))
ae <- readr::read_csv(paste0(astub,"_eigenvalues.csv"),show_col_types=FALSE)
ze <- readr::read_csv(file.path(td,"SnowyOwl_PCA_ANGSD_pseudohaploid_Z_eigenvalues.csv"),show_col_types=FALSE)
we <- readr::read_csv(file.path(td,"SnowyOwl_PCA_ANGSD_pseudohaploid_W_females_only_eigenvalues.csv"),show_col_types=FALSE)
pve <- function(x,pc) x$pve[match(pc,x$pc)]
s4 <- data.frame(dataset=c("GATK autosomes, filtered","ANGSD pseudohaploid autosomes","ANGSD pseudohaploid Z","ANGSD pseudohaploid W"),
 input_type=c("diploid SNP genotypes",rep("pseudohaploid alleles",3)),
 sample_filter=c(rep("50 individuals retained after QC and relatedness filtering; one representative retained per closely related group.",2),"66 non-contaminated individuals","44 genetically female individuals"),
 site_filter=c("Existing Run A markers: MAC >=3 and LD pruning in 66 individuals","Biallelic sites; MAC >=3 recalculated in 50 individuals",rep("Biallelic sites; MAC >=1",2)),
 n_samples=c(50,50,66,44),n_sites_or_snps=c(n_snps,ae$retained_biallelic_sites[1],ze$retained_biallelic_sites[1],we$retained_biallelic_sites[1]),
 pc1_variance_percent=round(c(pve(ge,"PC1"),pve(ae,"PC1"),pve(ze,"PC1"),pve(we,"PC1")),3),
 pc2_variance_percent=round(c(pve(ge,"PC2"),pve(ae,"PC2"),pve(ze,"PC2"),pve(we,"PC2")),3),
 recommended_use=c("Figure S2A","Figure S2B","Figure S2C; diagnostic","Figure S2D; diagnostic"),
 notes=c(if(is.na(denominator)) "Total-variance percentages unavailable: only 20 eigenvalues exported; do not normalize by their partial sum." else "Percentages use full eigenvalue spectrum.",rep("Percentages use full eigenvalue spectrum.",3)))
s4 <- s4[, !names(s4) %in% c("recommended_use", "notes"), drop = FALSE]
for (dir in c("manuscript/supporting_tables","manuscript/supporting_tables/curated")) {
 readr::write_csv(s4,file.path(dir,"Table_S4_PCA_dataset_summaries.csv"),na="")
 f <- file.path(dir,"Table_S2_sequencing_and_genotype_QC.csv")
 if(file.exists(f)) {
  x <- readr::read_csv(f,show_col_types=FALSE)
  x$include_pca_filtered <- ifelse(x$sample_short_id %in% meta$sample_short_id[meta$include_pca_filtered],"T","F")
  # read_csv infers existing T/F columns as logical; restore compact labels.
  inclusion_cols <- grep("^include_", names(x), value = TRUE)
  x[inclusion_cols] <- lapply(x[inclusion_cols], function(v) {
    v <- as.character(v)
    v[v == "TRUE"] <- "T"
    v[v == "FALSE"] <- "F"
    v
  })
  x$sex_field[is.na(x$sex_field) | trimws(x$sex_field) == ""] <- "NA"
  readr::write_csv(x,f,na="")
 }
}
online <- "manuscript/online_data_tables/Online_Data_Table_1_full_sample_metadata_analysis_inclusion.csv"
x <- readr::read_csv(online,show_col_types=FALSE);x$include_pca_filtered <- x$sample_id %in% keep
readr::write_csv(x,online)
source("scripts/02c_pca_supporting_figure.R")
message("Updated curated50 PCA figure, metadata, Tables S2/S4 and Online Data Table 1.")
print(s4)
