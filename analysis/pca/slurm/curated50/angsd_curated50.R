# Base R only; same recoding, centering and covariance method as 02_pca.R.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
samples <- readLines(args[2]); keep <- readLines(args[3]); out <- args[4]
stopifnot(length(keep) == 50, !anyDuplicated(keep), !anyDuplicated(samples),
          all(keep %in% samples))
haplo <- read.table(args[1], header = FALSE, colClasses = "character",
                    comment.char = "", quote = "", check.names = FALSE)
stopifnot(ncol(haplo) == length(samples) + 3)
# Preserve the original matrix column order and write it with the outputs.
selected <- which(samples %in% keep)
samples <- samples[selected]
geno_char <- as.matrix(haplo[, selected + 3, drop = FALSE])
recode <- function(x) {
  x[x %in% c("N", "0", ".", "")] <- NA_character_
  counts <- sort(table(x, useNA = "no"), decreasing = TRUE)
  z <- rep(NA_real_, length(x))
  if (length(counts) == 2) {
    z[which(x == names(counts)[1])] <- 0
    z[which(x == names(counts)[2])] <- 1
  }
  z
}
geno <- t(apply(geno_char, 1, recode))
called <- rowSums(!is.na(geno)); ac <- rowSums(geno == 1, na.rm = TRUE)
valid <- called >= 3 & pmin(ac, called - ac) >= 3
stopifnot(sum(valid) > 1)
centered <- sweep(geno[valid, , drop = FALSE], 1,
                  rowMeans(geno[valid, , drop = FALSE], na.rm = TRUE), "-")
centered[is.na(centered)] <- 0
covariance <- crossprod(centered) / (nrow(centered) - 1)
eig <- eigen(covariance, symmetric = TRUE)
stopifnot(all(is.finite(eig$values)), sum(eig$values) > 0,
          identical(sort(samples), sort(keep)))
scores <- data.frame(sample_id = samples, eig$vectors[, 1:10, drop = FALSE])
names(scores)[-1] <- paste0("PC", 1:10)
eigenvalues <- data.frame(pc = paste0("PC", seq_along(eig$values)),
  eigenvalue = eig$values, pve = eig$values / sum(eig$values) * 100,
  dataset = "ANGSD pseudohaploid autosomes curated50", retained_samples = 50,
  min_mac = 3, retained_biallelic_sites = sum(valid))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stub <- file.path(out, "SnowyOwl_PCA_ANGSD_pseudohaploid_autosomes_curated50_mac3")
write.csv(scores, paste0(stub, "_scores.csv"), row.names = FALSE)
write.csv(eigenvalues, paste0(stub, "_eigenvalues.csv"), row.names = FALSE)
writeLines(samples, file.path(out, "actual_samples.txt"))
write.table(haplo[valid, 1:3], file.path(out, "retained_sites.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
message("Retained ", length(samples), " samples and ", sum(valid), " sites.")
