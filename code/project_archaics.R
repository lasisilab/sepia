#!/usr/bin/env Rscript
# SEPIA — modern PCA + archaic projection from a plink --export A dosage matrix.
# Builds the PCA on the MODERN samples only, then projects the archaics onto it
# (no-mean-imputation style: a missing archaic genotype is set to the modern mean, so it
# contributes 0 to the centred score and doesn't drag the sample toward a wrong centroid).
# Usage: project_archaics.R <dose.raw> <archaic_names.txt> <out_prefix>
args <- commandArgs(TRUE)
raw <- read.table(args[1], header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("NA"))
arc_names <- readLines(args[2])
iid <- raw$IID
G <- as.matrix(raw[, -(1:6), drop = FALSE])           # .raw: FID IID PAT MAT SEX PHENOTYPE, then dosages
is_arc <- iid %in% arc_names
M <- G[!is_arc, , drop = FALSE]; A <- G[is_arc, , drop = FALSE]
mod_iid <- iid[!is_arc]; arc_iid <- iid[is_arc]
cat(sprintf("loaded %d modern + %d archaic; %d SNP columns\n", nrow(M), nrow(A), ncol(G)))

mu <- colMeans(M, na.rm = TRUE)
for (j in seq_len(ncol(M))) { na <- is.na(M[, j]); if (any(na)) M[na, j] <- mu[j] }
sdv <- apply(M, 2, sd); keep <- is.finite(sdv) & sdv > 0
M <- M[, keep, drop = FALSE]; A <- A[, keep, drop = FALSE]; mu <- mu[keep]
cat(sprintf("kept %d polymorphic SNPs\n", sum(keep)))
for (j in seq_len(ncol(A))) { na <- is.na(A[, j]); if (any(na)) A[na, j] <- mu[j] }

pr   <- prcomp(M, center = TRUE, scale. = TRUE)
proj <- predict(pr, A)
ve   <- round(100 * (pr$sdev^2 / sum(pr$sdev^2))[1:2], 1)
write.table(data.frame(IID = mod_iid, PC1 = pr$x[, 1], PC2 = pr$x[, 2]),
            paste0(args[3], "_modern.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(IID = arc_iid, PC1 = proj[, 1], PC2 = proj[, 2]),
            paste0(args[3], "_archaic.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: PC1 %.1f%%, PC2 %.1f%%; wrote _modern.tsv (%d) + _archaic.tsv (%d)\n",
            args[3], ve[1], ve[2], nrow(M), nrow(A)))
