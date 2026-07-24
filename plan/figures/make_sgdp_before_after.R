#!/usr/bin/env Rscript
# =====================================================================================
# SEPIA — SGDP whole-genome PCA: Lily's thesis reference vs the SEPIA rebuild.
#
# BEFORE: data/sgdp.wg.pca.eigenvec        — Lily's 15-sample SGDP WG PCA (thesis Fig 7 basis)
# AFTER : output/pca_wg/sgdp166_pca.eigenvec — SEPIA's 159-sample SGDP WG PCA (Sanity check 1)
#
# Both are healthy PCAs (the WG panel was never hit by the hg38/hg19 build bug); the point
# of this panel is the REFERENCE: 15 samples can't define global structure, 159 do. Colouring
# by SGDP region shows the SEPIA rebuild recovers the textbook continental layout that 15
# scattered samples only hint at. (The pigmentation-panel before/after — the bug fix itself —
# is the separate figure produced from the B2a run.)
#
# Run from the repo root:  Rscript plan/figures/make_sgdp_before_after.R
# =====================================================================================
suppressPackageStartupMessages({ library(ggplot2) })

repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
if (!dir.exists(file.path(repo, "data"))) repo <- getwd()          # fallback: run from repo root
setwd(repo)
source("plan/figures/_pca_style.R")

read_eig <- function(path, iid_col) {
  d <- read.table(path, header = FALSE, skip = 1, stringsAsFactors = FALSE)
  data.frame(IID = d[[iid_col]], PC1 = d[[iid_col + 1]], PC2 = d[[iid_col + 2]], stringsAsFactors = FALSE)
}
before <- read_eig("data/sgdp.wg.pca.eigenvec", 1)               # header: #IID PC1..
after  <- read_eig("output/pca_wg/sgdp166_pca.eigenvec", 2)      # header: #FID IID PC1..

meta <- read.table("data/sgdp_metadata.tsv", header = TRUE, sep = "\t",
                   quote = "", stringsAsFactors = FALSE)
reg  <- setNames(meta$region, meta$IID)
before$region <- reg[before$IID]; after$region <- reg[after$IID]
before$region[is.na(before$region)] <- "unknown"
after$region[is.na(after$region)]   <- "unknown"

# orient PC1 so Africa sits on the right in both panels (PCA sign is arbitrary)
orient <- function(df) {
  afr <- df$region == "Africa"
  if (sum(afr) > 0 && mean(df$PC1[afr], na.rm = TRUE) < 0) df$PC1 <- -df$PC1
  df
}
before <- orient(before); after <- orient(after)

before$panel <- sprintf("BEFORE — Lily's thesis reference (%d SGDP samples)", nrow(before))
after$panel  <- sprintf("AFTER — SEPIA rebuild (%d SGDP samples)", nrow(after))
all <- rbind(before, after)
all$panel <- factor(all$panel, levels = c(before$panel[1], after$panel[1]))

all$region <- region_factor(all$region)

p <- ggplot(all, aes(PC1, PC2, colour = region, shape = region)) +
  geom_point(size = 2.6, stroke = 0.9, alpha = 0.95) +
  facet_wrap(~panel, scales = "free", nrow = 1) +
  region_scales("SGDP region") +
  labs(title = "SGDP whole-genome PCA — the reference behind the pigmentation analysis",
       subtitle = "Same analysis, more samples: 15 scattered genomes → 159 recover the textbook continental structure (regions coded by colour AND shape)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "#dee2e6", fill = NA),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, colour = "#555"),
        legend.position = "right")

out <- "plan/figures/sgdp_pca_before_after.png"
ggsave(out, p, width = 11, height = 4.9, dpi = 150, bg = "white")
cat("wrote", out, "\n")
cat(sprintf("before: %d samples (%d with region)\n", nrow(before), sum(before$region != "unknown")))
cat(sprintf("after : %d samples (%d with region)\n", nrow(after),  sum(after$region  != "unknown")))
