#!/usr/bin/env Rscript
# =====================================================================================
# SEPIA — Pigmentation-panel PCA: the bug fix itself (A2), before vs after.
#
# BEFORE (committed, buggy — panel read at hg38 coords on hg19 data):
#   data/pigmentation_snps_pca.eigenvec   — 15 SGDP moderns, rank-1 (PC1 only; PC2-10 ~ 3e-18)
#   data/ancient.projected.pig.sscore     — every archaic projected to ONE point (0.0336)
# AFTER (B2a, panel read at CORRECT hg19 coords):
#   output/pca_pig/pig_panel_pca.eigenvec + .eigenval  — JOINT PCA of 159 SGDP + 6 high-cov archaics
#
# Unlike the SGDP WG figure, this IS the bug: the left panel is degenerate (the whole
# "archaics cluster with Africans / small number of loci" reading was a one-SNP artifact);
# the right panel is the same 222 SNPs read at the right coordinates, so real variation returns.
#
# Run from repo root AFTER B2a lands:  Rscript plan/figures/make_pig_before_after.R
# =====================================================================================
suppressPackageStartupMessages({ library(ggplot2) })
setwd(if (dir.exists("data")) "." else stop("run from repo root"))
source("plan/figures/_pca_style.R")

AFTER_EV <- "output/pca_pig/pig_panel_pca.eigenvec"
if (!file.exists(AFTER_EV)) stop("B2a output not present yet: ", AFTER_EV)

ARCHAICS <- c("Altai","Vindija33.19","Denisova3","Denisova25","Chagyrskaya8","Mez1")
SHORT <- c(Altai="Altai", Vindija33.19="Vindija", Denisova3="Den3",
           Denisova25="Den25", Chagyrskaya8="Chag", Mez1="Mez1")
meta <- read.table("data/sgdp_metadata.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
reg  <- setNames(meta$region, meta$IID)

# ---- BEFORE: Lily's degenerate modern PCA + collapsed archaic projection -------------
bm <- read.table("data/pigmentation_snps_pca.eigenvec", header = FALSE, skip = 1, stringsAsFactors = FALSE)
before_mod <- data.frame(IID = bm[[2]], PC1 = bm[[3]], PC2 = bm[[4]],
                         region = reg[bm[[2]]], type = "modern", label = NA, stringsAsFactors = FALSE)
ba <- read.csv("data/ancient.projected.pig.sscore", check.names = FALSE, stringsAsFactors = FALSE)
before_arc <- data.frame(IID = ba$IID, PC1 = ba$PC1_AVG, PC2 = ba$PC2_AVG,
                         region = "Archaic", type = "archaic", label = ba$IID, stringsAsFactors = FALSE)
before <- rbind(before_mod, before_arc)
before$panel <- "BEFORE — hg38 coords on hg19 data (rank-1)"

# ---- AFTER: B2a joint PCA (moderns by region + archaics highlighted) -----------------
am <- read.table(AFTER_EV, header = FALSE, skip = 1, stringsAsFactors = FALSE)
after <- data.frame(IID = am[[2]], PC1 = am[[3]], PC2 = am[[4]], stringsAsFactors = FALSE)
after$type   <- ifelse(after$IID %in% ARCHAICS, "archaic", "modern")
after$region <- ifelse(after$type == "archaic", "Archaic", reg[after$IID])
after$label  <- ifelse(after$type == "archaic", SHORT[after$IID], NA)
ev <- scan("output/pca_pig/pig_panel_pca.eigenval", quiet = TRUE)
after$panel <- sprintf("AFTER — correct hg19 coords (PC1 %.0f%%, PC2 %.0f%%)",
                       100*ev[1]/sum(ev), 100*ev[2]/sum(ev))

all <- rbind(before, after)
all$region[is.na(all$region)] <- "unknown"
all$panel <- factor(all$panel, levels = c(before$panel[1], after$panel[1]))
all$panel <- factor(all$panel, levels = c(before$panel[1], after$panel[1]))

mod <- subset(all, type == "modern"); arc <- subset(all, type == "archaic")
mod$region <- region_factor(ifelse(is.na(mod$region), "unknown", mod$region))
p <- ggplot(mapping = aes(PC1, PC2)) +
  geom_point(data = mod, aes(colour = region, shape = region), size = 2.5, stroke = .9, alpha = .95) +
  geom_point(data = arc, shape = 25, size = 3.3, fill = "#212529", colour = "white", stroke = .4) +
  geom_text(data = subset(arc, grepl("AFTER", panel)), aes(label = label),
            size = 2.4, vjust = -0.9, hjust = 0.5, colour = "#212529") +
  facet_wrap(~panel, scales = "free", nrow = 1) +
  region_scales("SGDP region") +
  labs(title = "Pigmentation-panel PCA — the build-mismatch bug (A2), before and after",
       subtitle = "Left: the panel read at hg38 positions on hg19 data — 221/222 SNPs monomorphic, everyone collapses (archaics ◆ on one point). Right: the same panel at the correct hg19 positions.") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "#dee2e6", fill = NA),
        strip.text = element_text(face = "bold", size = 10.5), plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "#555"), legend.position = "right")

out <- "plan/figures/pig_pca_before_after.png"
ggsave(out, p, width = 11.5, height = 5, dpi = 150, bg = "white")
cat("wrote", out, "\n")
cat("AFTER eigenvalues:", paste(round(ev,3), collapse=" "), "\n")
cat(sprintf("AFTER: %d modern + %d archaic\n", sum(after$type=="modern"), sum(after$type=="archaic")))
