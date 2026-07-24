#!/usr/bin/env Rscript
# =====================================================================================
# SEPIA — does the archaic + modern distribution differ in PIGMENTATION space vs the
# WHOLE GENOME? Modern PCA with archaics PROJECTED on, for WG and pigmentation, before/after.
#
#   WG  before : data/sgdp.wg.pca.eigenvec (15 mod) + data/ancient.projected.sgdp.sscore (arc)
#   pig before : data/pigmentation_snps_pca.eigenvec (15 mod) + data/ancient.projected.pig.sscore
#   WG  after  : output/projections/wg_after_{modern,archaic}.tsv   (159 mod + projected arc)
#   pig after  : output/projections/pig_after_{modern,archaic}.tsv
# The "after" files are produced by the cluster projection job; until they exist, only the
# BEFORE row is drawn. Regions coded by colour AND shape (shared _pca_style.R). Run from repo root.
# =====================================================================================
suppressPackageStartupMessages(library(ggplot2))
setwd(if (dir.exists("data")) "." else stop("run from repo root"))
source("plan/figures/_pca_style.R")

meta <- read.delim("data/sgdp_metadata.tsv", stringsAsFactors = FALSE)
reg  <- setNames(meta$region, meta$IID)
# high-cov archaics we label (the rest are low-cov and pile at the coverage attractor)
HICOV <- c("Chagyrskaya8","Denisova3","Denisova25","Denisova5","Vi33.19",
           "Vindija33.19","Altai","Mez1")
SHORT <- c(Chagyrskaya8="Chag", Denisova3="Den3", Denisova25="Den25", Denisova5="Den5(Altai)",
           Vi33.19="Vindija", Vindija33.19="Vindija", Altai="Altai", Mez1="Mez1")

# --- readers -------------------------------------------------------------------------
read_eigenvec <- function(path, iid_col) {   # modern PCA sample scores
  d <- read.table(path, header = FALSE, skip = 1, stringsAsFactors = FALSE)
  data.frame(IID = d[[iid_col]], PC1 = d[[iid_col + 1]], PC2 = d[[iid_col + 2]],
             type = "modern", stringsAsFactors = FALSE)
}
read_sscore <- function(path) {               # archaic projected scores (plink --score .sscore, CSV)
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  data.frame(IID = d$IID, PC1 = d$PC1_AVG, PC2 = d$PC2_AVG, type = "archaic", stringsAsFactors = FALSE)
}
read_tsv_coords <- function(path, type) {     # our after-projection coords
  d <- read.delim(path, stringsAsFactors = FALSE)
  data.frame(IID = d$IID, PC1 = d$PC1, PC2 = d$PC2, type = type, stringsAsFactors = FALSE)
}

assemble <- function(mod, arc, panel) {
  d <- rbind(mod, arc)
  d$region <- ifelse(d$type == "archaic", NA, reg[d$IID])
  d$region[is.na(d$region) & d$type == "modern"] <- "unknown"
  d$label  <- ifelse(d$type == "archaic" & d$IID %in% HICOV, SHORT[d$IID], NA)
  afr <- d$region == "Africa" & !is.na(d$region)         # orient Africa to the right
  if (sum(afr) && mean(d$PC1[afr]) < 0) d$PC1 <- -d$PC1
  d$panel <- panel
  d
}

rows <- list()
rows[["wg_b"]]  <- assemble(read_eigenvec("data/sgdp.wg.pca.eigenvec", 1),
                            read_sscore("data/ancient.projected.sgdp.sscore"),
                            "WHOLE GENOME — before (Lily, 15 modern)")
rows[["pig_b"]] <- assemble(read_eigenvec("data/pigmentation_snps_pca.eigenvec", 2),
                            read_sscore("data/ancient.projected.pig.sscore"),
                            "PIGMENTATION — before (Lily, 15 modern)")

have_after <- all(file.exists(c("output/projections/wg_after_modern.tsv",
                                "output/projections/wg_after_archaic.tsv",
                                "output/projections/pig_after_modern.tsv",
                                "output/projections/pig_after_archaic.tsv")))
if (have_after) {
  rows[["wg_a"]]  <- assemble(read_tsv_coords("output/projections/wg_after_modern.tsv", "modern"),
                              read_tsv_coords("output/projections/wg_after_archaic.tsv", "archaic"),
                              "WHOLE GENOME — after (SEPIA, 159 modern)")
  rows[["pig_a"]] <- assemble(read_tsv_coords("output/projections/pig_after_modern.tsv", "modern"),
                              read_tsv_coords("output/projections/pig_after_archaic.tsv", "archaic"),
                              "PIGMENTATION — after (SEPIA, 159 modern)")
  order <- c(rows[["wg_b"]]$panel[1], rows[["wg_a"]]$panel[1], rows[["pig_b"]]$panel[1], rows[["pig_a"]]$panel[1])
  ncol <- 2; h <- 8
} else {
  order <- c(rows[["wg_b"]]$panel[1], rows[["pig_b"]]$panel[1]); ncol <- 2; h <- 4.4
  cat("NOTE: after-projection coords not found — drawing BEFORE row only.\n")
}

all <- do.call(rbind, rows)
all$panel <- factor(all$panel, levels = order)
mod <- subset(all, type == "modern"); mod$region <- region_factor(mod$region)
arc <- subset(all, type == "archaic")

p <- ggplot(mapping = aes(PC1, PC2)) +
  geom_point(data = mod, aes(colour = region, shape = region), size = 2.2, stroke = 0.9, alpha = 0.9) +
  geom_point(data = arc, shape = 25, size = 2.9, fill = "#212529", colour = "white", stroke = 0.4) +
  geom_text(data = subset(arc, !is.na(label)), aes(label = label), size = 2.1, vjust = -0.9, colour = "#212529") +
  facet_wrap(~panel, scales = "free", ncol = ncol) +
  region_scales() +
  labs(title = "Where do samples fall — pigmentation space vs the whole genome?",
       subtitle = "Modern PCA with archaics (▽) projected on. Regions by colour AND shape. PCA sign flipped so Africa is on the right.") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "#dee2e6", fill = NA),
        strip.text = element_text(face = "bold", size = 9.5), plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "#555"), legend.position = "right")

ggsave("plan/figures/wg_vs_pig_projection.png", p, width = 12, height = h, dpi = 150, bg = "white")
cat("wrote plan/figures/wg_vs_pig_projection.png  (", if (have_after) "2x2 before+after" else "before only", ")\n")
