#!/usr/bin/env Rscript
# =====================================================================================
# SEPIA B3 — wider pigmentation PCA figure: 222-panel vs the extended set (all vs high-conf).
# Three joint PCAs of the same 159 SGDP + 5 high-cov archaics, at correct hg19 coordinates:
#   (1) 222-panel (B2a)            output/pca_pig/pig_panel_pca.*
#   (2) extended, all ~923 SNPs    output/pca_wider/pig_wider_pca.*
#   (3) extended, high-confidence  output/pca_wider/pig_wider_highconf_pca.*
# Regions get a distinct shape as well as colour (shared _pca_style.R). Run from repo root.
# =====================================================================================
suppressPackageStartupMessages(library(ggplot2))
setwd(if (dir.exists("data")) "." else stop("run from repo root"))
source("plan/figures/_pca_style.R")

ARCHAICS <- c("Altai","Vindija33.19","Denisova3","Denisova25","Chagyrskaya8","Mez1")
SHORT <- c(Altai="Altai", Vindija33.19="Vindija", Denisova3="Den3",
           Denisova25="Den25", Chagyrskaya8="Chag", Mez1="Mez1")
meta <- read.delim("data/sgdp_metadata.tsv", stringsAsFactors = FALSE)
reg  <- setNames(meta$region, meta$IID)

panels <- list(
  list(f = "output/pca_pig/pig_panel_pca",           lab = "222-panel · %d SNPs"),
  list(f = "output/pca_wider/pig_wider_pca",          lab = "Extended (all) · %d SNPs"),
  list(f = "output/pca_wider/pig_wider_highconf_pca", lab = "Extended (high-conf) · %d SNPs")
)
# per-panel SNP counts (from the .bim written alongside, else the run notes / eigenval length)
nsnp <- c(160, 923, 350)   # 222-panel biallelic (B2a), wider all, wider high-conf

acc <- list()
for (i in seq_along(panels)) {
  p <- panels[[i]]
  ev <- read.table(paste0(p$f, ".eigenvec"), header = FALSE, skip = 1, stringsAsFactors = FALSE)
  val <- scan(paste0(p$f, ".eigenval"), quiet = TRUE)
  d <- data.frame(IID = ev[[2]], PC1 = ev[[3]], PC2 = ev[[4]], stringsAsFactors = FALSE)
  d$type   <- ifelse(d$IID %in% ARCHAICS, "archaic", "modern")
  d$region <- ifelse(d$type == "archaic", NA, reg[d$IID])
  d$region[is.na(d$region) & d$type == "modern"] <- "unknown"
  d$label  <- ifelse(d$type == "archaic", SHORT[d$IID], NA)
  afr <- d$region == "Africa" & !is.na(d$region)          # orient Africa to the right
  if (sum(afr) && mean(d$PC1[afr]) < 0) d$PC1 <- -d$PC1
  d$panel <- sprintf(paste0(p$lab, " (PC1 %d%%)"), nsnp[i], round(100*val[1]/sum(val)))
  acc[[i]] <- d
}
all <- do.call(rbind, acc)
all$panel <- factor(all$panel, levels = sapply(acc, function(d) d$panel[1]))
mod <- subset(all, type == "modern"); mod$region <- region_factor(mod$region)
arc <- subset(all, type == "archaic")

p <- ggplot(mapping = aes(PC1, PC2)) +
  geom_point(data = mod, aes(colour = region, shape = region), size = 2.3, stroke = 0.9, alpha = 0.9) +
  geom_point(data = arc, shape = 25, size = 3, fill = "#212529", colour = "white", stroke = 0.4) +
  geom_text(data = arc, aes(label = label), size = 2.2, vjust = -0.9, colour = "#212529") +
  facet_wrap(~panel, scales = "free", nrow = 1) +
  region_scales() +
  labs(title = "Pigmentation PCA — Lily's 222-panel vs the extended SNP set",
       subtitle = "Same 159 SGDP + 5 high-cov archaics (◣), at correct hg19 coordinates. Regions coded by colour AND shape.") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "#dee2e6", fill = NA),
        strip.text = element_text(face = "bold", size = 9.5), plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, colour = "#555"), legend.position = "right")

ggsave("plan/figures/wider_pca.png", p, width = 13, height = 4.6, dpi = 150, bg = "white")
cat("wrote plan/figures/wider_pca.png\n")
