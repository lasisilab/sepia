#!/usr/bin/env Rscript
# =====================================================================================
# SEPIA — where do the archaics fall in THREE spaces (all: modern PCA + archaics projected)?
#   whole genome  : output/projections/wg_after_{modern,archaic}.tsv
#   pigment SNPs  : output/projections/pig_after_{modern,archaic}.tsv        (ascertained panel)
#   pigment genes : output/pca_generegion/generegion_after_{modern,archaic}.tsv (whole gene regions)
# If gene-region and SNP-panel agree, the pigmentation signal isn't a panel artifact.
# Regions coloured AND shaped; archaics ggrepel-labelled, high-cov black / low-cov grey. Repo root.
# =====================================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })
setwd(if (dir.exists("data")) "." else stop("run from repo root"))
source("plan/figures/_pca_style.R")
meta <- read.delim("data/sgdp_metadata.tsv", stringsAsFactors = FALSE)
reg  <- setNames(meta$region, meta$IID)
HICOV <- c("Chagyrskaya8","Denisova3","Denisova25","Denisova5","Vi33.19","Vindija33.19","Altai","Mez1")
SHORT <- c(Chagyrskaya8="Chag", Denisova3="Den3", Denisova25="Den25", Denisova5="Den5/Altai",
           Vi33.19="Vindija", Vindija33.19="Vindija", Altai="Altai", Mez1="Mez1",
           Denisova11="Den11", Goyet="Goyet", Hohlenstein_Stadel="HST", Les_Cottes="LesC",
           Mezmaiskaya2="Mez2", El_Sidron="ElSidron", Spy="Spy", Vindija87="Vi87", L4741="L4741", SP4903="SP4903")
rd <- function(p) { d <- read.delim(p, stringsAsFactors = FALSE); data.frame(IID=d$IID, PC1=d$PC1, PC2=d$PC2) }
panel <- function(mdir, mfile, afile, lab) {
  m <- rd(file.path(mdir, mfile)); m$type <- "modern"
  a <- rd(file.path(mdir, afile)); a$type <- "archaic"
  d <- rbind(m, a); d <- d[is.finite(d$PC1) & is.finite(d$PC2), ]
  d$region <- ifelse(d$type=="archaic", NA, reg[d$IID]); d$region[is.na(d$region) & d$type=="modern"] <- "unknown"
  d$label  <- ifelse(d$type=="archaic", ifelse(d$IID %in% names(SHORT), SHORT[d$IID], d$IID), NA)
  d$hicov  <- d$type=="archaic" & d$IID %in% HICOV
  afr <- d$region=="Africa" & !is.na(d$region); if (sum(afr) && mean(d$PC1[afr])<0) d$PC1 <- -d$PC1
  d$panel <- lab; d
}
need <- c("output/projections/wg_after_modern.tsv","output/pca_generegion/generegion_after_modern.tsv")
if (!all(file.exists(need))) stop("missing coords: ", paste(need[!file.exists(need)], collapse=", "))

rows <- rbind(
  panel("output/projections", "wg_after_modern.tsv",  "wg_after_archaic.tsv",  "WHOLE GENOME (~480k SNPs)"),
  panel("output/projections", "pig_after_modern.tsv", "pig_after_archaic.tsv", "PIGMENTATION SNPs (ascertained panel)"),
  panel("output/pca_generegion", "generegion_after_modern.tsv", "generegion_after_archaic.tsv", "PIGMENTATION GENES (whole regions)"))
rows$panel <- factor(rows$panel, levels=c("WHOLE GENOME (~480k SNPs)","PIGMENTATION SNPs (ascertained panel)","PIGMENTATION GENES (whole regions)"))
mod <- subset(rows, type=="modern"); mod$region <- region_factor(mod$region)
arc <- subset(rows, type=="archaic" & hicov)   # high-cov only (low-cov/excluded scatter as artifacts, e.g. El Sidron)

p <- ggplot(mapping=aes(PC1,PC2)) +
  geom_point(data=mod, aes(colour=region, shape=region), size=2.1, stroke=0.9, alpha=0.9) +
  geom_point(data=subset(arc,!hicov), shape=25, size=2.1, fill="#adb5bd", colour="white", stroke=0.3) +
  geom_point(data=subset(arc, hicov), shape=25, size=3.0, fill="#212529", colour="white", stroke=0.4) +
  geom_text_repel(data=arc, aes(label=label), size=2.1, colour="#343a40",
                  max.overlaps=Inf, min.segment.length=0, segment.size=0.2, box.padding=0.35, seed=7) +
  facet_wrap(~panel, scales="free", nrow=1) + region_scales() +
  labs(title="Where do the archaics fall — genome-wide vs pigmentation SNPs vs pigmentation genes?",
       subtitle="Modern PCA (159 SGDP) + the 5 high-coverage archaics projected (▼). Africa oriented right. (Low-cov/excluded archaics omitted.)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.border=element_rect(colour="#dee2e6", fill=NA),
        strip.text=element_text(face="bold", size=9), plot.title=element_text(face="bold", size=13),
        plot.subtitle=element_text(size=9, colour="#555"), legend.position="right")
ggsave("plan/figures/three_way_projection.png", p, width=13.5, height=5, dpi=150, bg="white")
cat("wrote plan/figures/three_way_projection.png\n")
