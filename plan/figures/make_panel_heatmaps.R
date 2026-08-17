#!/usr/bin/env Rscript
# =====================================================================================
# SEPIA — corrected-coordinate pigmentation-panel heatmaps (Lily's #1/#2/#3, repaired).
# Substrate: output/panel_matrix/panel_hg19_{genotypes,depth}.tsv — 159 SGDP + 5 high-cov
# archaics at the CORRECT hg19 panel coordinates (code/export_panel_matrix.sh).
#   (1) ATGC variability heatmap — homozygotes coloured by base, hets purple, missing black.
#   (2) depth heatmap — log10 depth, faceted by chromosome ("break x-axis by chromosome").
# Run from repo root:  Rscript plan/figures/make_panel_heatmaps.R
# =====================================================================================
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
setwd(if (dir.exists("data")) "." else stop("run from repo root"))
DIR  <- "output/panel_matrix"
ARCH <- c("Vindija33.19","Mez1","Denisova3","Altai","Denisova25")

meta <- fread("data/sgdp_metadata.tsv"); reg <- setNames(meta$region, meta$IID)
region_of <- function(s) ifelse(s %in% ARCH, "ARCHAIC", ifelse(s %in% names(reg), reg[s], "unknown"))
REGION_ORDER <- c("Africa","WestEurasia","SouthAsia","CentralAsiaSiberia","EastAsia","America","Oceania","ARCHAIC")
CHR_ORDER <- c(as.character(1:22),"X","Y")

geno <- fread(file.path(DIR,"panel_hg19_genotypes.tsv"))
dep  <- fread(file.path(DIR,"panel_hg19_depth.tsv"))
samp <- setdiff(names(geno), c("chrom","pos","ref","alt"))

# ---- sample ordering: by region (deep-to-shallow), archaics last ----
sdf <- data.table(sample=samp, region=region_of(samp))
sdf[, region := factor(region, levels=REGION_ORDER)]; setorder(sdf, region, sample)
samp_ord <- sdf$sample
# ---- locus ordering: chrom then pos ----
geno[, chrom := factor(as.character(chrom), levels=CHR_ORDER)]; setorder(geno, chrom, pos)
geno[, locus := factor(paste0(chrom,":",pos), levels=paste0(chrom,":",pos))]

# ======================= (1) ATGC variability heatmap =======================
gl <- melt(geno, id.vars=c("chrom","pos","ref","alt","locus"), measure.vars=samp,
           variable.name="sample", value.name="gt")
classify <- function(gt){ a <- tstrsplit(gt, "[/|]"); a1 <- a[[1]]; a2 <- a[[2]]
  ifelse(a1==".", "missing", ifelse(a1==a2, a1, "het")) }
gl[, call := factor(classify(gt), levels=c("A","C","G","T","het","missing"))]
gl[, sample := factor(sample, levels=rev(samp_ord))]
PAL <- c(A="#33a02c", C="#1f78b4", G="#ff7f00", T="#e31a1c", het="#6a3d9a", missing="#000000")

p1 <- ggplot(gl, aes(locus, sample, fill=call)) + geom_tile() +
  scale_fill_manual(values=PAL, drop=FALSE, name="genotype") +
  labs(title="Pigmentation panel — genotype variability at CORRECT hg19 coords",
       subtitle="161 panel loci x 164 samples (159 SGDP + 5 high-cov archaics). Homozygotes coloured by base; het=purple; missing=black.",
       x="panel locus (chrom:pos, ordered)", y=NULL) +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=3.2),
        panel.grid=element_blank(), plot.title=element_text(face="bold", size=12),
        plot.subtitle=element_text(size=8, colour="#555"), legend.position="right")
ggsave("plan/figures/panel_atgc_heatmap.png", p1, width=13, height=11, dpi=140, bg="white")
cat("wrote plan/figures/panel_atgc_heatmap.png\n")

# ======================= (2) depth heatmap, faceted by chromosome =======================
dep[, chrom := factor(as.character(chrom), levels=CHR_ORDER)]; setorder(dep, chrom, pos)
dep[, locus := factor(paste0(chrom,":",pos), levels=paste0(chrom,":",pos))]
dl <- melt(dep, id.vars=c("chrom","pos","locus"), measure.vars=samp, variable.name="sample", value.name="dp")
dl[dp==".", dp := NA]; dl[, dp := as.numeric(dp)]
dl[, sample := factor(sample, levels=rev(samp_ord))]; dl[, logdp := log10(dp+1)]

p2 <- ggplot(dl, aes(locus, sample, fill=logdp)) + geom_tile() +
  scale_fill_viridis_c(option="magma", na.value="#000000", name="log10(DP+1)") +
  facet_grid(~chrom, scales="free_x", space="free_x", switch="x") +
  labs(title="Pigmentation panel — read depth at CORRECT hg19 coords, by chromosome",
       subtitle="Per-sample DP at 161 panel loci. Missing call = black. (Corrects the thesis figure: right coordinates, moderns+archaics.)",
       x="chromosome", y=NULL) +
  theme_minimal(base_size=9) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=3.2),
        panel.grid=element_blank(), panel.spacing=unit(1,"pt"), strip.text=element_text(size=6),
        plot.title=element_text(face="bold", size=12), plot.subtitle=element_text(size=8, colour="#555"),
        legend.position="right")
ggsave("plan/figures/panel_depth_heatmap.png", p2, width=14, height=11, dpi=140, bg="white")
cat("wrote plan/figures/panel_depth_heatmap.png\n")

cat(sprintf("samples=%d loci=%d | missing genotypes=%.1f%%\n", length(samp), nrow(geno), 100*mean(gl$call=="missing")))
