#!/usr/bin/env Rscript
# =====================================================================================
# SEPIA — corrected-coordinate pigmentation-panel heatmaps (Lily's #1/#2/#3, repaired).
# Substrate: output/panel_matrix/panel_hg19_{genotypes,depth}.tsv — 159 SGDP + 5 high-cov
# archaics at the CORRECT hg19 panel coordinates (code/export_panel_matrix.sh).
#
# Each ROW = one individual genome; rows grouped by region (159 SGDP) with the 5 high-cov
# archaics as their own block (their names are the only y labels shown). The x-axis is the
# panel loci ordered by chrom:pos, with white chromosome boundaries and labels that carry the
# SNP count per chromosome — so the uneven chromosome spans read as "more panel SNPs here"
# (chr5 = SLC45A2 region, chr15 = HERC2/OCA2/SLC24A5 cluster), not chromosome length.
#   (1) ATGC variability — homozygotes coloured by base, hets purple, missing black.
#   (2) read depth — log10(DP+1), missing black.
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

# ---- sample order: region (deep->shallow), archaics last ----
sdf <- data.table(sample=samp, region=factor(region_of(samp), levels=REGION_ORDER))
setorder(sdf, region, sample); samp_ord <- sdf$sample

# ---- locus order + x index + per-chromosome annotation (shared by both figures) ----
loci <- unique(geno[, .(chrom=factor(as.character(chrom), levels=CHR_ORDER), pos)])
setorder(loci, chrom, pos); loci[, xidx := .I]; loci[, key := paste0(chrom,":",pos)]
chrsum <- loci[, .(n=.N, mid=mean(xidx), maxx=max(xidx)), by=chrom]; setorder(chrsum, chrom)
bnds    <- head(chrsum$maxx, -1) + 0.5
xbreaks <- chrsum$mid
xlabs   <- sprintf("%s (%d)", as.character(chrsum$chrom), chrsum$n)
ylab_archaics <- function(v) ifelse(v %in% ARCH, v, "")   # only the 5 archaics get y labels

prep <- function(dt, valuecol, idv) {
  m <- melt(dt, id.vars=idv, measure.vars=samp, variable.name="sample", value.name=valuecol)
  m[, key := paste0(chrom,":",pos)]
  m <- merge(m, loci[, .(key, xidx)], by="key")
  m[, sample := factor(sample, levels=rev(samp_ord))]
  m[, region := factor(region_of(as.character(sample)), levels=REGION_ORDER)]
  m[]
}

base_theme <- theme_minimal(base_size=9) + theme(
  panel.grid=element_blank(), panel.spacing.y=unit(2,"pt"),
  strip.text.y.left=element_text(angle=0, hjust=1, size=7, face="bold", colour="#333"),
  axis.text.x=element_text(size=6, angle=90, vjust=0.5, hjust=1),
  axis.text.y=element_text(size=5, colour="#6a040f"), axis.ticks=element_blank(),
  plot.title=element_text(face="bold", size=13), plot.subtitle=element_text(size=8.5, colour="#555"),
  plot.caption=element_text(size=7.5, colour="#777", hjust=0), legend.position="right")

SUBTITLE <- "Each row = one individual genome. Rows grouped by region (159 SGDP); the 5 high-coverage archaics are the bottom block (only their names are labelled)."
CAPTION  <- "x-axis: 161 panel loci ordered by chrom:pos; white lines = chromosome boundaries; chromosome label shows the SNP count (chr5=SLC45A2, chr15=HERC2/OCA2/SLC24A5). Width reflects SNP count, not chromosome length."

xscale <- scale_x_continuous(breaks=xbreaks, labels=xlabs, expand=c(0,0))
yscale <- scale_y_discrete(labels=ylab_archaics)
facets <- facet_grid(rows=vars(region), scales="free_y", space="free_y", switch="y")
vlines <- geom_vline(xintercept=bnds, colour="white", linewidth=0.35)

# ======================= (1) ATGC variability =======================
gl <- prep(geno, "gt", c("chrom","pos","ref","alt"))
classify <- function(gt){ a <- tstrsplit(gt, "[/|]"); ifelse(a[[1]]==".", "missing", ifelse(a[[1]]==a[[2]], a[[1]], "het")) }
gl[, call := factor(classify(gt), levels=c("A","C","G","T","het","missing"))]
PAL <- c(A="#33a02c", C="#1f78b4", G="#ff7f00", T="#e31a1c", het="#6a3d9a", missing="#000000")
p1 <- ggplot(gl, aes(xidx, sample, fill=call)) + geom_tile() + vlines +
  scale_fill_manual(values=PAL, drop=FALSE, name="genotype") + xscale + yscale + facets +
  labs(title="Pigmentation panel — genotype variability (correct hg19 coordinates)",
       subtitle=SUBTITLE, caption=CAPTION, x="chromosome (panel-SNP count)", y=NULL) + base_theme
ggsave("plan/figures/panel_atgc_heatmap.png", p1, width=13, height=12, dpi=140, bg="white")
cat("wrote plan/figures/panel_atgc_heatmap.png\n")

# ======================= (2) read depth =======================
dl <- prep(dep, "dp", c("chrom","pos"))
dl[dp==".", dp := NA]; dl[, dp := as.numeric(dp)]; dl[, logdp := log10(dp+1)]
p2 <- ggplot(dl, aes(xidx, sample, fill=logdp)) + geom_tile() + vlines +
  scale_fill_viridis_c(option="magma", na.value="#000000", name="log10(DP+1)") + xscale + yscale + facets +
  labs(title="Pigmentation panel — read depth (correct hg19 coordinates)",
       subtitle=SUBTITLE, caption=CAPTION, x="chromosome (panel-SNP count)", y=NULL) + base_theme
ggsave("plan/figures/panel_depth_heatmap.png", p2, width=13, height=12, dpi=140, bg="white")
cat("wrote plan/figures/panel_depth_heatmap.png\n")

cat(sprintf("samples=%d loci=%d | missing genotypes=%.1f%%\n", length(samp), nrow(loci), 100*mean(gl$call=="missing")))
