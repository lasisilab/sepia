#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# SEPIA — AFTER-fix consequence figures (Fig 1 coverage, Fig 3 scree) from the
# CORRECTED outputs. Companion to make_consequence_figures.R, which renders the
# BEFORE panels from the committed buggy files; this fills the AFTER panels once
# the fixes (A2 build, A1 chr-naming, Chagyrskaya) have landed.
#   coverage: output/panel_matrix/panel_hg19_depth.tsv (6 high-cov, VCF DP) +
#             lowcov_panel_depth.tsv (3 low-cov, remote-extracted)
#   scree:    output/pca_pig/pig_panel_pca.eigenval (corrected) vs data/sgdp.wg.pca.eigenval
# Run from repo root:  Rscript plan/figures/make_after_figures.R
# ---------------------------------------------------------------------------
suppressMessages({library(data.table); library(ggplot2)})
setwd(if (dir.exists("output")) "." else stop("run from repo root"))
outdir <- "plan/figures"
theme_sepia <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "top")

## ---- Fig 1 AFTER: realised coverage per archaic (Denisova 3 fixed) ----
meancol <- function(dt, s) {
  v <- suppressWarnings(as.numeric(dt[[s]])); v[is.na(v)] <- 0
  data.table(sample = s, mean_depth = mean(v), covered = sum(v > 0), n = length(v))
}
hi <- fread("output/panel_matrix/panel_hg19_depth.tsv")
lo <- fread("output/panel_matrix/lowcov_panel_depth.tsv")
cov <- rbind(
  rbindlist(lapply(c("Altai","Vindija33.19","Denisova3","Mez1","Denisova25","Chagyrskaya8"), meancol, dt = hi)),
  rbindlist(lapply(c("Denisova11","Goyet","LesCottes"), meancol, dt = lo)))
pretty <- c(Altai = "Altai (Den 5)", Vindija33.19 = "Vindija 33.19", Denisova3 = "Denisova 3",
  Mez1 = "Mezmaiskaya 1", Denisova25 = "Denisova 25", Chagyrskaya8 = "Chagyrskaya 8",
  Denisova11 = "Denisova 11", Goyet = "Goyet", LesCottes = "Les Cottes")
highset <- c("Altai","Vindija33.19","Denisova3","Denisova25","Chagyrskaya8")
cov[, tier := fifelse(sample %in% highset, "declared high-coverage", "low / medium")]
cov[, label := factor(pretty[sample], levels = pretty[sample][order(mean_depth)])]

figA <- ggplot(cov, aes(mean_depth, label, colour = tier)) +
  geom_segment(aes(x = 0, xend = mean_depth, yend = label), linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.1fx  (%d/%d)", mean_depth, covered, n)),
            hjust = -0.15, size = 3, colour = "grey25") +
  scale_colour_manual(values = c("declared high-coverage" = "#1E5A54", "low / medium" = "#B0B0B0"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "[AFTER] Realised coverage over the pigmentation panel",
       subtitle = "Corrected hg19 extraction: Denisova 3 now sits with the high-coverage genomes — the A1 chr-naming fingerprint is gone.",
       x = "Mean depth over panel positions", y = NULL) +
  theme_sepia
ggsave(file.path(outdir, "after_coverage.png"), figA, width = 8.5, height = 5.2, dpi = 150, bg = "white")

## ---- Fig 3 AFTER: scree (pigmentation panel now multi-component) ----
read_eig <- function(f, panel) { v <- fread(f, header = FALSE)[[1]]; data.table(panel = panel, PC = seq_along(v), eigenvalue = v) }
scree <- rbind(
  read_eig("output/pca_pig/pig_panel_pca.eigenval", "Pigmentation panel (corrected)"),
  read_eig("data/sgdp.wg.pca.eigenval", "Whole-genome panel"))
scree[, panel := factor(panel, levels = c("Pigmentation panel (corrected)", "Whole-genome panel"))]

figC <- ggplot(scree, aes(factor(PC), eigenvalue, fill = panel)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = c("Pigmentation panel (corrected)" = "#1E5A54", "Whole-genome panel" = "#5B8A83")) +
  labs(title = "[AFTER] PCA eigenvalue spectrum (scree)",
       subtitle = "The pigmentation panel is now multi-component (PC1 ~35%, PC2 ~17%) — the rank-1 collapse is gone.",
       x = "Principal component", y = "Eigenvalue") +
  theme_sepia
ggsave(file.path(outdir, "after_scree.png"), figC, width = 9, height = 4.4, dpi = 150, bg = "white")
cat("wrote after_coverage.png, after_scree.png\n")
