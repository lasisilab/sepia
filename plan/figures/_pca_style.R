# Shared PCA styling for SEPIA figures — regions get a distinct SHAPE as well as colour
# (redundant encoding, so regions are distinguishable even when colours are close), plus a
# consistent archaic marker. Sourced by make_sgdp_before_after.R, make_pig_before_after.R,
# make_wider_pca.R, and sgdp_reference_pca.qmd.
SEPIA_REGIONS <- c("Africa","WestEurasia","SouthAsia","CentralAsiaSiberia",
                   "EastAsia","Oceania","America","unknown")
SEPIA_COLORS  <- c(Africa="#E4572E", WestEurasia="#4C6EF5", SouthAsia="#9C36B5",
                   CentralAsiaSiberia="#0CA678", EastAsia="#F08C00", Oceania="#1098AD",
                   America="#E64980", unknown="#adb5bd")
# distinct shapes: ● ▲ ■ ◆ ✳ ⊠ ✕ ○  (all take `colour`, so colour+shape agree)
SEPIA_SHAPES  <- c(Africa=16, WestEurasia=17, SouthAsia=15, CentralAsiaSiberia=18,
                   EastAsia=8, Oceania=7, America=4, unknown=1)

# ggplot layers: map colour+shape to region with the palettes above.
region_scales <- function(name = "SGDP region") {
  list(
    ggplot2::scale_colour_manual(values = SEPIA_COLORS, name = name, drop = TRUE),
    ggplot2::scale_shape_manual(values = SEPIA_SHAPES, name = name, drop = TRUE)
  )
}
# put a region column onto the canonical level order (drops unused levels in the legend)
region_factor <- function(x) factor(x, levels = SEPIA_REGIONS[SEPIA_REGIONS %in% x])
