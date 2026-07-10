#!/usr/bin/env Rscript
# ==============================================================================
# 11_strobe_flowchart.R — Figure 1: STROBE-IPD Flow Diagram (Lancet style)
# ==============================================================================
# Creates a publication-quality STROBE participant flow diagram using grid graphics
# v2.0 — Lancet Public Health formatting (unified font size, softened palette,
#         simplified bottom layer boxes, increased height for readability)
# ==============================================================================

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
})

cat(strrep("=", 70), "\n")
cat("Figure 1: STROBE-IPD Participant Flow Diagram (Lancet style)\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# FLOW DIAGRAM NUMBERS (verified from data pipeline)
# ==============================================================================

n_cohorts       <- 5
n_total_persons <- "174,695"
n_total_waves   <- "526,455"

n_excl_adl_tot  <- "28,134"
n_excl_adl      <- "26,462"
n_excl_adl_na   <- "1,672"

n_eligible      <- "146,561"
n_excl_nofu     <- "44,555"

n_analysis      <- "102,006"
n_events_dfs    <- "26,758"
n_events_mort   <- "2,644"

n_ic_env        <- "88,409"
n_ic_only       <- "101,114"
n_env_only      <- "89,139"

# ==============================================================================
# DRAW FLOWCHART
# ==============================================================================

draw_flowchart <- function() {

  # --- Unified font sizes (Lancet style) ---
  cex_main      <- 0.58
  cex_exclusion <- 0.52
  cex_layer     <- 0.48

  # --- Softened color palette ---
  fill_source    <- "#E8F0FE"   # Light blue
  fill_main      <- "white"
  fill_analysis  <- "#E8F5E9"   # Light green
  fill_exclusion <- "#FFF8F0"   # Very pale orange (softened from #FFF3E0)
  fill_cohort    <- "#F3E5F5"   # Light purple
  fill_layers    <- c("#E3F2FD", "#FFF8E1", "#FCE4EC", "#E8EAF6")

  # --- Helper functions ---
  draw_box <- function(x, y, w, h, label, fill = "white", border = "grey40",
                       cex = cex_main, font = 1, line_height = 1.2) {
    grid.roundrect(x = x, y = y, width = w, height = h,
                   r = unit(3, "mm"),
                   gp = gpar(fill = fill, col = border, lwd = 1))
    lines <- strsplit(label, "\n")[[1]]
    n <- length(lines)
    for (i in seq_along(lines)) {
      offset <- (n/2 - i + 0.5) * line_height * cex * 0.012
      grid.text(lines[i], x = x, y = y + offset,
                gp = gpar(cex = cex, font = font, col = "grey10",
                          fontfamily = "sans"))
    }
  }

  draw_arrow <- function(x0, y0, x1, y1) {
    grid.lines(x = c(x0, x1), y = c(y0, y1),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               gp = gpar(fill = "grey30", col = "grey30", lwd = 1))
  }

  draw_arrow_right <- function(x0, y0, x1, y1) {
    grid.lines(x = c(x0, x0, x1), y = c(y0, y1, y1),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               gp = gpar(fill = "grey30", col = "grey30", lwd = 1))
  }

  # --- Layout coordinates ---
  cx <- 0.38
  rx <- 0.74
  bw <- 0.42
  ew <- 0.36

  y1 <- 0.92
  y2 <- 0.77
  y3 <- 0.62
  y4 <- 0.44
  y5 <- 0.20

  ye1 <- (y2 + y3) / 2
  ye2 <- (y3 + y4) / 2

  # --- Box 1: Data Sources ---
  label1 <- paste0(
    "5 International Ageing Cohorts\n",
    "526,455 person-wave observations\n",
    "174,695 unique participants\n",
    "HRS (28,696) | ELSA (14,709) | CHARLS (15,264)\n",
    "SHARE (102,570) | MHAS (13,456)"
  )
  draw_box(cx, y1, bw + 0.06, 0.12, label1, fill = fill_source, cex = cex_main)

  draw_arrow(cx, y1 - 0.06, cx, y2 + 0.04)

  # --- Box 2: Baseline Sample ---
  label2 <- paste0(
    "Baseline sample\n",
    "(first wave per participant)\n",
    "N = 174,695"
  )
  draw_box(cx, y2, bw, 0.065, label2, fill = fill_main, cex = cex_main)

  draw_arrow(cx, y2 - 0.033, cx, y3 + 0.04)

  # --- Exclusion Box 1: ADL ---
  label_e1 <- paste0(
    "Excluded (n = 28,134)\n",
    "ADL \u2265 1 at baseline: 26,462\n",
    "ADL missing: 1,672"
  )
  draw_box(rx, ye1, ew, 0.065, label_e1, fill = fill_exclusion, cex = cex_exclusion)
  draw_arrow_right(cx + bw/2, ye1, rx - ew/2, ye1)

  # --- Box 3: Eligible ---
  label3 <- paste0(
    "Eligible participants\n",
    "(community-dwelling, age \u226560, baseline ADL = 0)\n",
    "N = 146,561"
  )
  draw_box(cx, y3, bw, 0.065, label3, fill = fill_main, cex = cex_main)

  draw_arrow(cx, y3 - 0.033, cx, y4 + 0.055)

  # --- Exclusion Box 2: No follow-up ---
  label_e2 <- paste0(
    "Excluded (n = 44,555)\n",
    "No follow-up wave: 44,550\n",
    "Zero follow-up time: 5"
  )
  draw_box(rx, ye2, ew, 0.065, label_e2, fill = fill_exclusion, cex = cex_exclusion)
  draw_arrow_right(cx + bw/2, ye2, rx - ew/2, ye2)

  # --- Box 4: Analysis Cohort ---
  label4 <- paste0(
    "Analysis cohort\n",
    "N = 102,006 participants\n",
    "DFS events: 26,758 (26.2%)\n",
    "Mortality events: 2,644 (2.6%)\n",
    "Median follow-up: 6.0 years"
  )
  draw_box(cx, y4, bw + 0.02, 0.10, label4, fill = fill_analysis, cex = cex_main)

  # --- Per-cohort breakdown ---
  cohort_label <- paste0(
    "Per-cohort:\n",
    "HRS: 20,550 | ELSA: 9,834\n",
    "CHARLS: 10,056 | SHARE: 53,869\n",
    "MHAS: 7,697"
  )
  draw_box(rx, y4, ew, 0.075, cohort_label, fill = fill_cohort, cex = cex_exclusion)
  grid.lines(x = c(cx + (bw+0.02)/2, rx - ew/2), y = c(y4, y4),
             gp = gpar(col = "grey30", lwd = 1, lty = 2))

  draw_arrow(cx, y4 - 0.05, cx, y5 + 0.085)

  # --- Box 5: Four Analysis Layers (simplified text) ---
  layer_w <- 0.21
  layer_h <- 0.14
  layer_y <- y5
  layer_x <- c(0.12, 0.35, 0.58, 0.81)

  layer_labels <- c(
    "Layer 1\nTwo-stage RE-MA\nN = 88,409 (interaction)",
    "Layer 2\nTrajectory (LCGA)\n4 classes \u00d7 5 cohorts",
    "Layer 3\nMSM / IPTW\nTime-varying adjustment",
    "Layer 4\n20 Sensitivity analyses\n19/20 significant"
  )

  for (i in 1:4) {
    draw_box(layer_x[i], layer_y, layer_w, layer_h,
             layer_labels[i], fill = fill_layers[i], cex = cex_layer)
  }

  for (i in 1:4) {
    grid.lines(x = c(cx, layer_x[i]),
               y = c(y5 + 0.085, layer_y + layer_h/2),
               gp = gpar(col = "grey50", lwd = 0.8, lty = 1))
  }

  # No in-figure title (Lancet: title goes in figure legend)
}

# ==============================================================================
# SAVE
# ==============================================================================

cat("Generating Figure 1...\n")

# PDF (preferred for Lancet submission)
pdf("figures/fig1_strobe_flow.pdf", width = 10, height = 10)
grid.newpage()
draw_flowchart()
dev.off()
cat("  Saved: figures/fig1_strobe_flow.pdf\n")

# TIFF (backup raster format)
tiff("figures/fig1_strobe_flow.tiff", width = 10, height = 10,
     units = "in", res = 300, compression = "lzw")
grid.newpage()
draw_flowchart()
dev.off()
cat("  Saved: figures/fig1_strobe_flow.tiff\n")

cat("\nDone!\n")
