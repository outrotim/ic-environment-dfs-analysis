#!/usr/bin/env Rscript
# ==============================================================================
# 11_strobe_flowchart.R — Figure 1: STROBE-IPD Flow Diagram
# ==============================================================================
# Creates a publication-quality STROBE participant flow diagram using grid graphics
# Numbers derived from data pipeline (01–06)
# ==============================================================================

suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
})

cat(strrep("=", 70), "\n")
cat("Figure 1: STROBE-IPD Participant Flow Diagram\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# FLOW DIAGRAM NUMBERS (verified from data pipeline)
# ==============================================================================

# Stage 1: Initial data
n_cohorts       <- 5
cohort_detail   <- c(
  "HRS (USA): 28,696",
  "ELSA (UK): 14,709",
  "CHARLS (China): 15,264",
  "SHARE (Europe): 102,570",
  "MHAS (Mexico): 13,456"
)
n_total_persons <- "174,695"
n_total_waves   <- "526,455"

# Stage 2: Exclusion — ADL >= 1 at baseline
n_excl_adl      <- "26,462"
n_excl_adl_na   <- "1,672"
n_excl_adl_tot  <- "28,134"

# Stage 3: Eligible (ADL = 0)
n_eligible      <- "146,561"

# Stage 4: Exclusion — no follow-up
n_excl_nofu     <- "44,555"

# Stage 5: Analysis cohort
n_analysis      <- "102,006"
n_events_dfs    <- "26,758"
n_events_mort   <- "2,644"

# Stage 6: Subsets for specific analyses
n_ic_env        <- "88,409"   # Both IC + Env available (interaction model)
n_ic_only       <- "101,114"  # IC available
n_env_only      <- "89,139"   # Env available

# Per-cohort analysis sample
cohort_analysis <- c(
  "HRS: 20,550",
  "ELSA: 9,834",
  "CHARLS: 10,056",
  "SHARE: 53,869",
  "MHAS: 7,697"
)

# ==============================================================================
# DRAW FLOWCHART
# ==============================================================================

draw_flowchart <- function() {

  # --- Helper functions ---
  draw_box <- function(x, y, w, h, label, fill = "white", border = "grey30",
                       cex = 0.65, font = 1, line_height = 1.2) {
    grid.roundrect(x = x, y = y, width = w, height = h,
                   r = unit(2, "mm"),
                   gp = gpar(fill = fill, col = border, lwd = 1.2))
    lines <- strsplit(label, "\n")[[1]]
    n <- length(lines)
    for (i in seq_along(lines)) {
      offset <- (n/2 - i + 0.5) * line_height * cex * 0.012
      grid.text(lines[i], x = x, y = y + offset,
                gp = gpar(cex = cex, font = font, col = "grey10"))
    }
  }

  draw_arrow <- function(x0, y0, x1, y1) {
    grid.lines(x = c(x0, x1), y = c(y0, y1),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               gp = gpar(fill = "grey30", col = "grey30", lwd = 1))
  }

  draw_arrow_right <- function(x0, y0, x1, y1) {
    # L-shaped arrow: down then right
    grid.lines(x = c(x0, x0, x1), y = c(y0, y1, y1),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               gp = gpar(fill = "grey30", col = "grey30", lwd = 1))
  }

  # --- Layout coordinates ---
  cx <- 0.38   # center x of main flow
  rx <- 0.74   # x of exclusion boxes (right side)
  bw <- 0.42   # main box width
  ew <- 0.36   # exclusion box width

  # Y positions (top to bottom)
  y1 <- 0.92   # Data sources
  y2 <- 0.77   # Baseline persons
  y3 <- 0.62   # Eligible
  y4 <- 0.44   # Analysis cohort
  y5 <- 0.22   # Analysis layers

  # Exclusion box y-positions (between main boxes)
  ye1 <- (y2 + y3) / 2    # ADL exclusion
  ye2 <- (y3 + y4) / 2    # No follow-up exclusion

  # --- Box 1: Data Sources ---
  label1 <- paste0(
    "5 International Aging Cohorts\n",
    "526,455 person-wave observations\n",
    "174,695 unique participants\n",
    "HRS (28,696) | ELSA (14,709) | CHARLS (15,264)\n",
    "SHARE (102,570) | MHAS (13,456)"
  )
  draw_box(cx, y1, bw + 0.06, 0.12, label1, fill = "#E8F0FE", cex = 0.58)

  # Arrow 1→2
  draw_arrow(cx, y1 - 0.06, cx, y2 + 0.04)

  # --- Box 2: Baseline Sample ---
  label2 <- paste0(
    "Baseline sample\n",
    "(first wave per participant)\n",
    "N = 174,695"
  )
  draw_box(cx, y2, bw, 0.065, label2, cex = 0.6)

  # Arrow 2→3
  draw_arrow(cx, y2 - 0.033, cx, y3 + 0.04)

  # --- Exclusion Box 1: ADL ---
  label_e1 <- paste0(
    "Excluded (n = 28,134)\n",
    "- ADL >= 1 at baseline: 26,462\n",
    "- ADL missing: 1,672"
  )
  draw_box(rx, ye1, ew, 0.065, label_e1, fill = "#FFF3E0", cex = 0.55)
  draw_arrow_right(cx + bw/2, ye1, rx - ew/2, ye1)

  # --- Box 3: Eligible ---
  label3 <- paste0(
    "Eligible participants\n",
    "(community-dwelling, age >=60, baseline ADL = 0)\n",
    "N = 146,561"
  )
  draw_box(cx, y3, bw, 0.065, label3, cex = 0.6)

  # Arrow 3→4
  draw_arrow(cx, y3 - 0.033, cx, y4 + 0.055)

  # --- Exclusion Box 2: No follow-up ---
  label_e2 <- paste0(
    "Excluded (n = 44,555)\n",
    "- No follow-up wave: 44,550\n",
    "- Zero follow-up time: 5"
  )
  draw_box(rx, ye2, ew, 0.065, label_e2, fill = "#FFF3E0", cex = 0.55)
  draw_arrow_right(cx + bw/2, ye2, rx - ew/2, ye2)

  # --- Box 4: Analysis Cohort ---
  label4 <- paste0(
    "Analysis cohort\n",
    "N = 102,006 participants\n",
    "DFS events: 26,758 (26.2%)\n",
    "Mortality events: 2,644 (2.6%)\n",
    "Median follow-up: 6.0 years"
  )
  draw_box(cx, y4, bw + 0.02, 0.10, label4, fill = "#E8F5E9", cex = 0.6, font = 1)

  # --- Per-cohort breakdown ---
  cohort_label <- paste0(
    "Per-cohort:\n",
    "HRS: 20,550 | ELSA: 9,834\n",
    "CHARLS: 10,056 | SHARE: 53,869\n",
    "MHAS: 7,697"
  )
  draw_box(rx, y4, ew, 0.075, cohort_label, fill = "#F3E5F5", cex = 0.55)
  grid.lines(x = c(cx + (bw+0.02)/2, rx - ew/2), y = c(y4, y4),
             gp = gpar(col = "grey30", lwd = 1, lty = 2))

  # Arrow 4→5
  draw_arrow(cx, y4 - 0.05, cx, y5 + 0.085)

  # --- Box 5: Four Analysis Layers ---
  # Draw 4 sub-boxes side by side
  layer_w <- 0.21
  layer_h <- 0.14
  layer_y <- y5
  layer_x <- c(0.12, 0.35, 0.58, 0.81)
  layer_fills <- c("#E3F2FD", "#FFF8E1", "#FCE4EC", "#E8EAF6")

  layer_labels <- c(
    paste0("Layer 1: IPD-MA\n(Two-stage Cox + RE-MA)\n\n",
           "IC main: N = 101,114\n",
           "Env main: N = 89,139\n",
           "ICxEnv: N = 88,409"),
    paste0("Layer 2: Trajectory\n(LCGA, K=4)\n\n",
           "IC trajectories\n",
           "Env trajectories\n",
           "4-group DFS"),
    paste0("Layer 3: Causal\n(MSM / IPTW)\n\n",
           "Sustained Env effect\n",
           "ICxEnv causal\n",
           "IC-tertile stratified"),
    paste0("Layer 4: Sensitivity\n(17 pre-specified)\n\n",
           "5 domains (A-E)\n",
           "18/19 significant\n",
           "E-value = 1.47")
  )

  for (i in 1:4) {
    draw_box(layer_x[i], layer_y, layer_w, layer_h,
             layer_labels[i], fill = layer_fills[i], cex = 0.48)
  }

  # Arrows from center to each layer
  for (i in 1:4) {
    grid.lines(x = c(cx, layer_x[i]),
               y = c(y5 + 0.085, layer_y + layer_h/2),
               gp = gpar(col = "grey50", lwd = 0.8, lty = 1))
  }

  # --- Title ---
  grid.text("Figure 1. STROBE-IPD Participant Flow Diagram",
            x = 0.5, y = 0.99,
            gp = gpar(cex = 0.85, font = 2, col = "grey10"))
}

# ==============================================================================
# SAVE
# ==============================================================================

cat("Generating Figure 1...\n")

# PDF
pdf("figures/fig1_strobe_flow.pdf", width = 11, height = 9)
grid.newpage()
draw_flowchart()
dev.off()
cat("  Saved: figures/fig1_strobe_flow.pdf\n")

# PNG
png("figures/fig1_strobe_flow.png", width = 11, height = 9, units = "in", res = 300)
grid.newpage()
draw_flowchart()
dev.off()
cat("  Saved: figures/fig1_strobe_flow.png\n")

cat("\nDone!\n")
