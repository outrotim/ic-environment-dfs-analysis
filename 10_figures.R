#!/usr/bin/env Rscript
# 10_figures.R — Publication Figures
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# Figures:
#   Fig 1: Study flow diagram (STROBE) — text-based, finalize in manuscript
#   Fig 2: Forest plot — IPD-MA main results (IC, Env, IC×Env per cohort + pooled)
#   Fig 3: Kaplan-Meier DFS by IC×Env 4-group
#   Fig 4: Sensitivity analysis forest plot (IC×Env interaction across all analyses)
#   Fig 5: MSM causal effect — Env→DFS stratified by IC tertile
#   eFig 1: LCGA trajectory classes (IC and Env mean trajectories)
#   eFig 2: IPTW weight distribution diagnostics
#
# Inputs:  results/*.csv, data/*.rds
# Outputs: figures/*.pdf, figures/*.png

library(tidyverse)
library(survival)
library(metafor)
library(gridExtra)

set.seed(2024)
t0 <- Sys.time()

# Create output directory
dir.create("figures", showWarnings = FALSE)

# Theme for publication
theme_pub <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# ============================================================
cat("=== Fig 2: IPD-MA Forest Plot ===\n")
# ============================================================

# Load stage 1 + 2 results from 06_ipd_ma.R
stage1 <- read_csv("results/stage1_estimates.csv", show_col_types = FALSE)
stage2 <- read_csv("results/stage2_pooled.csv", show_col_types = FALSE)

# Map variable names to display terms
# stage1 uses 'variable' and 'model'; stage2 uses 'variable' and 'model'
# Key variables: ic_cfa, env_formative, ic_cfa_c:env_form_c (from interaction model M3)
key_vars <- c("ic_cfa", "env_formative", "ic_cfa_c:env_form_c")
var_labels <- c("ic_cfa" = "Intrinsic Capacity",
                "env_formative" = "Environment Support",
                "ic_cfa_c:env_form_c" = "IC × Environment")

# Build forest data from actual column names
forest_data <- tibble()

for (vv in key_vars) {
  # Per-cohort estimates from stage1
  s1 <- stage1 %>% filter(variable == vv)
  if (nrow(s1) == 0) next

  for (i in seq_len(nrow(s1))) {
    forest_data <- bind_rows(forest_data, tibble(
      term = vv, label = s1$cohort[i], type = "cohort",
      hr = s1$hr[i], ci_lo = s1$hr_lower[i], ci_hi = s1$hr_upper[i],
      n = s1$n[i], weight = NA_real_
    ))
  }

  # Pooled from stage2
  s2 <- stage2 %>% filter(variable == vv)
  if (nrow(s2) > 0) {
    forest_data <- bind_rows(forest_data, tibble(
      term = vv, label = "Pooled (RE)", type = "pooled",
      hr = s2$pooled_hr[1], ci_lo = s2$hr_lower[1], ci_hi = s2$hr_upper[1],
      n = sum(s1$n), weight = NA_real_
    ))
  }
}

# Term labels
term_labels <- c(
  "ic_cfa_c" = "Intrinsic Capacity",
  "env_form_c" = "Environment Support",
  "ic_cfa_c:env_form_c" = "IC × Environment"
)

forest_data <- forest_data %>%
  mutate(
    term_label = var_labels[term],
    label = factor(label, levels = rev(c("CHARLS", "ELSA", "HRS", "MHAS", "SHARE", "Pooled (RE)"))),
    ci_label = sprintf("%.2f (%.2f–%.2f)", hr, ci_lo, ci_hi)
  )

p2 <- ggplot(forest_data, aes(x = hr, y = label, color = type, shape = type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2) +
  facet_wrap(~ term_label, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c("cohort" = "grey30", "pooled" = "firebrick"), guide = "none") +
  scale_shape_manual(values = c("cohort" = 16, "pooled" = 18), guide = "none") +
  labs(x = "Hazard Ratio (95% CI)", y = NULL,
       title = "Figure 2. Two-Stage IPD Meta-Analysis: IC, Environment, and IC×Env Interaction") +
  theme_pub +
  theme(strip.text = element_text(size = 10))

ggsave("figures/fig2_forest_ipdma.pdf", p2, width = 8, height = 9)
ggsave("figures/fig2_forest_ipdma.png", p2, width = 8, height = 9, dpi = 300)
cat("  Saved fig2_forest_ipdma.pdf/png\n")

# ============================================================
cat("\n=== Fig 3: KM Curves by IC×Env 4-group ===\n")
# ============================================================

surv_base <- readRDS("data/survival_cohort.rds")

# Create 4-group (IC tertile × Env group) — use ic_cat + env_cat
surv_km <- surv_base %>%
  filter(!is.na(ic_cat) & !is.na(env_cat) & time_dfs > 0) %>%
  mutate(
    ic_env_4grp = case_when(
      ic_cat == "Low" & env_cat == "Low"  ~ "Low IC / Low Env",
      ic_cat == "Low" & env_cat == "High" ~ "Low IC / High Env",
      ic_cat != "Low" & env_cat == "Low"  ~ "High IC / Low Env",
      TRUE                                 ~ "High IC / High Env"
    ),
    ic_env_4grp = factor(ic_env_4grp, levels = c(
      "High IC / High Env", "High IC / Low Env",
      "Low IC / High Env", "Low IC / Low Env"
    ))
  )

km_fit <- survfit(Surv(time_dfs, event_dfs) ~ ic_env_4grp, data = surv_km)

# Extract KM data manually from survfit object
n_strata <- length(km_fit$strata)
strata_names <- gsub("ic_env_4grp=", "", names(km_fit$strata))
strata_sizes <- km_fit$strata
strata_vec <- rep(strata_names, strata_sizes)

km_df <- tibble(
  time      = km_fit$time,
  estimate  = km_fit$surv,
  conf.low  = km_fit$lower,
  conf.high = km_fit$upper,
  strata    = strata_vec
)

colors_4grp <- c(
  "High IC / High Env" = "#2166AC",
  "High IC / Low Env"  = "#92C5DE",
  "Low IC / High Env"  = "#F4A582",
  "Low IC / Low Env"   = "#B2182B"
)

p3 <- ggplot(km_df, aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_step(linewidth = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  scale_color_manual(values = colors_4grp, name = "IC × Environment") +
  scale_fill_manual(values = colors_4grp, guide = "none") +
  scale_x_continuous(breaks = seq(0, 20, 2), limits = c(0, 18)) +
  scale_y_continuous(labels = scales::percent, limits = c(0.3, 1)) +
  labs(x = "Time (years)", y = "Disability-Free Survival",
       title = "Figure 3. Disability-Free Survival by IC–Environment Group") +
  theme_pub +
  theme(legend.position = c(0.25, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80"))

ggsave("figures/fig3_km_4group.pdf", p3, width = 8, height = 6)
ggsave("figures/fig3_km_4group.png", p3, width = 8, height = 6, dpi = 300)
cat("  Saved fig3_km_4group.pdf/png\n")

# ============================================================
cat("\n=== Fig 4: Sensitivity Analysis Forest Plot ===\n")
# ============================================================

sens <- read_csv("results/sensitivity_summary.csv", show_col_types = FALSE)

# Clean labels
sens_plot <- sens %>%
  filter(!is.na(hr)) %>%
  mutate(
    label = case_when(
      analysis == "S16_main_reference"  ~ "Main analysis (reference)",
      analysis == "S1_ic_sum"           ~ "S1: IC sum score",
      analysis == "S2_env_4indicator"   ~ "S2: Env 4-indicator",
      analysis == "S4_uncentered_interaction" ~ "S4: Uncentered interaction",
      analysis == "S5_mortality_only"   ~ "S5: Mortality only",
      analysis == "S6_disability_only"  ~ "S6: Disability only",
      analysis == "S7_adl_ge2_dfs"      ~ "S7: ADL ≥ 2 threshold",
      analysis == "S8_age_ge65"         ~ "S8: Age ≥ 65",
      analysis == "S9_age_ge70"         ~ "S9: Age ≥ 70",
      analysis == "S10_female_only"     ~ "S10: Female only",
      analysis == "S11_male_only"       ~ "S11: Male only",
      analysis == "S12_adj_education"   ~ "S12: + Education",
      analysis == "S13_adj_marital"     ~ "S13: + Marital status",
      analysis == "S17_ge3_waves"       ~ "S17: ≥ 3 waves only",
      grepl("LOO_excl", analysis)       ~ gsub("S15_LOO_excl_", "S15: Excl. ", analysis),
      TRUE ~ analysis
    ),
    sig = ifelse(p < 0.05, "Significant", "Not significant"),
    is_main = ifelse(analysis == "S16_main_reference", "main", "sensitivity")
  )

# Order
sens_plot$label <- factor(sens_plot$label, levels = rev(sens_plot$label))

p4 <- ggplot(sens_plot, aes(x = hr, y = label, color = sig, shape = is_main)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = sens_plot$hr[sens_plot$is_main == "main"][1],
             linetype = "dotted", color = "firebrick", alpha = 0.5) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25) +
  scale_color_manual(values = c("Significant" = "steelblue", "Not significant" = "grey60"),
                     name = "p < 0.05") +
  scale_shape_manual(values = c("main" = 18, "sensitivity" = 16), guide = "none") +
  labs(x = "Hazard Ratio (95% CI) — IC × Environment Interaction",
       y = NULL,
       title = "Figure 4. Sensitivity Analyses for IC × Environment Interaction") +
  theme_pub

ggsave("figures/fig4_sensitivity_forest.pdf", p4, width = 9, height = 8)
ggsave("figures/fig4_sensitivity_forest.png", p4, width = 9, height = 8, dpi = 300)
cat("  Saved fig4_sensitivity_forest.pdf/png\n")

# ============================================================
cat("\n=== Fig 5: MSM Env effect by IC tertile ===\n")
# ============================================================

msm_strat <- read_csv("results/msm_stratified.csv", show_col_types = FALSE)
msm_main <- read_csv("results/msm_main.csv", show_col_types = FALSE) %>%
  filter(analysis == "MSM_main_pooled_strat" & term == "cum_env_prop")

msm_fig <- bind_rows(
  msm_strat %>%
    filter(term == "cum_env_prop") %>%
    select(stratum, estimate, conf.low, conf.high, p.value, n, n_events) %>%
    mutate(label = case_when(
      stratum == "Low_IC"  ~ "Low IC tertile",
      stratum == "Mid_IC"  ~ "Mid IC tertile",
      stratum == "High_IC" ~ "High IC tertile"
    )),
  msm_main %>%
    mutate(label = "Overall (IPTW-weighted)", stratum = "Overall",
           n = NA_integer_, n_events = NA_integer_) %>%
    select(stratum, estimate, conf.low, conf.high, p.value, n, n_events, label)
)

msm_fig$label <- factor(msm_fig$label,
  levels = rev(c("Overall (IPTW-weighted)", "Low IC tertile",
                 "Mid IC tertile", "High IC tertile")))

p5 <- ggplot(msm_fig, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = stratum), size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = stratum),
                 height = 0.25) +
  scale_color_manual(values = c(
    "Overall" = "black", "Low_IC" = "#B2182B",
    "Mid_IC" = "#F4A582", "High_IC" = "#2166AC"
  ), guide = "none") +
  labs(x = "Hazard Ratio (95% CI) — Causal Effect of Sustained High Environment",
       y = NULL,
       title = "Figure 5. MSM: Causal Effect of Environment on DFS by Baseline IC",
       subtitle = "IPTW-weighted Cox models, stratified by cohort") +
  theme_pub

ggsave("figures/fig5_msm_stratified.pdf", p5, width = 8, height = 4)
ggsave("figures/fig5_msm_stratified.png", p5, width = 8, height = 4, dpi = 300)
cat("  Saved fig5_msm_stratified.pdf/png\n")

# ============================================================
cat("\n=== eFig 1: LCGA Trajectory Classes ===\n")
# ============================================================

traj_classes <- read_csv("results/lcga_class_trajectories.csv", show_col_types = FALSE)

# CSV format: class_label, time, mean_y, sd_y, n, cohort, variable
# variable = "ic_cfa" or "env_formative"
if (nrow(traj_classes) > 0 && "variable" %in% names(traj_classes)) {
  # Split by domain
  ic_traj  <- traj_classes %>% filter(variable == "ic_cfa")
  env_traj <- traj_classes %>% filter(variable == "env_formative")

  plot_traj <- function(dat, title_text, y_lab) {
    ggplot(dat, aes(x = time, y = mean_y, color = class_label, group = class_label)) +
      geom_line(linewidth = 0.8) +
      geom_ribbon(aes(ymin = mean_y - sd_y / sqrt(pmax(n, 1)),
                      ymax = mean_y + sd_y / sqrt(pmax(n, 1)),
                      fill = class_label),
                  alpha = 0.15, colour = NA) +
      geom_point(size = 1.5) +
      facet_wrap(~ cohort, scales = "free", nrow = 1) +
      labs(x = "Wave (time)", y = y_lab, title = title_text,
           color = "Trajectory Class", fill = "Trajectory Class") +
      theme_pub +
      theme(legend.position = "bottom")
  }

  if (nrow(ic_traj) > 0) {
    pe1a <- plot_traj(ic_traj, "eFigure 1A. LCGA Trajectory Classes — Intrinsic Capacity", "IC Score (CFA)")
    ggsave("figures/efig1a_lcga_ic.pdf", pe1a, width = 14, height = 4)
    ggsave("figures/efig1a_lcga_ic.png", pe1a, width = 14, height = 4, dpi = 300)
    cat("  Saved efig1a_lcga_ic.pdf/png\n")
  }

  if (nrow(env_traj) > 0) {
    pe1b <- plot_traj(env_traj, "eFigure 1B. LCGA Trajectory Classes — Environment Support", "Env Score (Formative)")
    ggsave("figures/efig1b_lcga_env.pdf", pe1b, width = 14, height = 4)
    ggsave("figures/efig1b_lcga_env.png", pe1b, width = 14, height = 4, dpi = 300)
    cat("  Saved efig1b_lcga_env.pdf/png\n")
  }
} else {
  cat("  Trajectory class data not in expected format — skipping eFig 1\n")
}

# ============================================================
cat("\n=== eFig 2: IPTW Weight Diagnostics ===\n")
# ============================================================

wt_diag <- read_csv("results/msm_weights_diagnostics.csv", show_col_types = FALSE)

pe2 <- ggplot(wt_diag, aes(x = cohort, y = wt_mean)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = wt_mean - wt_sd, ymax = wt_mean + wt_sd),
                width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "Cohort", y = "Mean IPTW (stabilized, truncated)",
       title = "eFigure 2. IPTW Weight Distribution by Cohort",
       subtitle = "Error bars = ±1 SD; dashed line = ideal mean of 1.0") +
  theme_pub

ggsave("figures/efig2_iptw_diagnostics.pdf", pe2, width = 6, height = 4)
ggsave("figures/efig2_iptw_diagnostics.png", pe2, width = 6, height = 4, dpi = 300)
cat("  Saved efig2_iptw_diagnostics.pdf/png\n")

# ============================================================
cat("\n=== Summary ===\n")
# ============================================================

figs <- list.files("figures", pattern = "\\.(pdf|png)$")
cat(sprintf("Generated %d figure files:\n", length(figs)))
for (f in sort(figs)) cat(sprintf("  figures/%s\n", f))

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat(sprintf("\n=== FIGURES COMPLETE (%.1f min) ===\n", elapsed))
