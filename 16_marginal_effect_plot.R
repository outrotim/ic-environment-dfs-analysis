#!/usr/bin/env Rscript
# 16_marginal_effect_plot.R — Marginal Effect Plot: IC × Environment Interaction
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# Purpose: Visualize how the protective effect of environment varies across
#          the IC continuum. This is the "interaction made visible" figure.
#
# Method:
#   From M4: coxph(Surv(time_dfs, event_dfs) ~ ic_cfa_c + env_form_c +
#                  ic_env_cont + age + female + chronic_n, strata(cohort))
#   The marginal HR of Env at a given IC level x is:
#     log(HR_env | IC = x) = beta_env + beta_int * x
#     SE = sqrt(Var(beta_env) + x^2 * Var(beta_int) + 2*x*Cov(beta_env, beta_int))
#
# Outputs:
#   figures/fig3_marginal_effect.pdf   — NEW Figure 3
#   figures/fig3_marginal_effect.tiff
#   results/marginal_effect_data.csv

library(tidyverse)
library(survival)
library(metafor)
library(patchwork)

set.seed(2024)
cat("=== Marginal Effect Plot: Environment HR across IC continuum ===\n")

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# --- Lancet theme ---
base_family <- "sans"
lancet_colors <- c(
  blue = "#00468B", red = "#AD002A", green = "#42B540",
  purple = "#925E9F", orange = "#ED7D31", grey = "#7F7F7F",
  ltblue = "#0099B4", ltgrey = "#ADB6B6"
)
theme_lancet <- theme_minimal(base_size = 10, base_family = base_family) +
  theme(
    plot.title = element_blank(), plot.subtitle = element_blank(),
    axis.line.x = element_line(color = "grey30", linewidth = 0.4),
    axis.line.y = element_line(color = "grey30", linewidth = 0.4),
    axis.ticks = element_line(color = "grey30", linewidth = 0.3),
    axis.text = element_text(size = 8, color = "grey20"),
    axis.title = element_text(size = 9, color = "grey10"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    legend.background = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )
fig_width_2col <- 7.2
fig_dpi <- 300

# ==============================================================================
# STEP 1: Per-cohort Cox M4 models → extract beta_env, beta_int, vcov
# ==============================================================================

surv_base <- readRDS("data/survival_cohort.rds")

cohorts <- unique(surv_base$cohort)
cat("Cohorts:", paste(cohorts, collapse = ", "), "\n")

# Store per-cohort results
cohort_results <- list()

for (coh in cohorts) {
  dat <- surv_base %>%
    filter(cohort == coh & !is.na(ic_cfa_c) & !is.na(env_form_c) &
           !is.na(ic_env_cont) & time_dfs > 0)

  if (nrow(dat) < 100) {
    cat("  Skipping", coh, "— too few observations\n")
    next
  }

  fit <- coxph(Surv(time_dfs, event_dfs) ~ ic_cfa_c + env_form_c +
               ic_env_cont + age + female + chronic_n,
               data = dat)

  coefs <- coef(fit)
  vc <- vcov(fit)

  # Extract indices for env_form_c and ic_env_cont
  idx_env <- which(names(coefs) == "env_form_c")
  idx_int <- which(names(coefs) == "ic_env_cont")

  beta_env <- coefs[idx_env]
  beta_int <- coefs[idx_int]
  var_env  <- vc[idx_env, idx_env]
  var_int  <- vc[idx_int, idx_int]
  cov_ei   <- vc[idx_env, idx_int]

  cohort_results[[coh]] <- list(
    cohort = coh, n = nrow(dat),
    beta_env = beta_env, beta_int = beta_int,
    var_env = var_env, var_int = var_int, cov_ei = cov_ei
  )

  cat(sprintf("  %s: N=%d, beta_env=%.4f, beta_int=%.4f\n",
              coh, nrow(dat), beta_env, beta_int))
}

# ==============================================================================
# STEP 2: Pool using random-effects meta-analysis (two-stage approach)
# ==============================================================================

# We need to pool beta_env, beta_int, and their covariance
# Use univariate pooling for simplicity (consistent with main analysis)
cr <- bind_rows(lapply(cohort_results, function(x) {
  tibble(cohort = x$cohort, n = x$n,
         beta_env = x$beta_env, se_env = sqrt(x$var_env),
         beta_int = x$beta_int, se_int = sqrt(x$var_int),
         cov_ei = x$cov_ei)
}))

cat("\nPooling with REML random-effects...\n")

# Pool beta_env
pool_env <- rma(yi = cr$beta_env, sei = cr$se_env, method = "REML")
# Pool beta_int
pool_int <- rma(yi = cr$beta_int, sei = cr$se_int, method = "REML")

# For the covariance between pooled beta_env and pooled beta_int,
# use weighted average of within-study covariances (approximate)
w <- 1 / (cr$se_env^2 + pool_env$tau2)  # weights from env pooling
pooled_cov_ei <- sum(w * cr$cov_ei) / sum(w)

cat(sprintf("  Pooled beta_env: %.4f (SE=%.4f)\n", pool_env$beta, pool_env$se))
cat(sprintf("  Pooled beta_int: %.4f (SE=%.4f)\n", pool_int$beta, pool_int$se))
cat(sprintf("  Pooled cov(env,int): %.6f\n", pooled_cov_ei))

# ==============================================================================
# STEP 3: Compute marginal HR of Env across IC continuum
# ==============================================================================

# IC range: from -2 to +2 SD (covers virtually all data)
ic_grid <- seq(-2, 2, by = 0.05)

# Pooled marginal effect
pooled_marg <- tibble(
  ic = ic_grid,
  log_hr = as.numeric(pool_env$beta) + as.numeric(pool_int$beta) * ic_grid,
  se = sqrt(as.numeric(pool_env$se)^2 +
            ic_grid^2 * as.numeric(pool_int$se)^2 +
            2 * ic_grid * pooled_cov_ei),
  hr = exp(log_hr),
  hr_lo = exp(log_hr - 1.96 * se),
  hr_hi = exp(log_hr + 1.96 * se),
  source = "Pooled (RE)"
)

# Per-cohort marginal effects
cohort_marg <- bind_rows(lapply(cohort_results, function(x) {
  tibble(
    ic = ic_grid,
    log_hr = x$beta_env + x$beta_int * ic_grid,
    se = sqrt(x$var_env + ic_grid^2 * x$var_int + 2 * ic_grid * x$cov_ei),
    hr = exp(log_hr),
    hr_lo = exp(log_hr - 1.96 * se),
    hr_hi = exp(log_hr + 1.96 * se),
    source = x$cohort
  )
}))

# IC distribution for rug/histogram
ic_dist <- surv_base %>%
  filter(!is.na(ic_cfa_c)) %>%
  pull(ic_cfa_c)

cat(sprintf("  IC range in data: [%.2f, %.2f], mean=%.2f, SD=%.2f\n",
            min(ic_dist), max(ic_dist), mean(ic_dist), sd(ic_dist)))

# Save data
write_csv(bind_rows(pooled_marg, cohort_marg), "results/marginal_effect_data.csv")
cat("  Saved results/marginal_effect_data.csv\n")

# ==============================================================================
# STEP 4: Plot
# ==============================================================================

# Cohort colors (muted)
cohort_colors <- c(
  "CHARLS" = "#7FABD3", "ELSA" = "#B8D4A0", "HRS" = "#E8B88A",
  "MHAS" = "#D4A5C9", "SHARE" = "#A8C8C8"
)

# Panel A: Marginal effect plot
pA <- ggplot() +
  # Reference line at HR = 1
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  # Per-cohort lines (thin, semi-transparent)
  geom_line(data = cohort_marg,
            aes(x = ic, y = hr, color = source),
            linewidth = 0.4, alpha = 0.5) +
  # Pooled CI ribbon
  geom_ribbon(data = pooled_marg,
              aes(x = ic, ymin = hr_lo, ymax = hr_hi),
              fill = lancet_colors["blue"], alpha = 0.15) +
  # Pooled line
  geom_line(data = pooled_marg,
            aes(x = ic, y = hr),
            color = lancet_colors["blue"], linewidth = 1) +
  # Annotations: key IC values
  geom_point(data = pooled_marg %>% filter(ic %in% c(-1, 0, 1)),
             aes(x = ic, y = hr),
             shape = 21, size = 2.5, fill = lancet_colors["blue"],
             color = "white", stroke = 0.8) +
  geom_text(data = pooled_marg %>% filter(ic %in% c(-1, 0, 1)) %>%
              mutate(label = sprintf("%.2f", hr)),
            aes(x = ic, y = hr, label = label),
            vjust = -1.5, size = 2.8, family = base_family,
            color = lancet_colors["blue"]) +
  # Scales
  scale_color_manual(values = cohort_colors, name = "Cohort") +
  scale_x_continuous(breaks = seq(-2, 2, 0.5),
                     labels = function(x) ifelse(x == 0, "Mean", as.character(x))) +
  scale_y_continuous(breaks = seq(0.5, 1.2, 0.1)) +
  coord_cartesian(ylim = c(0.55, 1.15), xlim = c(-2, 2)) +
  labs(x = "Intrinsic capacity (SD units from cohort mean)",
       y = "HR per 1-SD increase in environment support") +
  theme_lancet +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7)
  )

# Panel B: IC distribution (density rug at bottom)
ic_dens <- tibble(ic = ic_dist) %>% filter(ic >= -2 & ic <= 2)

pB <- ggplot(ic_dens, aes(x = ic)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = lancet_colors["ltgrey"], color = NA) +
  scale_x_continuous(breaks = seq(-2, 2, 0.5), limits = c(-2, 2),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Density") +
  theme_lancet +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    plot.margin = margin(0, 15, 0, 10)
  )

# Combine: main plot (top) + density (bottom)
p_combined <- pA / pB + plot_layout(heights = c(5, 1))

ggsave("figures/fig3_marginal_effect.pdf", p_combined,
       width = fig_width_2col, height = 5.5, device = cairo_pdf)
ggsave("figures/fig3_marginal_effect.tiff", p_combined,
       width = fig_width_2col, height = 5.5, dpi = fig_dpi, compression = "lzw")
cat("\n  Saved figures/fig3_marginal_effect.pdf/tiff\n")

# ==============================================================================
# STEP 5: Key summary statistics
# ==============================================================================

cat("\n=== Key values for manuscript ===\n")
# HR of Env at IC = -1 SD (low IC)
v_lo <- pooled_marg %>% filter(abs(ic - (-1)) < 0.01)
cat(sprintf("  IC = -1 SD: HR = %.2f (%.2f-%.2f)\n", v_lo$hr, v_lo$hr_lo, v_lo$hr_hi))

# HR of Env at IC = 0 (mean)
v_mean <- pooled_marg %>% filter(abs(ic) < 0.01)
cat(sprintf("  IC = mean:  HR = %.2f (%.2f-%.2f)\n", v_mean$hr, v_mean$hr_lo, v_mean$hr_hi))

# HR of Env at IC = +1 SD (high IC)
v_hi <- pooled_marg %>% filter(abs(ic - 1) < 0.01)
cat(sprintf("  IC = +1 SD: HR = %.2f (%.2f-%.2f)\n", v_hi$hr, v_hi$hr_lo, v_hi$hr_hi))

# At what IC level does HR cross 1.0? (environment no longer protective)
crossing <- pooled_marg %>% filter(hr >= 1) %>% slice(1)
if (nrow(crossing) > 0) {
  cat(sprintf("  HR crosses 1.0 at IC = %.2f SD\n", crossing$ic))
} else {
  cat("  HR remains < 1.0 across entire IC range\n")
}

cat("\n=== MARGINAL EFFECT PLOT COMPLETE ===\n")
