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
#   figures/fig3_marginal_effect.pdf
#   figures/fig3_marginal_effect.tiff
#   results/marginal_effect_data.csv

library(tidyverse)
library(survival)
library(metafor)
library(patchwork)

set.seed(2024)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Lancet-style colors
lancet_colors <- c(
  blue = "#00468B", red = "#AD002A", green = "#42B540",
  purple = "#925E9F", orange = "#ED7D31", grey = "#7F7F7F",
  ltblue = "#0099B4", ltgrey = "#ADB6B6"
)
theme_lancet <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# Per-cohort Cox M4 + extract coefficients
surv_base <- readRDS("data/survival_cohort.rds")
cohorts <- unique(surv_base$cohort)

cohort_results <- list()
for (coh in cohorts) {
  dat <- surv_base %>%
    filter(cohort == coh & !is.na(ic_cfa_c) & !is.na(env_form_c) &
             !is.na(ic_env_cont) & time_dfs > 0)
  if (nrow(dat) < 100) next
  fit <- coxph(Surv(time_dfs, event_dfs) ~ ic_cfa_c + env_form_c +
                 ic_env_cont + age + female + chronic_n, data = dat)
  coefs <- coef(fit); vc <- vcov(fit)
  idx_env <- which(names(coefs) == "env_form_c")
  idx_int <- which(names(coefs) == "ic_env_cont")
  cohort_results[[coh]] <- list(
    cohort = coh, n = nrow(dat),
    beta_env = coefs[idx_env], beta_int = coefs[idx_int],
    var_env = vc[idx_env, idx_env], var_int = vc[idx_int, idx_int],
    cov_ei = vc[idx_env, idx_int]
  )
}

# REML pool beta_env and beta_int separately
cr <- bind_rows(lapply(cohort_results, function(x) {
  tibble(cohort = x$cohort, n = x$n,
         beta_env = x$beta_env, se_env = sqrt(x$var_env),
         beta_int = x$beta_int, se_int = sqrt(x$var_int),
         cov_ei = x$cov_ei)
}))
pool_env <- rma(yi = cr$beta_env, sei = cr$se_env, method = "REML")
pool_int <- rma(yi = cr$beta_int, sei = cr$se_int, method = "REML")
w <- 1 / (cr$se_env^2 + pool_env$tau2)
pooled_cov_ei <- sum(w * cr$cov_ei) / sum(w)

# Marginal effect across IC continuum via delta method
ic_grid <- seq(-2, 2, by = 0.05)
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
write_csv(bind_rows(pooled_marg, cohort_marg), "results/marginal_effect_data.csv")

# Plot
cohort_colors <- c("CHARLS" = "#7FABD3", "ELSA" = "#B8D4A0", "HRS" = "#E8B88A",
                   "MHAS" = "#D4A5C9", "SHARE" = "#A8C8C8")
pA <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_line(data = cohort_marg,
            aes(x = ic, y = hr, color = source), linewidth = 0.4, alpha = 0.5) +
  geom_ribbon(data = pooled_marg,
              aes(x = ic, ymin = hr_lo, ymax = hr_hi),
              fill = lancet_colors["blue"], alpha = 0.15) +
  geom_line(data = pooled_marg,
            aes(x = ic, y = hr),
            color = lancet_colors["blue"], linewidth = 1) +
  scale_color_manual(values = cohort_colors, name = "Cohort") +
  coord_cartesian(ylim = c(0.55, 1.15), xlim = c(-2, 2)) +
  labs(x = "Intrinsic capacity (SD units from cohort mean)",
       y = "HR per 1-SD increase in environment support") +
  theme_lancet

ic_dens <- tibble(ic = surv_base$ic_cfa_c[!is.na(surv_base$ic_cfa_c)]) %>%
  filter(ic >= -2 & ic <= 2)
pB <- ggplot(ic_dens, aes(x = ic)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = lancet_colors["ltgrey"], color = NA) +
  labs(x = NULL, y = "Density") +
  theme_lancet

p_combined <- pA / pB + plot_layout(heights = c(5, 1))
ggsave("figures/fig3_marginal_effect.pdf", p_combined,
       width = 7.2, height = 5.5, device = cairo_pdf)
ggsave("figures/fig3_marginal_effect.tiff", p_combined,
       width = 7.2, height = 5.5, dpi = 300, compression = "lzw")

# Summary at key IC levels
for (target in c(-1, 0, 1)) {
  v <- pooled_marg %>% filter(abs(ic - target) < 0.01)
  cat(sprintf("  IC = %+d SD: HR = %.2f (%.2f-%.2f)\n",
              target, v$hr, v$hr_lo, v$hr_hi))
}
