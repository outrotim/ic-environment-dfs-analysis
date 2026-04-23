#!/usr/bin/env Rscript
# 15_objective_env.R — Objective Environment Proxy Sensitivity Analysis
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# The primary environment index uses 3 indicators:
#   1. Social participation (z_env_soc) — behavioural, potentially endogenous
#   2. Financial capacity (z_env_fin) — reversed income quintile, structural
#   3. Living arrangement (z_env_liv) — household composition, structural
#
# This sensitivity analysis creates an "objective proxy" environment index
# using only the two structural/material indicators (financial + living),
# excluding social participation (the most behavioural component).
#
# Rationale: Reviewers may argue social participation is an outcome rather
# than an environmental exposure. This analysis tests whether the IC × Env
# interaction holds with a purely structural environment measure.
#
# Framework: Two-stage (cohort-specific Cox → REML random-effects pooling)
#
# Inputs:  data/survival_cohort.rds
# Outputs: results/objective_env_results.csv
#          results/objective_env_summary.txt

library(tidyverse)
library(survival)
library(metafor)

set.seed(2024)
t0 <- Sys.time()

# ============================================================
cat("=== Objective Environment Proxy: Sensitivity Analysis ===\n\n")
# ============================================================

# --- Load data ---
surv_base <- readRDS("data/survival_cohort.rds")
cat(sprintf("Survival cohort: %s persons\n", format(nrow(surv_base), big.mark = ",")))

# --- Check availability of individual env components ---
has_fin <- "z_env_fin" %in% names(surv_base)
has_liv <- "z_env_liv" %in% names(surv_base)
has_soc <- "z_env_soc" %in% names(surv_base)

cat(sprintf("  z_env_fin available: %s (%.1f%% non-NA)\n", has_fin,
            if (has_fin) mean(!is.na(surv_base$z_env_fin)) * 100 else 0))
cat(sprintf("  z_env_liv available: %s (%.1f%% non-NA)\n", has_liv,
            if (has_liv) mean(!is.na(surv_base$z_env_liv)) * 100 else 0))
cat(sprintf("  z_env_soc available: %s (%.1f%% non-NA)\n", has_soc,
            if (has_soc) mean(!is.na(surv_base$z_env_soc)) * 100 else 0))

if (!has_fin || !has_liv) {
  stop("Required env components (z_env_fin, z_env_liv) not found in survival_cohort.")
}

# --- Construct objective environment proxy ---
# 2-indicator formative composite: mean of financial capacity + living arrangement
# Require both indicators to be non-missing for this strict proxy
surv_obj <- surv_base %>%
  mutate(
    env_objective = ifelse(
      !is.na(z_env_fin) & !is.na(z_env_liv),
      (z_env_fin + z_env_liv) / 2,
      NA_real_
    )
  )

cat(sprintf("\nObjective env proxy: %.1f%% available (both fin + liv required)\n",
            mean(!is.na(surv_obj$env_objective)) * 100))

# Center and create interaction
surv_obj <- surv_obj %>%
  group_by(cohort) %>%
  mutate(
    env_obj_c = scale(env_objective, center = TRUE, scale = TRUE)[, 1],
    ic_cfa_c = scale(ic_cfa, center = TRUE, scale = TRUE)[, 1]
  ) %>%
  ungroup() %>%
  mutate(ic_env_obj = ic_cfa_c * env_obj_c)

# Filter to complete cases
surv_obj <- surv_obj %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_obj_c))

cat(sprintf("Analytic sample (IC + Obj Env complete): %s\n",
            format(nrow(surv_obj), big.mark = ",")))

# Per-cohort summary
cat("\nPer-cohort breakdown:\n")
surv_obj %>%
  group_by(cohort) %>%
  summarise(
    n = n(),
    events = sum(event_dfs),
    event_rate = round(mean(event_dfs) * 100, 1),
    env_obj_mean = round(mean(env_objective, na.rm = TRUE), 3),
    env_obj_sd = round(sd(env_objective, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print(n = 10)

# --- Correlation with primary env index ---
if ("env_formative" %in% names(surv_obj)) {
  r_overall <- cor(surv_obj$env_objective, surv_obj$env_formative,
                   use = "complete.obs")
  cat(sprintf("\nCorrelation with primary env index: r = %.3f\n", r_overall))

  cat("Per-cohort correlations:\n")
  surv_obj %>%
    filter(!is.na(env_formative)) %>%
    group_by(cohort) %>%
    summarise(
      r = round(cor(env_objective, env_formative, use = "complete.obs"), 3),
      .groups = "drop"
    ) %>%
    print(n = 10)
}

# ============================================================
cat("\n=== Two-Stage Analysis: IC × Objective Environment ===\n")
# ============================================================

# --- Stage 1: Cohort-specific Cox models ---
stage1 <- tibble()
cohorts <- sort(unique(surv_obj$cohort))

for (coh in cohorts) {
  d <- surv_obj %>% filter(cohort == coh)
  n_events <- sum(d$event_dfs)

  if (nrow(d) < 100 || n_events < 10) {
    cat(sprintf("  %s: SKIPPED (n=%d, events=%d)\n", coh, nrow(d), n_events))
    next
  }

  fit <- tryCatch(
    coxph(Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_obj_c + age + female + chronic_n,
          data = d),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    cat(sprintf("  %s: MODEL FAILED\n", coh))
    next
  }

  res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
           n = nrow(d), n_events = n_events)
  stage1 <- bind_rows(stage1, res)

  # Print interaction term
  int_row <- res %>% filter(grepl(":", term))
  if (nrow(int_row) > 0) {
    cat(sprintf("  %s: n=%d, events=%d, IC×ObjEnv HR=%.3f (%.3f–%.3f), p=%.4f\n",
                coh, nrow(d), n_events,
                int_row$estimate, int_row$conf.low, int_row$conf.high, int_row$p.value))
  }
}

# --- Stage 2: Random-effects meta-analysis ---
cat("\n--- Pooled IC × Objective Environment Interaction ---\n")
int_obj <- stage1 %>% filter(grepl(":", term))

results <- tibble()

if (nrow(int_obj) >= 2) {
  ma_obj <- rma(yi = log_hr, sei = se_log_hr, data = int_obj, method = "REML")

  hr <- exp(as.numeric(ma_obj$beta))
  ci_lo <- exp(ma_obj$ci.lb)
  ci_hi <- exp(ma_obj$ci.ub)
  p_val <- ma_obj$pval
  i2 <- ma_obj$I2

  cat(sprintf("  Pooled HR = %.3f (%.3f–%.3f), p = %.4f, I² = %.1f%%\n",
              hr, ci_lo, ci_hi, p_val, i2))
  cat(sprintf("  K = %d cohorts, N = %s, Events = %s\n",
              nrow(int_obj), format(sum(int_obj$n), big.mark = ","),
              format(sum(int_obj$n_events), big.mark = ",")))

  results <- bind_rows(results, tibble(
    analysis = "Objective_env_proxy",
    term = "IC × Env_objective",
    hr = hr, ci_lo = ci_lo, ci_hi = ci_hi,
    p = p_val, i2 = i2,
    n_cohorts = nrow(int_obj),
    n_total = sum(int_obj$n),
    n_events = sum(int_obj$n_events),
    note = "2-indicator env (financial + living), excl. social participation"
  ))

  # Also pool the main effects for comparison
  for (var_name in c("ic_cfa_c", "env_obj_c")) {
    var_sub <- stage1 %>% filter(term == var_name)
    if (nrow(var_sub) >= 2) {
      ma_var <- rma(yi = log_hr, sei = se_log_hr, data = var_sub, method = "REML")
      cat(sprintf("  %s: HR = %.3f (%.3f–%.3f), p = %.2e, I² = %.1f%%\n",
                  var_name,
                  exp(as.numeric(ma_var$beta)), exp(ma_var$ci.lb), exp(ma_var$ci.ub),
                  ma_var$pval, ma_var$I2))

      results <- bind_rows(results, tibble(
        analysis = "Objective_env_proxy",
        term = var_name,
        hr = exp(as.numeric(ma_var$beta)),
        ci_lo = exp(ma_var$ci.lb), ci_hi = exp(ma_var$ci.ub),
        p = ma_var$pval, i2 = ma_var$I2,
        n_cohorts = nrow(var_sub),
        n_total = sum(var_sub$n),
        n_events = sum(var_sub$n_events),
        note = "Main effect from objective env model"
      ))
    }
  }
} else {
  cat("  INSUFFICIENT DATA for pooling\n")
}

# ============================================================
cat("\n=== Comparison: Primary vs Objective Environment ===\n")
# ============================================================

# Load primary results for comparison
primary_file <- "results/ipdma_results.csv"
if (file.exists(primary_file)) {
  primary <- read_csv(primary_file, show_col_types = FALSE)
  int_primary <- primary %>%
    filter(grepl("interaction|ic_env", variable, ignore.case = TRUE) |
           grepl(":", variable))
  if (nrow(int_primary) > 0) {
    cat("Primary interaction (3-indicator env):\n")
    print(int_primary %>% select(variable, hr, ci_lo, ci_hi, p, i2), n = 5)
  }
}

cat("\nObjective env proxy interaction (2-indicator, no social participation):\n")
print(results %>% filter(grepl("×|:", term)) %>%
        select(term, hr, ci_lo, ci_hi, p, i2), n = 5)

# --- Save ---
write_csv(results, "results/objective_env_results.csv")
write_csv(stage1, "results/objective_env_stage1.csv")

# Summary text file
sink("results/objective_env_summary.txt")
cat("=== OBJECTIVE ENVIRONMENT PROXY: SENSITIVITY ANALYSIS ===\n\n")
cat("Purpose: Test whether IC × Env interaction holds with purely structural\n")
cat("environment measure (financial capacity + living arrangement), excluding\n")
cat("social participation (the most behavioural component).\n\n")
cat("Env proxy: mean(z_env_fin, z_env_liv), both required non-missing\n")
cat("Model: Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_obj_c + age + female + chronic_n\n")
cat("Pooling: REML random-effects meta-analysis across cohorts\n\n")
print(results, n = 10, width = 120)
cat("\nInterpretation:\n")
cat("  If the interaction remains significant with the objective proxy,\n")
cat("  the IC×Env synergy is not driven solely by social participation\n")
cat("  (which could reflect reverse causation from functional ability).\n")
sink()

elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat(sprintf("\n=== OBJECTIVE ENV ANALYSIS COMPLETE (%.1f sec) ===\n", elapsed))
cat("Output files:\n")
cat("  results/objective_env_results.csv\n")
cat("  results/objective_env_stage1.csv\n")
cat("  results/objective_env_summary.txt\n")
