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

# --- Load data ---
surv_base <- readRDS("data/survival_cohort.rds")

# --- Construct objective environment proxy ---
surv_obj <- surv_base %>%
  mutate(
    env_objective = ifelse(
      !is.na(z_env_fin) & !is.na(z_env_liv),
      (z_env_fin + z_env_liv) / 2,
      NA_real_
    )
  )

# Center and create interaction
surv_obj <- surv_obj %>%
  group_by(cohort) %>%
  mutate(
    env_obj_c = scale(env_objective, center = TRUE, scale = TRUE)[, 1],
    ic_cfa_c = scale(ic_cfa, center = TRUE, scale = TRUE)[, 1]
  ) %>%
  ungroup() %>%
  mutate(ic_env_obj = ic_cfa_c * env_obj_c) %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_obj_c))

# Stage 1: Cohort-specific Cox
stage1 <- tibble()
cohorts <- sort(unique(surv_obj$cohort))
for (coh in cohorts) {
  d <- surv_obj %>% filter(cohort == coh)
  n_events <- sum(d$event_dfs)
  if (nrow(d) < 100 || n_events < 10) next
  fit <- tryCatch(
    coxph(Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_obj_c + age + female + chronic_n,
          data = d),
    error = function(e) NULL)
  if (is.null(fit)) next
  res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
           n = nrow(d), n_events = n_events)
  stage1 <- bind_rows(stage1, res)
}

# Stage 2: REML pooling
int_obj <- stage1 %>% filter(grepl(":", term))
ma_obj <- rma(yi = log_hr, sei = se_log_hr, data = int_obj, method = "REML")

results <- tibble(
  analysis = "Objective_env_proxy",
  term = "IC × Env_objective",
  hr = exp(as.numeric(ma_obj$beta)),
  ci_lo = exp(ma_obj$ci.lb),
  ci_hi = exp(ma_obj$ci.ub),
  p = ma_obj$pval,
  i2 = ma_obj$I2,
  n_cohorts = nrow(int_obj),
  n_total = sum(int_obj$n),
  n_events = sum(int_obj$n_events),
  note = "2-indicator env (financial + living), excl. social participation"
)
print(results)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write_csv(results, "results/objective_env_results.csv")
write_csv(stage1, "results/objective_env_stage1.csv")

cat(sprintf("Completed in %.1f sec\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
