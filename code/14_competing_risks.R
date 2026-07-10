#!/usr/bin/env Rscript
# 14_competing_risks.R — Cause-Specific Hazard Models for Competing Risks
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# Addresses the competing risks between disability and death:
#   - CSH-Disability: disability as event, death censored
#   - CSH-Death: death as event, disability censored
#
# Both models include IC × Environment interaction.
# Two-stage framework: cohort-specific Cox → REML random-effects pooling.
#
# Inputs:  data/survival_cohort.rds
# Outputs: results/competing_risks_results.csv
#          results/competing_risks_summary.txt

library(tidyverse)
library(survival)
library(metafor)

set.seed(2024)
t0 <- Sys.time()

# ============================================================
cat("=== Competing Risks: Cause-Specific Hazard Models ===\n\n")
# ============================================================

# --- Load data ---
surv_base <- readRDS("data/survival_cohort.rds")
cat(sprintf("Survival cohort: %s persons\n", format(nrow(surv_base), big.mark = ",")))

# --- Construct competing risks variables ---
# t_disab, t_death, t_censor already exist in survival_cohort
# We need:
#   cs_time  = min(t_disab, t_death, t_censor)
#   cs_status: 0=censored, 1=disability (first), 2=death (without prior disability)

surv_cr <- surv_base %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_form_c)) %>%
  mutate(
    cs_time = pmin(t_disab, t_death, t_censor, na.rm = TRUE),
    cs_status = case_when(
      cs_time == t_disab & cs_time < Inf ~ 1L,    # disability occurred first
      cs_time == t_death & cs_time < Inf ~ 2L,     # death without prior disability
      TRUE ~ 0L                                      # censored
    ),
    cs_time = pmax(cs_time, 0.001)  # avoid zero times
  )

cat(sprintf("Analytic sample (IC + Env complete): %s\n", format(nrow(surv_cr), big.mark = ",")))
cat(sprintf("  Disability events (cause 1): %d\n", sum(surv_cr$cs_status == 1)))
cat(sprintf("  Death events (cause 2):      %d\n", sum(surv_cr$cs_status == 2)))
cat(sprintf("  Censored:                    %d\n", sum(surv_cr$cs_status == 0)))

# Per-cohort breakdown
cat("\nPer-cohort breakdown:\n")
surv_cr %>%
  group_by(cohort) %>%
  summarise(
    n = n(),
    disability = sum(cs_status == 1),
    death = sum(cs_status == 2),
    censored = sum(cs_status == 0),
    .groups = "drop"
  ) %>%
  print(n = 10)

# ============================================================
cat("\n=== Model 1: CSH for Disability (death censored) ===\n")
# ============================================================
# Cause-specific hazard: model disability as event, treat death as censoring.
# This gives the instantaneous risk of disability given survival to that point.

stage1_disab <- tibble()
cohorts <- sort(unique(surv_cr$cohort))

for (coh in cohorts) {
  d <- surv_cr %>% filter(cohort == coh)

  # For CSH-disability: event = (cs_status == 1), censor everything else
  d <- d %>% mutate(event_disab_cs = as.integer(cs_status == 1))

  n_events <- sum(d$event_disab_cs)
  if (nrow(d) < 100 || n_events < 10) {
    cat(sprintf("  %s: SKIPPED (n=%d, events=%d)\n", coh, nrow(d), n_events))
    next
  }

  fit <- tryCatch(
    coxph(Surv(cs_time, event_disab_cs) ~ ic_cfa_c * env_form_c + age + female + chronic_n,
          data = d),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    cat(sprintf("  %s: MODEL FAILED\n", coh))
    next
  }

  res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
           n = nrow(d), n_events = n_events, model = "CSH_disability")
  stage1_disab <- bind_rows(stage1_disab, res)

  # Print interaction term
  int_row <- res %>% filter(grepl(":", term))
  if (nrow(int_row) > 0) {
    cat(sprintf("  %s: n=%d, events=%d, IC×Env HR=%.3f (%.3f–%.3f), p=%.4f\n",
                coh, nrow(d), n_events,
                int_row$estimate, int_row$conf.low, int_row$conf.high, int_row$p.value))
  }
}

# Pool CSH-disability interaction
cat("\n--- Pooled CSH-Disability (IC×Env interaction) ---\n")
int_disab <- stage1_disab %>% filter(grepl(":", term))

if (nrow(int_disab) >= 2) {
  ma_disab <- rma(yi = log_hr, sei = se_log_hr, data = int_disab, method = "REML")
  cat(sprintf("  Pooled HR = %.3f (%.3f–%.3f), p = %.4f, I² = %.1f%%\n",
              exp(as.numeric(ma_disab$beta)), exp(ma_disab$ci.lb), exp(ma_disab$ci.ub),
              ma_disab$pval, ma_disab$I2))
  cat(sprintf("  K = %d cohorts, N = %s, Events = %s\n",
              nrow(int_disab), format(sum(int_disab$n), big.mark = ","),
              format(sum(int_disab$n_events), big.mark = ",")))
} else {
  cat("  INSUFFICIENT DATA for pooling\n")
}

# ============================================================
cat("\n=== Model 2: CSH for Death (disability censored) ===\n")
# ============================================================
# Cause-specific hazard: model death as event, treat disability as censoring.

stage1_death <- tibble()

for (coh in cohorts) {
  d <- surv_cr %>% filter(cohort == coh)

  # For CSH-death: event = (cs_status == 2), censor everything else
  d <- d %>% mutate(event_death_cs = as.integer(cs_status == 2))

  n_events <- sum(d$event_death_cs)
  if (nrow(d) < 100 || n_events < 10) {
    cat(sprintf("  %s: SKIPPED (n=%d, events=%d)\n", coh, nrow(d), n_events))
    next
  }

  fit <- tryCatch(
    coxph(Surv(cs_time, event_death_cs) ~ ic_cfa_c * env_form_c + age + female + chronic_n,
          data = d),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    cat(sprintf("  %s: MODEL FAILED\n", coh))
    next
  }

  res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
           n = nrow(d), n_events = n_events, model = "CSH_death")
  stage1_death <- bind_rows(stage1_death, res)

  int_row <- res %>% filter(grepl(":", term))
  if (nrow(int_row) > 0) {
    cat(sprintf("  %s: n=%d, events=%d, IC×Env HR=%.3f (%.3f–%.3f), p=%.4f\n",
                coh, nrow(d), n_events,
                int_row$estimate, int_row$conf.low, int_row$conf.high, int_row$p.value))
  }
}

# Pool CSH-death interaction
cat("\n--- Pooled CSH-Death (IC×Env interaction) ---\n")
int_death <- stage1_death %>% filter(grepl(":", term))

if (nrow(int_death) >= 2) {
  ma_death <- rma(yi = log_hr, sei = se_log_hr, data = int_death, method = "REML")
  cat(sprintf("  Pooled HR = %.3f (%.3f–%.3f), p = %.4f, I² = %.1f%%\n",
              exp(as.numeric(ma_death$beta)), exp(ma_death$ci.lb), exp(ma_death$ci.ub),
              ma_death$pval, ma_death$I2))
  cat(sprintf("  K = %d cohorts, N = %s, Events = %s\n",
              nrow(int_death), format(sum(int_death$n), big.mark = ","),
              format(sum(int_death$n_events), big.mark = ",")))
} else {
  cat("  INSUFFICIENT DATA for pooling\n")
}

# ============================================================
cat("\n=== Summary Comparison ===\n")
# ============================================================

summary_rows <- tibble()

# Row 1: Primary DFS (for reference — read from survival_cohort)
# Re-run the primary interaction model for comparison
stage1_dfs <- tibble()
for (coh in cohorts) {
  d <- surv_cr %>% filter(cohort == coh)
  fit <- tryCatch(
    coxph(Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + female + chronic_n,
          data = d),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
             n = nrow(d), n_events = sum(d$event_dfs))
    stage1_dfs <- bind_rows(stage1_dfs, res)
  }
}

int_dfs <- stage1_dfs %>% filter(grepl(":", term))
if (nrow(int_dfs) >= 2) {
  ma_dfs <- rma(yi = log_hr, sei = se_log_hr, data = int_dfs, method = "REML")
  summary_rows <- bind_rows(summary_rows, tibble(
    outcome = "DFS (primary composite)",
    hr = exp(as.numeric(ma_dfs$beta)),
    ci_lo = exp(ma_dfs$ci.lb), ci_hi = exp(ma_dfs$ci.ub),
    p = ma_dfs$pval, i2 = ma_dfs$I2,
    n = sum(int_dfs$n), events = sum(int_dfs$n_events)
  ))
}

# Row 2: CSH-Disability
if (exists("ma_disab")) {
  summary_rows <- bind_rows(summary_rows, tibble(
    outcome = "CSH-Disability (death censored)",
    hr = exp(as.numeric(ma_disab$beta)),
    ci_lo = exp(ma_disab$ci.lb), ci_hi = exp(ma_disab$ci.ub),
    p = ma_disab$pval, i2 = ma_disab$I2,
    n = sum(int_disab$n), events = sum(int_disab$n_events)
  ))
}

# Row 3: CSH-Death
if (exists("ma_death")) {
  summary_rows <- bind_rows(summary_rows, tibble(
    outcome = "CSH-Death (disability censored)",
    hr = exp(as.numeric(ma_death$beta)),
    ci_lo = exp(ma_death$ci.lb), ci_hi = exp(ma_death$ci.ub),
    p = ma_death$pval, i2 = ma_death$I2,
    n = sum(int_death$n), events = sum(int_death$n_events)
  ))
}

cat("\n")
cat(sprintf("%-40s %8s %18s %8s %6s %10s %10s\n",
            "Outcome", "HR", "95% CI", "p", "I²", "N", "Events"))
cat(paste(rep("-", 105), collapse = ""), "\n")
for (i in seq_len(nrow(summary_rows))) {
  r <- summary_rows[i, ]
  cat(sprintf("%-40s %8.3f %8.3f–%-8.3f %8.4f %5.1f%% %10s %10s\n",
              r$outcome, r$hr, r$ci_lo, r$ci_hi, r$p, r$i2,
              format(r$n, big.mark = ","), format(r$events, big.mark = ",")))
}

# Save
write_csv(summary_rows, "results/competing_risks_results.csv")

# Save full stage-1 results for transparency
all_stage1 <- bind_rows(stage1_disab, stage1_death) %>%
  select(model, cohort, term, estimate, conf.low, conf.high, p.value,
         log_hr, se_log_hr, n, n_events)
write_csv(all_stage1, "results/competing_risks_stage1.csv")

# Summary text file
sink("results/competing_risks_summary.txt")
cat("=== COMPETING RISKS ANALYSIS: CAUSE-SPECIFIC HAZARD MODELS ===\n\n")
cat("Two cause-specific Cox models decompose the composite DFS endpoint:\n")
cat("  CSH-Disability: disability as event, death treated as censoring\n")
cat("  CSH-Death: death as event, disability treated as censoring\n\n")
cat("Both models: Surv(cs_time, event) ~ ic_cfa_c * env_form_c + age + female + chronic_n\n")
cat("Pooling: REML random-effects meta-analysis across 5 cohorts\n\n")
print(summary_rows, n = 10, width = 120)
cat("\nInterpretation:\n")
cat("  If CSH-Disability interaction is significant but CSH-Death is not,\n")
cat("  the IC×Env synergy operates specifically through disability prevention,\n")
cat("  consistent with the primary analysis finding.\n")
sink()

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat(sprintf("\n=== COMPETING RISKS ANALYSIS COMPLETE (%.1f min) ===\n", elapsed))
cat("Output files:\n")
cat("  results/competing_risks_results.csv\n")
cat("  results/competing_risks_stage1.csv\n")
cat("  results/competing_risks_summary.txt\n")
