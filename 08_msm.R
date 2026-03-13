#!/usr/bin/env Rscript
# 08_msm.R — Marginal Structural Model with IPTW
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# Purpose: Estimate the causal effect of time-varying Env support on DFS,
#          adjusting for time-varying confounding by IC (which is both
#          affected by prior Env and affects future Env).
#
# Key causal question:
#   "If we could intervene to maintain high Env support, what would be the
#    effect on DFS — and does this effect differ by IC level?"
#
# Pipeline:
#   Step 1: Prepare person-period dataset
#   Step 2: Construct IPTW (stabilized) per cohort
#   Step 3: Weighted MSM for Env → DFS (overall)
#   Step 4: MSM with IC×Env interaction (causal interaction)
#   Step 5: Stratified MSM by baseline IC tertile
#   Step 6: Weight diagnostics (truncation, balance)
#   Step 7: Summary + per-cohort → RE-MA
#
# Inputs:  data/analytic_final.rds
# Outputs: results/msm_weights_diagnostics.csv, results/msm_main.csv,
#          results/msm_interaction.csv, results/msm_stratified.csv

library(tidyverse)
library(survival)
library(WeightIt)
library(cobalt)
library(metafor)

set.seed(2024)
t0 <- Sys.time()

# ============================================================
cat("=== Step 1: Prepare person-period dataset ===\n")
# ============================================================

pooled <- readRDS("data/analytic_final.rds")

# Need: id, cohort, wave (ordered), time-varying IC/Env/chronic/ADL, baseline covariates
# Keep persons with ≥2 waves and valid IC + Env

pp <- pooled %>%
  filter(!is.na(ic_cfa) & !is.na(env_formative)) %>%
  group_by(cohort, id) %>%
  filter(n() >= 2) %>%
  arrange(wave) %>%
  mutate(
    period   = row_number(),          # sequential period index
    time_yr  = interview_year - first(interview_year),  # years from baseline
    # Binary treatment: high Env (above cohort-wave median)
    env_high = NA_integer_   # will be filled per cohort-wave
  ) %>%
  ungroup()

# Compute env_high as above cohort-wave median
pp <- pp %>%
  group_by(cohort, wave) %>%
  mutate(env_high = as.integer(env_formative >= median(env_formative, na.rm = TRUE))) %>%
  ungroup()

# Create lagged variables (previous wave values)
pp <- pp %>%
  group_by(cohort, id) %>%
  mutate(
    env_high_lag   = lag(env_high),
    ic_cfa_lag     = lag(ic_cfa),
    chronic_n_lag  = lag(chronic_n),
    adl_total_lag  = lag(adl_total)
  ) %>%
  ungroup()

# Baseline covariates (from first period)
baseline_cov <- pp %>%
  filter(period == 1) %>%
  select(cohort, id,
         age_bl   = age,
         female_bl = female,
         ic_bl    = ic_cfa,
         env_bl   = env_formative,
         chronic_bl = chronic_n) %>%
  mutate(
    ic_tertile = case_when(
      is.na(ic_bl)                                 ~ NA_character_,
      ic_bl <= quantile(ic_bl, 1/3, na.rm = TRUE)  ~ "Low_IC",
      ic_bl <= quantile(ic_bl, 2/3, na.rm = TRUE)  ~ "Mid_IC",
      TRUE                                          ~ "High_IC"
    )
  )

# Compute ic_tertile per cohort (not pooled)
baseline_cov <- pp %>%
  filter(period == 1) %>%
  select(cohort, id, age_bl = age, female_bl = female,
         ic_bl = ic_cfa, env_bl = env_formative, chronic_bl = chronic_n) %>%
  group_by(cohort) %>%
  mutate(
    ic_tertile = case_when(
      is.na(ic_bl) ~ NA_character_,
      ic_bl <= quantile(ic_bl, 1/3, na.rm = TRUE) ~ "Low_IC",
      ic_bl <= quantile(ic_bl, 2/3, na.rm = TRUE) ~ "Mid_IC",
      TRUE ~ "High_IC"
    )
  ) %>%
  ungroup()

pp <- pp %>% left_join(baseline_cov, by = c("cohort", "id"))

# Construct DFS event per person (same logic as 06_ipd_ma.R)
# First DFS event: first ADL ≥ 1 (after baseline) or death
pp_baseline <- pp %>% filter(period == 1)

# Find first post-baseline disability
disability_events <- pp %>%
  filter(period > 1 & !is.na(adl_total) & adl_total >= 1) %>%
  group_by(cohort, id) %>%
  slice_min(period, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cohort, id, disab_period = period, disab_year = interview_year)

# Person-level survival
person_surv <- pp %>%
  group_by(cohort, id) %>%
  summarise(
    baseline_year     = first(interview_year),
    last_year         = last(interview_year),
    n_periods         = n(),
    died              = any(person_died == 1, na.rm = TRUE),
    death_year        = first(death_year[!is.na(death_year)]),
    .groups           = "drop"
  ) %>%
  left_join(disability_events, by = c("cohort", "id")) %>%
  left_join(baseline_cov, by = c("cohort", "id")) %>%
  mutate(
    # Exclude baseline ADL ≥ 1
    bl_adl = pp$adl_total[match(paste(cohort, id), paste(pp$cohort, pp$id))],
    # Time to disability
    t_disab = ifelse(!is.na(disab_year), disab_year - baseline_year, Inf),
    # Time to death
    t_death = ifelse(died & !is.na(death_year) & death_year > baseline_year,
                     death_year - baseline_year, Inf),
    # Censor time
    t_censor = last_year - baseline_year,
    # DFS
    time_dfs_raw = pmin(t_disab, t_death, t_censor, na.rm = TRUE),
    event_dfs = as.integer(time_dfs_raw < Inf &
                             (time_dfs_raw == t_disab | time_dfs_raw == t_death)),
    time_dfs = ifelse(event_dfs == 1, time_dfs_raw, t_censor)
  ) %>%
  filter(!is.na(bl_adl) & bl_adl == 0 & time_dfs > 0)

cat(sprintf("Person-period: %d rows, %d persons\n", nrow(pp), n_distinct(pp$uid)))
cat(sprintf("MSM survival cohort: %d persons (baseline ADL=0, FU>0)\n", nrow(person_surv)))

# ============================================================
cat("\n=== Step 2: Construct IPTW (stabilized) ===\n")
# ============================================================

# For each person at each post-baseline period:
#   Denominator: P(env_high_t | env_high_{t-1}, ic_cfa_t, chronic_n_t, age, female)
#   Numerator:   P(env_high_t | env_high_{t-1})
#   sw_t = num / denom
#   SW = cumulative product of sw_t

# Work on post-baseline periods only (period > 1)
pp_post <- pp %>%
  filter(period > 1 & !is.na(env_high) & !is.na(env_high_lag)) %>%
  # Remove persons excluded from survival analysis
  semi_join(person_surv, by = c("cohort", "id"))

cat(sprintf("Post-baseline person-periods for IPTW: %d rows\n", nrow(pp_post)))

# Fill missing confounders with LOCF within person
pp_post <- pp_post %>%
  group_by(cohort, id) %>%
  fill(chronic_n, adl_total, .direction = "down") %>%
  ungroup() %>%
  mutate(
    chronic_n  = replace_na(chronic_n, 0),
    adl_total  = replace_na(adl_total, 0)
  )

# Fit denominator and numerator models PER COHORT
compute_iptw <- function(dat, cohort_name) {
  cat(sprintf("  [%s] n=%d person-periods...\n", cohort_name, nrow(dat)))

  # Denominator model (full)
  denom_fit <- glm(env_high ~ env_high_lag + ic_cfa + chronic_n + age + female,
                   family = binomial, data = dat)
  p_denom <- predict(denom_fit, type = "response")

  # Numerator model (stabilized — only prior treatment)
  num_fit <- glm(env_high ~ env_high_lag,
                 family = binomial, data = dat)
  p_num <- predict(num_fit, type = "response")

  # Stabilized weight for this period
  dat$sw <- ifelse(dat$env_high == 1, p_num / p_denom, (1 - p_num) / (1 - p_denom))

  # Cumulative product within person
  dat <- dat %>%
    group_by(cohort, id) %>%
    mutate(cum_sw = cumprod(sw)) %>%
    ungroup()

  # Truncate at 1st and 99th percentile
  q01 <- quantile(dat$cum_sw, 0.01, na.rm = TRUE)
  q99 <- quantile(dat$cum_sw, 0.99, na.rm = TRUE)
  dat$cum_sw_trunc <- pmin(pmax(dat$cum_sw, q01), q99)

  cat(sprintf("    Weights: mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
              mean(dat$cum_sw_trunc), median(dat$cum_sw_trunc),
              min(dat$cum_sw_trunc), max(dat$cum_sw_trunc)))
  cat(sprintf("    Truncated at [%.3f, %.3f]\n", q01, q99))

  dat
}

# Apply per cohort
pp_weighted <- tibble()
for (coh in sort(unique(pp_post$cohort))) {
  dat_coh <- pp_post %>% filter(cohort == coh)
  if (nrow(dat_coh) < 100) next
  weighted_coh <- tryCatch(compute_iptw(dat_coh, coh), error = function(e) {
    cat(sprintf("    FAILED: %s\n", e$message)); NULL
  })
  if (!is.null(weighted_coh)) pp_weighted <- bind_rows(pp_weighted, weighted_coh)
}

# Get final weight per person (last period's cumulative weight)
final_weights <- pp_weighted %>%
  group_by(cohort, id) %>%
  slice_max(period, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cohort, id, iptw = cum_sw_trunc)

# Merge weights into survival data
msm_surv <- person_surv %>%
  left_join(final_weights, by = c("cohort", "id")) %>%
  filter(!is.na(iptw))

cat(sprintf("\nMSM analysis dataset: %d persons with valid IPTW\n", nrow(msm_surv)))

# Weight diagnostics
wt_diag <- msm_surv %>%
  group_by(cohort) %>%
  summarise(
    n = n(),
    wt_mean = round(mean(iptw), 3),
    wt_sd   = round(sd(iptw), 3),
    wt_min  = round(min(iptw), 3),
    wt_max  = round(max(iptw), 3),
    .groups = "drop"
  )
cat("\nWeight diagnostics:\n")
print(wt_diag)
write_csv(wt_diag, "results/msm_weights_diagnostics.csv")

# ============================================================
cat("\n=== Step 3: Weighted MSM — Env → DFS (overall) ===\n")
# ============================================================

# Compute cumulative Env exposure: proportion of periods with high Env
cum_env <- pp %>%
  filter(!is.na(env_high)) %>%
  group_by(cohort, id) %>%
  summarise(
    cum_env_prop = mean(env_high),
    mean_env     = mean(env_formative, na.rm = TRUE),
    .groups = "drop"
  )

msm_surv <- msm_surv %>%
  left_join(cum_env, by = c("cohort", "id"))

# --- Model A: Cumulative Env exposure → DFS (IPTW-weighted, stratified by cohort) ---
msm_A <- coxph(Surv(time_dfs, event_dfs) ~ cum_env_prop + age_bl + female_bl + chronic_bl +
                  strata(cohort),
                data = msm_surv, weights = iptw, robust = TRUE)

cat("Model A: Cumulative Env → DFS (IPTW-weighted)\n")
msm_A_res <- broom::tidy(msm_A, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(!grepl("strata", term))
print(msm_A_res %>% select(term, estimate, conf.low, conf.high, p.value, robust.se))

# --- Model B: Unweighted (for comparison) ---
msm_B <- coxph(Surv(time_dfs, event_dfs) ~ cum_env_prop + age_bl + female_bl + chronic_bl +
                  strata(cohort),
                data = msm_surv, robust = TRUE)

cat("\nModel B: Same model UNWEIGHTED (conventional)\n")
msm_B_res <- broom::tidy(msm_B, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(!grepl("strata", term))
print(msm_B_res %>% select(term, estimate, conf.low, conf.high, p.value))

# ============================================================
cat("\n=== Step 4: MSM with IC×Env interaction ===\n")
# ============================================================

# Add baseline IC + interaction
msm_C <- coxph(Surv(time_dfs, event_dfs) ~
                  cum_env_prop + ic_bl + cum_env_prop:ic_bl +
                  age_bl + female_bl + chronic_bl + strata(cohort),
                data = msm_surv, weights = iptw, robust = TRUE)

cat("Model C: Env + IC + IC×Env interaction (IPTW-weighted)\n")
msm_C_res <- broom::tidy(msm_C, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(!grepl("strata", term))
print(msm_C_res %>% select(term, estimate, conf.low, conf.high, p.value))

# ============================================================
cat("\n=== Step 5: Stratified MSM by baseline IC tertile ===\n")
# ============================================================

strat_results <- tibble()

for (tert in c("Low_IC", "Mid_IC", "High_IC")) {
  dat_tert <- msm_surv %>% filter(ic_tertile == tert)
  if (nrow(dat_tert) < 100) next

  fit <- tryCatch(
    coxph(Surv(time_dfs, event_dfs) ~ cum_env_prop + age_bl + female_bl + chronic_bl +
            strata(cohort),
          data = dat_tert, weights = iptw, robust = TRUE),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == "cum_env_prop") %>%
      mutate(stratum = tert, n = nrow(dat_tert),
             n_events = sum(dat_tert$event_dfs))
    strat_results <- bind_rows(strat_results, res)
  }
}

cat("Env → DFS stratified by baseline IC tertile (IPTW-weighted):\n")
print(strat_results %>% select(stratum, estimate, conf.low, conf.high, p.value, n, n_events))

# Interaction p-value from Model C
int_p <- msm_C_res %>% filter(grepl(":", term)) %>% pull(p.value)
cat(sprintf("\nInteraction p-value (Model C): %.4f\n", int_p))

# ============================================================
cat("\n=== Step 6: Per-cohort MSM → RE-MA ===\n")
# ============================================================

# Per-cohort weighted Cox for main effect of Env
stage1_msm <- tibble()

for (coh in sort(unique(msm_surv$cohort))) {
  dat_coh <- msm_surv %>% filter(cohort == coh)
  if (nrow(dat_coh) < 50) next

  # Main effect
  fit_main <- tryCatch(
    coxph(Surv(time_dfs, event_dfs) ~ cum_env_prop + age_bl + female_bl + chronic_bl,
          data = dat_coh, weights = iptw, robust = TRUE),
    error = function(e) NULL
  )

  if (!is.null(fit_main)) {
    res <- broom::tidy(fit_main, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term == "cum_env_prop") %>%
      mutate(cohort = coh, model = "main",
             log_hr = log(estimate), se_log_hr = robust.se,
             n = nrow(dat_coh), n_events = sum(dat_coh$event_dfs))
    stage1_msm <- bind_rows(stage1_msm, res)
  }

  # With interaction
  fit_int <- tryCatch(
    coxph(Surv(time_dfs, event_dfs) ~ cum_env_prop + ic_bl + cum_env_prop:ic_bl +
            age_bl + female_bl + chronic_bl,
          data = dat_coh, weights = iptw, robust = TRUE),
    error = function(e) NULL
  )

  if (!is.null(fit_int)) {
    res_int <- broom::tidy(fit_int, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(grepl("cum_env_prop", term)) %>%
      mutate(cohort = coh, model = "interaction",
             log_hr = log(estimate), se_log_hr = robust.se,
             n = nrow(dat_coh))
    stage1_msm <- bind_rows(stage1_msm, res_int)
  }
}

cat("Stage 1 per-cohort MSM results:\n")
print(stage1_msm %>% select(cohort, model, term, estimate, conf.low, conf.high, p.value))

# Stage 2: pool main effect of Env
main_sub <- stage1_msm %>% filter(model == "main")
if (nrow(main_sub) >= 2) {
  ma_main <- rma(yi = log_hr, sei = se_log_hr, data = main_sub, method = "REML")
  cat(sprintf("\nPooled MSM Env → DFS: HR=%.3f (%.3f–%.3f), p=%.4f, I²=%.1f%%\n",
              exp(ma_main$beta), exp(ma_main$ci.lb), exp(ma_main$ci.ub),
              ma_main$pval, ma_main$I2))
}

# Pool interaction term
int_sub <- stage1_msm %>% filter(model == "interaction" & grepl(":", term))
if (nrow(int_sub) >= 2) {
  ma_int <- rma(yi = log_hr, sei = se_log_hr, data = int_sub, method = "REML")
  cat(sprintf("Pooled MSM IC×Env interaction: HR=%.3f (%.3f–%.3f), p=%.4f, I²=%.1f%%\n",
              exp(ma_int$beta), exp(ma_int$ci.lb), exp(ma_int$ci.ub),
              ma_int$pval, ma_int$I2))
}

# ============================================================
cat("\n=== Step 7: Summary ===\n")
# ============================================================

# Compile main results
all_results <- bind_rows(
  msm_A_res %>% filter(term == "cum_env_prop") %>%
    mutate(analysis = "MSM_main_pooled_strat"),
  msm_C_res %>% mutate(analysis = "MSM_interaction_pooled_strat"),
  strat_results %>% mutate(analysis = paste0("MSM_stratified_", stratum))
)

write_csv(all_results, "results/msm_main.csv")
write_csv(stage1_msm, "results/msm_percohort.csv")
write_csv(strat_results, "results/msm_stratified.csv")

cat("\n--- KEY MSM FINDINGS ---\n")
cat("1. Causal effect of sustained high Env on DFS (IPTW-weighted):\n")
env_hr <- msm_A_res %>% filter(term == "cum_env_prop")
cat(sprintf("   HR = %.3f (%.3f–%.3f), p = %.4f\n",
            env_hr$estimate, env_hr$conf.low, env_hr$conf.high, env_hr$p.value))

cat("2. IC × Env causal interaction:\n")
int_hr <- msm_C_res %>% filter(grepl(":", term))
if (nrow(int_hr) > 0) {
  cat(sprintf("   HR = %.3f (%.3f–%.3f), p = %.4f\n",
              int_hr$estimate, int_hr$conf.low, int_hr$conf.high, int_hr$p.value))
}

cat("3. Env effect stratified by baseline IC:\n")
for (i in seq_len(nrow(strat_results))) {
  r <- strat_results[i, ]
  cat(sprintf("   %s: HR = %.3f (%.3f–%.3f), p = %.4f\n",
              r$stratum, r$estimate, r$conf.low, r$conf.high, r$p.value))
}

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat(sprintf("\n=== MSM ANALYSIS COMPLETE (%.1f min) ===\n", elapsed))
cat("Output files:\n")
cat("  results/msm_weights_diagnostics.csv\n")
cat("  results/msm_main.csv\n")
cat("  results/msm_percohort.csv\n")
cat("  results/msm_stratified.csv\n")
