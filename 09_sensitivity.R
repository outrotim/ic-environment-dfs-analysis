#!/usr/bin/env Rscript
# 09_sensitivity.R — Sensitivity Analyses (17 pre-specified)
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# Organized into 5 domains:
#   A. Exposure measurement robustness (S1–S4)
#   B. Outcome definition (S5–S7)
#   C. Population scope (S8–S11)
#   D. Statistical model variants (S12–S15)
#   E. Bias assessment (S16–S17)
#
# All analyses replicate the core two-stage IPD-MA (06_ipd_ma.R) framework:
#   Stage 1: cohort-specific Cox → Stage 2: RE-MA pooling
#
# Primary target: IC×Env interaction effect on DFS
#
# Inputs:  data/analytic_final.rds, data/survival_cohort.rds
# Outputs: results/sensitivity_all.csv, results/sensitivity_summary.csv

library(tidyverse)
library(survival)
library(metafor)

set.seed(2024)
t0 <- Sys.time()

# ============================================================
cat("=== Loading data ===\n")
# ============================================================

pooled <- readRDS("data/analytic_final.rds")
surv_base <- readRDS("data/survival_cohort.rds")

cat(sprintf("Pooled panel: %d rows, %d persons\n", nrow(pooled), n_distinct(paste(pooled$cohort, pooled$id))))
cat(sprintf("Survival cohort (base): %d persons\n", nrow(surv_base)))

# ============================================================
# Helper: Two-stage IPD-MA (reusable for each sensitivity)
# ============================================================

run_two_stage <- function(dat, formula_str, label) {
  # Extract event column from formula (e.g., "Surv(time_dfs, event_dfs)" → "event_dfs")
  event_col <- sub(".*Surv\\([^,]+,\\s*([^)]+)\\).*", "\\1", formula_str)
  event_col <- trimws(event_col)

  # Stage 1: per-cohort Cox
  stage1 <- tibble()
  for (coh in sort(unique(dat$cohort))) {
    d <- dat %>% filter(cohort == coh)
    if (nrow(d) < 100 || sum(d[[event_col]], na.rm = TRUE) < 10) next

    fit <- tryCatch(
      coxph(as.formula(formula_str), data = d),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
             n = nrow(d), n_events = sum(d[[event_col]], na.rm = TRUE))
    stage1 <- bind_rows(stage1, res)
  }

  if (nrow(stage1) == 0) {
    return(tibble(analysis = label, term = NA_character_, hr = NA, ci_lo = NA,
                  ci_hi = NA, p = NA, i2 = NA, n_cohorts = 0, note = "No valid models"))
  }

  # Stage 2: pool interaction term if exists, else main term
  target_terms <- stage1 %>% distinct(term) %>% pull(term)
  int_term <- target_terms[grepl(":", target_terms)]

  results <- tibble()

  for (tt in target_terms) {
    sub <- stage1 %>% filter(term == tt)
    if (nrow(sub) < 2) {
      results <- bind_rows(results, tibble(
        analysis = label, term = tt, hr = sub$estimate[1],
        ci_lo = sub$conf.low[1], ci_hi = sub$conf.high[1],
        p = sub$p.value[1], i2 = NA, n_cohorts = 1, note = "single cohort"
      ))
      next
    }

    ma <- tryCatch(rma(yi = log_hr, sei = se_log_hr, data = sub, method = "REML"),
                   error = function(e) NULL)
    if (is.null(ma)) next

    results <- bind_rows(results, tibble(
      analysis = label, term = tt,
      hr = exp(as.numeric(ma$beta)),
      ci_lo = exp(ma$ci.lb), ci_hi = exp(ma$ci.ub),
      p = ma$pval, i2 = ma$I2,
      n_cohorts = nrow(sub),
      note = sprintf("N=%d, events=%d", sum(sub$n), sum(sub$n_events))
    ))
  }

  results
}

# Collector
all_sens <- tibble()

# ============================================================
cat("\n=== Domain A: Exposure measurement robustness ===\n")
# ============================================================

# --- S1: IC sum score instead of CFA ---
cat("S1: IC sum score (vs CFA)...\n")
surv_s1 <- surv_base %>%
  left_join(
    pooled %>% filter(!is.na(ic_sum)) %>%
      group_by(cohort, id) %>% slice_min(wave, n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(cohort, id, ic_sum_bl = ic_sum),
    by = c("cohort", "id")
  ) %>%
  filter(!is.na(ic_sum_bl)) %>%
  mutate(ic_sum_c = ic_sum_bl - mean(ic_sum_bl, na.rm = TRUE))

s1 <- run_two_stage(surv_s1,
  "Surv(time_dfs, event_dfs) ~ ic_sum_c * env_form_c + age + female + chronic_n",
  "S1_ic_sum")
all_sens <- bind_rows(all_sens, s1)
cat(sprintf("  Done: %d rows\n", nrow(s1)))

# --- S2: Env 4-indicator (add social contact) ---
cat("S2: Env 4-indicator (with social contact)...\n")
surv_s2 <- surv_base %>%
  left_join(
    pooled %>% filter(!is.na(env_formative_4)) %>%
      group_by(cohort, id) %>% slice_min(wave, n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(cohort, id, env4_bl = env_formative_4),
    by = c("cohort", "id")
  ) %>%
  filter(!is.na(env4_bl)) %>%
  mutate(env4_c = env4_bl - mean(env4_bl, na.rm = TRUE))

s2 <- run_two_stage(surv_s2,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env4_c + age + female + chronic_n",
  "S2_env_4indicator")
all_sens <- bind_rows(all_sens, s2)
cat(sprintf("  Done: %d rows\n", nrow(s2)))

# --- S3: Categorical IC×Env (6-group) ---
cat("S3: Categorical IC×Env 6-group...\n")
surv_s3 <- surv_base %>%
  left_join(
    pooled %>% filter(!is.na(ic_env_6grp)) %>%
      group_by(cohort, id) %>% slice_min(wave, n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(cohort, id, grp6 = ic_env_6grp),
    by = c("cohort", "id")
  ) %>%
  filter(!is.na(grp6))

# Reference: High IC + High Env
if ("High_IC/High_Env" %in% surv_s3$grp6) {
  surv_s3$grp6 <- relevel(factor(surv_s3$grp6), ref = "High_IC/High_Env")
} else {
  lvls <- sort(unique(surv_s3$grp6))
  surv_s3$grp6 <- factor(surv_s3$grp6, levels = lvls)
}

s3 <- run_two_stage(surv_s3,
  "Surv(time_dfs, event_dfs) ~ grp6 + age + female + chronic_n",
  "S3_categorical_6grp")
all_sens <- bind_rows(all_sens, s3)
cat(sprintf("  Done: %d rows\n", nrow(s3)))

# --- S4: Continuous IC×Env (no centering) ---
cat("S4: Uncentered continuous interaction...\n")
surv_s4 <- surv_base %>% filter(!is.na(ic_cfa_c) & !is.na(env_form_c))
# Use raw values
surv_s4 <- surv_s4 %>%
  left_join(
    pooled %>% filter(!is.na(ic_cfa) & !is.na(env_formative)) %>%
      group_by(cohort, id) %>% slice_min(wave, n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(cohort, id, ic_raw = ic_cfa, env_raw = env_formative),
    by = c("cohort", "id")
  ) %>%
  filter(!is.na(ic_raw) & !is.na(env_raw))

s4 <- run_two_stage(surv_s4,
  "Surv(time_dfs, event_dfs) ~ ic_raw * env_raw + age + female + chronic_n",
  "S4_uncentered_interaction")
all_sens <- bind_rows(all_sens, s4)
cat(sprintf("  Done: %d rows\n", nrow(s4)))

# ============================================================
cat("\n=== Domain B: Outcome definition ===\n")
# ============================================================

# --- S5: All-cause mortality only ---
cat("S5: Mortality only (not DFS)...\n")
surv_s5 <- surv_base %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_form_c) & !is.na(time_mort) & time_mort > 0)

s5 <- run_two_stage(surv_s5,
  "Surv(time_mort, event_mort) ~ ic_cfa_c * env_form_c + age + female + chronic_n",
  "S5_mortality_only")
all_sens <- bind_rows(all_sens, s5)
cat(sprintf("  Done: %d rows\n", nrow(s5)))

# --- S6: ADL disability only (no death) ---
cat("S6: ADL disability only...\n")
# Use disability_year from survival_cohort + censor at last_interview_year
surv_s6 <- surv_base %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_form_c)) %>%
  mutate(
    event_disab = as.integer(!is.na(disability_year) & disability_year > baseline_year),
    time_disab  = ifelse(event_disab == 1, disability_year - baseline_year,
                         last_interview_year - baseline_year),
    time_disab  = pmax(time_disab, 0.01)
  )

s6 <- run_two_stage(surv_s6,
  "Surv(time_disab, event_disab) ~ ic_cfa_c * env_form_c + age + female + chronic_n",
  "S6_disability_only")
all_sens <- bind_rows(all_sens, s6)
cat(sprintf("  Done: %d rows\n", nrow(s6)))

# --- S7: ADL ≥ 2 threshold (stricter disability) ---
cat("S7: Stricter disability (ADL ≥ 2)...\n")
disab2_events <- pooled %>%
  group_by(cohort, id) %>%
  arrange(wave) %>%
  mutate(period = row_number()) %>%
  filter(period > 1 & !is.na(adl_total) & adl_total >= 2) %>%
  slice_min(period, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cohort, id, disab2_year = interview_year)

surv_s7 <- surv_base %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_form_c)) %>%
  left_join(disab2_events, by = c("cohort", "id")) %>%
  mutate(
    t_disab2   = ifelse(!is.na(disab2_year), disab2_year - baseline_year, Inf),
    t_death    = ifelse(person_died == 1 & !is.na(death_year) & death_year > baseline_year,
                        death_year - baseline_year, Inf),
    t_censor   = last_interview_year - baseline_year,
    time_dfs2  = pmin(t_disab2, t_death, t_censor, na.rm = TRUE),
    event_dfs2 = as.integer(time_dfs2 < Inf & (time_dfs2 == t_disab2 | time_dfs2 == t_death)),
    time_dfs2  = ifelse(time_dfs2 == Inf | time_dfs2 <= 0, t_censor, time_dfs2),
    time_dfs2  = pmax(time_dfs2, 0.01)
  )

s7 <- run_two_stage(surv_s7,
  "Surv(time_dfs2, event_dfs2) ~ ic_cfa_c * env_form_c + age + female + chronic_n",
  "S7_adl_ge2_dfs")
all_sens <- bind_rows(all_sens, s7)
cat(sprintf("  Done: %d rows\n", nrow(s7)))

# ============================================================
cat("\n=== Domain C: Population scope ===\n")
# ============================================================

# --- S8: Age ≥ 65 only ---
cat("S8: Age ≥ 65 only...\n")
surv_s8 <- surv_base %>% filter(age >= 65 & !is.na(ic_cfa_c) & !is.na(env_form_c))
s8 <- run_two_stage(surv_s8,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + female + chronic_n",
  "S8_age_ge65")
all_sens <- bind_rows(all_sens, s8)
cat(sprintf("  Done: n=%d\n", nrow(surv_s8)))

# --- S9: Age ≥ 70 only ---
cat("S9: Age ≥ 70 only...\n")
surv_s9 <- surv_base %>% filter(age >= 70 & !is.na(ic_cfa_c) & !is.na(env_form_c))
s9 <- run_two_stage(surv_s9,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + female + chronic_n",
  "S9_age_ge70")
all_sens <- bind_rows(all_sens, s9)
cat(sprintf("  Done: n=%d\n", nrow(surv_s9)))

# --- S10: Female only ---
cat("S10: Female only...\n")
surv_s10 <- surv_base %>% filter(female == 1 & !is.na(ic_cfa_c) & !is.na(env_form_c))
s10 <- run_two_stage(surv_s10,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + chronic_n",
  "S10_female_only")
all_sens <- bind_rows(all_sens, s10)
cat(sprintf("  Done: n=%d\n", nrow(surv_s10)))

# --- S11: Male only ---
cat("S11: Male only...\n")
surv_s11 <- surv_base %>% filter(female == 0 & !is.na(ic_cfa_c) & !is.na(env_form_c))
s11 <- run_two_stage(surv_s11,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + chronic_n",
  "S11_male_only")
all_sens <- bind_rows(all_sens, s11)
cat(sprintf("  Done: n=%d\n", nrow(surv_s11)))

# ============================================================
cat("\n=== Domain D: Statistical model variants ===\n")
# ============================================================

# --- S12: Adjustment for education ---
cat("S12: Adjusted for education...\n")
surv_s12 <- surv_base %>%
  left_join(
    pooled %>% filter(!is.na(edu_years)) %>%
      group_by(cohort, id) %>% slice_min(wave, n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(cohort, id, edu_years_bl = edu_years),
    by = c("cohort", "id")
  ) %>%
  filter(!is.na(edu_years_bl) & !is.na(ic_cfa_c) & !is.na(env_form_c))

s12 <- run_two_stage(surv_s12,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + female + chronic_n + edu_years_bl",
  "S12_adj_education")
all_sens <- bind_rows(all_sens, s12)
cat(sprintf("  Done: n=%d\n", nrow(surv_s12)))

# --- S13: Adjustment for marital status ---
cat("S13: Adjusted for marital status...\n")
surv_s13 <- surv_base %>%
  left_join(
    pooled %>% filter(!is.na(marital)) %>%
      group_by(cohort, id) %>% slice_min(wave, n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(cohort, id, marital_bl = marital),
    by = c("cohort", "id")
  ) %>%
  filter(!is.na(marital_bl) & !is.na(ic_cfa_c) & !is.na(env_form_c))

s13 <- run_two_stage(surv_s13,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + female + chronic_n + factor(marital_bl)",
  "S13_adj_marital")
all_sens <- bind_rows(all_sens, s13)
cat(sprintf("  Done: n=%d\n", nrow(surv_s13)))

# --- S14: Fine-Gray competing risk (death as competing event for disability) ---
cat("S14: Fine-Gray competing risk for disability...\n")

# Recode status: 0=censored, 1=disability, 2=death (without disability)
surv_s14 <- surv_base %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_form_c)) %>%
  mutate(
    t_disab_s14   = ifelse(!is.na(disability_year) & disability_year > baseline_year,
                           disability_year - baseline_year, Inf),
    t_death_s14   = ifelse(person_died == 1 & !is.na(death_year) & death_year > baseline_year,
                           death_year - baseline_year, Inf),
    t_censor_s14  = last_interview_year - baseline_year,
    cr_time = pmin(t_disab_s14, t_death_s14, t_censor_s14, na.rm = TRUE),
    cr_status = case_when(
      cr_time == t_disab_s14 & cr_time < Inf ~ 1L,  # disability
      cr_time == t_death_s14 & cr_time < Inf ~ 2L,   # death (competing)
      TRUE ~ 0L                                        # censored
    ),
    cr_time = pmax(cr_time, 0.01)
  )

# Per-cohort Fine-Gray → pool
stage1_cr <- tibble()
for (coh in sort(unique(surv_s14$cohort))) {
  d <- surv_s14 %>% filter(cohort == coh)
  if (nrow(d) < 100 || sum(d$cr_status == 1) < 10) next

  fg <- tryCatch({
    fg_fit <- coxph(Surv(cr_time, factor(cr_status)) ~ ic_cfa_c * env_form_c +
                      age + female + chronic_n,
                    data = d, id = id)
    # Extract from multistate cox (Fine-Gray-like in survival package)
    broom::tidy(fg_fit, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(grepl("1:ic_cfa_c:env_form_c", term) | grepl("^1:", term)) %>%
      mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
             n = nrow(d), n_events = sum(d$cr_status == 1))
  }, error = function(e) NULL)

  if (!is.null(fg) && nrow(fg) > 0) stage1_cr <- bind_rows(stage1_cr, fg)
}

if (nrow(stage1_cr) > 0) {
  # Pool interaction
  int_sub <- stage1_cr %>% filter(grepl(":", term))
  if (nrow(int_sub) >= 2) {
    ma_cr <- tryCatch(rma(yi = log_hr, sei = se_log_hr, data = int_sub, method = "REML"),
                      error = function(e) NULL)
    if (!is.null(ma_cr)) {
      s14 <- tibble(analysis = "S14_competing_risk", term = "IC×Env_CR",
                    hr = exp(as.numeric(ma_cr$beta)), ci_lo = exp(ma_cr$ci.lb),
                    ci_hi = exp(ma_cr$ci.ub), p = ma_cr$pval, i2 = ma_cr$I2,
                    n_cohorts = nrow(int_sub),
                    note = sprintf("N=%d", sum(int_sub$n)))
      all_sens <- bind_rows(all_sens, s14)
      cat(sprintf("  Done: pooled HR=%.3f\n", exp(as.numeric(ma_cr$beta))))
    }
  }
} else {
  # Fallback: just note it
  all_sens <- bind_rows(all_sens, tibble(
    analysis = "S14_competing_risk", term = "IC×Env_CR",
    hr = NA, ci_lo = NA, ci_hi = NA, p = NA, i2 = NA,
    n_cohorts = 0, note = "Fine-Gray multistate failed"))
  cat("  Fine-Gray: no valid models (will note in manuscript)\n")
}

# --- S15: Stratified by cohort (leave-one-out) ---
cat("S15: Leave-one-out cohort sensitivity...\n")
cohorts <- sort(unique(surv_base$cohort))
s15_results <- tibble()

for (leave_out in cohorts) {
  dat_loo <- surv_base %>%
    filter(cohort != leave_out & !is.na(ic_cfa_c) & !is.na(env_form_c))

  res <- run_two_stage(dat_loo,
    "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + female + chronic_n",
    paste0("S15_LOO_excl_", leave_out))

  int_res <- res %>% filter(grepl(":", term))
  if (nrow(int_res) > 0) s15_results <- bind_rows(s15_results, int_res)
}

all_sens <- bind_rows(all_sens, s15_results)
cat(sprintf("  Done: %d LOO analyses\n", nrow(s15_results)))

# ============================================================
cat("\n=== Domain E: Bias assessment ===\n")
# ============================================================

# --- S16: E-value calculation for unmeasured confounding ---
cat("S16: E-value for unmeasured confounding...\n")
# Get main analysis interaction HR
main_int <- all_sens %>% filter(analysis == "S1_ic_sum") %>% slice(0) # placeholder

# Use base analysis from survival_cohort (the main 06 result)
# Recompute quickly for reference
surv_main <- surv_base %>% filter(!is.na(ic_cfa_c) & !is.na(env_form_c))
main_res <- run_two_stage(surv_main,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + female + chronic_n",
  "S16_main_reference")
main_int_row <- main_res %>% filter(grepl(":", term))

if (nrow(main_int_row) > 0) {
  hr <- main_int_row$hr[1]
  ci_lo <- main_int_row$ci_lo[1]

  # E-value formula: HR + sqrt(HR * (HR - 1))
  # For protective (HR < 1): use 1/HR
  hr_for_eval <- ifelse(hr < 1, 1/hr, hr)
  e_value <- hr_for_eval + sqrt(hr_for_eval * (hr_for_eval - 1))

  ci_for_eval <- ifelse(ci_lo < 1 & hr < 1, 1/ci_lo,
                         ifelse(ci_lo > 1, ci_lo, 1))
  e_value_ci <- ifelse(ci_for_eval > 1,
                        ci_for_eval + sqrt(ci_for_eval * (ci_for_eval - 1)),
                        1)

  s16 <- tibble(analysis = "S16_e_value", term = "IC×Env interaction",
                hr = hr, ci_lo = ci_lo, ci_hi = main_int_row$ci_hi[1],
                p = main_int_row$p[1],
                i2 = main_int_row$i2[1], n_cohorts = main_int_row$n_cohorts[1],
                note = sprintf("E-value=%.2f (CI=%.2f)", e_value, e_value_ci))
  all_sens <- bind_rows(all_sens, s16)
  cat(sprintf("  E-value: %.2f (CI bound: %.2f)\n", e_value, e_value_ci))
} else {
  cat("  Could not compute E-value (no interaction term)\n")
}

# Also add the main reference result
all_sens <- bind_rows(all_sens, main_res)

# --- S17: Excluding persons with < 3 waves (follow-up adequacy) ---
cat("S17: Excluding persons with < 3 waves...\n")
n_waves <- pooled %>%
  filter(!is.na(ic_cfa) & !is.na(env_formative)) %>%
  group_by(cohort, id) %>%
  summarise(n_waves = n(), .groups = "drop")

surv_s17 <- surv_base %>%
  inner_join(n_waves %>% filter(n_waves >= 3), by = c("cohort", "id")) %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_form_c))

s17 <- run_two_stage(surv_s17,
  "Surv(time_dfs, event_dfs) ~ ic_cfa_c * env_form_c + age + female + chronic_n",
  "S17_ge3_waves")
all_sens <- bind_rows(all_sens, s17)
cat(sprintf("  Done: n=%d (of %d base)\n", nrow(surv_s17), nrow(surv_base)))

# ============================================================
cat("\n=== Final Summary ===\n")
# ============================================================

write_csv(all_sens, "results/sensitivity_all.csv")

# Summary table: interaction term only
int_summary <- all_sens %>%
  filter(grepl(":", term) | grepl("E-value|_CR", analysis)) %>%
  select(analysis, term, hr, ci_lo, ci_hi, p, i2, n_cohorts, note) %>%
  arrange(analysis)

write_csv(int_summary, "results/sensitivity_summary.csv")

cat("\n--- SENSITIVITY ANALYSIS SUMMARY (IC×Env interaction) ---\n")
cat(sprintf("%-35s %8s %18s %8s %6s\n", "Analysis", "HR", "95% CI", "p", "I²"))
cat(paste(rep("-", 85), collapse = ""), "\n")

for (i in seq_len(nrow(int_summary))) {
  r <- int_summary[i, ]
  ci_str <- sprintf("%.3f–%.3f", r$ci_lo, r$ci_hi)
  cat(sprintf("%-35s %8.3f %18s %8.4f %5.1f%%\n",
              r$analysis, r$hr, ci_str, r$p, ifelse(is.na(r$i2), 0, r$i2)))
}

# Count how many sensitivity analyses show significant interaction
n_sig <- sum(int_summary$p < 0.05, na.rm = TRUE)
n_total <- nrow(int_summary)
cat(sprintf("\nSignificant (p<0.05): %d / %d (%.0f%%)\n", n_sig, n_total, 100*n_sig/n_total))

# Direction consistency
n_protective <- sum(int_summary$hr < 1, na.rm = TRUE)
cat(sprintf("Protective direction (HR<1): %d / %d (%.0f%%)\n",
            n_protective, n_total, 100*n_protective/n_total))

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat(sprintf("\n=== SENSITIVITY ANALYSES COMPLETE (%.1f min) ===\n", elapsed))
cat("Output files:\n")
cat("  results/sensitivity_all.csv\n")
cat("  results/sensitivity_summary.csv\n")
