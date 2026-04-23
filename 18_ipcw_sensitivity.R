# ==============================================================================
# 18_ipcw_sensitivity.R
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# Purpose: Quantify and adjust for potential selection bias from the exclusion
#          of 44,555 participants who had no post-baseline follow-up wave
#          (loss-to-follow-up, LTFU) in the primary analytic cohort (n=102,006).
#
# Rationale:
#   - Main analysis excluded: 174,695 baseline -> 146,561 (drop 26,462 ADL>=1
#     + 1,672 ADL NA) -> 102,006 analytic cohort (drop 44,555 zero-follow-up).
#   - The 44,555 LTFU represent ~30% of the eligible baseline (ADL=0) population.
#   - If LTFU is non-random with respect to IC/Env, estimated interaction HR
#     could be biased.
#
# Approach:
#   (1) Baseline-characteristics comparison between LTFU-excluded (n=44,555)
#       and analytic (n=102,006) groups, with standardized mean differences.
#   (2) Inverse-probability-of-censoring weighting (IPCW): logistic model
#       predicting P(retained | baseline covariates), weights applied to
#       the analytic cohort and the M4 Cox re-fitted with those weights.
#   (3) Compare IPCW-weighted pooled HR to the unweighted primary estimate.
#
# Input:   data/analytic_final.rds  (from 04_env_index.R, includes all baseline)
#          data/survival_cohort.rds (from 06_ipd_ma.R,   102,006 analytic)
#
# Output:  results/ipcw_baseline_compare.csv  — SMDs between LTFU vs analyzed
#          results/ipcw_weights_diagnostics.csv  — weight distribution per cohort
#          results/ipcw_stage1.csv  — per-cohort M4 Cox with IPCW
#          results/ipcw_stage2_pooled.csv  — REML-pooled IPCW-weighted HR
#
# Dependencies: tidyverse, survival, metafor
# Author:  Study 17 Team
# Date:    2026-04-22
# Version: 1.0
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(metafor)
})

SCRIPT_DIR <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) here::here()  # set working directory to the repo root before running
)
DATA_DIR    <- file.path(SCRIPT_DIR, "data")
RESULTS_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# STEP 1: Reconstruct the eligible baseline (n=146,561) and LTFU flag
# ============================================================================

step1_baseline_flag <- function() {
  cat(strrep("=", 70), "\n"); cat("Step 1: Reconstructing eligible baseline + LTFU flag\n"); cat(strrep("=", 70), "\n\n")

  pooled  <- readRDS(file.path(DATA_DIR, "analytic_final.rds"))
  surv    <- readRDS(file.path(DATA_DIR, "survival_cohort.rds"))

  # Baseline = first wave per person per cohort, with ADL present
  baseline <- pooled %>%
    group_by(cohort, id) %>% slice_min(wave, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    rename(baseline_wave = wave, baseline_year = interview_year)

  cat(sprintf("  Baseline rows: %s\n", format(nrow(baseline), big.mark=",")))

  eligible <- baseline %>% filter(!is.na(adl_total) & adl_total == 0)
  cat(sprintf("  Eligible (ADL=0 at baseline): %s\n", format(nrow(eligible), big.mark=",")))

  # LTFU flag: eligible but NOT in survival_cohort
  surv_keys <- surv %>% select(cohort, id) %>% distinct() %>% mutate(in_analysis = 1L)
  eligible <- eligible %>%
    left_join(surv_keys, by = c("cohort", "id")) %>%
    mutate(in_analysis = ifelse(is.na(in_analysis), 0L, 1L))

  n_in  <- sum(eligible$in_analysis == 1)
  n_out <- sum(eligible$in_analysis == 0)
  cat(sprintf("  Analytic (retained): %s\n", format(n_in, big.mark=",")))
  cat(sprintf("  LTFU (excluded):     %s\n", format(n_out, big.mark=",")))
  cat(sprintf("  LTFU rate: %.1f%%\n", 100 * n_out / nrow(eligible)))

  eligible
}

# ============================================================================
# STEP 2: Baseline-characteristics comparison (LTFU vs analytic)
# ============================================================================

smd <- function(x, g) {
  # Standardized mean difference for continuous variable x between groups g=0/1
  m1 <- mean(x[g == 1], na.rm = TRUE); m0 <- mean(x[g == 0], na.rm = TRUE)
  s1 <- var(x[g == 1], na.rm = TRUE);  s0 <- var(x[g == 0], na.rm = TRUE)
  sp <- sqrt((s1 + s0) / 2)
  if (!is.finite(sp) || sp == 0) NA_real_ else (m1 - m0) / sp
}

step2_smd <- function(eligible) {
  cat("\n", strrep("=", 70), "\n"); cat("Step 2: Baseline comparison (LTFU vs Analytic)\n"); cat(strrep("=", 70), "\n\n")

  vars_cont <- c("age", "chronic_n", "ic_cfa", "env_formative", "edu_years")
  vars_bin  <- c("female")

  out <- list()
  for (v in c(vars_cont, vars_bin)) {
    if (!v %in% names(eligible)) next
    vv <- eligible[[v]]
    g  <- eligible$in_analysis
    m1 <- mean(vv[g == 1], na.rm = TRUE); sd1 <- sd(vv[g == 1], na.rm = TRUE)
    m0 <- mean(vv[g == 0], na.rm = TRUE); sd0 <- sd(vv[g == 0], na.rm = TRUE)
    n1 <- sum(!is.na(vv[g == 1])); n0 <- sum(!is.na(vv[g == 0]))
    out[[v]] <- tibble(
      variable = v,
      mean_analytic = m1, sd_analytic = sd1, n_analytic = n1,
      mean_ltfu     = m0, sd_ltfu     = sd0, n_ltfu     = n0,
      smd = smd(vv, g),
      abs_smd = abs(smd(vv, g))
    )
  }

  # Also stratify by cohort
  cohort_smd <- eligible %>%
    group_by(cohort) %>%
    summarise(
      n_analytic = sum(in_analysis == 1),
      n_ltfu     = sum(in_analysis == 0),
      pct_ltfu   = round(100 * n_ltfu / (n_analytic + n_ltfu), 1),
      age_smd    = smd(age, in_analysis),
      ic_smd     = smd(ic_cfa, in_analysis),
      env_smd    = smd(env_formative, in_analysis),
      chr_smd    = smd(chronic_n, in_analysis),
      .groups = "drop"
    )

  overall <- bind_rows(out)
  write.csv(overall, file.path(RESULTS_DIR, "ipcw_baseline_compare.csv"), row.names = FALSE)
  write.csv(cohort_smd, file.path(RESULTS_DIR, "ipcw_ltfu_by_cohort.csv"), row.names = FALSE)

  cat("  --- Overall baseline SMDs (|SMD| > 0.1 = material imbalance) ---\n")
  print(overall %>% select(variable, mean_analytic, mean_ltfu, smd, abs_smd) %>%
          mutate(across(where(is.numeric), ~ round(., 3))))
  cat("\n  --- LTFU rate by cohort ---\n")
  print(cohort_smd)

  list(overall = overall, by_cohort = cohort_smd)
}

# ============================================================================
# STEP 3: Fit IPCW weights (logistic: P(retained | baseline covariates))
# ============================================================================

step3_ipcw_weights <- function(eligible) {
  cat("\n", strrep("=", 70), "\n"); cat("Step 3: IPCW weight estimation\n"); cat(strrep("=", 70), "\n\n")

  # Use baseline variables available across all cohorts.
  # Avoid variables that themselves depend on retention (e.g. future IC/Env).
  model_data <- eligible %>%
    filter(!is.na(age) & !is.na(female) & !is.na(chronic_n)) %>%
    mutate(
      age_c    = (age - 60) / 10,             # centered, decadal
      chronic  = as.numeric(chronic_n),
      retained = in_analysis
    )

  # Stratified by cohort (each cohort may have different LTFU mechanism)
  cohorts <- sort(unique(model_data$cohort))
  weight_rows <- list()
  diag_rows   <- list()

  for (coh in cohorts) {
    d <- model_data %>% filter(cohort == coh)
    # Simple model; using ic_cfa and env_formative if available adds information
    # but risks conditioning on post-baseline (survivor) pattern; use only age, female, chronic.
    fit <- glm(retained ~ age_c + female + chronic, data = d, family = binomial())
    d$p_ret   <- predict(fit, newdata = d, type = "response")
    d$p_ret   <- pmax(pmin(d$p_ret, 0.999), 0.01)
    d$w_ipcw  <- 1 / d$p_ret
    # Stabilized weights: numerator = marginal P(retained in cohort)
    p_marg <- mean(d$retained)
    d$sw   <- p_marg / d$p_ret
    # Truncation at 1st/99th percentile of sw for analytic participants
    d_analytic <- d %>% filter(retained == 1)
    lo <- quantile(d_analytic$sw, 0.01, na.rm = TRUE)
    hi <- quantile(d_analytic$sw, 0.99, na.rm = TRUE)
    d$sw_tr <- pmin(pmax(d$sw, lo), hi)

    diag_rows[[coh]] <- tibble(
      cohort = coh,
      n = nrow(d),
      p_ret_mean = mean(d$p_ret), p_ret_min = min(d$p_ret), p_ret_max = max(d$p_ret),
      sw_mean = mean(d_analytic$sw), sw_min = min(d_analytic$sw), sw_max = max(d_analytic$sw),
      sw_tr_mean = mean(d_analytic$sw_tr), sw_tr_max = max(d_analytic$sw_tr)
    )
    weight_rows[[coh]] <- d %>% filter(retained == 1) %>% select(cohort, id, sw, sw_tr)
  }

  weights_df <- bind_rows(weight_rows)
  diag_df    <- bind_rows(diag_rows)

  write.csv(diag_df, file.path(RESULTS_DIR, "ipcw_weights_diagnostics.csv"), row.names = FALSE)

  cat("  --- IPCW weight diagnostics ---\n")
  print(diag_df %>% mutate(across(where(is.numeric), ~ round(., 3))))

  weights_df
}

# ============================================================================
# STEP 4: IPCW-weighted M4 Cox + REML pooled
# ============================================================================

step4_ipcw_cox <- function(weights_df) {
  cat("\n", strrep("=", 70), "\n"); cat("Step 4: IPCW-weighted M4 Cox + pooled\n"); cat(strrep("=", 70), "\n\n")

  surv <- readRDS(file.path(DATA_DIR, "survival_cohort.rds"))
  surv <- surv %>%
    inner_join(weights_df, by = c("cohort", "id"))

  cat(sprintf("  Survival cohort with IPCW weights: %s\n", format(nrow(surv), big.mark=",")))

  cohorts <- sort(unique(surv$cohort))
  rows <- list()
  for (coh in cohorts) {
    d <- surv %>%
      filter(cohort == coh &
             !is.na(ic_cfa) & !is.na(env_formative) & !is.na(ic_env_cont) &
             !is.na(age) & !is.na(female) & !is.na(chronic_n))
    if (nrow(d) < 100 || sum(d$event_dfs) < 20) next

    fit <- coxph(Surv(time_dfs, event_dfs) ~ ic_cfa + env_formative + ic_env_cont +
                   age + female + chronic_n,
                 data = d, weights = sw_tr, robust = TRUE)
    s <- summary(fit)$coefficients %>% as.data.frame() %>%
      tibble::rownames_to_column("variable")
    names(s)[names(s) == "coef"]    <- "log_hr"
    # With weights=, coxph returns "robust se" in column "robust se"
    se_col <- intersect(c("robust se", "se(coef)"), names(s))[1]
    names(s)[names(s) == se_col]    <- "se_log_hr"
    names(s)[names(s) == "exp(coef)"] <- "hr"
    p_col <- intersect(c("Pr(>|z|)"), names(s))[1]
    names(s)[names(s) == p_col]     <- "p_value"

    s <- s %>% transmute(
      cohort = coh, model = "M4_IPCW",
      variable, log_hr, se_log_hr, hr, p_value,
      hr_lower = exp(log_hr - 1.96 * se_log_hr),
      hr_upper = exp(log_hr + 1.96 * se_log_hr),
      n = nrow(d), n_events = sum(d$event_dfs)
    )
    rows[[coh]] <- s
  }
  stage1 <- bind_rows(rows)
  write.csv(stage1, file.path(RESULTS_DIR, "ipcw_stage1.csv"), row.names = FALSE)

  int <- stage1 %>% filter(variable == "ic_env_cont")
  cat("\n  --- Per-cohort IPCW-weighted IC x Env HR ---\n")
  print(int %>% select(cohort, hr, hr_lower, hr_upper, p_value) %>%
          mutate(across(where(is.numeric), ~ round(., 3))))

  # REML pooled
  if (nrow(int) >= 2) {
    ma <- rma(yi = log_hr, sei = se_log_hr, data = int, method = "REML")
    pooled <- tibble(
      outcome = "DFS", model = "M4_IPCW_pooled",
      hr = exp(ma$b[1,1]),
      hr_lower = exp(ma$ci.lb), hr_upper = exp(ma$ci.ub),
      p_value = ma$pval, i2 = ma$I2, tau2 = ma$tau2, k = ma$k
    )
    cat(sprintf("\n  *** IPCW-weighted pooled IC x Env HR = %.3f (%.3f-%.3f), p = %.4f, I2 = %.1f%% ***\n",
                pooled$hr, pooled$hr_lower, pooled$hr_upper, pooled$p_value, pooled$i2))
    write.csv(pooled, file.path(RESULTS_DIR, "ipcw_stage2_pooled.csv"), row.names = FALSE)
  }

  stage1
}

# ============================================================================
# MAIN
# ============================================================================

main <- function() {
  t0 <- Sys.time()
  elig <- step1_baseline_flag()
  smd_res <- step2_smd(elig)
  w_df    <- step3_ipcw_weights(elig)
  step4_ipcw_cox(w_df)
  cat(sprintf("\nTotal runtime: %.1f min\n",
              as.numeric(Sys.time() - t0, units="mins")))
}

if (sys.nframe() == 0) main()
