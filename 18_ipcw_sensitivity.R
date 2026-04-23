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
# Dependencies: tidyverse, survival, metafor, here
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

# --- SMD helper ---
smd <- function(x, g) {
  m1 <- mean(x[g == 1], na.rm = TRUE); m0 <- mean(x[g == 0], na.rm = TRUE)
  s1 <- var(x[g == 1], na.rm = TRUE);  s0 <- var(x[g == 0], na.rm = TRUE)
  sp <- sqrt((s1 + s0) / 2)
  if (!is.finite(sp) || sp == 0) NA_real_ else (m1 - m0) / sp
}

# STEP 1: Reconstruct eligible baseline + LTFU flag
step1_baseline_flag <- function() {
  pooled  <- readRDS(file.path(DATA_DIR, "analytic_final.rds"))
  surv    <- readRDS(file.path(DATA_DIR, "survival_cohort.rds"))
  baseline <- pooled %>%
    group_by(cohort, id) %>% slice_min(wave, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    rename(baseline_wave = wave, baseline_year = interview_year)
  eligible <- baseline %>% filter(!is.na(adl_total) & adl_total == 0)
  surv_keys <- surv %>% select(cohort, id) %>% distinct() %>% mutate(in_analysis = 1L)
  eligible <- eligible %>%
    left_join(surv_keys, by = c("cohort", "id")) %>%
    mutate(in_analysis = ifelse(is.na(in_analysis), 0L, 1L))
  cat(sprintf("  Analytic (retained): %s / LTFU (excluded): %s / LTFU rate: %.1f%%\n",
              format(sum(eligible$in_analysis == 1), big.mark=","),
              format(sum(eligible$in_analysis == 0), big.mark=","),
              100 * sum(eligible$in_analysis == 0) / nrow(eligible)))
  eligible
}

# STEP 2: Baseline comparison (LTFU vs analytic)
step2_smd <- function(eligible) {
  vars_cont <- c("age", "chronic_n", "ic_cfa", "env_formative", "edu_years")
  vars_bin  <- c("female")
  out <- list()
  for (v in c(vars_cont, vars_bin)) {
    if (!v %in% names(eligible)) next
    vv <- eligible[[v]]; g <- eligible$in_analysis
    out[[v]] <- tibble(
      variable = v,
      mean_analytic = mean(vv[g == 1], na.rm = TRUE),
      sd_analytic = sd(vv[g == 1], na.rm = TRUE),
      n_analytic = sum(!is.na(vv[g == 1])),
      mean_ltfu = mean(vv[g == 0], na.rm = TRUE),
      sd_ltfu = sd(vv[g == 0], na.rm = TRUE),
      n_ltfu = sum(!is.na(vv[g == 0])),
      smd = smd(vv, g), abs_smd = abs(smd(vv, g))
    )
  }
  cohort_smd <- eligible %>%
    group_by(cohort) %>%
    summarise(
      n_analytic = sum(in_analysis == 1),
      n_ltfu = sum(in_analysis == 0),
      pct_ltfu = round(100 * n_ltfu / (n_analytic + n_ltfu), 1),
      age_smd = smd(age, in_analysis),
      ic_smd = smd(ic_cfa, in_analysis),
      env_smd = smd(env_formative, in_analysis),
      chr_smd = smd(chronic_n, in_analysis),
      .groups = "drop")
  overall <- bind_rows(out)
  write.csv(overall, file.path(RESULTS_DIR, "ipcw_baseline_compare.csv"), row.names = FALSE)
  write.csv(cohort_smd, file.path(RESULTS_DIR, "ipcw_ltfu_by_cohort.csv"), row.names = FALSE)
  list(overall = overall, by_cohort = cohort_smd)
}

# STEP 3: IPCW weights (cohort-stratified logistic)
step3_ipcw_weights <- function(eligible) {
  model_data <- eligible %>%
    filter(!is.na(age) & !is.na(female) & !is.na(chronic_n)) %>%
    mutate(age_c = (age - 60) / 10,
           chronic = as.numeric(chronic_n),
           retained = in_analysis)
  cohorts <- sort(unique(model_data$cohort))
  weight_rows <- list(); diag_rows <- list()
  for (coh in cohorts) {
    d <- model_data %>% filter(cohort == coh)
    fit <- glm(retained ~ age_c + female + chronic, data = d, family = binomial())
    d$p_ret <- predict(fit, newdata = d, type = "response")
    d$p_ret <- pmax(pmin(d$p_ret, 0.999), 0.01)
    d$w_ipcw <- 1 / d$p_ret
    p_marg <- mean(d$retained)
    d$sw <- p_marg / d$p_ret
    d_analytic <- d %>% filter(retained == 1)
    lo <- quantile(d_analytic$sw, 0.01, na.rm = TRUE)
    hi <- quantile(d_analytic$sw, 0.99, na.rm = TRUE)
    d$sw_tr <- pmin(pmax(d$sw, lo), hi)
    diag_rows[[coh]] <- tibble(
      cohort = coh, n = nrow(d),
      p_ret_mean = mean(d$p_ret), p_ret_min = min(d$p_ret), p_ret_max = max(d$p_ret),
      sw_mean = mean(d_analytic$sw), sw_max = max(d_analytic$sw),
      sw_tr_max = max(d_analytic$sw_tr)
    )
    weight_rows[[coh]] <- d %>% filter(retained == 1) %>%
      select(cohort, id, sw, sw_tr)
  }
  weights_df <- bind_rows(weight_rows)
  diag_df <- bind_rows(diag_rows)
  write.csv(diag_df, file.path(RESULTS_DIR, "ipcw_weights_diagnostics.csv"), row.names = FALSE)
  weights_df
}

# STEP 4: IPCW-weighted M4 Cox + REML pooled
step4_ipcw_cox <- function(weights_df) {
  surv <- readRDS(file.path(DATA_DIR, "survival_cohort.rds")) %>%
    inner_join(weights_df, by = c("cohort", "id"))
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
    names(s)[names(s) == "coef"] <- "log_hr"
    se_col <- intersect(c("robust se", "se(coef)"), names(s))[1]
    names(s)[names(s) == se_col] <- "se_log_hr"
    names(s)[names(s) == "exp(coef)"] <- "hr"
    p_col <- intersect(c("Pr(>|z|)"), names(s))[1]
    names(s)[names(s) == p_col] <- "p_value"
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
  if (nrow(int) >= 2) {
    ma <- rma(yi = log_hr, sei = se_log_hr, data = int, method = "REML")
    pooled <- tibble(
      outcome = "DFS", model = "M4_IPCW_pooled",
      hr = exp(ma$b[1,1]),
      hr_lower = exp(ma$ci.lb), hr_upper = exp(ma$ci.ub),
      p_value = ma$pval, i2 = ma$I2, tau2 = ma$tau2, k = ma$k
    )
    cat(sprintf("  IPCW-weighted pooled IC x Env HR = %.3f (%.3f-%.3f), p = %.4f\n",
                pooled$hr, pooled$hr_lower, pooled$hr_upper, pooled$p_value))
    write.csv(pooled, file.path(RESULTS_DIR, "ipcw_stage2_pooled.csv"), row.names = FALSE)
  }
  stage1
}

main <- function() {
  t0 <- Sys.time()
  elig <- step1_baseline_flag()
  step2_smd(elig)
  w_df <- step3_ipcw_weights(elig)
  step4_ipcw_cox(w_df)
  cat(sprintf("\nTotal runtime: %.1f min\n",
              as.numeric(Sys.time() - t0, units="mins")))
}

if (sys.nframe() == 0) main()
