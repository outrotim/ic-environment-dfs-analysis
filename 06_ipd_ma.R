# ==============================================================================
# 06_ipd_ma.R
# Study 17: IC-Environment Fit and Disability-Free Survival in Older Adults
#
# Purpose: Two-stage IPD Meta-Analysis
#   Stage 1: Cohort-specific Cox proportional hazards models
#   Stage 2: Random-effects meta-analysis (DerSimonian-Laird / REML)
#
# Input:   data/analytic_final.rds (from 04_env_index.R)
# Output:  results/stage1_estimates.csv  — cohort-specific HRs
#          results/stage2_pooled.csv     — pooled HRs with I²
#          results/ipd_ma_summary.csv    — analysis summary
#          data/survival_cohort.rds      — constructed survival dataset
#
# Dependencies: tidyverse, survival, metafor
# Author:  Study 17 Team
# Date:    2026-03-15
# Version: 1.0
# ==============================================================================

library(tidyverse)
library(survival)

if (!requireNamespace("metafor", quietly = TRUE)) {
  install.packages("metafor", repos = "https://cloud.r-project.org")
}
library(metafor)

SCRIPT_DIR <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) {
    here::here()
  }
)
RESULTS_DIR <- file.path(SCRIPT_DIR, "results")
DATA_DIR    <- file.path(SCRIPT_DIR, "data")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# STEP 1: Construct survival dataset from panel data
# ==============================================================================

step1_survival <- function() {
  cat("=" %>% strrep(70), "\n")
  cat("Step 1: Constructing survival dataset\n")
  cat("=" %>% strrep(70), "\n\n")

  pooled <- readRDS(file.path(DATA_DIR, "analytic_final.rds"))
  cat("  Full dataset:", format(nrow(pooled), big.mark = ","), "rows,",
      n_distinct(paste(pooled$cohort, pooled$id)), "persons\n")

  # --- 1a: Get baseline (first wave per person per cohort) ---
  baseline <- pooled %>%
    group_by(cohort, id) %>%
    slice_min(wave, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(baseline_wave = wave, baseline_year = interview_year)

  cat("  Baseline sample:", format(nrow(baseline), big.mark = ","), "\n")

  # --- 1b: Exclude baseline ADL >= 1 (not at risk for incident disability) ---
  n_adl_ge1 <- sum(baseline$adl_total >= 1, na.rm = TRUE)
  n_adl_na  <- sum(is.na(baseline$adl_total))

  eligible <- baseline %>%
    filter(!is.na(adl_total) & adl_total == 0)

  cat("  Excluded: ADL >= 1 at baseline =", format(n_adl_ge1, big.mark = ","),
      ", ADL NA =", n_adl_na, "\n")
  cat("  Eligible (ADL = 0 at baseline):", format(nrow(eligible), big.mark = ","), "\n")

  # --- 1c: Find first post-baseline ADL disability ---
  disability_onset <- pooled %>%
    inner_join(eligible %>% select(cohort, id, baseline_wave),
               by = c("cohort", "id")) %>%
    filter(wave > baseline_wave & !is.na(adl_total) & adl_total >= 1) %>%
    group_by(cohort, id) %>%
    slice_min(wave, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(cohort, id, disability_year = interview_year)

  cat("  Incident ADL disability detected:", format(nrow(disability_onset), big.mark = ","), "persons\n")

  # --- 1d: Find last interview year per person ---
  last_obs <- pooled %>%
    semi_join(eligible, by = c("cohort", "id")) %>%
    group_by(cohort, id) %>%
    summarise(last_interview_year = max(interview_year, na.rm = TRUE),
              n_waves = n(),
              .groups = "drop")

  # --- 1e: Construct survival variables ---
  surv <- eligible %>%
    left_join(disability_onset, by = c("cohort", "id")) %>%
    left_join(last_obs, by = c("cohort", "id")) %>%
    mutate(
      # Time to each event type (NA → Inf for min-comparison)
      t_disab = ifelse(!is.na(disability_year), disability_year - baseline_year, Inf),
      t_death = ifelse(person_died == 1 & !is.na(death_year) & death_year > baseline_year,
                       death_year - baseline_year, Inf),
      t_censor = last_interview_year - baseline_year,

      # === DFS outcome: first of disability, death, or censoring ===
      time_dfs_raw = pmin(t_disab, t_death, t_censor, na.rm = TRUE),
      event_dfs = as.integer(time_dfs_raw < Inf &
                               (time_dfs_raw == t_disab | time_dfs_raw == t_death)),
      # Ensure censored persons use actual follow-up time
      time_dfs = ifelse(event_dfs == 1, time_dfs_raw, t_censor),

      # === Mortality outcome ===
      time_mort_raw = pmin(t_death, t_censor, na.rm = TRUE),
      event_mort = as.integer(t_death < Inf & time_mort_raw == t_death),
      time_mort = ifelse(event_mort == 1, time_mort_raw, t_censor)
    )

  # --- 1f: Exclude persons with zero follow-up ---
  n_zero <- sum(surv$time_dfs <= 0 | is.na(surv$time_dfs))
  surv_clean <- surv %>%
    filter(time_dfs > 0 & !is.na(time_dfs))

  cat("  Excluded: zero/NA follow-up =", n_zero, "\n")
  cat("  Final analytic cohort:", format(nrow(surv_clean), big.mark = ","), "\n\n")

  # --- 1g: Summary ---
  cat("  DFS events:", sum(surv_clean$event_dfs),
      "(", round(mean(surv_clean$event_dfs)*100, 1), "%)\n")
  cat("    - Disability first:", sum(surv_clean$event_dfs == 1 & surv_clean$t_disab <= surv_clean$t_death), "\n")
  cat("    - Death first:", sum(surv_clean$event_dfs == 1 & surv_clean$t_death < surv_clean$t_disab), "\n")
  cat("  Mortality events:", sum(surv_clean$event_mort),
      "(", round(mean(surv_clean$event_mort)*100, 1), "%)\n")
  cat("  Median follow-up:", round(median(surv_clean$time_dfs), 1), "years\n")

  # Per-cohort summary
  cat("\n  Per-cohort breakdown:\n")
  cohort_sum <- surv_clean %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      dfs_events = sum(event_dfs),
      dfs_rate = round(mean(event_dfs)*100, 1),
      mort_events = sum(event_mort),
      mort_rate = round(mean(event_mort)*100, 1),
      median_fu = round(median(time_dfs), 1),
      .groups = "drop"
    )
  print(cohort_sum)

  # --- 1h: Ensure exposure completeness ---
  cat("\n  Exposure completeness in analytic cohort:\n")
  cat("    ic_cfa:", round(mean(!is.na(surv_clean$ic_cfa))*100, 1), "%\n")
  cat("    env_formative:", round(mean(!is.na(surv_clean$env_formative))*100, 1), "%\n")
  cat("    ic_env_cont:", round(mean(!is.na(surv_clean$ic_env_cont))*100, 1), "%\n")
  cat("    ic_env_6grp:", round(mean(!is.na(surv_clean$ic_env_6grp))*100, 1), "%\n")

  # Save
  saveRDS(surv_clean, file.path(DATA_DIR, "survival_cohort.rds"))
  cat("\n  Saved: data/survival_cohort.rds\n")

  surv_clean
}

# ==============================================================================
# STEP 2: Stage 1 — Cohort-specific Cox models
# ==============================================================================

step2_stage1 <- function(surv) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 2: Stage 1 — Cohort-specific Cox models\n")
  cat("=" %>% strrep(70), "\n\n")

  cohorts <- sort(unique(surv$cohort))
  results_list <- list()

  # Model definitions (primary covariates: age + female + chronic_n)
  # CHARLS lacks edu_years, so primary model excludes it for all cohorts
  models <- list(
    M1_IC     = event_dfs ~ ic_cfa + age + female + chronic_n,
    M2_Env    = event_dfs ~ env_formative + age + female + chronic_n,
    M3_Both   = event_dfs ~ ic_cfa + env_formative + age + female + chronic_n,
    M4_Interaction = event_dfs ~ ic_cfa + env_formative + ic_env_cont + age + female + chronic_n,
    M5_6group = event_dfs ~ ic_env_6grp + age + female + chronic_n
  )

  for (coh in cohorts) {
    cat("  --- ", coh, " ---\n")
    dat <- surv %>%
      filter(cohort == coh) %>%
      # Complete cases for key variables
      filter(!is.na(ic_cfa) & !is.na(env_formative) &
             !is.na(ic_env_cont) & !is.na(ic_env_6grp) &
             !is.na(age) & !is.na(female) & !is.na(chronic_n))

    cat("    N =", nrow(dat), ", events =", sum(dat$event_dfs), "\n")

    if (nrow(dat) < 50 || sum(dat$event_dfs) < 20) {
      cat("    SKIPPED: insufficient sample/events\n")
      next
    }

    for (model_name in names(models)) {
      tryCatch({
        # Build Surv object
        formula_rhs <- as.character(models[[model_name]])[3]
        full_formula <- as.formula(paste0("Surv(time_dfs, event_dfs) ~ ", formula_rhs))

        fit <- coxph(full_formula, data = dat)

        # Extract coefficients
        coef_df <- as.data.frame(summary(fit)$coefficients)
        coef_df$variable <- rownames(coef_df)
        coef_df$model <- model_name
        coef_df$cohort <- coh
        coef_df$n <- nrow(dat)
        coef_df$n_events <- sum(dat$event_dfs)

        # Rename columns for consistency
        names(coef_df)[names(coef_df) == "coef"] <- "log_hr"
        names(coef_df)[names(coef_df) == "se(coef)"] <- "se_log_hr"
        names(coef_df)[names(coef_df) == "exp(coef)"] <- "hr"
        names(coef_df)[names(coef_df) == "Pr(>|z|)"] <- "p_value"

        # Add 95% CI
        coef_df$hr_lower <- exp(coef_df$log_hr - 1.96 * coef_df$se_log_hr)
        coef_df$hr_upper <- exp(coef_df$log_hr + 1.96 * coef_df$se_log_hr)

        results_list[[paste(coh, model_name, sep = "_")]] <- coef_df

      }, error = function(e) {
        cat("    WARNING:", model_name, "failed:", conditionMessage(e), "\n")
      })
    }
  }

  # Combine all results
  stage1 <- bind_rows(results_list)

  # Clean up variable names for 6-group model
  stage1 <- stage1 %>%
    mutate(variable = str_replace(variable, "ic_env_6grp", ""))

  # Save Stage 1 results
  write_csv(stage1, file.path(RESULTS_DIR, "stage1_estimates.csv"))
  cat("\n  Stage 1 total estimates:", nrow(stage1), "\n")
  cat("  Saved: results/stage1_estimates.csv\n")

  # Print key results
  cat("\n  Key exposure effects (DFS):\n")
  key_vars <- c("ic_cfa", "env_formative", "ic_env_cont")
  for (v in key_vars) {
    sub <- stage1 %>% filter(variable == v)
    if (nrow(sub) > 0) {
      cat("    ", v, ":\n")
      for (i in seq_len(nrow(sub))) {
        cat("      ", sub$cohort[i], ": HR =", round(sub$hr[i], 3),
            "(", round(sub$hr_lower[i], 3), "-", round(sub$hr_upper[i], 3), ")",
            "p =", format(sub$p_value[i], digits = 3), "\n")
      }
    }
  }

  stage1
}

# ==============================================================================
# STEP 3: Stage 2 — Random-effects meta-analysis
# ==============================================================================

step3_stage2 <- function(stage1) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 3: Stage 2 — Random-effects meta-analysis\n")
  cat("=" %>% strrep(70), "\n\n")

  # Identify exposure variables to pool (from models M1-M5)
  # Pool each variable within each model

  pool_targets <- stage1 %>%
    filter(!(variable %in% c("age", "female", "chronic_n"))) %>%
    distinct(model, variable) %>%
    arrange(model, variable)

  cat("  Pooling", nrow(pool_targets), "variable-model combinations\n\n")

  pooled_list <- list()

  for (i in seq_len(nrow(pool_targets))) {
    mod <- pool_targets$model[i]
    var <- pool_targets$variable[i]

    sub <- stage1 %>%
      filter(model == mod & variable == var & !is.na(log_hr) & !is.na(se_log_hr))

    if (nrow(sub) < 2) {
      cat("    SKIP:", mod, "/", var, "— fewer than 2 cohorts\n")
      next
    }

    tryCatch({
      # Random-effects meta-analysis (REML)
      re_fit <- rma(yi = log_hr, sei = se_log_hr, data = sub, method = "REML")

      pooled_list[[paste(mod, var, sep = "___")]] <- tibble(
        model        = mod,
        variable     = var,
        k            = re_fit$k,
        pooled_loghr = as.numeric(re_fit$beta),
        pooled_se    = as.numeric(re_fit$se),
        pooled_hr    = exp(as.numeric(re_fit$beta)),
        hr_lower     = exp(as.numeric(re_fit$ci.lb)),
        hr_upper     = exp(as.numeric(re_fit$ci.ub)),
        p_value      = as.numeric(re_fit$pval),
        tau2         = re_fit$tau2,
        I2           = re_fit$I2,
        Q            = re_fit$QE,
        Q_p          = re_fit$QEp
      )

    }, error = function(e) {
      cat("    WARNING:", mod, "/", var, "meta-analysis failed:", conditionMessage(e), "\n")
    })
  }

  pooled <- bind_rows(pooled_list)

  # Save
  write_csv(pooled, file.path(RESULTS_DIR, "stage2_pooled.csv"))
  cat("  Pooled estimates:", nrow(pooled), "\n")
  cat("  Saved: results/stage2_pooled.csv\n\n")

  # === Print key results ===
  cat("  ┌─────────────────────────────────────────────────────────────┐\n")
  cat("  │                     KEY RESULTS (DFS)                      │\n")
  cat("  └─────────────────────────────────────────────────────────────┘\n\n")

  # M1: IC main effect
  m1_ic <- pooled %>% filter(model == "M1_IC" & variable == "ic_cfa")
  if (nrow(m1_ic) > 0) {
    cat("  M1 — IC main effect:\n")
    cat("    HR =", round(m1_ic$pooled_hr, 3),
        "(", round(m1_ic$hr_lower, 3), "-", round(m1_ic$hr_upper, 3), ")",
        "p =", format(m1_ic$p_value, digits = 3),
        "I² =", round(m1_ic$I2, 1), "%\n")
    cat("    → Per 1-SD increase in IC: ", round((1-m1_ic$pooled_hr)*100, 1),
        "% ", ifelse(m1_ic$pooled_hr < 1, "lower", "higher"), " DFS risk\n\n")
  }

  # M2: Env main effect
  m2_env <- pooled %>% filter(model == "M2_Env" & variable == "env_formative")
  if (nrow(m2_env) > 0) {
    cat("  M2 — Env main effect:\n")
    cat("    HR =", round(m2_env$pooled_hr, 3),
        "(", round(m2_env$hr_lower, 3), "-", round(m2_env$hr_upper, 3), ")",
        "p =", format(m2_env$p_value, digits = 3),
        "I² =", round(m2_env$I2, 1), "%\n\n")
  }

  # M4: Interaction
  m4_int <- pooled %>% filter(model == "M4_Interaction" & variable == "ic_env_cont")
  if (nrow(m4_int) > 0) {
    cat("  M4 — IC×Env interaction:\n")
    cat("    HR =", round(m4_int$pooled_hr, 3),
        "(", round(m4_int$hr_lower, 3), "-", round(m4_int$hr_upper, 3), ")",
        "p =", format(m4_int$p_value, digits = 3),
        "I² =", round(m4_int$I2, 1), "%\n")
    if (m4_int$p_value < 0.05) {
      cat("    ★ SIGNIFICANT multiplicative interaction\n\n")
    } else {
      cat("    → Non-significant on multiplicative scale (check additive)\n\n")
    }
  }

  # M5: 6-group
  cat("  M5 — IC×Env 6-group (ref: High IC / High Env):\n")
  m5 <- pooled %>% filter(model == "M5_6group") %>% arrange(variable)
  if (nrow(m5) > 0) {
    for (j in seq_len(nrow(m5))) {
      sig <- ifelse(m5$p_value[j] < 0.001, "***",
             ifelse(m5$p_value[j] < 0.01, "**",
             ifelse(m5$p_value[j] < 0.05, "*", "")))
      cat("    ", sprintf("%-25s", m5$variable[j]),
          "HR =", sprintf("%.3f", m5$pooled_hr[j]),
          "(", sprintf("%.3f", m5$hr_lower[j]), "-", sprintf("%.3f", m5$hr_upper[j]), ")",
          sig, "\n")
    }
    cat("\n")
  }

  pooled
}

# ==============================================================================
# STEP 4: Secondary outcomes — Mortality only
# ==============================================================================

step4_mortality <- function(surv) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 4: Secondary outcome — Mortality\n")
  cat("=" %>% strrep(70), "\n\n")

  cohorts <- sort(unique(surv$cohort))
  results_list <- list()

  models <- list(
    M1_IC_mort     = event_mort ~ ic_cfa + age + female + chronic_n,
    M3_Both_mort   = event_mort ~ ic_cfa + env_formative + age + female + chronic_n,
    M4_Int_mort    = event_mort ~ ic_cfa + env_formative + ic_env_cont + age + female + chronic_n,
    M5_6grp_mort   = event_mort ~ ic_env_6grp + age + female + chronic_n
  )

  for (coh in cohorts) {
    dat <- surv %>%
      filter(cohort == coh) %>%
      filter(!is.na(ic_cfa) & !is.na(env_formative) &
             !is.na(ic_env_cont) & !is.na(ic_env_6grp) &
             !is.na(age) & !is.na(female) & !is.na(chronic_n))

    if (nrow(dat) < 50 || sum(dat$event_mort) < 20) next

    for (model_name in names(models)) {
      tryCatch({
        formula_rhs <- as.character(models[[model_name]])[3]
        full_formula <- as.formula(paste0("Surv(time_mort, event_mort) ~ ", formula_rhs))
        fit <- coxph(full_formula, data = dat)

        coef_df <- as.data.frame(summary(fit)$coefficients)
        coef_df$variable <- rownames(coef_df)
        coef_df$model <- model_name
        coef_df$cohort <- coh
        coef_df$n <- nrow(dat)
        coef_df$n_events <- sum(dat$event_mort)
        names(coef_df)[names(coef_df) == "coef"] <- "log_hr"
        names(coef_df)[names(coef_df) == "se(coef)"] <- "se_log_hr"
        names(coef_df)[names(coef_df) == "exp(coef)"] <- "hr"
        names(coef_df)[names(coef_df) == "Pr(>|z|)"] <- "p_value"
        coef_df$hr_lower <- exp(coef_df$log_hr - 1.96 * coef_df$se_log_hr)
        coef_df$hr_upper <- exp(coef_df$log_hr + 1.96 * coef_df$se_log_hr)

        results_list[[paste(coh, model_name, sep = "_")]] <- coef_df
      }, error = function(e) NULL)
    }
  }

  stage1_mort <- bind_rows(results_list) %>%
    mutate(variable = str_replace(variable, "ic_env_6grp", ""))

  # Pool mortality results
  pool_targets <- stage1_mort %>%
    filter(!(variable %in% c("age", "female", "chronic_n"))) %>%
    distinct(model, variable)

  pooled_mort <- list()
  for (i in seq_len(nrow(pool_targets))) {
    mod <- pool_targets$model[i]
    var <- pool_targets$variable[i]
    sub <- stage1_mort %>% filter(model == mod & variable == var)
    if (nrow(sub) < 2) next
    tryCatch({
      re_fit <- rma(yi = log_hr, sei = se_log_hr, data = sub, method = "REML")
      pooled_mort[[paste(mod, var)]] <- tibble(
        model = mod, variable = var, k = re_fit$k,
        pooled_hr = exp(as.numeric(re_fit$beta)),
        hr_lower = exp(as.numeric(re_fit$ci.lb)),
        hr_upper = exp(as.numeric(re_fit$ci.ub)),
        p_value = as.numeric(re_fit$pval),
        I2 = re_fit$I2
      )
    }, error = function(e) NULL)
  }

  pooled_mort_df <- bind_rows(pooled_mort)

  cat("  Key mortality results:\n")
  # IC effect on mortality
  m1m <- pooled_mort_df %>% filter(model == "M1_IC_mort" & variable == "ic_cfa")
  if (nrow(m1m) > 0) {
    cat("    IC → Mortality: HR =", round(m1m$pooled_hr, 3),
        "(", round(m1m$hr_lower, 3), "-", round(m1m$hr_upper, 3), ")\n")
  }
  # Interaction on mortality
  m4m <- pooled_mort_df %>% filter(model == "M4_Int_mort" & variable == "ic_env_cont")
  if (nrow(m4m) > 0) {
    cat("    IC×Env → Mortality: HR =", round(m4m$pooled_hr, 3),
        "(", round(m4m$hr_lower, 3), "-", round(m4m$hr_upper, 3), ")",
        "p =", format(m4m$p_value, digits = 3), "\n")
  }

  # Append to stage1 and pooled CSVs
  write_csv(stage1_mort, file.path(RESULTS_DIR, "stage1_mortality.csv"))
  write_csv(pooled_mort_df, file.path(RESULTS_DIR, "stage2_mortality.csv"))
  cat("  Saved: results/stage1_mortality.csv, results/stage2_mortality.csv\n")

  pooled_mort_df
}

# ==============================================================================
# STEP 5: Subgroup analyses (sex, age)
# ==============================================================================

step5_subgroups <- function(surv) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 5: Subgroup analyses\n")
  cat("=" %>% strrep(70), "\n\n")

  cohorts <- sort(unique(surv$cohort))
  subgroup_results <- list()

  # --- 5a: By sex ---
  for (sex_val in c(0, 1)) {
    sex_label <- ifelse(sex_val == 0, "Male", "Female")
    cat("  Subgroup: ", sex_label, "\n")

    for (coh in cohorts) {
      dat <- surv %>%
        filter(cohort == coh & female == sex_val) %>%
        filter(!is.na(ic_cfa) & !is.na(env_formative) &
               !is.na(ic_env_cont) & !is.na(age) & !is.na(chronic_n))

      if (nrow(dat) < 50 || sum(dat$event_dfs) < 20) next

      tryCatch({
        fit <- coxph(Surv(time_dfs, event_dfs) ~
                       ic_cfa + env_formative + ic_env_cont + age + chronic_n,
                     data = dat)
        coef_df <- as.data.frame(summary(fit)$coefficients)
        coef_df$variable <- rownames(coef_df)
        coef_df$cohort <- coh
        coef_df$subgroup <- sex_label
        coef_df$subgroup_type <- "sex"
        coef_df$n <- nrow(dat)
        coef_df$n_events <- sum(dat$event_dfs)
        names(coef_df)[names(coef_df) == "coef"] <- "log_hr"
        names(coef_df)[names(coef_df) == "se(coef)"] <- "se_log_hr"
        names(coef_df)[names(coef_df) == "exp(coef)"] <- "hr"
        names(coef_df)[names(coef_df) == "Pr(>|z|)"] <- "p_value"
        coef_df$hr_lower <- exp(coef_df$log_hr - 1.96 * coef_df$se_log_hr)
        coef_df$hr_upper <- exp(coef_df$log_hr + 1.96 * coef_df$se_log_hr)
        subgroup_results[[paste(sex_label, coh)]] <- coef_df
      }, error = function(e) NULL)
    }
  }

  # --- 5b: By age group ---
  for (age_grp in c("younger", "older")) {
    age_label <- ifelse(age_grp == "younger", "Age < 70", "Age >= 70")
    cat("  Subgroup: ", age_label, "\n")

    for (coh in cohorts) {
      dat <- surv %>%
        filter(cohort == coh) %>%
        filter(if (age_grp == "younger") age < 70 else age >= 70) %>%
        filter(!is.na(ic_cfa) & !is.na(env_formative) &
               !is.na(ic_env_cont) & !is.na(age) &
               !is.na(female) & !is.na(chronic_n))

      if (nrow(dat) < 50 || sum(dat$event_dfs) < 20) next

      tryCatch({
        fit <- coxph(Surv(time_dfs, event_dfs) ~
                       ic_cfa + env_formative + ic_env_cont + age + female + chronic_n,
                     data = dat)
        coef_df <- as.data.frame(summary(fit)$coefficients)
        coef_df$variable <- rownames(coef_df)
        coef_df$cohort <- coh
        coef_df$subgroup <- age_label
        coef_df$subgroup_type <- "age_group"
        coef_df$n <- nrow(dat)
        coef_df$n_events <- sum(dat$event_dfs)
        names(coef_df)[names(coef_df) == "coef"] <- "log_hr"
        names(coef_df)[names(coef_df) == "se(coef)"] <- "se_log_hr"
        names(coef_df)[names(coef_df) == "exp(coef)"] <- "hr"
        names(coef_df)[names(coef_df) == "Pr(>|z|)"] <- "p_value"
        coef_df$hr_lower <- exp(coef_df$log_hr - 1.96 * coef_df$se_log_hr)
        coef_df$hr_upper <- exp(coef_df$log_hr + 1.96 * coef_df$se_log_hr)
        subgroup_results[[paste(age_label, coh)]] <- coef_df
      }, error = function(e) NULL)
    }
  }

  subgroup_df <- bind_rows(subgroup_results)

  # Pool subgroup interaction terms
  cat("\n  Pooled IC×Env interaction by subgroup:\n")
  sub_pooled <- list()

  for (sg in unique(subgroup_df$subgroup)) {
    sub_data <- subgroup_df %>%
      filter(subgroup == sg & variable == "ic_env_cont")

    if (nrow(sub_data) < 2) next

    tryCatch({
      re <- rma(yi = log_hr, sei = se_log_hr, data = sub_data, method = "REML")
      sub_pooled[[sg]] <- tibble(
        subgroup = sg,
        pooled_hr = exp(as.numeric(re$beta)),
        hr_lower = exp(as.numeric(re$ci.lb)),
        hr_upper = exp(as.numeric(re$ci.ub)),
        p_value = as.numeric(re$pval),
        I2 = re$I2,
        k = re$k
      )
      cat("    ", sprintf("%-12s", sg), ": HR =", sprintf("%.3f", exp(as.numeric(re$beta))),
          "(", sprintf("%.3f", exp(as.numeric(re$ci.lb))), "-",
          sprintf("%.3f", exp(as.numeric(re$ci.ub))), ")",
          "p =", format(as.numeric(re$pval), digits = 3), "\n")
    }, error = function(e) NULL)
  }

  sub_pooled_df <- bind_rows(sub_pooled)

  # Save
  write_csv(subgroup_df, file.path(RESULTS_DIR, "subgroup_stage1.csv"))
  write_csv(sub_pooled_df, file.path(RESULTS_DIR, "subgroup_pooled.csv"))
  cat("\n  Saved: results/subgroup_stage1.csv, results/subgroup_pooled.csv\n")

  sub_pooled_df
}

# ==============================================================================
# STEP 6: PH assumption check + RERI (additive interaction)
# ==============================================================================

step6_diagnostics <- function(surv) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 6: Diagnostics — PH check + RERI\n")
  cat("=" %>% strrep(70), "\n\n")

  # --- 6a: PH assumption test (Schoenfeld) for largest cohort (SHARE) ---
  cat("  6a: PH assumption test (SHARE, M4 model)\n")
  share_dat <- surv %>%
    filter(cohort == "SHARE") %>%
    filter(!is.na(ic_cfa) & !is.na(env_formative) &
           !is.na(ic_env_cont) & !is.na(age) &
           !is.na(female) & !is.na(chronic_n))

  tryCatch({
    fit_share <- coxph(Surv(time_dfs, event_dfs) ~
                         ic_cfa + env_formative + ic_env_cont +
                         age + female + chronic_n,
                       data = share_dat)
    ph_test <- cox.zph(fit_share)
    cat("    Schoenfeld residuals test:\n")
    print(ph_test$table)
    cat("\n")
  }, error = function(e) {
    cat("    PH test failed:", conditionMessage(e), "\n")
  })

  # --- 6b: RERI (Relative Excess Risk due to Interaction) ---
  # Additive interaction measure from the 6-group model
  cat("  6b: RERI (additive interaction) from 6-group model\n")
  cat("  Using M5 pooled estimates: RERI = HR_11 - HR_10 - HR_01 + 1\n")
  cat("  Where: HR_11 = Low IC / Low Env, HR_10 = Low IC / High Env,\n")
  cat("         HR_01 = High IC / Low Env (all vs High IC / High Env)\n\n")

  # Load pooled 6-group estimates
  pooled <- read_csv(file.path(RESULTS_DIR, "stage2_pooled.csv"), show_col_types = FALSE)
  m5 <- pooled %>% filter(model == "M5_6group")

  hr_low_low  <- m5 %>% filter(variable == "Low IC / Low Env") %>% pull(pooled_hr)
  hr_low_high <- m5 %>% filter(variable == "Low IC / High Env") %>% pull(pooled_hr)
  hr_high_low <- m5 %>% filter(variable == "High IC / Low Env") %>% pull(pooled_hr)

  if (length(hr_low_low) > 0 & length(hr_low_high) > 0 & length(hr_high_low) > 0) {
    reri <- hr_low_low - hr_low_high - hr_high_low + 1
    cat("    HR(Low IC/Low Env) =", round(hr_low_low, 3), "\n")
    cat("    HR(Low IC/High Env) =", round(hr_low_high, 3), "\n")
    cat("    HR(High IC/Low Env) =", round(hr_high_low, 3), "\n")
    cat("    RERI =", round(reri, 3), "\n")
    if (reri > 0) {
      cat("    → Positive RERI: SYNERGISTIC interaction on additive scale\n")
      cat("      Combined risk EXCEEDS sum of individual risks\n")
    } else {
      cat("    → Negative RERI: ANTAGONISTIC interaction on additive scale\n")
    }
  } else {
    cat("    Could not compute RERI — missing group estimates\n")
  }

  cat("\n")
}

# ==============================================================================
# STEP 7: Summary
# ==============================================================================

step7_summary <- function(surv, pooled_dfs, pooled_mort, pooled_sub) {
  cat("=" %>% strrep(70), "\n")
  cat("IPD META-ANALYSIS COMPLETE\n")
  cat("=" %>% strrep(70), "\n\n")

  cat("  Analytic cohort:      ", format(nrow(surv), big.mark = ","), "persons\n")
  cat("  DFS events:           ", sum(surv$event_dfs), "\n")
  cat("  Mortality events:     ", sum(surv$event_mort), "\n")
  cat("  Median follow-up:     ", round(median(surv$time_dfs), 1), "years\n\n")

  # Compile summary table
  summary_rows <- list()

  # IC main effect
  ic_dfs <- pooled_dfs %>% filter(model == "M1_IC" & variable == "ic_cfa")
  if (nrow(ic_dfs) > 0) {
    summary_rows[["IC_DFS"]] <- tibble(
      outcome = "DFS", exposure = "IC (per 1-SD)", model = "M1",
      hr = ic_dfs$pooled_hr, hr_lower = ic_dfs$hr_lower, hr_upper = ic_dfs$hr_upper,
      p_value = ic_dfs$p_value, I2 = ic_dfs$I2
    )
  }

  # Env main effect
  env_dfs <- pooled_dfs %>% filter(model == "M2_Env" & variable == "env_formative")
  if (nrow(env_dfs) > 0) {
    summary_rows[["Env_DFS"]] <- tibble(
      outcome = "DFS", exposure = "Env (per 1-SD)", model = "M2",
      hr = env_dfs$pooled_hr, hr_lower = env_dfs$hr_lower, hr_upper = env_dfs$hr_upper,
      p_value = env_dfs$p_value, I2 = env_dfs$I2
    )
  }

  # Interaction
  int_dfs <- pooled_dfs %>% filter(model == "M4_Interaction" & variable == "ic_env_cont")
  if (nrow(int_dfs) > 0) {
    summary_rows[["ICxEnv_DFS"]] <- tibble(
      outcome = "DFS", exposure = "IC×Env interaction", model = "M4",
      hr = int_dfs$pooled_hr, hr_lower = int_dfs$hr_lower, hr_upper = int_dfs$hr_upper,
      p_value = int_dfs$p_value, I2 = int_dfs$I2
    )
  }

  # Mortality interaction
  int_mort <- pooled_mort %>% filter(model == "M4_Int_mort" & variable == "ic_env_cont")
  if (nrow(int_mort) > 0) {
    summary_rows[["ICxEnv_Mort"]] <- tibble(
      outcome = "Mortality", exposure = "IC×Env interaction", model = "M4",
      hr = int_mort$pooled_hr, hr_lower = int_mort$hr_lower, hr_upper = int_mort$hr_upper,
      p_value = int_mort$p_value, I2 = int_mort$I2
    )
  }

  summary_df <- bind_rows(summary_rows)
  write_csv(summary_df, file.path(RESULTS_DIR, "ipd_ma_summary.csv"))

  cat("  Summary table:\n")
  for (r in seq_len(nrow(summary_df))) {
    cat("    ", sprintf("%-20s", summary_df$exposure[r]),
        sprintf("%-10s", summary_df$outcome[r]),
        "HR =", sprintf("%.3f", summary_df$hr[r]),
        "(", sprintf("%.3f", summary_df$hr_lower[r]), "-",
        sprintf("%.3f", summary_df$hr_upper[r]), ")",
        "p =", format(summary_df$p_value[r], digits = 3),
        "I² =", sprintf("%.1f%%", summary_df$I2[r]), "\n")
  }

  cat("\n  Output files:\n")
  cat("    data/survival_cohort.rds\n")
  cat("    results/stage1_estimates.csv\n")
  cat("    results/stage2_pooled.csv\n")
  cat("    results/stage1_mortality.csv\n")
  cat("    results/stage2_mortality.csv\n")
  cat("    results/subgroup_stage1.csv\n")
  cat("    results/subgroup_pooled.csv\n")
  cat("    results/ipd_ma_summary.csv\n")
  cat("=" %>% strrep(70), "\n")
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

surv       <- step1_survival()
stage1     <- step2_stage1(surv)
pooled_dfs <- step3_stage2(stage1)
pooled_mort <- step4_mortality(surv)
pooled_sub <- step5_subgroups(surv)
step6_diagnostics(surv)
step7_summary(surv, pooled_dfs, pooled_mort, pooled_sub)
