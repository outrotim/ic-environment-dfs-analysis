#!/usr/bin/env Rscript
# 07_trajectory.R — LCGA Trajectory Analysis for IC and Env
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# Purpose: Identify latent trajectory classes of IC & Env over time,
#          cross-classify into 4 groups, and link to DFS outcomes.
#
# Pipeline:
#   Step 1: Prepare trajectory dataset (≥3 obs per person)
#   Step 2: IC trajectory LCGA per cohort (model selection + full fit)
#   Step 3: Env trajectory LCGA per cohort
#   Step 4: Dichotomize + cross-classify → 4-group (IC↓Env↓, IC↓Env→, IC→Env↓, IC→Env→)
#   Step 5: Link 4-group trajectories to DFS (Cox per-cohort → RE-MA)
#   Step 6: Sensitivity — empirical OLS slope approach
#   Step 7: Output summary
#
# Inputs:  data/analytic_final.rds, data/survival_cohort.rds
# Outputs: data/trajectory_classes.rds, results/lcga_*_fit.csv,
#          results/trajectory_dfs_*.csv, results/trajectory_summary.csv

library(tidyverse)
library(lcmm)
library(survival)
library(metafor)

set.seed(2024)
t0 <- Sys.time()

# ============================================================
cat("=== Step 1: Prepare trajectory dataset ===\n")
# ============================================================

pooled <- readRDS("data/analytic_final.rds")
surv   <- readRDS("data/survival_cohort.rds")

# Filter: ≥3 valid IC+Env observations per person
traj_data <- pooled %>%
  filter(!is.na(ic_cfa) & !is.na(env_formative)) %>%
  group_by(cohort, id) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  # Time variable: years from first observation
  group_by(cohort, id) %>%
  mutate(time = interview_year - min(interview_year, na.rm = TRUE)) %>%
  ungroup() %>%
  # Numeric UID for lcmm (must be integer)
  mutate(uid = as.integer(factor(paste(cohort, id, sep = "___"))))

uid_map <- traj_data %>% distinct(uid, cohort, id)

cat(sprintf("Trajectory dataset: %d persons, %d observations\n",
            n_distinct(traj_data$uid), nrow(traj_data)))
cat("Per cohort:\n")
traj_data %>%
  group_by(cohort) %>%
  summarise(n_persons = n_distinct(uid), n_obs = n(),
            med_waves = median(table(uid)), .groups = "drop") %>%
  print()

# ============================================================
cat("\n=== Step 2: IC Trajectory LCGA ===\n")
# ============================================================

#' Fit LCGA models for K = 1..max_k, select best K by BIC,
#' refit on full data, and return class assignments.
#'
#' Classes are labelled C1 (lowest mean) to Ck (highest mean).
fit_lcga_cohort <- function(dat, outcome_var, cohort_name,
                            max_k = 4, subsample_n = 3000) {

  cat(sprintf("  [%s] %s — n=%d persons\n", cohort_name, outcome_var,
              n_distinct(dat$uid)))

  # --- subsample for model selection ---
  uids_all <- unique(dat$uid)
  if (length(uids_all) > subsample_n) {
    sub_uids <- sample(uids_all, subsample_n)
    sub_dat  <- dat %>% filter(uid %in% sub_uids)
    cat(sprintf("    Subsampled to %d for model selection\n", subsample_n))
  } else {
    sub_dat <- dat
  }

  sub_dat$y <- sub_dat[[outcome_var]]

  # K = 1 (baseline)
  m1 <- tryCatch(
    hlme(y ~ time, subject = "uid", ng = 1, data = sub_dat, verbose = FALSE),
    error = function(e) { cat("    K=1 FAILED\n"); NULL }
  )
  if (is.null(m1)) return(NULL)

  models <- list(`1` = m1)
  cat(sprintf("    K=1  BIC=%.1f\n", m1$BIC))

  # K = 2 .. max_k

  for (k in 2:max_k) {
    mk <- tryCatch(
      hlme(y ~ time, mixture = ~ time, subject = "uid", ng = k,
           data = sub_dat, B = m1, verbose = FALSE, maxiter = 500),
      error = function(e) NULL
    )
    if (!is.null(mk) && mk$conv == 1) {
      models[[as.character(k)]] <- mk
      cat(sprintf("    K=%d  BIC=%.1f  (converged)\n", k, mk$BIC))
    } else {
      cat(sprintf("    K=%d  FAILED or not converged\n", k))
    }
  }

  # --- Model comparison ---
  fit_stats <- tibble(
    cohort   = cohort_name,
    variable = outcome_var,
    K        = as.integer(names(models)),
    BIC      = sapply(models, \(m) m$BIC),
    loglik   = sapply(models, \(m) m$loglik),
    npm      = sapply(models, \(m) length(m$best))
  ) %>% mutate(delta_BIC = BIC - min(BIC))

  best_k <- fit_stats %>% slice_min(BIC, n = 1) %>% pull(K)
  cat(sprintf("    ★ Best K = %d\n", best_k))

  # --- Refit best K on FULL data ---
  full_dat   <- dat
  full_dat$y <- full_dat[[outcome_var]]

  m1_full <- hlme(y ~ time, subject = "uid", ng = 1,
                  data = full_dat, verbose = FALSE)

  if (best_k == 1) {
    best_model <- m1_full
  } else {
    best_model <- tryCatch(
      hlme(y ~ time, mixture = ~ time, subject = "uid", ng = best_k,
           data = full_dat, B = m1_full, verbose = FALSE, maxiter = 500),
      error = function(e) {
        cat("    Full-data refit FAILED; falling back to K=1\n")
        m1_full
      }
    )
    # If didn't converge, fall back
    if (best_model$conv != 1) {
      cat("    Full-data refit did NOT converge; falling back to K=1\n")
      best_model <- m1_full
      best_k <- 1L
    }
  }

  # --- Extract class assignments ---
  classes <- best_model$pprob %>%
    select(uid, class) %>%
    as_tibble()

  # Relabel C1 = lowest mean, Ck = highest mean
  if (best_k > 1) {
    class_means <- dat %>%
      inner_join(classes, by = "uid") %>%
      group_by(class) %>%
      summarise(mean_val = mean(.data[[outcome_var]], na.rm = TRUE),
                .groups = "drop") %>%
      arrange(mean_val)

    relabel <- setNames(paste0("C", seq_len(best_k)), class_means$class)
    classes$class_label <- relabel[as.character(classes$class)]
  } else {
    classes$class_label <- "C1"
  }

  # --- Class trajectory summary (mean ± SD per time point) ---
  class_traj <- dat %>%
    inner_join(classes, by = "uid") %>%
    group_by(class_label, time) %>%
    summarise(mean_y = mean(.data[[outcome_var]], na.rm = TRUE),
              sd_y   = sd(.data[[outcome_var]], na.rm = TRUE),
              n      = n(), .groups = "drop") %>%
    mutate(cohort = cohort_name, variable = outcome_var)

  list(
    fit_stats  = fit_stats,
    best_k     = best_k,
    classes    = classes,
    class_traj = class_traj
  )
}

# --- Run IC LCGA per cohort ---
cohorts <- sort(unique(traj_data$cohort))
ic_results  <- list()
ic_fit_all  <- tibble()
ic_traj_all <- tibble()

for (coh in cohorts) {
  dat_coh <- traj_data %>% filter(cohort == coh)
  res <- fit_lcga_cohort(dat_coh, "ic_cfa", coh, max_k = 4, subsample_n = 3000)
  if (!is.null(res)) {
    ic_results[[coh]] <- res
    ic_fit_all  <- bind_rows(ic_fit_all, res$fit_stats)
    ic_traj_all <- bind_rows(ic_traj_all, res$class_traj)
  }
}

cat("\nIC LCGA model selection summary:\n")
print(ic_fit_all)

# ============================================================
cat("\n=== Step 3: Env Trajectory LCGA ===\n")
# ============================================================

env_results  <- list()
env_fit_all  <- tibble()
env_traj_all <- tibble()

for (coh in cohorts) {
  dat_coh <- traj_data %>% filter(cohort == coh)
  res <- fit_lcga_cohort(dat_coh, "env_formative", coh, max_k = 4, subsample_n = 3000)
  if (!is.null(res)) {
    env_results[[coh]] <- res
    env_fit_all  <- bind_rows(env_fit_all, res$fit_stats)
    env_traj_all <- bind_rows(env_traj_all, res$class_traj)
  }
}

cat("\nEnv LCGA model selection summary:\n")
print(env_fit_all)

# ============================================================
cat("\n=== Step 4: Dichotomize + Cross-classify ===\n")
# ============================================================

# Strategy: dichotomize IC and Env trajectory classes into
# "Declining" (C1 = lowest mean class) vs "Stable/High" (all other classes).
# This yields 4 interpretable groups uniform across cohorts:
#   IC↓Env↓  IC↓Env→  IC→Env↓  IC→Env→

trajectory_classes <- tibble()

for (coh in cohorts) {
  ic_cls <- ic_results[[coh]]$classes %>%
    mutate(
      ic_declining = ifelse(class_label == "C1", "IC_decline", "IC_stable")
    ) %>%
    select(uid, ic_class = class_label, ic_declining)

  env_cls <- env_results[[coh]]$classes %>%
    mutate(
      env_declining = ifelse(class_label == "C1", "Env_decline", "Env_stable")
    ) %>%
    select(uid, env_class = class_label, env_declining)

  merged <- inner_join(ic_cls, env_cls, by = "uid") %>%
    mutate(
      cohort = coh,
      traj_4grp = paste0(ic_declining, "/", env_declining)
    )

  trajectory_classes <- bind_rows(trajectory_classes, merged)
}

# Merge back original IDs
trajectory_classes <- trajectory_classes %>%
  left_join(uid_map, by = c("uid", "cohort"))

cat("4-group trajectory distribution:\n")
print(trajectory_classes %>%
        count(traj_4grp) %>%
        mutate(pct = round(n / sum(n) * 100, 1)) %>%
        arrange(desc(n)))

cat("\nPer cohort:\n")
print(trajectory_classes %>%
        count(cohort, traj_4grp) %>%
        group_by(cohort) %>%
        mutate(pct = round(n / sum(n) * 100, 1)) %>%
        ungroup())

saveRDS(trajectory_classes, "data/trajectory_classes.rds")
cat("Saved: data/trajectory_classes.rds\n")

# ============================================================
cat("\n=== Step 5: Trajectory 4-group → DFS ===\n")
# ============================================================

# Merge with survival data
traj_surv <- surv %>%
  inner_join(
    trajectory_classes %>% select(cohort, id, ic_declining, env_declining, traj_4grp),
    by = c("cohort", "id")
  )

cat(sprintf("Trajectory-survival overlap: %d persons (%.1f%% of survival cohort)\n",
            nrow(traj_surv), nrow(traj_surv) / nrow(surv) * 100))

# Reference: IC_stable/Env_stable (best trajectory)
traj_surv$traj_4grp <- relevel(factor(traj_surv$traj_4grp),
                                ref = "IC_stable/Env_stable")

# --- Stage 1: Per-cohort Cox ---
stage1_traj <- tibble()

for (coh in cohorts) {
  dat_coh <- traj_surv %>% filter(cohort == coh)
  if (nrow(dat_coh) < 50 || length(unique(dat_coh$traj_4grp)) < 2) next

  fit <- tryCatch(
    coxph(Surv(time_dfs, event_dfs) ~ traj_4grp + age + female + chronic_n,
          data = dat_coh),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    coefs <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(grepl("traj_4grp", term)) %>%
      mutate(cohort = coh, n = nrow(dat_coh),
             log_hr = log(estimate), se_log_hr = std.error)
    stage1_traj <- bind_rows(stage1_traj, coefs)
  }
}

cat("\nStage 1 results (trajectory 4-group → DFS):\n")
stage1_traj %>%
  select(cohort, term, estimate, conf.low, conf.high, p.value) %>%
  print(n = 30)

# --- Stage 2: RE-MA for each trajectory contrast ---
stage2_traj <- tibble()
contrasts <- unique(stage1_traj$term)

for (ct in contrasts) {
  sub <- stage1_traj %>% filter(term == ct)
  if (nrow(sub) < 2) next

  ma <- tryCatch(
    rma(yi = log_hr, sei = se_log_hr, data = sub, method = "REML"),
    error = function(e) NULL
  )
  if (!is.null(ma)) {
    stage2_traj <- bind_rows(stage2_traj, tibble(
      term     = ct,
      pooled_hr = exp(ma$beta[1,1]),
      ci_lower  = exp(ma$ci.lb),
      ci_upper  = exp(ma$ci.ub),
      p_value   = ma$pval[1],
      I2        = ma$I2,
      k         = ma$k
    ))
  }
}

cat("\nStage 2 pooled trajectory → DFS:\n")
print(stage2_traj)

# --- Key contrast: environment buffering of IC decline ---
# IC_decline/Env_stable vs IC_decline/Env_decline
cat("\n--- Key buffering contrast ---\n")
buff_data <- traj_surv %>%
  filter(ic_declining == "IC_decline") %>%
  mutate(env_buffer = ifelse(env_declining == "Env_stable",
                             "Buffered", "Unbuffered"))

buff_stage1 <- tibble()
for (coh in cohorts) {
  dat_coh <- buff_data %>% filter(cohort == coh)
  if (nrow(dat_coh) < 30 || length(unique(dat_coh$env_buffer)) < 2) next

  dat_coh$env_buffer <- relevel(factor(dat_coh$env_buffer), ref = "Unbuffered")
  fit <- tryCatch(
    coxph(Surv(time_dfs, event_dfs) ~ env_buffer + age + female + chronic_n,
          data = dat_coh),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    coef <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(grepl("env_buffer", term)) %>%
      mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error)
    buff_stage1 <- bind_rows(buff_stage1, coef)
  }
}

if (nrow(buff_stage1) >= 2) {
  ma_buff <- rma(yi = log_hr, sei = se_log_hr, data = buff_stage1, method = "REML")
  cat(sprintf("Buffered (IC↓ + Env stable) vs Unbuffered (IC↓ + Env↓):\n"))
  cat(sprintf("  Pooled HR = %.3f (%.3f–%.3f), p = %.4f, I² = %.1f%%\n",
              exp(ma_buff$beta), exp(ma_buff$ci.lb), exp(ma_buff$ci.ub),
              ma_buff$pval, ma_buff$I2))
}

write_csv(stage1_traj, "results/trajectory_dfs_stage1.csv")
write_csv(stage2_traj, "results/trajectory_dfs_stage2.csv")
cat("Saved: results/trajectory_dfs_stage1.csv, results/trajectory_dfs_stage2.csv\n")

# ============================================================
cat("\n=== Step 6: Sensitivity — Empirical OLS slope ===\n")
# ============================================================

# Per-person OLS slope for IC and Env (fast, no latent class)
slopes <- traj_data %>%
  group_by(cohort, id, uid) %>%
  summarise(
    n_obs        = n(),
    ic_slope     = coef(lm(ic_cfa ~ time))[2],
    env_slope    = coef(lm(env_formative ~ time))[2],
    .groups      = "drop"
  )

# Dichotomize by median slope (per cohort)
slopes <- slopes %>%
  group_by(cohort) %>%
  mutate(
    ic_trend  = ifelse(ic_slope  < median(ic_slope,  na.rm = TRUE),
                       "IC_decline", "IC_stable"),
    env_trend = ifelse(env_slope < median(env_slope, na.rm = TRUE),
                       "Env_decline", "Env_stable"),
    slope_4grp = paste0(ic_trend, "/", env_trend)
  ) %>%
  ungroup()

cat("Empirical slope 4-group distribution:\n")
slopes %>% count(slope_4grp) %>% mutate(pct = round(n / sum(n) * 100, 1)) %>% print()

# Merge with survival
slope_surv <- surv %>%
  inner_join(slopes %>% select(cohort, id, ic_trend, env_trend, slope_4grp),
             by = c("cohort", "id"))

cat(sprintf("Slope-survival overlap: %d persons\n", nrow(slope_surv)))

# Pooled Cox (stratified by cohort)
slope_surv$slope_4grp <- relevel(factor(slope_surv$slope_4grp),
                                  ref = "IC_stable/Env_stable")

slope_cox <- coxph(Surv(time_dfs, event_dfs) ~ slope_4grp + age + female + chronic_n +
                     strata(cohort),
                   data = slope_surv)

cat("\nEmpirical slope 4-group → DFS (stratified Cox):\n")
slope_res <- broom::tidy(slope_cox, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(grepl("slope_4grp", term))
print(slope_res %>% select(term, estimate, conf.low, conf.high, p.value))

# Buffering contrast in slope approach
slope_buff <- slope_surv %>%
  filter(ic_trend == "IC_decline") %>%
  mutate(env_buffer = ifelse(env_trend == "Env_stable", "Buffered", "Unbuffered"))

slope_buff$env_buffer <- relevel(factor(slope_buff$env_buffer), ref = "Unbuffered")
buff_cox <- coxph(Surv(time_dfs, event_dfs) ~ env_buffer + age + female + chronic_n +
                    strata(cohort),
                  data = slope_buff)

buff_slope_res <- broom::tidy(buff_cox, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(grepl("env_buffer", term))

cat(sprintf("\nSlope sensitivity — Buffered vs Unbuffered (IC decline subgroup):\n"))
cat(sprintf("  HR = %.3f (%.3f–%.3f), p = %.4f\n",
            buff_slope_res$estimate, buff_slope_res$conf.low,
            buff_slope_res$conf.high, buff_slope_res$p.value))

write_csv(slope_res, "results/slope_sensitivity_dfs.csv")
cat("Saved: results/slope_sensitivity_dfs.csv\n")

# ============================================================
cat("\n=== Step 7: Summary ===\n")
# ============================================================

# Compile overall summary
summary_lines <- c(
  sprintf("IC LCGA best K: %s",
          paste(sapply(cohorts, \(c) sprintf("%s=K%d", c, ic_results[[c]]$best_k)),
                collapse = ", ")),
  sprintf("Env LCGA best K: %s",
          paste(sapply(cohorts, \(c) sprintf("%s=K%d", c, env_results[[c]]$best_k)),
                collapse = ", ")),
  sprintf("Trajectory-survival overlap: %d / %d (%.1f%%)",
          nrow(traj_surv), nrow(surv), nrow(traj_surv) / nrow(surv) * 100),
  sprintf("Slope-survival overlap: %d / %d (%.1f%%)",
          nrow(slope_surv), nrow(surv), nrow(slope_surv) / nrow(surv) * 100)
)

cat("\n--- SUMMARY ---\n")
for (s in summary_lines) cat(s, "\n")

if (nrow(stage2_traj) > 0) {
  cat("\nPooled trajectory 4-group → DFS (Stage 2):\n")
  for (i in seq_len(nrow(stage2_traj))) {
    r <- stage2_traj[i, ]
    cat(sprintf("  %s: HR=%.3f (%.3f–%.3f), p=%.4f, I²=%.1f%%\n",
                r$term, r$pooled_hr, r$ci_lower, r$ci_upper, r$p_value, r$I2))
  }
}

# Save fit stats and trajectory data for figures
write_csv(ic_fit_all, "results/lcga_ic_fit.csv")
write_csv(env_fit_all, "results/lcga_env_fit.csv")
write_csv(bind_rows(ic_traj_all, env_traj_all), "results/lcga_class_trajectories.csv")

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat(sprintf("\n=== TRAJECTORY ANALYSIS COMPLETE (%.1f min) ===\n", elapsed))
cat("Output files:\n")
cat("  data/trajectory_classes.rds\n")
cat("  results/lcga_ic_fit.csv\n")
cat("  results/lcga_env_fit.csv\n")
cat("  results/lcga_class_trajectories.csv\n")
cat("  results/trajectory_dfs_stage1.csv\n")
cat("  results/trajectory_dfs_stage2.csv\n")
cat("  results/slope_sensitivity_dfs.csv\n")
