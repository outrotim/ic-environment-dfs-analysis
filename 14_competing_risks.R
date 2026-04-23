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

# --- Load data ---
surv_base <- readRDS("data/survival_cohort.rds")
cat(sprintf("Survival cohort: %s persons\n", format(nrow(surv_base), big.mark = ",")))

# --- Construct competing risks variables ---
surv_cr <- surv_base %>%
  filter(!is.na(ic_cfa_c) & !is.na(env_form_c)) %>%
  mutate(
    cs_time = pmin(t_disab, t_death, t_censor, na.rm = TRUE),
    cs_status = case_when(
      cs_time == t_disab & cs_time < Inf ~ 1L,
      cs_time == t_death & cs_time < Inf ~ 2L,
      TRUE ~ 0L
    ),
    cs_time = pmax(cs_time, 0.001)
  )

# Model 1: CSH for Disability (death censored)
stage1_disab <- tibble()
cohorts <- sort(unique(surv_cr$cohort))
for (coh in cohorts) {
  d <- surv_cr %>% filter(cohort == coh) %>%
    mutate(event_disab_cs = as.integer(cs_status == 1))
  n_events <- sum(d$event_disab_cs)
  if (nrow(d) < 100 || n_events < 10) next
  fit <- tryCatch(
    coxph(Surv(cs_time, event_disab_cs) ~ ic_cfa_c * env_form_c + age + female + chronic_n,
          data = d),
    error = function(e) NULL)
  if (is.null(fit)) next
  res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
           n = nrow(d), n_events = n_events, model = "CSH_disability")
  stage1_disab <- bind_rows(stage1_disab, res)
}
int_disab <- stage1_disab %>% filter(grepl(":", term))
ma_disab <- rma(yi = log_hr, sei = se_log_hr, data = int_disab, method = "REML")

# Model 2: CSH for Death (disability censored)
stage1_death <- tibble()
for (coh in cohorts) {
  d <- surv_cr %>% filter(cohort == coh) %>%
    mutate(event_death_cs = as.integer(cs_status == 2))
  n_events <- sum(d$event_death_cs)
  if (nrow(d) < 100 || n_events < 10) next
  fit <- tryCatch(
    coxph(Surv(cs_time, event_death_cs) ~ ic_cfa_c * env_form_c + age + female + chronic_n,
          data = d),
    error = function(e) NULL)
  if (is.null(fit)) next
  res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(cohort = coh, log_hr = log(estimate), se_log_hr = std.error,
           n = nrow(d), n_events = n_events, model = "CSH_death")
  stage1_death <- bind_rows(stage1_death, res)
}
int_death <- stage1_death %>% filter(grepl(":", term))
ma_death <- rma(yi = log_hr, sei = se_log_hr, data = int_death, method = "REML")

# Summary
summary_rows <- tibble(
  outcome = c("CSH-Disability (death censored)", "CSH-Death (disability censored)"),
  hr = c(exp(as.numeric(ma_disab$beta)), exp(as.numeric(ma_death$beta))),
  ci_lo = c(exp(ma_disab$ci.lb), exp(ma_death$ci.lb)),
  ci_hi = c(exp(ma_disab$ci.ub), exp(ma_death$ci.ub)),
  p = c(ma_disab$pval, ma_death$pval),
  i2 = c(ma_disab$I2, ma_death$I2),
  n = c(sum(int_disab$n), sum(int_death$n)),
  events = c(sum(int_disab$n_events), sum(int_death$n_events))
)
print(summary_rows)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write_csv(summary_rows, "results/competing_risks_results.csv")
all_stage1 <- bind_rows(stage1_disab, stage1_death) %>%
  select(model, cohort, term, estimate, conf.low, conf.high, p.value,
         log_hr, se_log_hr, n, n_events)
write_csv(all_stage1, "results/competing_risks_stage1.csv")

cat(sprintf("Completed in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
