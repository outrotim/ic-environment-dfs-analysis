#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
test_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else normalizePath("tests/test_consistency_followup.R")
study_dir <- normalizePath(file.path(dirname(test_file), ".."))
audit_script <- paste(
  readLines(file.path(study_dir, "12_consistency_check.R"), warn = FALSE),
  collapse = "\n"
)

pooled <- readRDS(file.path(study_dir, "data", "analytic_final.rds"))
authoritative <- readRDS(file.path(study_dir, "data", "survival_cohort.rds"))

baseline <- pooled %>%
  group_by(cohort, id) %>%
  slice_min(wave, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(baseline_wave = wave, baseline_year = interview_year)

eligible_baseline <- baseline %>%
  filter(!is.na(adl_total) & adl_total == 0)

disability_onset <- pooled %>%
  inner_join(
    eligible_baseline %>% select(cohort, id, baseline_wave),
    by = c("cohort", "id")
  ) %>%
  filter(wave > baseline_wave & !is.na(adl_total) & adl_total >= 1) %>%
  group_by(cohort, id) %>%
  slice_min(wave, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cohort, id, disability_year = interview_year)

last_observation <- pooled %>%
  group_by(cohort, id) %>%
  summarise(last_interview_year = max(interview_year, na.rm = TRUE), .groups = "drop")

eligible <- eligible_baseline %>%
  left_join(disability_onset, by = c("cohort", "id")) %>%
  left_join(last_observation, by = c("cohort", "id")) %>%
  mutate(
    t_disab = ifelse(!is.na(disability_year), disability_year - baseline_year, Inf),
    t_death = ifelse(
      person_died == 1 & !is.na(death_year) & death_year > baseline_year,
      death_year - baseline_year,
      Inf
    ),
    t_censor = last_interview_year - baseline_year,
    time_dfs_raw = pmin(t_disab, t_death, t_censor, na.rm = TRUE),
    event_dfs = as.integer(
      time_dfs_raw < Inf & (time_dfs_raw == t_disab | time_dfs_raw == t_death)
    ),
    time_dfs = ifelse(event_dfs == 1, time_dfs_raw, t_censor)
  )

analytic <- eligible %>% filter(!is.na(time_dfs) & time_dfs > 0)
expected_by_cohort <- authoritative %>% count(cohort, name = "n_authoritative")
actual_by_cohort <- analytic %>% count(cohort, name = "n_rederived")
comparison <- full_join(actual_by_cohort, expected_by_cohort, by = "cohort")

stopifnot(
  "Audit script must reproduce DFS analysis time rather than use row count" =
    grepl("time_dfs_raw", audit_script, fixed = TRUE) &&
      !grepl("filter(n_waves > 1)", audit_script, fixed = TRUE),
  "Zero/NA DFS-time exclusions must equal the flow diagram" =
    sum(is.na(eligible$time_dfs) | eligible$time_dfs <= 0) == 44555L,
  "DFS-time cascade must reproduce N=102,006" = nrow(analytic) == 102006L,
  "Re-derived cohort counts must match survival_cohort.rds" =
    all(comparison$n_rederived == comparison$n_authoritative)
)

cat("PASS: consistency audit reproduces DFS analysis time and N=102,006.\n")
