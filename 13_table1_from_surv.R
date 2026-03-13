#!/usr/bin/env Rscript
# ==============================================================================
# 13_table1_from_surv.R
# Study 17: Regenerate Table 1 from survival_cohort.rds (N=102,006)
#
# Purpose:
#   The analytic sample used in all IPD-MA and MSM analyses is defined by
#   survival_cohort.rds (N=102,006). This file is already one row per person
#   (no wave column). This script generates the definitive Table 1.
#
# Output:
#   results/table1_surv_by_cohort.csv    (Panel A, by cohort)
#   results/table1_surv_by_icenv.csv     (Panel B, by IC×Env 6-group)
# ==============================================================================

library(tidyverse)

SCRIPT_DIR <- here::here()
DATA_DIR    <- file.path(SCRIPT_DIR, "data")
RESULTS_DIR <- file.path(SCRIPT_DIR, "results")

cat("=" |> strrep(70), "\n")
cat("TABLE 1 REGENERATION FROM survival_cohort.rds\n")
cat("=" |> strrep(70), "\n\n")

# ==============================================================================
# STEP 1: Load survival cohort (already one row per person)
# ==============================================================================

surv <- readRDS(file.path(DATA_DIR, "survival_cohort.rds"))
cat("Loaded survival_cohort.rds:", format(nrow(surv), big.mark = ","), "rows\n")
cat("Unique IDs:", n_distinct(surv$id), "\n\n")

# Per-cohort breakdown
cat("Per-cohort N:\n")
surv %>% count(cohort) %>% print()

# ==============================================================================
# STEP 2: Derive binary variables for Table 1
# ==============================================================================

surv <- surv %>%
  mutate(
    multimorbid = as.integer(chronic_n >= 2),
    iadl_any    = as.integer(iadl_total >= 1)
  )

# ==============================================================================
# STEP 3: Table-building functions
# ==============================================================================

fmt_ms <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

fmt_np <- function(x) {
  x <- x[!is.na(x)]
  n_pos <- sum(x == 1 | x == TRUE)
  sprintf("%s (%.1f)", format(n_pos, big.mark = ","), n_pos / length(x) * 100)
}

build_table <- function(data, group_var) {
  groups <- data %>% pull(!!sym(group_var)) %>% unique() %>% sort()

  vars_cont <- c("age", "chronic_n", "imrc", "dlrc", "depression",
                  "grip_max", "srh_r", "ic_cfa", "env_formative",
                  "social_participation", "n_waves")
  labels_cont <- c("Age, years, mean (SD)", "Chronic diseases, mean (SD)",
                    "Immediate recall (0-10)", "Delayed recall (0-10)",
                    "Depression score", "Grip strength, kg",
                    "Self-rated health (1-5)", "IC CFA score, mean (SD)",
                    "Env formative index, mean (SD)",
                    "Social participation, count", "Follow-up waves, mean (SD)")

  vars_bin <- c("female", "married", "lives_alone", "multimorbid",
                "iadl_any", "event_dfs", "person_died")
  labels_bin <- c("Female, n (%)", "Married/partnered, n (%)",
                  "Living alone, n (%)", "Multimorbidity (>=2), n (%)",
                  "Any IADL disability, n (%)",
                  "Any ADL disability, n (%)",
                  "Died (ever), n (%)")

  # N row
  overall_n <- format(nrow(data), big.mark = ",")
  grp_ns <- data %>% count(!!sym(group_var)) %>%
    mutate(n_fmt = format(n, big.mark = ","))

  result <- tibble(Variable = "N", Overall = overall_n)
  for (g in groups) {
    result[[as.character(g)]] <- grp_ns %>%
      filter(!!sym(group_var) == g) %>% pull(n_fmt)
  }

  # Continuous variables
  for (i in seq_along(vars_cont)) {
    v <- vars_cont[i]
    if (!v %in% names(data)) next
    row <- tibble(Variable = labels_cont[i], Overall = fmt_ms(data[[v]]))
    for (g in groups) {
      sub <- data %>% filter(!!sym(group_var) == g)
      row[[as.character(g)]] <- fmt_ms(sub[[v]])
    }
    result <- bind_rows(result, row)
  }

  # Binary variables
  for (i in seq_along(vars_bin)) {
    v <- vars_bin[i]
    if (!v %in% names(data)) next
    row <- tibble(Variable = labels_bin[i], Overall = fmt_np(data[[v]]))
    for (g in groups) {
      sub <- data %>% filter(!!sym(group_var) == g)
      row[[as.character(g)]] <- fmt_np(sub[[v]])
    }
    result <- bind_rows(result, row)
  }

  result
}

# ==============================================================================
# STEP 4: Generate Panel A (by cohort)
# ==============================================================================

cat("\n--- PANEL A: By Cohort ---\n\n")
table1a <- build_table(surv, "cohort")
print(table1a, n = 30, width = 200)
write_csv(table1a, file.path(RESULTS_DIR, "table1_surv_by_cohort.csv"))
cat("\nSaved: results/table1_surv_by_cohort.csv\n")

# ==============================================================================
# STEP 5: Generate Panel B (by IC×Env 6-group)
# ==============================================================================

cat("\n--- PANEL B: By IC×Env 6-group ---\n\n")
surv_icenv <- surv %>% filter(!is.na(ic_env_6grp))
cat("IC×Env 6-group available:", format(nrow(surv_icenv), big.mark = ","), "\n\n")

table1b <- build_table(surv_icenv, "ic_env_6grp")
print(table1b, n = 30, width = 200)
write_csv(table1b, file.path(RESULTS_DIR, "table1_surv_by_icenv.csv"))
cat("\nSaved: results/table1_surv_by_icenv.csv\n")

# ==============================================================================
# STEP 6: Summary statistics
# ==============================================================================

cat("\n--- KEY SUMMARY STATISTICS ---\n\n")
n_total <- nrow(surv)
cat("Total analytic N:", format(n_total, big.mark = ","), "\n")
cat("Mean age:", round(mean(surv$age, na.rm = TRUE), 1),
    "(SD", round(sd(surv$age, na.rm = TRUE), 1), ")\n")
cat("Female %:", round(mean(surv$female, na.rm = TRUE) * 100, 1), "\n")

# Events
dfs_events <- sum(surv$event_dfs == 1, na.rm = TRUE)
deaths <- sum(surv$event_mort == 1, na.rm = TRUE)
cat("DFS events:", format(dfs_events, big.mark = ","),
    sprintf("(%.1f%%)\n", dfs_events / n_total * 100))
cat("Deaths:", format(deaths, big.mark = ","),
    sprintf("(%.1f%%)\n", deaths / n_total * 100))

# IC + Env availability
n_both <- sum(!is.na(surv$ic_cfa) & !is.na(surv$env_formative))
cat("IC+Env both available:", format(n_both, big.mark = ","), "\n")

cat("\n✓ Table 1 regeneration complete.\n")
