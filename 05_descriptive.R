# ==============================================================================
# 05_descriptive.R
# Study 17: IC-Environment Fit and Disability-Free Survival in Older Adults
#
# Purpose: Generate descriptive statistics and Table 1
#   (1) Select baseline sample (first observation per person)
#   (2) Table 1a: Baseline characteristics by cohort (5 cohorts + overall)
#   (3) Table 1b: Baseline characteristics by IC×Env 6-group
#   (4) Supplementary: variable availability, follow-up summary
#   (5) Export publication-ready CSV tables
#
# Input:   data/analytic_final.rds (from 04_env_index.R)
# Output:  results/table1_by_cohort.csv
#          results/table1_by_icenv_group.csv
#          results/baseline_summary.csv
#          results/variable_availability.csv
#
# Dependencies: tidyverse
# Author:  Study 17 Team
# Date:    2026-03-14
# Version: 1.0
# ==============================================================================

library(tidyverse)

# --- Configuration ------------------------------------------------------------

SCRIPT_DIR <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) {
    here::here()
  }
)
DATA_DIR    <- file.path(SCRIPT_DIR, "data")
RESULTS_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Helper Functions
# ==============================================================================

# Format: mean (SD) or median [IQR]
fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

fmt_median_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  q <- quantile(x, c(0.5, 0.25, 0.75))
  sprintf("%.1f [%.1f, %.1f]", q[1], q[2], q[3])
}

# Format: n (%)
fmt_n_pct <- function(x, condition = NULL) {
  if (!is.null(condition)) x <- x[!is.na(x)]
  n_total <- length(x)
  if (n_total == 0) return("—")
  if (is.null(condition)) {
    # Binary 0/1 variable
    n_pos <- sum(x == 1, na.rm = TRUE)
  } else {
    n_pos <- sum(condition, na.rm = TRUE)
    n_total <- sum(!is.na(condition))
  }
  sprintf("%s (%.1f)", format(n_pos, big.mark = ","), n_pos / n_total * 100)
}

# Compute descriptive row for a continuous variable
desc_continuous <- function(data, var, label, groups, fmt = "mean_sd") {
  func <- if (fmt == "mean_sd") fmt_mean_sd else fmt_median_iqr

  overall <- func(data[[var]])

  by_group <- data %>%
    group_by(across(all_of(groups))) %>%
    summarise(val = func(.data[[var]]), .groups = "drop") %>%
    pivot_wider(names_from = all_of(groups), values_from = val)

  bind_cols(
    tibble(Variable = label, Overall = overall),
    by_group
  )
}

# Compute descriptive row for a binary variable
desc_binary <- function(data, var, label, groups) {
  overall_n <- sum(data[[var]] == 1, na.rm = TRUE)
  overall_total <- sum(!is.na(data[[var]]))
  overall <- sprintf("%s (%.1f)", format(overall_n, big.mark = ","),
                     overall_n / overall_total * 100)

  by_group <- data %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(across(all_of(groups))) %>%
    summarise(
      val = {
        n_pos <- sum(.data[[var]] == 1, na.rm = TRUE)
        n_tot <- n()
        sprintf("%s (%.1f)", format(n_pos, big.mark = ","), n_pos / n_tot * 100)
      },
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = all_of(groups), values_from = val)

  bind_cols(
    tibble(Variable = label, Overall = overall),
    by_group
  )
}

# Compute N row
desc_n <- function(data, groups) {
  overall <- format(nrow(data), big.mark = ",")

  by_group <- data %>%
    group_by(across(all_of(groups))) %>%
    summarise(val = format(n(), big.mark = ","), .groups = "drop") %>%
    pivot_wider(names_from = all_of(groups), values_from = val)

  bind_cols(
    tibble(Variable = "N", Overall = overall),
    by_group
  )
}

# ==============================================================================
# STEP 1: Load Data and Select Baseline Sample
# ==============================================================================

step1_baseline <- function() {
  cat("=" %>% strrep(70), "\n")
  cat("Step 1: Loading data and selecting baseline sample\n")
  cat("=" %>% strrep(70), "\n\n")

  pooled <- readRDS(file.path(DATA_DIR, "analytic_final.rds"))
  cat("  Full dataset:", format(nrow(pooled), big.mark = ","), "rows,",
      n_distinct(pooled$id), "persons,", ncol(pooled), "columns\n")

  # Baseline: first observation per person
  baseline <- pooled %>%
    group_by(id) %>%
    slice_min(wave, n = 1, with_ties = FALSE) %>%
    ungroup()

  cat("  Baseline sample:", format(nrow(baseline), big.mark = ","), "persons\n")
  cat("  Cohort distribution:\n")
  baseline %>%
    count(cohort) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    print()

  # Create derived variables for Table 1
  baseline <- baseline %>%
    mutate(
      # ADL/IADL disability (binary)
      adl_any  = as.integer(adl_total >= 1),
      iadl_any = as.integer(iadl_total >= 1),

      # Age groups
      age_grp = cut(age, breaks = c(59, 64, 69, 74, 79, 84, 120),
                    labels = c("60-64", "65-69", "70-74",
                               "75-79", "80-84", "85+"),
                    right = TRUE),

      # Married/partnered (binary)
      married = {
        m <- tolower(as.character(marital))
        as.integer(m %in% c("1", "5", "7",
                            "1.married", "1.married, spouse present",
                            "married", "partnered", "1.married/partnered",
                            "5.registered partnership"))
      },

      # Living alone (binary)
      lives_alone = as.integer(living_arrangement == 0),

      # Chronic disease burden
      multimorbid = as.integer(chronic_n >= 2),

      # Number of follow-up waves
      n_waves_person = NA_integer_  # will compute from full data
    )

  # Compute waves per person from full dataset
  waves_per_person <- pooled %>%
    group_by(id) %>%
    summarise(n_waves = n(), .groups = "drop")

  baseline <- baseline %>%
    left_join(waves_per_person, by = "id")

  cat("\n  Waves per person: mean =", round(mean(baseline$n_waves), 1),
      ", median =", median(baseline$n_waves), "\n")

  # Follow-up summary
  cat("\n  Mortality:\n")
  cat("    Deaths:", sum(baseline$person_died, na.rm = TRUE),
      sprintf("(%.1f%%)\n", mean(baseline$person_died, na.rm = TRUE) * 100))
  cat("    ADL any:", sum(baseline$adl_any, na.rm = TRUE),
      sprintf("(%.1f%%)\n", mean(baseline$adl_any, na.rm = TRUE) * 100))
  cat("    IADL any:", sum(baseline$iadl_any, na.rm = TRUE),
      sprintf("(%.1f%%)\n", mean(baseline$iadl_any, na.rm = TRUE) * 100))

  cat("\n")
  list(pooled = pooled, baseline = baseline)
}

# ==============================================================================
# STEP 2: Table 1 — By Cohort
# ==============================================================================

step2_table1_cohort <- function(baseline) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 2: Table 1 — Baseline characteristics by cohort\n")
  cat("=" %>% strrep(70), "\n\n")

  grp <- "cohort"

  rows <- bind_rows(
    # --- N ---
    desc_n(baseline, grp),

    # --- Demographics ---
    tibble(Variable = "Demographics", Overall = "", .rows = 1) %>%
      bind_cols(tibble()),
    desc_continuous(baseline, "age", "  Age, years, mean (SD)", grp),
    desc_binary(baseline, "female", "  Female, n (%)", grp),
    desc_continuous(baseline, "edu_years", "  Education, years, mean (SD)", grp),
    desc_binary(baseline, "married", "  Married/partnered, n (%)", grp),
    desc_binary(baseline, "lives_alone", "  Living alone, n (%)", grp),
    desc_continuous(baseline, "chronic_n", "  Chronic diseases, mean (SD)", grp),
    desc_binary(baseline, "multimorbid", "  Multimorbidity (>=2), n (%)", grp),

    # --- IC indicators (raw) ---
    tibble(Variable = "IC Indicators", Overall = "", .rows = 1) %>%
      bind_cols(tibble()),
    desc_continuous(baseline, "imrc", "  Immediate recall (0-10)", grp),
    desc_continuous(baseline, "dlrc", "  Delayed recall (0-10)", grp),
    desc_continuous(baseline, "depression", "  Depression score", grp),
    desc_continuous(baseline, "grip_max", "  Grip strength, kg", grp),
    desc_continuous(baseline, "srh_r", "  Self-rated health (1-5, higher=better)", grp),

    # --- IC composites ---
    tibble(Variable = "IC Composite Scores", Overall = "", .rows = 1) %>%
      bind_cols(tibble()),
    desc_continuous(baseline, "ic_cfa", "  IC CFA score, mean (SD)", grp),
    desc_continuous(baseline, "ic_sum", "  IC sum score, mean (SD)", grp),

    # --- Env indicators ---
    tibble(Variable = "Environment Indicators", Overall = "", .rows = 1) %>%
      bind_cols(tibble()),
    desc_continuous(baseline, "social_participation",
                    "  Social participation, count", grp),
    desc_continuous(baseline, "env_formative",
                    "  Env formative index, mean (SD)", grp),

    # --- IC × Env ---
    tibble(Variable = "IC × Env Classification", Overall = "", .rows = 1) %>%
      bind_cols(tibble()),

    # --- Outcomes ---
    tibble(Variable = "Outcomes", Overall = "", .rows = 1) %>%
      bind_cols(tibble()),
    desc_binary(baseline, "adl_any", "  Any ADL disability, n (%)", grp),
    desc_binary(baseline, "iadl_any", "  Any IADL disability, n (%)", grp),
    desc_binary(baseline, "person_died", "  Died (ever), n (%)", grp),
    desc_continuous(baseline, "n_waves", "  Follow-up waves, mean (SD)", grp)
  )

  # Fill NA columns for section headers
  for (col in names(rows)) {
    rows[[col]][is.na(rows[[col]])] <- ""
  }

  cat("  Table 1 by Cohort:\n")
  print(rows, n = 50, width = 200)

  # Save
  out_csv <- file.path(RESULTS_DIR, "table1_by_cohort.csv")
  write_csv(rows, out_csv)
  cat("\n  Saved:", out_csv, "\n\n")

  rows
}

# ==============================================================================
# STEP 3: Table 1 — By IC×Env 6-Group
# ==============================================================================

step3_table1_icenv <- function(baseline) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 3: Table 1 — Baseline characteristics by IC×Env group\n")
  cat("=" %>% strrep(70), "\n\n")

  # Filter to those with 6-group assignment
  data <- baseline %>% filter(!is.na(ic_env_6grp))
  cat("  Persons with IC×Env group:", format(nrow(data), big.mark = ","), "\n\n")

  grp <- "ic_env_6grp"

  rows <- bind_rows(
    desc_n(data, grp),

    # Demographics
    desc_continuous(data, "age", "Age, years", grp),
    desc_binary(data, "female", "Female, n (%)", grp),
    desc_continuous(data, "edu_years", "Education, years", grp),
    desc_binary(data, "married", "Married/partnered, n (%)", grp),
    desc_continuous(data, "chronic_n", "Chronic diseases", grp),

    # IC & Env scores
    desc_continuous(data, "ic_cfa", "IC CFA score", grp),
    desc_continuous(data, "env_formative", "Env formative index", grp),

    # Outcomes
    desc_binary(data, "adl_any", "Any ADL disability, n (%)", grp),
    desc_binary(data, "iadl_any", "Any IADL disability, n (%)", grp),
    desc_binary(data, "person_died", "Died (ever), n (%)", grp),
    desc_continuous(data, "n_waves", "Follow-up waves", grp)
  )

  for (col in names(rows)) {
    rows[[col]][is.na(rows[[col]])] <- ""
  }

  cat("  Table 1 by IC×Env Group:\n")
  print(rows, n = 30, width = 250)

  out_csv <- file.path(RESULTS_DIR, "table1_by_icenv_group.csv")
  write_csv(rows, out_csv)
  cat("\n  Saved:", out_csv, "\n\n")

  rows
}

# ==============================================================================
# STEP 4: Variable Availability Matrix
# ==============================================================================

step4_availability <- function(baseline) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 4: Variable availability matrix\n")
  cat("=" %>% strrep(70), "\n\n")

  key_vars <- c(
    # Demographics
    "age", "female", "edu_years", "marital",
    # IC raw indicators
    "imrc", "dlrc", "ser7", "depression", "vision_r", "hearing_r",
    "grip_max", "gait_speed", "bmi", "srh",
    # IC composites
    "ic_cfa", "ic_sum",
    # Env raw
    "social_participation", "financial_strain", "living_arrangement",
    "social_contact",
    # Env composites
    "env_formative",
    # Interaction
    "ic_env_6grp",
    # Outcomes
    "adl_total", "iadl_total", "person_died", "death_year",
    # Other
    "chronic_n", "loneliness", "country"
  )

  avail_vars <- intersect(key_vars, names(baseline))

  avail <- baseline %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      across(all_of(avail_vars),
             ~ round(mean(!is.na(.)) * 100, 1),
             .names = "pct_{.col}"),
      .groups = "drop"
    )

  # Transpose for readability
  avail_long <- avail %>%
    pivot_longer(
      cols = starts_with("pct_"),
      names_to = "variable",
      values_to = "pct_available"
    ) %>%
    mutate(variable = str_remove(variable, "^pct_")) %>%
    pivot_wider(
      names_from = cohort,
      values_from = pct_available,
      id_cols = variable
    )

  # Add overall column
  overall_pct <- baseline %>%
    summarise(
      across(all_of(avail_vars),
             ~ round(mean(!is.na(.)) * 100, 1))
    ) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "Overall")

  avail_long <- avail_long %>%
    left_join(overall_pct, by = "variable") %>%
    relocate(Overall, .after = variable)

  cat("  Variable availability (% non-missing at baseline):\n")
  print(avail_long, n = 40, width = 120)

  out_csv <- file.path(RESULTS_DIR, "variable_availability.csv")
  write_csv(avail_long, out_csv)
  cat("\n  Saved:", out_csv, "\n\n")

  avail_long
}

# ==============================================================================
# STEP 5: Outcome & Follow-Up Summary
# ==============================================================================

step5_followup <- function(pooled, baseline) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 5: Outcome and follow-up summary\n")
  cat("=" %>% strrep(70), "\n\n")

  # --- Waves per person by cohort ---
  cat("  5a: Follow-up waves per person:\n")
  fu <- pooled %>%
    group_by(cohort, id) %>%
    summarise(n_waves = n(), .groups = "drop") %>%
    group_by(cohort) %>%
    summarise(
      n_persons = n(),
      waves_mean = round(mean(n_waves), 1),
      waves_median = median(n_waves),
      waves_min = min(n_waves),
      waves_max = max(n_waves),
      .groups = "drop"
    )
  print(fu)

  # --- Mortality by cohort ---
  cat("\n  5b: Mortality by cohort:\n")
  mort <- baseline %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      n_died = sum(person_died, na.rm = TRUE),
      pct_died = round(mean(person_died, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    )
  print(mort)

  # --- ADL disability prevalence at baseline by cohort ---
  cat("\n  5c: Baseline ADL disability by cohort:\n")
  adl <- baseline %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      adl_any_pct = round(mean(adl_any, na.rm = TRUE) * 100, 1),
      iadl_any_pct = round(mean(iadl_any, na.rm = TRUE) * 100, 1),
      adl_mean = round(mean(adl_total, na.rm = TRUE), 2),
      .groups = "drop"
    )
  print(adl)

  # --- ADL disability by IC×Env group (key result preview) ---
  cat("\n  5d: Baseline ADL disability by IC×Env group:\n")
  adl_grp <- baseline %>%
    filter(!is.na(ic_env_6grp)) %>%
    group_by(ic_env_6grp) %>%
    summarise(
      n = n(),
      adl_any_pct = round(mean(adl_any, na.rm = TRUE) * 100, 1),
      iadl_any_pct = round(mean(iadl_any, na.rm = TRUE) * 100, 1),
      mort_pct = round(mean(person_died, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    )
  print(adl_grp)
  cat("  Expected: Low IC / Low Env should have highest disability & mortality.\n")

  # --- Mortality by IC×Env group (the main research question) ---
  cat("\n  5e: Mortality by IC×Env group (preview of interaction):\n")
  cat("  Comparing Low Env vs High Env WITHIN each IC level:\n")
  mort_interaction <- baseline %>%
    filter(!is.na(ic_cat) & !is.na(env_cat)) %>%
    group_by(ic_cat, env_cat) %>%
    summarise(
      n = n(),
      mort_pct = round(mean(person_died, na.rm = TRUE) * 100, 1),
      adl_pct = round(mean(adl_any, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(ic_cat, env_cat)
  print(mort_interaction)

  # Compute the "buffering effect" (difference in outcomes between Low/High Env)
  cat("\n  Environment buffering effect (Low Env - High Env):\n")
  buffering <- mort_interaction %>%
    select(ic_cat, env_cat, mort_pct, adl_pct) %>%
    pivot_wider(
      names_from = env_cat,
      values_from = c(mort_pct, adl_pct)
    ) %>%
    mutate(
      mort_diff = `mort_pct_Low` - `mort_pct_High`,
      adl_diff  = `adl_pct_Low` - `adl_pct_High`
    )
  print(buffering)
  cat("  If IC×Env interaction exists: buffering effect should be LARGER\n")
  cat("  for Low IC (environment compensates more for those with low capacity).\n")

  # Save summary
  summary_df <- bind_rows(
    mort %>% mutate(table = "mortality"),
    adl %>% mutate(table = "adl_disability")
  )
  out_csv <- file.path(RESULTS_DIR, "baseline_summary.csv")
  write_csv(summary_df, out_csv)
  cat("\n  Saved:", out_csv, "\n\n")

  invisible(list(mortality = mort, adl = adl, buffering = buffering))
}

# ==============================================================================
# STEP 6: Final Summary
# ==============================================================================

step6_summary <- function(baseline) {
  cat("=" %>% strrep(70), "\n")
  cat("DESCRIPTIVE ANALYSIS COMPLETE\n")
  cat("=" %>% strrep(70), "\n")
  cat("  Baseline N:        ", format(nrow(baseline), big.mark = ","), "\n")
  cat("  Cohorts:           ", paste(unique(baseline$cohort), collapse = ", "), "\n")
  cat("  Deaths:            ", sum(baseline$person_died, na.rm = TRUE),
      sprintf("(%.1f%%)\n", mean(baseline$person_died, na.rm = TRUE) * 100))
  cat("  ADL disability:    ", sum(baseline$adl_any, na.rm = TRUE),
      sprintf("(%.1f%%)\n", mean(baseline$adl_any, na.rm = TRUE) * 100))
  cat("  IC×Env classified: ",
      sum(!is.na(baseline$ic_env_6grp)),
      sprintf("(%.1f%%)\n", mean(!is.na(baseline$ic_env_6grp)) * 100))
  cat("\n  Output files:\n")
  cat("    results/table1_by_cohort.csv\n")
  cat("    results/table1_by_icenv_group.csv\n")
  cat("    results/variable_availability.csv\n")
  cat("    results/baseline_summary.csv\n")
  cat("=" %>% strrep(70), "\n")
}

# ==============================================================================
# MAIN
# ==============================================================================

main <- function() {
  data <- step1_baseline()
  pooled   <- data$pooled
  baseline <- data$baseline

  step2_table1_cohort(baseline)
  step3_table1_icenv(baseline)
  step4_availability(baseline)
  step5_followup(pooled, baseline)
  step6_summary(baseline)
}

if (interactive() || !exists(".main05_called")) {
  .main05_called <- TRUE
  main()
}

# ==============================================================================
# NOTES
# ==============================================================================
#
# 1. BASELINE DEFINITION:
#    First observed wave per person (slice_min(wave)). This preserves the
#    largest possible sample. Alternative: first wave with complete IC data.
#
# 2. PERSON_DIED:
#    This is a person-level indicator (ever died during follow-up).
#    It is time-invariant across waves. The actual event time is
#    death_year (available for ~21% who died).
#
# 3. TABLE FORMAT:
#    Mean (SD) for continuous variables.
#    N (%) for categorical/binary variables.
#    Section headers separate variable domains.
#
# 4. IC×ENV BUFFERING PREVIEW:
#    Step 5e provides a crude preview of the interaction hypothesis.
#    If environment buffers IC decline, then:
#    - Low Env vs High Env mortality gap should be LARGER for Low IC
#    - High Env should "protect" Low IC individuals more than High IC
#    This cross-tabulation is DESCRIPTIVE only; confounding is not controlled.
#    Formal testing occurs in 06_ipd_ma.R (Cox models with interaction term).
# ==============================================================================
