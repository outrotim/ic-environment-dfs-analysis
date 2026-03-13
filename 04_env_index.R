# ==============================================================================
# 04_env_index.R
# Study 17: IC-Environment Fit and Disability-Free Survival in Older Adults
#
# Purpose: Construct Environment Support Index (formative composite)
#   (1) Review existing environment indicators from 02_harmonization.R
#   (2) PCA diagnostics: dimensionality and variance structure
#   (3) Build formative composite: equal-weight z-score average
#   (4) Create primary & sensitivity versions (with/without social_contact)
#   (5) Update IC×Env interaction using ic_cfa (from 03) × env_formative
#   (6) Validation: distributions, correlations, cohort-specific checks
#
# Methodology: Formative (causal) composite indicator
#   - Env support is a formative construct: indicators CAUSE the concept
#   - Social participation, financial capacity, living arrangement, social contact
#     are distinct enablers, not reflective manifestations of a latent trait
#   - CFA (reflective model) is inappropriate; equal-weight sum is standard
#   - Justified by: Bollen & Lennox 1991; Diamantopoulos & Winklhofer 2001;
#     WHO ICOPE handbook (environment as multi-dimensional enabler)
#
# Input:   data/ic_cfa_scores.rds (from 03_ic_cfa.R, contains env vars from 02)
# Output:  data/analytic_final.rds    — final analytic dataset with env + interaction
#          results/env_validation.csv  — env composite descriptive statistics
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
# STEP 1: Load Data and Review Environment Indicators
# ==============================================================================

step1_review <- function() {
  cat("=" %>% strrep(70), "\n")
  cat("Step 1: Load data and review environment indicators\n")
  cat("=" %>% strrep(70), "\n\n")

  rds_path <- file.path(DATA_DIR, "ic_cfa_scores.rds")
  if (!file.exists(rds_path)) {
    stop("ic_cfa_scores.rds not found. Run 03_ic_cfa.R first.")
  }

  pooled <- readRDS(rds_path)
  cat("  Loaded:", format(nrow(pooled), big.mark = ","), "rows,",
      n_distinct(pooled$id), "persons,", ncol(pooled), "columns\n\n")

  # --- Identify available environment z-score columns ---
  env_z_vars <- c("z_soc_part", "z_financial", "z_living", "z_social_contact")
  env_raw_vars <- c("social_participation", "financial_strain",
                     "living_arrangement", "social_contact")

  cat("  Environment indicator coverage (all person-waves):\n")
  cov_all <- pooled %>%
    summarise(
      across(all_of(env_z_vars[env_z_vars %in% names(pooled)]),
             ~ round(mean(!is.na(.)) * 100, 1),
             .names = "pct_{.col}")
    )
  print(cov_all)

  cat("\n  Coverage by cohort:\n")
  cov_cohort <- pooled %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      across(all_of(env_z_vars[env_z_vars %in% names(pooled)]),
             ~ round(mean(!is.na(.)) * 100, 1),
             .names = "pct_{.col}"),
      .groups = "drop"
    )
  print(cov_cohort, width = 120)

  # --- Identify the social_contact gap ---
  cat("\n  CRITICAL GAP: z_social_contact availability:\n")
  pooled %>%
    group_by(cohort) %>%
    summarise(
      pct_social_contact = round(mean(!is.na(z_social_contact)) * 100, 1),
      .groups = "drop"
    ) %>%
    print()
  cat("  → Only HRS and MHAS have social contact data.\n")
  cat("  → Primary Env index will use 3 indicators (excl. social_contact).\n")
  cat("  → Sensitivity: 4-indicator version for HRS/MHAS.\n\n")

  pooled
}

# ==============================================================================
# STEP 2: PCA Diagnostics (Informational Only)
# ==============================================================================

step2_pca_diagnostics <- function(pooled) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 2: PCA diagnostics (dimensionality check)\n")
  cat("=" %>% strrep(70), "\n\n")

  # Use 3 core indicators available across all cohorts
  pca_vars <- c("z_soc_part", "z_financial", "z_living")
  avail_vars <- intersect(pca_vars, names(pooled))

  if (length(avail_vars) < 3) {
    cat("  WARNING: Fewer than 3 environment z-score variables available.\n")
    cat("  Skipping PCA diagnostics.\n")
    return(invisible(NULL))
  }

  # Complete cases for PCA
  pca_data <- pooled %>%
    select(cohort, all_of(avail_vars)) %>%
    drop_na()

  cat("  PCA sample (complete cases on 3 indicators):",
      format(nrow(pca_data), big.mark = ","), "\n\n")

  # Overall PCA
  pca_mat <- as.matrix(pca_data[, avail_vars])
  pca_res <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

  cat("  Overall PCA (3 core env indicators):\n")
  var_prop <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  var_cum  <- cumsum(var_prop)
  pca_summary <- tibble(
    PC = paste0("PC", seq_along(var_prop)),
    eigenvalue = round(pca_res$sdev^2, 3),
    var_pct = round(var_prop * 100, 1),
    cum_pct = round(var_cum * 100, 1)
  )
  print(pca_summary)

  cat("\n  Loadings:\n")
  loadings_df <- as_tibble(pca_res$rotation, rownames = "variable")
  print(loadings_df)

  # PCA by cohort
  cat("\n  PCA PC1 variance explained by cohort:\n")
  cohort_pca <- pca_data %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      pc1_var_pct = {
        mat <- as.matrix(pick(all_of(avail_vars)))
        pca <- prcomp(mat, center = TRUE, scale. = TRUE)
        round(pca$sdev[1]^2 / sum(pca$sdev^2) * 100, 1)
      },
      .groups = "drop"
    )
  print(cohort_pca)

  cat("\n  Interpretation:\n")
  if (pca_summary$var_pct[1] < 50) {
    cat("  → PC1 explains < 50% variance: indicators are multidimensional.\n")
    cat("  → This SUPPORTS the formative (not reflective) approach.\n")
    cat("  → Equal-weight composite is appropriate.\n\n")
  } else {
    cat("  → PC1 explains >= 50% variance: moderate unidimensionality.\n")
    cat("  → Either formative or reflective approach could work.\n")
    cat("  → Using formative for theoretical consistency.\n\n")
  }

  invisible(pca_res)
}

# ==============================================================================
# STEP 3: Construct Formative Environment Composite
# ==============================================================================

step3_env_composite <- function(pooled) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 3: Constructing formative Environment Support Index\n")
  cat("=" %>% strrep(70), "\n\n")

  # --- z-score helper (within cohort-wave) ---
  z_func <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - m) / s
  }

  # --- Re-standardize env indicators within cohort-wave ---
  # The z_* vars from 02_harmonization.R may have been computed before
  # some data cleaning in 03_ic_cfa.R. Re-compute for consistency.
  cat("  Re-standardizing environment indicators within cohort-wave...\n")

  # Check which raw variables exist for re-standardization
  has_soc_part <- "social_participation" %in% names(pooled)
  has_fin_strain <- "financial_strain" %in% names(pooled)
  has_living <- "living_arrangement" %in% names(pooled)
  has_contact <- "social_contact" %in% names(pooled)

  pooled <- pooled %>%
    group_by(cohort, wave) %>%
    mutate(
      # Social participation: count → z-score (higher = more participation)
      z_env_soc = if (has_soc_part) z_func(social_participation) else NA_real_,

      # Financial capacity: reverse strain (higher = better off)
      # financial_strain: 5=most strained → negate so higher=better
      z_env_fin = if (has_fin_strain) z_func(-financial_strain) else NA_real_,

      # Living support: ordinal recoded → z-score
      # 0=alone, 1=with spouse(best), 2=with others(middle)
      z_env_liv = if (has_living) {
        z_func(case_when(
          living_arrangement == 0 ~ 0,    # alone (least support)
          living_arrangement == 2 ~ 1,    # with others (some support)
          living_arrangement == 1 ~ 2,    # with spouse (most support)
          TRUE ~ NA_real_
        ))
      } else {
        NA_real_
      },

      # Social contact: count of contact types → z-score (higher = more)
      z_env_con = if (has_contact) z_func(social_contact) else NA_real_
    ) %>%
    ungroup()

  # --- PRIMARY: 3-indicator Environment Index (all 5 cohorts) ---
  cat("  Building primary Env index (3 indicators, all cohorts)...\n")

  env_3_cols <- c("z_env_soc", "z_env_fin", "z_env_liv")
  env_3_mat <- as.matrix(pooled[, env_3_cols])
  n_avail_3 <- rowSums(!is.na(env_3_mat))

  pooled$env_formative <- ifelse(
    n_avail_3 >= 2,   # Require at least 2/3 indicators
    rowMeans(env_3_mat, na.rm = TRUE),
    NA_real_
  )

  cat(sprintf("  env_formative (3-ind): %.1f%% available (≥2/3 rule)\n",
              mean(!is.na(pooled$env_formative)) * 100))

  # --- SENSITIVITY: 4-indicator version (with social contact) ---
  cat("  Building sensitivity Env index (4 indicators, with contact)...\n")

  env_4_cols <- c("z_env_soc", "z_env_fin", "z_env_liv", "z_env_con")
  env_4_mat <- as.matrix(pooled[, env_4_cols])
  n_avail_4 <- rowSums(!is.na(env_4_mat))

  pooled$env_formative_4 <- ifelse(
    n_avail_4 >= 2,   # Require at least 2/4
    rowMeans(env_4_mat, na.rm = TRUE),
    NA_real_
  )

  cat(sprintf("  env_formative_4 (4-ind): %.1f%% available (≥2/4 rule)\n",
              mean(!is.na(pooled$env_formative_4)) * 100))

  # --- Coverage by cohort ---
  cat("\n  Environment composite coverage by cohort:\n")
  pooled %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      env_form_3_pct = round(mean(!is.na(env_formative)) * 100, 1),
      env_form_4_pct = round(mean(!is.na(env_formative_4)) * 100, 1),
      n_ind_mean = round(mean(n_avail_3), 1),
      .groups = "drop"
    ) %>%
    print()

  # --- Correlation between primary and sensitivity versions ---
  valid_both <- pooled %>%
    filter(!is.na(env_formative) & !is.na(env_formative_4))
  if (nrow(valid_both) > 100) {
    r_env <- cor(valid_both$env_formative, valid_both$env_formative_4)
    cat(sprintf("\n  Correlation env_formative vs env_formative_4: r = %.4f (n=%s)\n",
                r_env, format(nrow(valid_both), big.mark = ",")))
  }

  # --- Correlation with old env_sum from 02 ---
  if ("env_sum" %in% names(pooled)) {
    valid_old <- pooled %>%
      filter(!is.na(env_formative) & !is.na(env_sum))
    if (nrow(valid_old) > 100) {
      r_old <- cor(valid_old$env_formative, valid_old$env_sum)
      cat(sprintf("  Correlation env_formative vs env_sum (02): r = %.4f (n=%s)\n",
                  r_old, format(nrow(valid_old), big.mark = ",")))
    }
  }

  cat("\n")
  pooled
}

# ==============================================================================
# STEP 4: Update IC × Env Interaction Variables
# ==============================================================================

step4_interaction <- function(pooled) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 4: Updating IC × Env interaction variables\n")
  cat("=" %>% strrep(70), "\n\n")

  # Verify ic_cfa exists (from 03_ic_cfa.R)
  if (!"ic_cfa" %in% names(pooled)) {
    warning("ic_cfa not found. Using ic_sum as fallback.")
    pooled$ic_cfa <- pooled$ic_sum
  }

  # --- 4a: Continuous interaction (mean-centered within cohort) ---
  cat("  4a: Continuous interaction (ic_cfa × env_formative)...\n")
  pooled <- pooled %>%
    group_by(cohort) %>%
    mutate(
      ic_cfa_c  = ic_cfa - mean(ic_cfa, na.rm = TRUE),
      env_form_c = env_formative - mean(env_formative, na.rm = TRUE),
      ic_env_cont = ic_cfa_c * env_form_c
    ) %>%
    ungroup()

  cat(sprintf("  ic_env_cont (continuous interaction): %.1f%% available\n",
              mean(!is.na(pooled$ic_env_cont)) * 100))

  # --- 4b: Categorical IC × Env groups ---
  cat("  4b: Categorical grouping (IC tertile × Env median split)...\n")
  pooled <- pooled %>%
    group_by(cohort) %>%
    mutate(
      # IC tertiles within cohort (using ic_cfa)
      ic_cat = case_when(
        is.na(ic_cfa) ~ NA_character_,
        ic_cfa <= quantile(ic_cfa, 1/3, na.rm = TRUE) ~ "Low",
        ic_cfa <= quantile(ic_cfa, 2/3, na.rm = TRUE) ~ "Medium",
        TRUE ~ "High"
      ),
      # Env median split within cohort
      env_cat = case_when(
        is.na(env_formative) ~ NA_character_,
        env_formative <= median(env_formative, na.rm = TRUE) ~ "Low",
        TRUE ~ "High"
      ),
      # 6-group cross-classification
      ic_env_6grp = case_when(
        is.na(ic_cat) | is.na(env_cat) ~ NA_character_,
        TRUE ~ paste0(ic_cat, " IC / ", env_cat, " Env")
      )
    ) %>%
    ungroup()

  # --- Set factor levels in meaningful order ---
  grp_levels <- c(
    "High IC / High Env", "High IC / Low Env",
    "Medium IC / High Env", "Medium IC / Low Env",
    "Low IC / High Env", "Low IC / Low Env"
  )
  pooled$ic_env_6grp <- factor(pooled$ic_env_6grp, levels = grp_levels)
  pooled$ic_cat <- factor(pooled$ic_cat, levels = c("High", "Medium", "Low"))
  pooled$env_cat <- factor(pooled$env_cat, levels = c("High", "Low"))

  # Report group distribution
  cat("\n  IC × Env 6-group distribution:\n")
  pooled %>%
    filter(!is.na(ic_env_6grp)) %>%
    count(ic_env_6grp, .drop = FALSE) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    print()

  cat(sprintf("\n  6-group coverage: %.1f%%\n",
              mean(!is.na(pooled$ic_env_6grp)) * 100))

  # Report by cohort
  cat("\n  6-group coverage by cohort:\n")
  pooled %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      grp_pct = round(mean(!is.na(ic_env_6grp)) * 100, 1),
      .groups = "drop"
    ) %>%
    print()

  cat("\n")
  pooled
}

# ==============================================================================
# STEP 5: Validation
# ==============================================================================

step5_validation <- function(pooled) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 5: Environment composite validation\n")
  cat("=" %>% strrep(70), "\n\n")

  # --- 5a: Descriptive statistics by cohort ---
  cat("  5a: env_formative descriptive statistics by cohort:\n")
  desc <- pooled %>%
    filter(!is.na(env_formative)) %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      mean = round(mean(env_formative), 3),
      sd   = round(sd(env_formative), 3),
      min  = round(min(env_formative), 3),
      p25  = round(quantile(env_formative, 0.25), 3),
      median = round(median(env_formative), 3),
      p75  = round(quantile(env_formative, 0.75), 3),
      max  = round(max(env_formative), 3),
      .groups = "drop"
    )
  print(desc)

  # --- 5b: IC-Env correlation by cohort ---
  cat("\n  5b: Correlation between ic_cfa and env_formative:\n")
  ic_env_cor <- pooled %>%
    filter(!is.na(ic_cfa) & !is.na(env_formative)) %>%
    group_by(cohort) %>%
    summarise(
      n = n(),
      r_ic_env = round(cor(ic_cfa, env_formative), 3),
      .groups = "drop"
    )
  print(ic_env_cor)

  overall_r <- pooled %>%
    filter(!is.na(ic_cfa) & !is.na(env_formative)) %>%
    {cor(.$ic_cfa, .$env_formative)}
  cat(sprintf("  Overall r(ic_cfa, env_formative) = %.3f\n", overall_r))

  # --- 5c: Env composite by known-group validity ---
  cat("\n  5c: Env composite by living_arrangement (known-group validity):\n")
  if ("living_arrangement" %in% names(pooled)) {
    kg <- pooled %>%
      filter(!is.na(env_formative) & !is.na(living_arrangement)) %>%
      group_by(living_arrangement) %>%
      summarise(
        n = n(),
        mean_env = round(mean(env_formative), 3),
        sd_env   = round(sd(env_formative), 3),
        .groups = "drop"
      ) %>%
      mutate(label = case_when(
        living_arrangement == 0 ~ "Alone",
        living_arrangement == 1 ~ "With spouse",
        living_arrangement == 2 ~ "With others"
      ))
    print(kg)
    cat("  Expected: With spouse > With others > Alone\n")
  }

  # --- 5d: Env composite by age group ---
  cat("\n  5d: Env composite by age group:\n")
  age_env <- pooled %>%
    filter(!is.na(env_formative) & !is.na(age)) %>%
    mutate(age_grp = cut(age, breaks = c(59, 64, 69, 74, 79, 84, 120),
                          labels = c("60-64", "65-69", "70-74",
                                     "75-79", "80-84", "85+"))) %>%
    filter(!is.na(age_grp)) %>%
    group_by(age_grp) %>%
    summarise(
      n = n(),
      mean_env = round(mean(env_formative), 3),
      sd_env   = round(sd(env_formative), 3),
      mean_ic  = round(mean(ic_cfa, na.rm = TRUE), 3),
      .groups = "drop"
    )
  print(age_env)
  cat("  Expected: Env generally declines with age (fewer social ties, income)\n")

  # --- 5e: Component inter-correlations ---
  cat("\n  5e: Inter-correlations among env components:\n")
  env_comp_vars <- c("z_env_soc", "z_env_fin", "z_env_liv")
  avail_comp <- intersect(env_comp_vars, names(pooled))

  if (length(avail_comp) >= 2) {
    comp_data <- pooled %>% select(all_of(avail_comp)) %>% drop_na()
    cor_mat <- round(cor(comp_data), 3)
    print(cor_mat)
    cat("  Note: Low inter-correlations support the formative approach.\n")
    cat("  (Formative indicators need not correlate; reflective indicators should.)\n")
  }

  # --- 5f: Event-rate validation (if died/disability available) ---
  cat("\n  5f: Mortality rate by IC × Env group:\n")
  if ("died" %in% names(pooled)) {
    mort_grp <- pooled %>%
      filter(!is.na(ic_env_6grp) & !is.na(died)) %>%
      group_by(ic_env_6grp) %>%
      summarise(
        n = n(),
        n_died = sum(died, na.rm = TRUE),
        mort_pct = round(mean(died, na.rm = TRUE) * 100, 1),
        .groups = "drop"
      )
    print(mort_grp)
    cat("  Expected gradient: Low IC / Low Env should have highest mortality.\n")
  } else {
    cat("  (died variable not available; skipping event-rate validation)\n")
  }

  # --- Save validation results ---
  validation <- list(
    descriptives = desc,
    ic_env_cor = ic_env_cor,
    overall_r = overall_r
  )

  val_csv <- file.path(RESULTS_DIR, "env_validation.csv")
  write_csv(desc, val_csv)
  cat("\n  Saved:", val_csv, "\n\n")

  invisible(validation)
}

# ==============================================================================
# STEP 6: Save Final Analytic Dataset
# ==============================================================================

step6_save <- function(pooled) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 6: Saving final analytic dataset\n")
  cat("=" %>% strrep(70), "\n\n")

  out_rds <- file.path(DATA_DIR, "analytic_final.rds")
  saveRDS(pooled, out_rds)
  cat("  Saved:", out_rds, "\n")

  # --- Final summary ---
  cat("\n", "=" %>% strrep(70), "\n")
  cat("ENVIRONMENT INDEX CONSTRUCTION COMPLETE\n")
  cat("=" %>% strrep(70), "\n")
  cat("  Rows:             ", format(nrow(pooled), big.mark = ","), "\n")
  cat("  Persons:          ", format(n_distinct(pooled$id), big.mark = ","), "\n")
  cat("  Total columns:    ", ncol(pooled), "\n")
  cat("  ic_cfa avail:     ", round(mean(!is.na(pooled$ic_cfa)) * 100, 1), "%\n")
  cat("  env_formative:    ", round(mean(!is.na(pooled$env_formative)) * 100, 1), "%\n")
  cat("  env_formative_4:  ", round(mean(!is.na(pooled$env_formative_4)) * 100, 1), "%\n")
  cat("  ic_env_cont:      ", round(mean(!is.na(pooled$ic_env_cont)) * 100, 1), "%\n")
  cat("  ic_env_6grp:      ", round(mean(!is.na(pooled$ic_env_6grp)) * 100, 1), "%\n")
  cat("=" %>% strrep(70), "\n")

  # --- Column inventory ---
  new_cols <- c("z_env_soc", "z_env_fin", "z_env_liv", "z_env_con",
                "env_formative", "env_formative_4",
                "ic_cfa_c", "env_form_c", "ic_env_cont",
                "ic_cat", "env_cat", "ic_env_6grp")
  cat("\n  New columns added by 04_env_index.R:\n")
  for (col in new_cols) {
    if (col %in% names(pooled)) {
      cat(sprintf("    %-20s  %s available\n", col,
                  sprintf("%.1f%%", mean(!is.na(pooled[[col]])) * 100)))
    }
  }

  invisible(pooled)
}

# ==============================================================================
# MAIN
# ==============================================================================

main <- function() {
  pooled <- step1_review()
  step2_pca_diagnostics(pooled)
  pooled <- step3_env_composite(pooled)
  pooled <- step4_interaction(pooled)
  step5_validation(pooled)
  step6_save(pooled)
}

if (interactive() || !exists(".main04_called")) {
  .main04_called <- TRUE
  main()
}

# ==============================================================================
# NOTES & DESIGN RATIONALE
# ==============================================================================
#
# 1. FORMATIVE vs REFLECTIVE:
#    Environment support is a formative construct. Its components
#    (social participation, financial capacity, living arrangement)
#    are conceptually distinct enablers that do not necessarily co-vary.
#    A person can have high social participation but live alone and have
#    low income. CFA (which assumes reflective measurement) would
#    incorrectly require these indicators to share a common cause.
#    References: Bollen & Lennox (1991, Psych Bull);
#    Diamantopoulos & Winklhofer (2001, JMR).
#
# 2. INDICATOR SELECTION:
#    Primary (3 indicators, all 5 cohorts):
#      - Social participation (count of activities, z-scored)
#      - Financial capacity (reversed income quintile, z-scored)
#      - Living support (ordinal: alone < others < spouse, z-scored)
#    Sensitivity (4 indicators, HRS/MHAS only for full version):
#      - Above 3 + social contact (count of contact types, z-scored)
#    Social contact excluded from primary due to >60% missing in
#    ELSA, CHARLS, SHARE. Sensitivity analysis assesses impact.
#
# 3. WEIGHTING:
#    Equal weights used. Justifications:
#    - No theoretical basis to weight one enabler over another
#    - WHO ICOPE framework treats all environment domains equally
#    - PCA-derived weights would be sample-specific and less generalizable
#    - Equal weights are robust across different populations
#    Reference: Wainer (1976, "Estimating coefficients in linear models:
#    it don't make no nevermind").
#
# 4. IC × ENV INTERACTION:
#    Two specifications for the core research question:
#    - Continuous: ic_cfa_c × env_form_c (mean-centered within cohort)
#    - Categorical: IC tertile × Env median → 6 groups
#    Centering prevents multicollinearity in Cox models.
#    Categorical groups enable non-parametric visualization of the
#    interaction pattern and clinical interpretation.
#
# 5. DATA FLOW:
#    01_data_load.R → analytic_pooled.rds (36 cols)
#    02_harmonization.R → harmonized_pooled.rds (79 cols)
#    03_ic_cfa.R → ic_cfa_scores.rds (85 cols)
#    04_env_index.R → analytic_final.rds (~97 cols)  ← THIS SCRIPT
#    05_descriptive.R → Table 1, figures
#    06_ipd_ma.R → Main analysis results
#
# 6. MISSING DATA STRATEGY:
#    env_formative requires ≥2/3 core indicators present.
#    When 1 indicator is missing, the composite is the mean of 2.
#    This is a person-mean imputation for the composite score.
#    Full multiple imputation is handled at the analysis stage (06_ipd_ma.R).
# ==============================================================================
