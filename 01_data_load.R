# ==============================================================================
# 01_data_load.R
# Study 17: IC-Environment Fit and Disability-Free Survival in Older Adults
#
# Purpose: Load 5 core cohorts (HRS, ELSA, CHARLS, SHARE, MHAS) from
#          pre-cleaned CSV files, select/rename core variables to a common
#          schema, filter eligible participants (≥60y, community-dwelling),
#          compute derived variables, and output a unified long-format dataset.
#
# Input:   Pre-cleaned CSV files at DATA_DIR (long format, one row per
#          person-wave). Headers formatted as "varname (中文描述)".
# Output:  data/analytic_pooled.rds  — unified long-format analytic dataset
#          data/cohort_summary.csv   — per-cohort variable availability summary
#
# Author:  Study 17 Team
# Date:    2026-03-11
# Version: 1.0
# ==============================================================================

# --- Dependencies -------------------------------------------------------------
library(tidyverse)

# --- Configuration ------------------------------------------------------------

DATA_DIR <- file.path(here::here(), "data", "raw")  # Place raw CSV files here
OUTPUT_DIR <- tryCatch(
  file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data"),
  error = function(e) {
    # Fallback: use script's own directory (reliable across environments)
    script_dir <- here::here()
    file.path(script_dir, "data")
  }
)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

AGE_MIN <- 60L   # inclusion criterion

# CSV files for each cohort
COHORT_FILES <- c(
  HRS   = "hrs.csv",
  ELSA  = "elsa.csv",
  CHARLS = "charls.csv",
  SHARE = "share.csv",
  MHAS  = "mhas.csv"
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Clean CSV column headers
#' CSV headers are formatted as "variable_name (中文描述)" — extract var name only
clean_headers <- function(df) {
  names(df) <- str_extract(names(df), "^[^\\s(]+")
  df
}

#' Safely extract a numeric column (returns NA if column doesn't exist)
safe_num <- function(df, varname) {
  if (varname %in% names(df)) {
    suppressWarnings(as.numeric(df[[varname]]))
  } else {
    rep(NA_real_, nrow(df))
  }
}

#' Safely extract a character column
safe_chr <- function(df, varname) {
  if (varname %in% names(df)) {
    as.character(df[[varname]])
  } else {
    rep(NA_character_, nrow(df))
  }
}

#' Detect female from ragender variable (handles multiple encodings)
#' Gateway to Global Aging Data uses 0=female, 1=male
#' RAND HRS convention: 1=male, 2=female
detect_female <- function(x) {
  x <- tolower(as.character(x))
  case_when(
    x %in% c("0", "2", "2.female", "female") ~ 1L,
    x %in% c("1", "1.male", "male")          ~ 0L,
    TRUE ~ NA_integer_
  )
}

#' Compute chronic condition count from binary disease indicators
#' Standard 8 conditions: hypertension, diabetes, cancer, lung disease,
#' heart disease, stroke, psychiatric problems, arthritis
compute_chronic_count <- function(df) {
  vars <- c("hibpe", "diabe", "cancre", "lunge", "hearte",
            "stroke", "psyche", "arthre")
  avail <- intersect(vars, names(df))
  if (length(avail) == 0) return(rep(NA_real_, nrow(df)))

  mat <- df %>%
    select(all_of(avail)) %>%
    mutate(across(everything(), ~ {
      x <- tolower(as.character(.))
      case_when(
        x %in% c("1", "1.yes", "yes") ~ 1L,
        x %in% c("0", "5", "0.no", "5.not applicable", "no") ~ 0L,
        TRUE ~ NA_integer_
      )
    }))

  # Count non-NA conditions; set to NA if all are missing
  counts <- rowSums(mat, na.rm = TRUE)
  all_na <- rowSums(!is.na(mat)) == 0
  counts[all_na] <- NA_real_
  counts
}

#' Compute maximum grip strength from left/right hand measurements
compute_grip_max <- function(df) {
  l <- safe_num(df, "lgrip")
  r <- safe_num(df, "rgrip")
  pmax(l, r, na.rm = TRUE)
}

#' Compute ADL total from 6 individual binary items
#' Items: walking across room, dressing, bathing, eating, getting in/out of bed, toileting
compute_adl6 <- function(df) {
  items <- c("walkra", "dressa", "batha", "eata", "beda", "toilta")
  avail <- intersect(items, names(df))
  if (length(avail) < 4) {
    warning("Only ", length(avail), " of 6 ADL items found; returning NA")
    return(rep(NA_real_, nrow(df)))
  }

  mat <- df %>%
    select(all_of(avail)) %>%
    mutate(across(everything(), ~ {
      x <- as.numeric(.)
      if_else(!is.na(x) & x > 0, 1L, 0L)
    }))

  counts <- rowSums(mat, na.rm = TRUE)
  all_na <- rowSums(!is.na(mat)) == 0
  counts[all_na] <- NA_real_
  counts
}

#' Compute IADL total from specified items
compute_iadl <- function(df, items) {
  avail <- intersect(items, names(df))
  if (length(avail) < 3) {
    warning("Only ", length(avail), " IADL items found; returning NA")
    return(rep(NA_real_, nrow(df)))
  }

  mat <- df %>%
    select(all_of(avail)) %>%
    mutate(across(everything(), ~ {
      x <- as.numeric(.)
      if_else(!is.na(x) & x > 0, 1L, 0L)
    }))

  counts <- rowSums(mat, na.rm = TRUE)
  all_na <- rowSums(!is.na(mat)) == 0
  counts[all_na] <- NA_real_
  counts
}

#' Count social participation activities (binary: participated or not)
count_social_activities <- function(df, items) {
  avail <- intersect(items, names(df))
  if (length(avail) == 0) return(rep(NA_real_, nrow(df)))

  mat <- df %>%
    select(all_of(avail)) %>%
    mutate(across(everything(), ~ {
      x <- as.numeric(.)
      if_else(!is.na(x) & x > 0, 1L, 0L)
    }))

  counts <- rowSums(mat, na.rm = TRUE)
  all_na <- rowSums(!is.na(mat)) == 0
  counts[all_na] <- NA_real_
  counts
}

# ==============================================================================
# COHORT-SPECIFIC LOADING FUNCTIONS
# ==============================================================================
# Each function reads the raw CSV, maps cohort-specific variable names to a
# common schema, and returns a tibble with standardized column names.
#
# Common schema columns:
#   cohort, id, wave, age, female, edu_years, marital, interview_year,
#   imrc, dlrc, ser7, depression, dep_scale, dep_max,
#   vision, hearing, grip_max, gait_speed, bmi, srh,
#   adl_total, iadl_total, iwstat, death_year,
#   chronic_n, loneliness, social_participation,
#   country (SHARE only)
# ==============================================================================

load_hrs <- function() {
  message("  Loading HRS...")
  raw <- read_csv(file.path(DATA_DIR, "hrs.csv"), show_col_types = FALSE) %>%
    clean_headers()

  tibble(
    cohort      = "HRS",
    id          = as.character(raw$hhidpn),
    wave        = as.integer(raw$wave),
    age         = safe_num(raw, "ragey_e"),
    female      = detect_female(raw$ragender),
    edu_years   = safe_num(raw, "raedyrs"),
    marital     = safe_chr(raw, "mstath"),
    interview_year = safe_num(raw, "iwendy"),
    # --- IC: Cognition ---
    imrc        = safe_num(raw, "imrc"),        # 0-10
    dlrc        = safe_num(raw, "dlrc"),        # 0-10
    ser7        = safe_num(raw, "ser7"),        # 0-5
    # --- IC: Psychological ---
    depression  = safe_num(raw, "cesd"),        # CES-D 8-item, 0-8
    dep_scale   = "CES-D 8",
    dep_max     = 8L,
    # --- IC: Sensory ---
    vision      = safe_num(raw, "sight"),       # 1-5 (1=excellent)
    hearing     = safe_num(raw, "hearing"),     # 1-5
    # --- IC: Locomotion ---
    grip_max    = compute_grip_max(raw),        # kg
    gait_speed  = safe_num(raw, "wspeed"),      # seconds (varies by distance)
    # --- IC: Vitality ---
    bmi         = safe_num(raw, "bmi"),
    srh         = safe_num(raw, "shlt"),        # 1-5 (1=excellent)
    # --- Outcomes ---
    adl_total   = safe_num(raw, "adl6a"),       # 0-6
    iadl_total  = safe_num(raw, "iadl5a"),      # 0-5
    iwstat      = safe_chr(raw, "iwstat"),
    death_year  = safe_num(raw, "radyear"),
    # --- Covariates ---
    chronic_n   = compute_chronic_count(raw),
    # --- Environment ---
    loneliness  = safe_num(raw, "lblonely3"),   # UCLA-3
    social_participation = NA_real_,            # TODO: map HRS social vars
    # --- Meta ---
    country     = "US"
  )
}

load_elsa <- function() {
  message("  Loading ELSA...")
  raw <- read_csv(file.path(DATA_DIR, "elsa.csv"), show_col_types = FALSE) %>%
    clean_headers()

  # ELSA has no interview year variable; derive from wave number
  # W1=2002, W2=2004, W3=2006, ..., W10=2020
  wave_year_map <- c(
    `1` = 2002L, `2` = 2004L, `3` = 2006L, `4` = 2008L, `5` = 2010L,
    `6` = 2012L, `7` = 2014L, `8` = 2016L, `9` = 2018L, `10` = 2020L
  )

  tibble(
    cohort      = "ELSA",
    id          = as.character(raw$idauniqc),
    wave        = as.integer(raw$wave),
    age         = safe_num(raw, "agey"),
    female      = detect_female(raw$ragender),
    edu_years   = safe_num(raw, "raedyrs_e"),
    marital     = safe_chr(raw, "mstath"),
    interview_year = wave_year_map[as.character(as.integer(raw$wave))],
    # --- IC: Cognition ---
    imrc        = safe_num(raw, "imrc"),
    dlrc        = safe_num(raw, "dlrc"),
    ser7        = safe_num(raw, "ser7"),
    # --- IC: Psychological ---
    depression  = safe_num(raw, "cesd"),        # CES-D 8-item
    dep_scale   = "CES-D 8",
    dep_max     = 8L,
    # --- IC: Sensory ---
    vision      = safe_num(raw, "sight"),
    hearing     = safe_num(raw, "hearing"),
    # --- IC: Locomotion ---
    grip_max    = {
      gs <- safe_num(raw, "gripsum")            # preferred: dominant hand max
      gm <- compute_grip_max(raw)               # fallback: max(left, right)
      if_else(is.na(gs), gm, gs)
    },
    gait_speed  = safe_num(raw, "wspeed"),
    # --- IC: Vitality ---
    bmi         = safe_num(raw, "mbmi"),         # measured BMI
    srh         = safe_num(raw, "shlt"),
    # --- Outcomes ---
    adl_total   = safe_num(raw, "adltot6"),
    iadl_total  = safe_num(raw, "iadltot2_e"),
    iwstat      = safe_chr(raw, "iwstat"),
    death_year  = safe_num(raw, "radyear"),
    # --- Covariates ---
    chronic_n   = compute_chronic_count(raw),
    # --- Environment ---
    loneliness  = safe_num(raw, "lnlys3"),      # UCLA-3 (3-item average)
    social_participation = count_social_activities(
      raw, paste0("group", 1:8)                  # 8 org types
    ),
    # --- Meta ---
    country     = "GB",
    # --- Exclusion ---
    nursing_home = safe_num(raw, "nhmliv")       # 1=nursing home
  )
}

load_charls <- function() {
  message("  Loading CHARLS...")
  raw <- read_csv(file.path(DATA_DIR, "charls.csv"), show_col_types = FALSE) %>%
    clean_headers()

  tibble(
    cohort      = "CHARLS",
    id          = as.character(raw$ID),
    wave        = as.integer(raw$wave),
    age         = safe_num(raw, "age"),
    female      = detect_female(raw$ragender),
    edu_years   = NA_real_,                     # raedyrs not in CSV; raeduc_c is categorical
    marital     = safe_chr(raw, "marry"),
    interview_year = safe_num(raw, "iwy"),
    # --- IC: Cognition ---
    imrc        = safe_num(raw, "imrc"),
    dlrc        = safe_num(raw, "dlrc"),
    ser7        = safe_num(raw, "ser7"),
    # --- IC: Psychological ---
    depression  = safe_num(raw, "cesd10"),       # CES-D 10-item, 0-30
    dep_scale   = "CES-D 10",
    dep_max     = 30L,
    # --- IC: Sensory ---
    vision      = safe_num(raw, "eyesight_distance"),  # distance vision
    hearing     = safe_num(raw, "hear"),
    # --- IC: Locomotion ---
    grip_max    = {
      gs <- safe_num(raw, "gripsum")
      gm <- compute_grip_max(raw)
      if_else(is.na(gs), gm, gs)
    },
    gait_speed  = safe_num(raw, "wspeed"),
    # --- IC: Vitality ---
    bmi         = safe_num(raw, "bmi"),
    srh         = safe_num(raw, "srh"),          # named 'srh' in CHARLS CSV
    # --- Outcomes ---
    adl_total   = safe_num(raw, "adlab_c"),
    iadl_total  = safe_num(raw, "iadl"),
    iwstat      = safe_chr(raw, "iwstat"),
    death_year  = safe_num(raw, "radyear"),
    # --- Covariates ---
    chronic_n   = {
      cn <- safe_num(raw, "chronic_num")         # pre-computed in CHARLS CSV
      cc <- compute_chronic_count(raw)            # fallback
      if_else(is.na(cn), cc, cn)
    },
    # --- Environment ---
    loneliness  = safe_num(raw, "flonel"),        # single CES-D item (1-4), NOT UCLA
    social_participation = count_social_activities(
      raw, paste0("social", 1:11)                 # 11 social activity types
    ),
    # --- Meta ---
    country     = "CN"
  )
}

load_share <- function() {
  message("  Loading SHARE...")
  raw <- read_csv(file.path(DATA_DIR, "share.csv"), show_col_types = FALSE) %>%
    clean_headers()

  tibble(
    cohort      = "SHARE",
    id          = as.character(raw$mergeid),
    wave        = as.integer(raw$wave),
    age         = safe_num(raw, "agey"),
    female      = detect_female(raw$ragender),
    edu_years   = safe_num(raw, "raedyrs"),
    marital     = safe_chr(raw, "mstath"),
    interview_year = safe_num(raw, "iwy"),
    # --- IC: Cognition ---
    imrc        = safe_num(raw, "imrc"),
    dlrc        = safe_num(raw, "dlrc"),
    ser7        = safe_num(raw, "ser7"),
    # --- IC: Psychological ---
    depression  = safe_num(raw, "eurod"),         # EURO-D 12-item, 0-12
    dep_scale   = "EURO-D 12",
    dep_max     = 12L,
    # --- IC: Sensory ---
    vision      = safe_num(raw, "dsight"),        # distance vision
    hearing     = safe_num(raw, "hearing"),
    # --- IC: Locomotion ---
    grip_max    = compute_grip_max(raw),
    gait_speed  = safe_num(raw, "wspeed"),        # only W1-2; sparse
    # --- IC: Vitality ---
    bmi         = safe_num(raw, "bmi"),
    srh         = safe_num(raw, "shlt"),
    # --- Outcomes ---
    adl_total   = compute_adl6(raw),              # no pre-computed total
    iadl_total  = compute_iadl(raw, c("phonea", "medsa", "moneya", "shopa", "mealsa")),
    iwstat      = safe_chr(raw, "iwstat"),
    death_year  = safe_num(raw, "radyear"),
    # --- Covariates ---
    chronic_n   = compute_chronic_count(raw),
    # --- Environment ---
    loneliness  = NA_real_,                       # no UCLA scale
    social_participation = NA_real_,              # TODO: map SHARE ac035_ vars
    # --- Meta ---
    country     = safe_chr(raw, "country")
  )
}

load_mhas <- function() {
  message("  Loading MHAS...")
  raw <- read_csv(file.path(DATA_DIR, "mhas.csv"), show_col_types = FALSE) %>%
    clean_headers()

  tibble(
    cohort      = "MHAS",
    id          = as.character(raw$rahhidnp),
    wave        = as.integer(raw$wave),
    age         = safe_num(raw, "agey"),
    female      = detect_female(raw$ragender),
    edu_years   = safe_num(raw, "raedyrs"),
    marital     = safe_chr(raw, "mstath"),
    interview_year = safe_num(raw, "iwy"),
    # --- IC: Cognition ---
    imrc        = safe_num(raw, "imrc8"),         # 8-word list (not 10)
    dlrc        = safe_num(raw, "dlrc8"),
    ser7        = safe_num(raw, "ser7"),
    # --- IC: Psychological ---
    depression  = safe_num(raw, "cesd_m"),        # CES-D 9-item
    dep_scale   = "CES-D 9",
    dep_max     = 9L,
    # --- IC: Sensory ---
    vision      = safe_num(raw, "sight"),
    hearing     = safe_num(raw, "hearing"),
    # --- IC: Locomotion ---
    grip_max    = compute_grip_max(raw),
    gait_speed  = safe_num(raw, "wspeed"),        # 3m walk
    # --- IC: Vitality ---
    bmi         = safe_num(raw, "bmi"),
    srh         = safe_num(raw, "shlt"),
    # --- Outcomes ---
    adl_total   = safe_num(raw, "adltot6"),
    iadl_total  = safe_num(raw, "iadlfour"),      # 4-item IADL
    iwstat      = safe_chr(raw, "iwstat"),
    death_year  = safe_num(raw, "radyear"),
    # --- Covariates ---
    chronic_n   = compute_chronic_count(raw),
    # --- Environment ---
    loneliness  = NA_real_,                       # not available
    social_participation = NA_real_,              # TODO: map MHAS social vars
    # --- Meta ---
    country     = "MX"
  )
}

# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

main <- function() {
  cat("=" %>% strrep(70), "\n")
  cat("Study 17: Loading 5 core cohorts\n")
  cat("Data source:", DATA_DIR, "\n")
  cat("=" %>% strrep(70), "\n\n")

  # --- Step 0: Verify data directory exists ---
  if (!dir.exists(DATA_DIR)) {
    stop("Data directory not found: ", DATA_DIR,
         "\nPlease ensure the external drive is mounted.")
  }

  # Verify all CSV files exist
  missing <- COHORT_FILES[!file.exists(file.path(DATA_DIR, COHORT_FILES))]
  if (length(missing) > 0) {
    stop("Missing CSV files: ", paste(missing, collapse = ", "))
  }

  # --- Step 1: Load each cohort ---
  cat("Step 1: Loading individual cohorts...\n")

  cohorts <- list(
    load_hrs(),
    load_elsa(),
    load_charls(),
    load_share(),
    load_mhas()
  )

  # --- Step 2: Harmonize columns before binding ---
  # Ensure all tibbles have the same columns (some cohort-specific columns
  # like 'nursing_home' or 'country' may not be in all)
  cat("\nStep 2: Harmonizing columns...\n")

  all_cols <- unique(unlist(lapply(cohorts, names)))

  # Determine the correct type for each column from the first tibble that has it
  col_types <- setNames(rep("character", length(all_cols)), all_cols)
  for (col in all_cols) {
    for (df in cohorts) {
      if (col %in% names(df)) {
        col_types[[col]] <- class(df[[col]])[1]
        break
      }
    }
  }

  cohorts <- lapply(cohorts, function(df) {
    for (col in setdiff(all_cols, names(df))) {
      if (col_types[[col]] %in% c("numeric", "double")) {
        df[[col]] <- NA_real_
      } else if (col_types[[col]] == "integer") {
        df[[col]] <- NA_integer_
      } else {
        df[[col]] <- NA_character_
      }
    }
    df
  })

  # --- Step 3: Stack all cohorts ---
  cat("Step 3: Stacking cohorts...\n")
  pooled <- bind_rows(cohorts)

  cat("  Total rows (all ages, all waves): ", format(nrow(pooled), big.mark = ","), "\n")
  cat("  Unique persons: ", format(n_distinct(pooled$id), big.mark = ","), "\n")
  cat("  By cohort:\n")
  pooled %>%
    group_by(cohort) %>%
    summarise(
      n_rows = n(),
      n_persons = n_distinct(id),
      waves = paste(sort(unique(wave)), collapse = ","),
      .groups = "drop"
    ) %>%
    print()

  # --- Step 4: Filter eligible participants ---
  cat("\nStep 4: Applying eligibility filters...\n")

  n_before <- nrow(pooled)

  # 4a. Age ≥ 60
  pooled <- pooled %>% filter(!is.na(age) & age >= AGE_MIN)
  n_after_age <- nrow(pooled)
  cat("  After age ≥ ", AGE_MIN, ": ", format(n_after_age, big.mark = ","),
      " (removed ", format(n_before - n_after_age, big.mark = ","), ")\n")

  # 4b. Exclude nursing home residents (ELSA has nhmliv; others: use proxy)
  if ("nursing_home" %in% names(pooled)) {
    n_nh <- sum(pooled$nursing_home == 1, na.rm = TRUE)
    pooled <- pooled %>% filter(is.na(nursing_home) | nursing_home != 1)
    cat("  Excluded nursing home (ELSA): ", n_nh, "\n")
  }

  # 4c. Exclude baseline severe ADL disability (≥2 ADL at first wave)
  # This is applied later per-person after identifying baseline wave

  cat("  Final eligible rows: ", format(nrow(pooled), big.mark = ","), "\n")

  # --- Step 5: Standardize depression scores ---
  cat("\nStep 5: Creating standardized depression indicator...\n")

  # Binary depression flag using established cutoffs:
  #   CES-D 8 ≥ 3, CES-D 9 ≥ 5, CES-D 10 ≥ 10, EURO-D ≥ 4
  pooled <- pooled %>%
    mutate(
      depressed = case_when(
        dep_scale == "CES-D 8"  & depression >= 3  ~ 1L,
        dep_scale == "CES-D 8"  & depression < 3   ~ 0L,
        dep_scale == "CES-D 9"  & depression >= 5  ~ 1L,
        dep_scale == "CES-D 9"  & depression < 5   ~ 0L,
        dep_scale == "CES-D 10" & depression >= 10 ~ 1L,
        dep_scale == "CES-D 10" & depression < 10  ~ 0L,
        dep_scale == "EURO-D 12" & depression >= 4  ~ 1L,
        dep_scale == "EURO-D 12" & depression < 4   ~ 0L,
        TRUE ~ NA_integer_
      ),
      # Percentile-normalized depression (within-cohort-wave)
      dep_pctile = NA_real_  # placeholder: computed in 02_harmonization.R
    )

  # --- Step 6: Create derived variables ---
  cat("Step 6: Computing derived variables...\n")

  pooled <- pooled %>%
    mutate(
      # BMI U-shape: distance from optimal range (18.5-25)
      bmi_risk = case_when(
        bmi < 18.5 ~ 18.5 - bmi,
        bmi > 30   ~ bmi - 30,
        bmi >= 18.5 & bmi <= 30 ~ 0,
        TRUE ~ NA_real_
      ),
      # Vision/hearing: ensure higher = better (reverse if needed)
      # Original coding: 1=excellent ... 5=poor → reverse so 5=excellent
      vision_r = 6L - as.integer(vision),
      hearing_r = 6L - as.integer(hearing),
      # Age group
      age_group = cut(age, breaks = c(60, 65, 70, 75, 80, 85, Inf),
                      right = FALSE,
                      labels = c("60-64", "65-69", "70-74", "75-79", "80-84", "85+")),
      # Death indicator: primarily from death_year (propagated across waves),
      # supplemented by iwstat coding (0=alive, 1=dead in cleaned CSV;
      # or "5"/"5.died this wave" in Harmonized .dta coding)
      person_died = as.integer(
        !is.na(death_year) |
        iwstat %in% c("1", "5", "5.died this wave", "5.die")
      )
    )

  # --- Step 7: Variable availability summary ---
  cat("\nStep 7: Variable availability summary...\n")

  ic_env_vars <- c("imrc", "dlrc", "ser7", "depression", "vision", "hearing",
                    "grip_max", "gait_speed", "bmi", "srh", "adl_total",
                    "iadl_total", "chronic_n", "loneliness", "social_participation")

  var_summary <- pooled %>%
    group_by(cohort) %>%
    summarise(
      n_total = n(),
      n_persons = n_distinct(id),
      across(all_of(ic_env_vars), ~ round(mean(!is.na(.)) * 100, 1),
             .names = "pct_{.col}"),
      .groups = "drop"
    )

  print(var_summary)

  # --- Step 8: Save outputs ---
  cat("\nStep 8: Saving outputs...\n")

  saveRDS(pooled, file.path(OUTPUT_DIR, "analytic_pooled.rds"))
  cat("  Saved: ", file.path(OUTPUT_DIR, "analytic_pooled.rds"), "\n")

  write_csv(var_summary, file.path(OUTPUT_DIR, "cohort_summary.csv"))
  cat("  Saved: ", file.path(OUTPUT_DIR, "cohort_summary.csv"), "\n")

  # --- Final summary ---
  cat("\n", "=" %>% strrep(70), "\n")
  cat("DONE. Unified dataset:\n")
  cat("  Rows:       ", format(nrow(pooled), big.mark = ","), "\n")
  cat("  Persons:    ", format(n_distinct(pooled$id), big.mark = ","), "\n")
  cat("  Cohorts:    ", paste(unique(pooled$cohort), collapse = ", "), "\n")
  cat("  Columns:    ", ncol(pooled), "\n")
  cat("  Age range:  ", min(pooled$age, na.rm = TRUE), "-",
      max(pooled$age, na.rm = TRUE), "\n")
  cat("  Person-waves of deceased: ", sum(pooled$person_died == 1, na.rm = TRUE), "\n")
  cat("  Unique deceased persons:  ", n_distinct(pooled$id[pooled$person_died == 1]), "\n")
  cat("=" %>% strrep(70), "\n")

  invisible(pooled)
}

# ==============================================================================
# RUN
# ==============================================================================

if (interactive() || !exists(".main_called")) {
  .main_called <- TRUE
  pooled <- main()
}

# ==============================================================================
# NOTES & KNOWN LIMITATIONS
# ==============================================================================
#
# 1. DEPRESSION SCALE HETEROGENEITY:
#    - HRS/ELSA: CES-D 8-item (0-8)
#    - CHARLS: CES-D 10-item (0-30)
#    - SHARE: EURO-D 12-item (0-12) — different construct emphasis
#    - MHAS: CES-D 9-item (0-9)
#    Binary thresholds applied (CES-D 8≥3, CES-D 9≥5, CES-D 10≥10, EURO-D≥4).
#    Percentile normalization and IRT equating planned for 02_harmonization.R.
#    Key ref: Courtin et al. 2015, Int J Methods Psychiatr Res.
#
# 2. COGNITION WORD LIST LENGTH:
#    - HRS/ELSA/CHARLS/SHARE: 10-word list
#    - MHAS: 8-word list (imrc8, dlrc8)
#    Scores need scaling (×10/8) or percentile normalization.
#
# 3. GAIT SPEED LIMITATIONS:
#    - SHARE: only Waves 1-2 (~11% of respondents); 2.5m distance
#    - HRS: available from 2006 onward
#    - MHAS: 3m distance
#    - ELSA: available most waves
#    - CHARLS: available all waves
#    Distance differences affect speed comparability.
#
# 4. LONELINESS MEASUREMENT:
#    - HRS: UCLA-3 (lblonely3)
#    - ELSA: UCLA-3 (lnlys3, average of 3 items)
#    - CHARLS: single CES-D item only (flonel, 1-4); NOT UCLA
#    - SHARE: not available
#    - MHAS: not available
#    → Main analysis uses 4-variable Env index (without loneliness).
#
# 5. ENVIRONMENT VARIABLES (TODO for next version):
#    Social participation partially mapped (ELSA group1-8, CHARLS social1-11).
#    HRS, SHARE, MHAS social participation, social support, living arrangement,
#    and economic strain variables need further mapping.
#
# 6. EDUCATION (CHARLS):
#    CHARLS CSV lacks continuous education years (raedyrs).
#    Only categorical raeduc_c available. Needs conversion or imputation.
#
# 7. SELF-RATED HEALTH CODING:
#    Most cohorts use 1-5 (1=excellent, 5=poor).
#    CHARLS uses 'srh' (not 'shlt'). Verify coding direction.
#
# 8. ADL/IADL ITEM COUNTS:
#    - Standard: 6 ADL, 5 IADL
#    - MHAS: 6 ADL, 4 IADL (iadlfour)
#    - SHARE: computed from individual items (no pre-computed total)
#    Item-level harmonization may be needed.
# ==============================================================================
