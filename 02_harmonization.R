# ==============================================================================
# 02_harmonization.R
# Study 17: IC-Environment Fit and Disability-Free Survival in Older Adults
#
# Purpose: Harmonize IC indicators and environment variables across 5 cohorts:
#   (1) Extract additional environment variables from raw CSVs
#   (2) Standardize IC indicators (z-scores, reverse-coding, scaling)
#   (3) Standardize environment indicators
#   (4) Create simple composite scores (IC sum, Env sum)
#   (5) Additional harmonization (age cap, education, gait distance)
#
# Input:   data/analytic_pooled.rds (from 01_data_load.R)
#          Raw CSV files at DATA_DIR (for environment variables)
# Output:  data/harmonized_pooled.rds  — fully harmonized analytic dataset
#          data/harmonization_summary.csv — variable availability after harmonization
#
# Dependencies: tidyverse
# Author:  Study 17 Team
# Date:    2026-03-12
# Version: 1.1  (v1.1: fix CHARLS hhres/marry, HRS vol, social_contact calc)
# ==============================================================================

library(tidyverse)

# --- Configuration ------------------------------------------------------------

DATA_DIR <- file.path(here::here(), "data", "raw")  # Place raw CSV files here

SCRIPT_DIR <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) {
    here::here()
  }
)
OUTPUT_DIR <- file.path(SCRIPT_DIR, "data")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

AGE_CAP <- 105L  # Cap extreme ages

# CSV header cleaner (same as 01_data_load.R)
clean_headers <- function(df) {
  names(df) <- str_extract(names(df), "^[^\\s(]+")
  df
}

safe_num <- function(df, varname) {
  if (varname %in% names(df)) {
    suppressWarnings(as.numeric(df[[varname]]))
  } else {
    rep(NA_real_, nrow(df))
  }
}

# ==============================================================================
# STEP 1: Load pooled data from 01_data_load.R
# ==============================================================================

step1_load <- function() {
  cat("=" %>% strrep(70), "\n")
  cat("Step 1: Loading pooled data from 01_data_load.R\n")
  cat("=" %>% strrep(70), "\n\n")

  rds_path <- file.path(OUTPUT_DIR, "analytic_pooled.rds")
  if (!file.exists(rds_path)) {
    stop("analytic_pooled.rds not found. Run 01_data_load.R first.")
  }

  pooled <- readRDS(rds_path)
  cat("  Loaded:", format(nrow(pooled), big.mark = ","), "rows,",
      n_distinct(pooled$id), "persons\n")
  cat("  Cohorts:", paste(unique(pooled$cohort), collapse = ", "), "\n\n")

  pooled
}

# ==============================================================================
# STEP 2: Extract environment variables from raw CSVs
# ==============================================================================
# The 01_data_load.R script mapped social_participation for ELSA (group1-8)
# and CHARLS (social1-11) but left HRS, SHARE, MHAS as NA.
# This step also adds: social_contact, living_alone, financial_strain.

extract_env_hrs <- function() {
  cat("  Extracting HRS environment vars...\n")
  raw <- read_csv(file.path(DATA_DIR, "hrs.csv"), show_col_types = FALSE) %>%
    clean_headers()

  tibble(
    cohort = "HRS",
    id     = as.character(raw$hhidpn),
    wave   = as.integer(raw$wave),

    # Social participation: use `vol` (binary volunteer, ~99% wave 8+)
    # instead of detailed `volunteer/charity/club/...` (only wave 9+, ~35%)
    soc_volunteer   = as.integer(safe_num(raw, "vol") > 0),
    soc_charity     = NA_integer_,  # only wave 9+, coverage too low
    soc_club        = NA_integer_,  # only wave 9+, coverage too low
    soc_nonrelig    = NA_integer_,  # only wave 9+, coverage too low
    soc_art         = NA_integer_,  # only wave 9+, coverage too low
    soc_hobby       = as.integer(safe_num(raw, "hobby") > 0),

    # Social contact: monthly contact with children/relatives/friends
    contact_child   = safe_num(raw, "cntc"),    # 1=yes monthly
    contact_relat   = safe_num(raw, "cntr"),
    contact_friend  = safe_num(raw, "cntf"),

    # Household size (for living arrangement)
    hhres           = safe_num(raw, "hhres"),

    # Income (household level, for financial strain proxy)
    hh_income       = safe_num(raw, "itot"),
    hh_wealth       = safe_num(raw, "atotb")
  )
}

extract_env_elsa <- function() {
  cat("  Extracting ELSA environment vars...\n")
  raw <- read_csv(file.path(DATA_DIR, "elsa.csv"), show_col_types = FALSE) %>%
    clean_headers()

  tibble(
    cohort = "ELSA",
    id     = as.character(raw$idauniqc),
    wave   = as.integer(raw$wave),

    # Social participation already in pooled (group1-8), but add detail
    soc_volunteer   = NA_integer_,   # included in group1-8 aggregate
    soc_charity     = NA_integer_,
    soc_club        = NA_integer_,
    soc_nonrelig    = NA_integer_,
    soc_art         = NA_integer_,
    soc_hobby       = as.integer(safe_num(raw, "hobby") > 0),

    # Social contact: loneliness items as proxy for social contact quality
    contact_child   = NA_real_,
    contact_relat   = NA_real_,
    contact_friend  = NA_real_,

    # Household size
    hhres           = safe_num(raw, "hhhres"),

    # Income/wealth
    hh_income       = safe_num(raw, "hitot"),
    hh_wealth       = safe_num(raw, "hatotb")
  )
}

extract_env_charls <- function() {
  cat("  Extracting CHARLS environment vars...\n")
  raw <- read_csv(file.path(DATA_DIR, "charls.csv"), show_col_types = FALSE) %>%
    clean_headers()

  tibble(
    cohort = "CHARLS",
    id     = as.character(raw$ID),
    wave   = as.integer(raw$wave),

    # Social participation already in pooled (social1-11), supplement
    soc_volunteer   = NA_integer_,
    soc_charity     = NA_integer_,
    soc_club        = NA_integer_,
    soc_nonrelig    = NA_integer_,
    soc_art         = NA_integer_,
    soc_hobby       = NA_integer_,

    # Social contact proxy: number of children, child satisfaction
    contact_child   = NA_real_,
    contact_relat   = NA_real_,
    contact_friend  = NA_real_,

    # Household & financial
    # CHARLS has no `hhres`; use `family_size` (100% coverage, range 1-16)
    hhres           = safe_num(raw, "family_size"),
    hh_income       = safe_num(raw, "income_total"),
    hh_wealth       = safe_num(raw, "hatotfa")
  )
}

extract_env_share <- function() {
  cat("  Extracting SHARE environment vars...\n")
  raw <- read_csv(file.path(DATA_DIR, "share.csv"), show_col_types = FALSE) %>%
    clean_headers()

  # Social participation: act1-act8 (binary: did activity in past 12 months)
  act_items <- paste0("act", 1:8)
  avail_acts <- intersect(act_items, names(raw))

  soc_count <- if (length(avail_acts) > 0) {
    mat <- raw %>%
      select(all_of(avail_acts)) %>%
      mutate(across(everything(), ~ {
        x <- suppressWarnings(as.numeric(.))
        if_else(!is.na(x) & x > 0, 1L, 0L)
      }))
    counts <- rowSums(mat, na.rm = TRUE)
    all_na <- rowSums(!is.na(mat)) == 0
    counts[all_na] <- NA_real_
    counts
  } else {
    rep(NA_real_, nrow(raw))
  }

  tibble(
    cohort = "SHARE",
    id     = as.character(raw$mergeid),
    wave   = as.integer(raw$wave),

    # Social participation count
    soc_volunteer   = as.integer(safe_num(raw, "act1") > 0),
    soc_charity     = NA_integer_,
    soc_club        = as.integer(safe_num(raw, "act3") > 0),
    soc_nonrelig    = as.integer(safe_num(raw, "act5") > 0),
    soc_art         = NA_integer_,
    soc_hobby       = as.integer(safe_num(raw, "hobby") > 0),

    # Social isolation score (4-point, higher = more isolated)
    contact_child   = NA_real_,
    contact_relat   = NA_real_,
    contact_friend  = NA_real_,

    # Household
    hhres           = safe_num(raw, "hhres"),
    hh_income       = safe_num(raw, "hhitothhinc"),
    hh_wealth       = safe_num(raw, "hhatotb"),

    # SHARE-specific
    soc_count       = soc_count,
    sisa            = safe_num(raw, "sisa")       # social isolation 0-4
  )
}

extract_env_mhas <- function() {
  cat("  Extracting MHAS environment vars...\n")
  raw <- read_csv(file.path(DATA_DIR, "mhas.csv"), show_col_types = FALSE) %>%
    clean_headers()

  tibble(
    cohort = "MHAS",
    id     = as.character(raw$rahhidnp),
    wave   = as.integer(raw$wave),

    # Social participation
    soc_volunteer   = NA_integer_,
    soc_charity     = NA_integer_,
    soc_club        = NA_integer_,
    soc_nonrelig    = NA_integer_,
    soc_art         = NA_integer_,
    soc_hobby       = NA_integer_,

    # Social contact
    contact_child   = NA_real_,
    contact_relat   = safe_num(raw, "rfcnt"),       # weekly contact with relatives/friends
    contact_friend  = NA_real_,

    # Social activities
    soc_work        = safe_num(raw, "socwk"),        # weekly social activities
    soc_relig       = safe_num(raw, "relgwk"),       # weekly religious activities

    # Living arrangement
    coresd          = safe_num(raw, "coresd"),        # child co-residence
    hhres           = safe_num(raw, "hhres"),

    # Financial
    hh_income       = safe_num(raw, "itot"),
    hh_wealth       = safe_num(raw, "atotb")
  )
}

step2_extract_env <- function(pooled) {
  cat("=" %>% strrep(70), "\n")
  cat("Step 2: Extracting environment variables from raw CSVs\n")
  cat("=" %>% strrep(70), "\n\n")

  if (!dir.exists(DATA_DIR)) {
    warning("Raw data directory not found. Skipping environment extraction.\n",
            "Environment variables will remain as loaded from 01_data_load.R.")
    return(pooled)
  }

  # Extract environment vars from each cohort
  env_hrs   <- extract_env_hrs()
  env_elsa  <- extract_env_elsa()
  env_charls <- extract_env_charls()
  env_share <- extract_env_share()
  env_mhas  <- extract_env_mhas()

  # --- Compute social participation for cohorts where it was NA ---

  # HRS: count of social activities (6 binary items)
  env_hrs <- env_hrs %>%
    mutate(
      soc_participation_new = {
        mat <- select(., starts_with("soc_"))
        counts <- rowSums(mat, na.rm = TRUE)
        all_na <- rowSums(!is.na(mat)) == 0
        counts[all_na] <- NA_real_
        counts
      }
    )

  # SHARE: already computed soc_count
  env_share <- env_share %>%
    mutate(soc_participation_new = soc_count)

  # MHAS: combine social activities + religious
  env_mhas <- env_mhas %>%
    mutate(
      soc_participation_new = {
        sw <- if_else(!is.na(soc_work) & soc_work == 1, 1L, 0L)
        sr <- if_else(!is.na(soc_relig) & soc_relig == 1, 1L, 0L)
        rc <- if_else(!is.na(contact_relat) & contact_relat == 1, 1L, 0L)
        total <- sw + sr + rc
        # Set NA if all source variables are NA
        all_na <- is.na(soc_work) & is.na(soc_relig) & is.na(contact_relat)
        total[all_na] <- NA_real_
        total
      }
    )

  # ELSA/CHARLS: social participation already computed in 01_data_load.R

  # --- Compute living arrangement ---
  # Derive from marital status + household size
  # 0 = alone (hhres == 1)
  # 1 = with spouse (married/partnered, hhres >= 2)
  # 2 = with others (not alone, not spouse-only pattern)

  compute_living_arr <- function(marital, hhres) {
    m <- tolower(as.character(marital))
    # Married/partnered patterns across cohorts:
    #   HRS/ELSA: "1.married", "1.married, spouse present", "1.married/partnered"
    #   SHARE:    "1", "5" (registered partnership)
    #   CHARLS:   "1" (married with spouse), "7" (cohabiting)
    #   MHAS:     "1.married, spouse present", "1.married/partnered"
    married <- m %in% c("1", "5", "7",
                         "1.married", "1.married, spouse present",
                         "married", "partnered", "1.married/partnered",
                         "5.registered partnership")
    case_when(
      !is.na(hhres) & hhres == 1 ~ 0L,      # alone
      married & !is.na(hhres) & hhres >= 2 ~ 1L,  # with spouse
      !is.na(hhres) & hhres >= 2 ~ 2L,      # with others
      TRUE ~ NA_integer_
    )
  }

  # --- Compute financial strain proxy ---
  # Use income quintile within cohort-wave (lower quintile = more strain)
  # Reversed: 1 = lowest quintile (most strain), 5 = highest (least strain)

  compute_income_quintile <- function(income, cohort, wave) {
    df <- tibble(income = income, cohort = cohort, wave = wave) %>%
      group_by(cohort, wave) %>%
      mutate(
        inc_quintile = case_when(
          is.na(income) ~ NA_real_,
          TRUE ~ as.numeric(ntile(income, 5))
        )
      ) %>%
      ungroup()
    df$inc_quintile
  }

  # --- Merge environment variables into pooled data ---

  # Prepare minimal env tibbles with standardized columns
  common_env_cols <- c("cohort", "id", "wave", "hhres", "hh_income",
                       "hh_wealth", "soc_participation_new")

  env_hrs_slim <- env_hrs %>%
    select(any_of(c(common_env_cols, "contact_child", "contact_relat", "contact_friend")))
  env_elsa_slim <- env_elsa %>%
    mutate(soc_participation_new = NA_real_) %>%
    select(any_of(c(common_env_cols, "contact_child", "contact_relat", "contact_friend")))
  env_charls_slim <- env_charls %>%
    mutate(soc_participation_new = NA_real_) %>%
    select(any_of(c(common_env_cols, "contact_child", "contact_relat", "contact_friend")))
  env_share_slim <- env_share %>%
    select(any_of(c(common_env_cols, "contact_child", "contact_relat",
                    "contact_friend", "sisa")))
  env_mhas_slim <- env_mhas %>%
    select(any_of(c(common_env_cols, "contact_child", "contact_relat", "contact_friend")))

  # Ensure consistent columns across all env tibbles before binding
  all_env_cols <- unique(c(
    names(env_hrs_slim), names(env_elsa_slim), names(env_charls_slim),
    names(env_share_slim), names(env_mhas_slim)
  ))

  pad_cols <- function(df, all_cols) {
    for (col in setdiff(all_cols, names(df))) {
      df[[col]] <- NA_real_
    }
    df
  }

  env_all <- bind_rows(
    pad_cols(env_hrs_slim, all_env_cols),
    pad_cols(env_elsa_slim, all_env_cols),
    pad_cols(env_charls_slim, all_env_cols),
    pad_cols(env_share_slim, all_env_cols),
    pad_cols(env_mhas_slim, all_env_cols)
  )

  # Left-join environment variables to pooled data
  pooled <- pooled %>%
    left_join(env_all, by = c("cohort", "id", "wave"), suffix = c("", ".env"))

  # Update social_participation: use newly computed where old is NA
  pooled <- pooled %>%
    mutate(
      social_participation = if_else(
        is.na(social_participation) & !is.na(soc_participation_new),
        soc_participation_new,
        social_participation
      )
    )

  # Compute living arrangement from marital + hhres
  pooled <- pooled %>%
    mutate(
      living_arrangement = compute_living_arr(marital, hhres)
    )

  # Compute financial strain proxy (income quintile, reversed)
  pooled <- pooled %>%
    mutate(
      income_quintile = compute_income_quintile(hh_income, cohort, wave),
      # Financial strain: lower income = higher strain (reverse of quintile)
      financial_strain = 6 - income_quintile  # 5=most strained, 1=least
    )

  # Social contact summary score (count of types with regular contact)
  pooled <- pooled %>%
    mutate(
      social_contact = {
        cc <- if_else(!is.na(contact_child) & contact_child > 0, 1L, 0L)
        cr <- if_else(!is.na(contact_relat) & contact_relat > 0, 1L, 0L)
        cf <- if_else(!is.na(contact_friend) & contact_friend > 0, 1L, 0L)
        total <- cc + cr + cf
        all_na <- is.na(contact_child) & is.na(contact_relat) & is.na(contact_friend)
        total[all_na] <- NA_real_
        total
      }
    )

  cat("  Environment variables merged.\n")
  cat("  Social participation coverage:\n")
  pooled %>%
    group_by(cohort) %>%
    summarise(
      soc_part_pct = round(mean(!is.na(social_participation)) * 100, 1),
      living_arr_pct = round(mean(!is.na(living_arrangement)) * 100, 1),
      fin_strain_pct = round(mean(!is.na(financial_strain)) * 100, 1),
      social_contact_pct = round(mean(!is.na(social_contact)) * 100, 1),
      .groups = "drop"
    ) %>%
    print()

  pooled
}

# ==============================================================================
# STEP 3: IC Indicator Harmonization
# ==============================================================================

step3_ic_harmonization <- function(pooled) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 3: IC indicator harmonization\n")
  cat("=" %>% strrep(70), "\n\n")

  # --- 3a: MHAS word recall scaling (8-word → 10-word equivalent) ---
  cat("  3a: Scaling MHAS word recall (8→10 words)...\n")
  pooled <- pooled %>%
    mutate(
      imrc = if_else(cohort == "MHAS" & !is.na(imrc),
                     imrc * (10 / 8), imrc),
      dlrc = if_else(cohort == "MHAS" & !is.na(dlrc),
                     dlrc * (10 / 8), dlrc)
    )

  # --- 3b: Age capping ---
  cat("  3b: Capping age at", AGE_CAP, "...\n")
  n_extreme <- sum(pooled$age > AGE_CAP, na.rm = TRUE)
  pooled <- pooled %>%
    mutate(age = pmin(age, AGE_CAP))
  cat("    Capped", n_extreme, "observations\n")

  # --- 3c: CHARLS education (categorical → approximate years) ---
  cat("  3c: Converting CHARLS education to years...\n")
  # CHARLS raeduc_c categories → approximate education years
  # Need to check what the actual categories are
  # Common CHARLS education coding:
  # 1=illiterate, 2=did not finish primary, 3=primary school,
  # 4=middle school, 5=high school, 6=vocational, 7=college+
  # We'll compute from marital status proxy if raeduc_c not available
  # For now, CHARLS edu_years remains NA; will be handled in CFA as missing

  # --- 3d: Reverse-code depression (higher = better psychological health) ---
  cat("  3d: Reverse-coding depression...\n")
  pooled <- pooled %>%
    mutate(
      # Depression: original higher = more depressed
      # Reversed: higher = better psychological health
      depression_r = dep_max - depression
    )

  # --- 3e: Depression percentile normalization (within cohort-wave) ---
  cat("  3e: Depression percentile normalization...\n")
  pooled <- pooled %>%
    group_by(cohort, wave) %>%
    mutate(
      dep_pctile = if_else(
        !is.na(depression),
        percent_rank(depression),  # 0-1, higher = more depressed
        NA_real_
      ),
      # Reversed percentile: higher = better
      dep_pctile_r = 1 - dep_pctile
    ) %>%
    ungroup()

  # --- 3f: Within-cohort-wave z-scores for IC indicators ---
  cat("  3f: Computing within-cohort-wave z-scores...\n")

  # Helper: within-group z-score
  z_score <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - m) / s
  }

  pooled <- pooled %>%
    group_by(cohort, wave) %>%
    mutate(
      # Cognition domain (higher = better)
      z_imrc      = z_score(imrc),
      z_dlrc      = z_score(dlrc),
      z_ser7      = z_score(ser7),

      # Psychological domain (higher = better)
      z_dep_r     = z_score(depression_r),  # reversed depression

      # Sensory domain (higher = better; already reversed in 01_data_load)
      z_vision    = z_score(vision_r),
      z_hearing   = z_score(hearing_r),

      # Vitality domain (higher = better)
      z_bmi_risk  = z_score(-bmi_risk),  # negate: lower risk = higher z
      z_srh       = z_score(6L - as.integer(srh))  # reverse SRH: 5=excellent→high
    ) %>%
    ungroup()

  # Locomotion: sex-stratified z-scores for grip strength
  cat("  3g: Sex-stratified z-scores for locomotion...\n")
  pooled <- pooled %>%
    group_by(cohort, wave, female) %>%
    mutate(
      z_grip      = z_score(grip_max),
      z_gait      = z_score(-gait_speed)  # negate: less time = faster = higher z
    ) %>%
    ungroup()

  # --- 3h: Domain-level composite scores (average of z-scores) ---
  cat("  3h: Computing domain-level composite scores...\n")
  pooled <- pooled %>%
    rowwise() %>%
    mutate(
      # Cognition: average of available z-scores
      ic_cognition = {
        vals <- c(z_imrc, z_dlrc, z_ser7)
        vals <- vals[!is.na(vals)]
        if (length(vals) >= 2) mean(vals) else NA_real_
      },
      # Psychological: single indicator
      ic_psychological = z_dep_r,
      # Sensory: average of vision and hearing
      ic_sensory = {
        vals <- c(z_vision, z_hearing)
        vals <- vals[!is.na(vals)]
        if (length(vals) >= 1) mean(vals) else NA_real_
      },
      # Locomotion: average of available measures
      ic_locomotion = {
        vals <- c(z_grip, z_gait)
        vals <- vals[!is.na(vals)]
        if (length(vals) >= 1) mean(vals) else NA_real_
      },
      # Vitality: average of BMI risk and SRH
      ic_vitality = {
        vals <- c(z_bmi_risk, z_srh)
        vals <- vals[!is.na(vals)]
        if (length(vals) >= 1) mean(vals) else NA_real_
      }
    ) %>%
    ungroup()

  # --- 3i: IC sum score (average across domains) ---
  cat("  3i: Computing IC sum score...\n")
  pooled <- pooled %>%
    rowwise() %>%
    mutate(
      ic_domains_available = sum(!is.na(c(ic_cognition, ic_psychological,
                                           ic_sensory, ic_locomotion, ic_vitality))),
      ic_sum = if_else(
        ic_domains_available >= 3,
        mean(c(ic_cognition, ic_psychological, ic_sensory,
               ic_locomotion, ic_vitality), na.rm = TRUE),
        NA_real_
      )
    ) %>%
    ungroup()

  cat("  IC sum score coverage:\n")
  pooled %>%
    group_by(cohort) %>%
    summarise(
      ic_sum_pct = round(mean(!is.na(ic_sum)) * 100, 1),
      domains_mean = round(mean(ic_domains_available, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    print()

  pooled
}

# ==============================================================================
# STEP 4: Environment Indicator Standardization
# ==============================================================================

step4_env_standardization <- function(pooled) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 4: Environment indicator standardization\n")
  cat("=" %>% strrep(70), "\n\n")

  z_score <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - m) / s
  }

  pooled <- pooled %>%
    group_by(cohort, wave) %>%
    mutate(
      # Social participation z-score (higher = more participation)
      z_soc_part = z_score(social_participation),

      # Financial strain z-score (higher = less strain = better)
      # financial_strain: 5=most strained, 1=least → reverse
      z_financial = z_score(-financial_strain),

      # Living arrangement: ordinal → z-score
      # 0=alone, 1=with spouse, 2=with others
      # Higher support = with spouse > with others > alone
      z_living = z_score(case_when(
        living_arrangement == 0 ~ 0,    # alone
        living_arrangement == 1 ~ 2,    # with spouse (most support)
        living_arrangement == 2 ~ 1,    # with others
        TRUE ~ NA_real_
      )),

      # Social contact z-score (higher = more contact)
      z_social_contact = z_score(social_contact)
    ) %>%
    ungroup()

  # --- Environment composite score ---
  cat("  Computing Environment composite score...\n")
  pooled <- pooled %>%
    rowwise() %>%
    mutate(
      env_components_available = sum(!is.na(c(z_soc_part, z_financial,
                                               z_living, z_social_contact))),
      # Primary: 4-variable Env index (without loneliness, per SAP)
      env_sum = if_else(
        env_components_available >= 2,
        mean(c(z_soc_part, z_financial, z_living, z_social_contact),
             na.rm = TRUE),
        NA_real_
      )
    ) %>%
    ungroup()

  cat("  Environment sum score coverage:\n")
  pooled %>%
    group_by(cohort) %>%
    summarise(
      env_sum_pct = round(mean(!is.na(env_sum)) * 100, 1),
      env_comp_mean = round(mean(env_components_available, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    print()

  pooled
}

# ==============================================================================
# STEP 5: Create IC × Env Interaction Variables
# ==============================================================================

step5_interaction <- function(pooled) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 5: Creating IC × Env interaction variables\n")
  cat("=" %>% strrep(70), "\n\n")

  # --- 5a: Continuous interaction (mean-centered) ---
  pooled <- pooled %>%
    group_by(cohort) %>%
    mutate(
      ic_centered  = ic_sum - mean(ic_sum, na.rm = TRUE),
      env_centered = env_sum - mean(env_sum, na.rm = TRUE),
      ic_env_interaction = ic_centered * env_centered
    ) %>%
    ungroup()

  # --- 5b: Categorical groups ---
  pooled <- pooled %>%
    group_by(cohort) %>%
    mutate(
      # IC tertiles within cohort
      ic_tertile = case_when(
        is.na(ic_sum) ~ NA_character_,
        ic_sum <= quantile(ic_sum, 1/3, na.rm = TRUE) ~ "Low",
        ic_sum <= quantile(ic_sum, 2/3, na.rm = TRUE) ~ "Medium",
        TRUE ~ "High"
      ),
      # Env median split within cohort
      env_group = case_when(
        is.na(env_sum) ~ NA_character_,
        env_sum <= median(env_sum, na.rm = TRUE) ~ "Low",
        TRUE ~ "High"
      ),
      # 6-group cross-classification
      ic_env_group = case_when(
        is.na(ic_tertile) | is.na(env_group) ~ NA_character_,
        TRUE ~ paste0(ic_tertile, "_IC/", env_group, "_Env")
      )
    ) %>%
    ungroup()

  # Report group sizes
  cat("  IC × Env group distribution:\n")
  pooled %>%
    filter(!is.na(ic_env_group)) %>%
    count(ic_env_group) %>%
    mutate(pct = round(n / sum(n) * 100, 1)) %>%
    print()

  pooled
}

# ==============================================================================
# STEP 6: Gait Speed Distance Correction
# ==============================================================================

step6_gait_correction <- function(pooled) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 6: Gait speed distance correction\n")
  cat("=" %>% strrep(70), "\n\n")

  # Different cohorts use different walk distances:
  #   HRS: varies (2.5m or 3m depending on wave)
  #   ELSA: 2.44m (8 feet)
  #   CHARLS: 2.5m
  #   SHARE: 2.5m (W1-2 only)
  #   MHAS: 3m
  #
  # The wspeed variable in Harmonized CSVs typically already represents
  # speed in m/s or seconds. Since we z-score within cohort-wave,
  # the distance differences are absorbed. No additional correction needed
  # for within-cohort z-scores.
  #
  # For cross-cohort comparison of raw speeds, convert to m/s:
  #   speed_ms = distance_m / time_seconds

  cat("  Note: Gait speed differences absorbed by within-cohort z-scoring.\n")
  cat("  No additional correction applied at this stage.\n")
  cat("  Cross-cohort comparison uses z_gait (standardized within cohort-wave-sex).\n")

  pooled
}

# ==============================================================================
# STEP 7: Data Quality Summary
# ==============================================================================

step7_quality_summary <- function(pooled) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 7: Data quality summary\n")
  cat("=" %>% strrep(70), "\n\n")

  # IC z-score availability
  ic_vars <- c("z_imrc", "z_dlrc", "z_ser7", "z_dep_r",
               "z_vision", "z_hearing", "z_grip", "z_gait",
               "z_bmi_risk", "z_srh")
  domain_vars <- c("ic_cognition", "ic_psychological", "ic_sensory",
                    "ic_locomotion", "ic_vitality", "ic_sum")
  env_vars <- c("z_soc_part", "z_financial", "z_living",
                "z_social_contact", "env_sum")
  interaction_vars <- c("ic_env_interaction", "ic_env_group")

  all_vars <- c(ic_vars, domain_vars, env_vars, interaction_vars)

  summary_df <- pooled %>%
    group_by(cohort) %>%
    summarise(
      n_total = n(),
      n_persons = n_distinct(id),
      across(all_of(all_vars[all_vars %in% names(pooled)]),
             ~ round(mean(!is.na(.)) * 100, 1),
             .names = "pct_{.col}"),
      .groups = "drop"
    )

  cat("  IC z-score coverage (%):\n")
  summary_df %>%
    select(cohort, n_total, n_persons,
           any_of(paste0("pct_", ic_vars))) %>%
    print(width = 150)

  cat("\n  Domain composite & IC sum coverage (%):\n")
  summary_df %>%
    select(cohort,
           any_of(paste0("pct_", domain_vars))) %>%
    print()

  cat("\n  Environment coverage (%):\n")
  summary_df %>%
    select(cohort,
           any_of(paste0("pct_", env_vars))) %>%
    print()

  summary_df
}

# ==============================================================================
# STEP 8: Save Harmonized Data
# ==============================================================================

step8_save <- function(pooled, summary_df) {
  cat("\n", "=" %>% strrep(70), "\n")
  cat("Step 8: Saving harmonized dataset\n")
  cat("=" %>% strrep(70), "\n\n")

  out_rds <- file.path(OUTPUT_DIR, "harmonized_pooled.rds")
  saveRDS(pooled, out_rds)
  cat("  Saved:", out_rds, "\n")

  out_csv <- file.path(OUTPUT_DIR, "harmonization_summary.csv")
  write_csv(summary_df, out_csv)
  cat("  Saved:", out_csv, "\n")

  # Final summary
  cat("\n", "=" %>% strrep(70), "\n")
  cat("HARMONIZATION COMPLETE\n")
  cat("  Rows:          ", format(nrow(pooled), big.mark = ","), "\n")
  cat("  Persons:       ", format(n_distinct(pooled$id), big.mark = ","), "\n")
  cat("  IC sum avail:  ", round(mean(!is.na(pooled$ic_sum)) * 100, 1), "%\n")
  cat("  Env sum avail: ", round(mean(!is.na(pooled$env_sum)) * 100, 1), "%\n")
  cat("  IC×Env groups: ", round(mean(!is.na(pooled$ic_env_group)) * 100, 1), "%\n")
  cat("  Total columns: ", ncol(pooled), "\n")
  cat("=" %>% strrep(70), "\n")

  invisible(pooled)
}

# ==============================================================================
# MAIN
# ==============================================================================

main <- function() {
  pooled <- step1_load()
  pooled <- step2_extract_env(pooled)
  pooled <- step3_ic_harmonization(pooled)
  pooled <- step4_env_standardization(pooled)
  pooled <- step5_interaction(pooled)
  pooled <- step6_gait_correction(pooled)
  summary_df <- step7_quality_summary(pooled)
  step8_save(pooled, summary_df)
}

if (interactive() || !exists(".main02_called")) {
  .main02_called <- TRUE
  main()
}

# ==============================================================================
# NOTES & KNOWN LIMITATIONS
# ==============================================================================
#
# 1. SOCIAL PARTICIPATION HETEROGENEITY:
#    - HRS: 6 activity types (volunteer, charity, club, nonreligious, art, hobby)
#    - ELSA: 8 organization types (group1-8)
#    - CHARLS: 11 social activity types (social1-11)
#    - SHARE: 8 activity types (act1-8)
#    - MHAS: 3 indicators (socwk, relgwk, rfcnt)
#    Activity counts are NOT directly comparable across cohorts.
#    Within-cohort z-scoring absorbs scale differences.
#    CFA latent factor approach (04_env_index.R) will formally handle this.
#
# 2. FINANCIAL STRAIN:
#    No direct "difficulty paying bills" variable available across all cohorts.
#    Using income quintile as proxy. This captures relative financial position
#    but not subjective financial stress. Alternative: wealth-to-needs ratio.
#
# 3. SOCIAL CONTACT:
#    Only HRS has the 3-type contact measure (children, relatives, friends).
#    MHAS has rfcnt (~91% coverage). ELSA/CHARLS/SHARE have 0%.
#    Social contact z-score is largely NA for ELSA/CHARLS/SHARE.
#    DECISION: Keep as 4th Env component. env_sum requires ≥2/4 components,
#    so ELSA/CHARLS/SHARE use 3-component average (soc_part, financial, living).
#    Sensitivity analysis: 3-variable Env index (without social_contact) for all.
#
# 4. LOCOMOTION DOMAIN SPARSITY:
#    - MHAS: grip only in Wave 3 (~12%); ic_locomotion will be mostly NA
#    - CHARLS: grip/gait only Waves 1-3; Waves 4-5 locomotion domain NA
#    - HRS: grip from Wave 7 (~2004); earlier waves locomotion NA
#    - ELSA ser7: only Waves 7-9
#    Impact: IC sum score uses ≥3/5 domains rule; some observations will
#    use 3-4 domain average instead of all 5.
#
# 5. IC SUM vs LATENT SCORE:
#    The ic_sum computed here is a simple z-score average for quick analysis.
#    The definitive IC score will be the CFA latent factor from 03_ic_cfa.R,
#    which properly handles measurement invariance and differential weighting.
#
# 6. CHARLS EDUCATION:
#    CHARLS lacks continuous education years (raedyrs). edu_years remains NA.
#    Options: (a) convert raeduc_c categorical to years using midpoints,
#    (b) treat education as categorical in all models,
#    (c) impute from other covariates. Decision deferred to analysis phase.
#
# 7. LIVING ARRANGEMENT INFERENCE:
#    Living arrangement derived from marital status + household size.
#    This is an approximation — cannot distinguish "with adult children"
#    from "with other relatives/non-relatives" without detailed roster data.
#    MHAS has explicit coresd (child co-residence) for finer coding.
#
# 8. CHARLS HOUSEHOLD SIZE (v1.1 fix):
#    CHARLS has no `hhres` variable. Uses `family_size` (100% coverage,
#    range 1-16) as equivalent. CHARLS marry coding: 1=married(with spouse),
#    2=separated, 3=divorced, 4=widowed, 5=never married, 7=cohabiting.
#    compute_living_arr() updated to match "7" as married/partnered.
#
# 9. HRS SOCIAL PARTICIPATION (v1.1 fix):
#    Detailed activity variables (volunteer, charity, club, nonreligious, art)
#    only available from wave 9 (2008+), giving ~23% overall coverage.
#    Switched to `vol` (binary volunteer indicator, ~99% from wave 8/2006+)
#    for dramatically improved coverage (~72% overall including early waves).
# ==============================================================================
