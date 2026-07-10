# ==============================================================================
# 17_klosa_replication.R  (v2.0 — reduced-indicator qualitative replication)
# Study 17: IC-Environment Fit and Disability-Free Survival
#
# IMPORTANT — DATA LIMITATION ACKNOWLEDGED
# ------------------------------------------------------------------------------
# Harmonized KLoSA E.2 (2006-2020) does NOT contain any cognitive indicators
# (no immediate recall, delayed recall, serial 7s, or MMSE-equivalent scores).
# It also lacks a non-religious "count of social activities" variable and a
# direct financial-strain item. This script therefore performs a REDUCED-
# INDICATOR QUALITATIVE REPLICATION rather than a methodologically identical
# replication of the main 5-cohort analysis.
#
# Reduced indicator set used here (documented in revision_notes/A2_klosa_limitation.md):
#   IC (5 indicators, no cognition):
#     - depression  (r{W}depresl  or  r{W}cesd10a/b)
#     - sight       (r{W}sighta,  1-5, reversed)
#     - hearing     (r{W}hearinga, 1-5, reversed)
#     - grip_max    (r{W}gripsum, kg)
#     - srh         (r{W}shlt,    1-5, reversed)
#   Env (3 indicators):
#     - financial   log(r{W}itot + 1),  z-scored within wave
#     - living      h{W}hhres,          z-scored within wave
#     - social_part r{W}socwk,          z-scored within wave
# CFA model: 2-factor (sensory + general_health). Single-group on KLoSA.
#
# Framing for manuscript:
#   "Qualitative external replication in KLoSA (Korea, n≈X,XXX) using the
#    reduced indicator set available in the Harmonized KLoSA E.2 release;
#    cognition was not replicated owing to absence of cognitive measures
#    in KLoSA E.2."
#
# Output:  results/klosa_cfa_fit.csv       — 2-factor Model B' fit on KLoSA
#          results/klosa_replication_m4.csv — KLoSA M4 Cox estimates
#          results/klosa_indicator_cov.csv  — per-indicator availability
#          data/klosa_analytic.rds           — KLoSA analytic file
#
# Dependencies: tidyverse, haven, lavaan, survival
# Author:  Study 17 Team
# Date:    2026-04-22
# Version: 2.0
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(haven)
  library(lavaan)
  library(survival)
})

SCRIPT_DIR <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) normalizePath(Sys.getenv("STUDY17_PROJECT_DIR", unset = getwd()), mustWork = FALSE)
)
DATA_DIR    <- file.path(SCRIPT_DIR, "data")
RESULTS_DIR <- file.path(SCRIPT_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

KLOSA_DTA <- Sys.getenv("STUDY17_KLOSA_DTA", unset = file.path(getwd(), "raw_data", "H_KLoSA_e2.dta"))

# ============================================================================
# Helpers
# ============================================================================

safe_col <- function(raw, col) {
  if (col %in% names(raw)) as.numeric(raw[[col]]) else rep(NA_real_, nrow(raw))
}
safe_col_chr <- function(raw, col) {
  if (col %in% names(raw)) as.character(raw[[col]]) else rep(NA_character_, nrow(raw))
}
# CESD: W1-4 used r{W}cesd10a, W5-8 used r{W}cesd10b
pick_cesd <- function(raw, w) {
  a <- safe_col(raw, paste0("r", w, "cesd10a"))
  b <- safe_col(raw, paste0("r", w, "cesd10b"))
  depresl <- safe_col(raw, paste0("r", w, "depresl"))
  # prefer a, then b, then depresl
  out <- rep(NA_real_, length(a))
  out <- ifelse(is.na(out) & !is.na(a), a, out)
  out <- ifelse(is.na(out) & !is.na(b), b, out)
  out <- ifelse(is.na(out) & !is.na(depresl), depresl, out)
  out
}

z_within <- function(x, grp) {
  ave(x, grp, FUN = function(v) {
    m <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(NA_real_, length(v)) else (v - m) / s
  })
}

# ============================================================================
# STEP 1: Load & reshape KLoSA wide → long (waves 1-8, 2006-2020)
# ============================================================================

step1_load <- function() {
  cat(strrep("=", 70), "\n"); cat("Step 1: Loading Harmonized KLoSA E.2\n"); cat(strrep("=", 70), "\n\n")

  if (!file.exists(KLOSA_DTA)) stop("KLoSA .dta not found at:\n  ", KLOSA_DTA)
  raw <- read_dta(KLOSA_DTA)
  cat(sprintf("  Rows %s, Cols %s\n",
              format(nrow(raw), big.mark=","),
              format(ncol(raw), big.mark=",")))

  # Invariants
  pid    <- if ("pid" %in% names(raw)) as.character(raw$pid) else as.character(raw$hhidpn)
  female <- as.integer(raw[["ragender"]] == 2)  # 1=M, 2=F
  edu_k  <- safe_col(raw, "raeduc_k")           # Korea-specific edu code (1-5)
  dy     <- safe_col(raw, "radyear")

  wave_vec <- 1:8
  panel <- lapply(wave_vec, function(w) {
    # Chronic count
    chr_cols <- c("hibpe","diabe","cancre","lunge","hearte","stroke","arthre")
    chr_mat  <- vapply(chr_cols, function(x) safe_col(raw, paste0("r", w, x)),
                       numeric(nrow(raw)))
    chr_any  <- rowSums(!is.na(chr_mat)) > 0
    chr_n    <- ifelse(chr_any, rowSums(chr_mat >= 1, na.rm = TRUE), NA_real_)

    # Financial proxy: log(total income + 1)
    itot_w <- safe_col(raw, paste0("r", w, "itot"))
    log_itot <- ifelse(!is.na(itot_w), log1p(pmax(itot_w, 0)), NA_real_)

    tibble(
      cohort      = "KLoSA",
      id          = pid,
      wave        = w,
      interview_year = safe_col(raw, paste0("r", w, "iwy")),
      age         = safe_col(raw, paste0("r", w, "agey")),
      female      = female,
      edu_k       = edu_k,
      marital     = safe_col_chr(raw, paste0("r", w, "mstath")),
      iwstat      = safe_col_chr(raw, paste0("r", w, "iwstat")),
      death_year  = dy,
      # --- IC (5 indicators, NO cognition) ---
      depression  = pick_cesd(raw, w),                  # 0-10 (higher=worse)
      sight       = safe_col(raw, paste0("r", w, "sighta")),      # 1-5
      hearing     = safe_col(raw, paste0("r", w, "hearinga")),    # 1-5
      grip_max    = safe_col(raw, paste0("r", w, "gripsum")),     # kg
      srh         = safe_col(raw, paste0("r", w, "shlt")),        # 1-5
      # --- Env (3 indicators) ---
      financial_raw = log_itot,                                   # log(income+1)
      living      = safe_col(raw, paste0("h", w, "hhres")),       # # residents
      social_part = safe_col(raw, paste0("r", w, "socwk")),       # community attend freq
      # --- Outcomes ---
      adl_total   = safe_col(raw, paste0("r", w, "adltotb_k")),   # Korean ADL total
      iadl_total  = safe_col(raw, paste0("r", w, "iadltotb_k")),  # IADL total
      chronic_n   = chr_n
    )
  })

  long <- bind_rows(panel) %>% filter(!is.na(age))
  cat(sprintf("  Long panel: %s rows, %s persons\n",
              format(nrow(long), big.mark=","),
              format(n_distinct(long$id), big.mark=",")))

  # Age eligibility
  long <- long %>% filter(age >= 60)
  cat(sprintf("  After age≥60: %s rows, %s persons\n",
              format(nrow(long), big.mark=","),
              format(n_distinct(long$id), big.mark=",")))

  # Indicator availability report
  av <- long %>% summarise(
    age = mean(!is.na(age)), female = mean(!is.na(female)),
    depression = mean(!is.na(depression)), sight = mean(!is.na(sight)),
    hearing = mean(!is.na(hearing)), grip = mean(!is.na(grip_max)),
    srh = mean(!is.na(srh)),
    financial = mean(!is.na(financial_raw)),
    living = mean(!is.na(living)), social_part = mean(!is.na(social_part)),
    adl = mean(!is.na(adl_total)), chronic = mean(!is.na(chronic_n))
  )
  cat("\n  Indicator availability (% of person-waves):\n")
  av_print <- round(as.numeric(av) * 100, 1); names(av_print) <- names(av)
  print(av_print)
  write.csv(
    tibble(indicator = names(av_print), pct_available = as.numeric(av_print)),
    file.path(RESULTS_DIR, "klosa_indicator_cov.csv"),
    row.names = FALSE
  )

  long
}

# ============================================================================
# STEP 2: z-scores + 2-factor CFA + IC/Env composites
# ============================================================================

step2_composites <- function(long) {
  cat("\n", strrep("=", 70), "\n"); cat("Step 2: z-scores and composites\n"); cat(strrep("=", 70), "\n\n")

  # Orient higher=better
  d <- long %>% mutate(
    .dep_h    = -depression,          # higher dep → worse
    .sight_h  = 6 - sight,            # 1=excellent → reverse
    .hearing_h= 6 - hearing,
    .grip_h   = grip_max,
    .srh_h    = 6 - srh,
    .fin_h    = financial_raw,        # log income: higher=better
    .liv_h    = living,               # more residents
    .soc_h    = social_part           # higher freq=better
  ) %>% mutate(
    z_dep     = z_within(.dep_h,      wave),
    z_sight   = z_within(.sight_h,    wave),
    z_hearing = z_within(.hearing_h,  wave),
    z_grip    = z_within(.grip_h,     wave),
    z_srh     = z_within(.srh_h,      wave),
    z_fin     = z_within(.fin_h,      wave),
    z_liv     = z_within(.liv_h,      wave),
    z_soc     = z_within(.soc_h,      wave)
  )

  # --- 2-factor CFA (Model B', reduced) ---
  cfa_model <- '
    sensory =~ z_sight + z_hearing
    general =~ z_dep + z_grip + z_srh
  '
  cat("  Fitting single-group 2-factor CFA Model B' on KLoSA...\n")
  fit <- tryCatch(
    cfa(cfa_model, data = d, missing = "fiml", estimator = "MLR"),
    error = function(e) { message("    CFA failed: ", conditionMessage(e)); NULL }
  )

  if (!is.null(fit)) {
    fm <- fitMeasures(fit, c("cfi","tli","rmsea","srmr","chisq","df","pvalue"))
    cat(sprintf("  KLoSA 2-factor fit: CFI=%.3f, TLI=%.3f, RMSEA=%.3f, SRMR=%.3f\n",
                fm["cfi"], fm["tli"], fm["rmsea"], fm["srmr"]))
    write.csv(
      tibble(cohort="KLoSA", model="Bprime_2factor",
             cfi=fm["cfi"], tli=fm["tli"], rmsea=fm["rmsea"], srmr=fm["srmr"],
             chisq=fm["chisq"], df=fm["df"], pvalue=fm["pvalue"]),
      file.path(RESULTS_DIR, "klosa_cfa_fit.csv"), row.names=FALSE
    )
    # |loading| weights
    std_load <- standardizedSolution(fit) %>%
      filter(op == "=~") %>% select(rhs, est.std) %>%
      group_by(rhs) %>% summarise(mean_loading = mean(est.std, na.rm=TRUE), .groups="drop")
    inds <- c("z_dep","z_sight","z_hearing","z_grip","z_srh")
    w <- setNames(rep(NA_real_, length(inds)), inds)
    for (ind in inds) {
      if (ind %in% std_load$rhs) w[ind] <- abs(std_load$mean_loading[std_load$rhs == ind])
    }
    w[is.na(w) | w < 0.05] <- 1.0  # fallback if indicator not in CFA or weight near 0
    cat("  |loading| weights (2-factor CFA, fallback=1 if weak):\n"); print(round(w, 3))
  } else {
    w <- c(z_dep=1, z_sight=1, z_hearing=1, z_grip=1, z_srh=1)
    cat("  CFA unavailable → equal weights\n")
  }

  # --- ic_cfa: |loading|-weighted composite, ≥3/5 rule ---
  z_cols <- names(w)
  z_mat  <- as.matrix(d[, z_cols])
  wval   <- sweep(z_mat, 2, w, "*")
  wavail <- sweep((!is.na(z_mat)) * 1, 2, w, "*")
  n_av   <- rowSums(!is.na(z_mat))
  d$ic_cfa <- ifelse(n_av >= 3,
                     rowSums(wval, na.rm=TRUE) / rowSums(wavail),
                     NA_real_)
  cat(sprintf("  ic_cfa coverage: %.1f%%\n", 100 * mean(!is.na(d$ic_cfa))))

  # --- env_formative: equal-weight z-avg, ≥2/3 rule ---
  env_mat <- as.matrix(d[, c("z_fin","z_liv","z_soc")])
  n_env   <- rowSums(!is.na(env_mat))
  d$env_formative <- ifelse(n_env >= 2, rowMeans(env_mat, na.rm=TRUE), NA_real_)
  cat(sprintf("  env_formative coverage: %.1f%%\n", 100 * mean(!is.na(d$env_formative))))

  d$ic_env_cont <- d$ic_cfa * d$env_formative

  d
}

# ============================================================================
# STEP 3: Construct survival dataset
# ============================================================================

step3_survival <- function(d) {
  cat("\n", strrep("=", 70), "\n"); cat("Step 3: Survival dataset\n"); cat(strrep("=", 70), "\n\n")

  # Baseline: first wave per person with non-missing adl_total
  baseline <- d %>%
    filter(!is.na(adl_total)) %>%
    group_by(id) %>% slice_min(wave, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    rename(baseline_wave = wave, baseline_year = interview_year)

  cat(sprintf("  Baseline (first wave with ADL) N = %s\n",
              format(nrow(baseline), big.mark=",")))

  eligible <- baseline %>% filter(adl_total == 0)
  cat(sprintf("  After ADL=0 filter: N = %s\n", format(nrow(eligible), big.mark=",")))

  disability <- d %>%
    inner_join(eligible %>% select(id, baseline_wave), by="id") %>%
    filter(wave > baseline_wave & !is.na(adl_total) & adl_total >= 1) %>%
    group_by(id) %>% slice_min(wave, n=1, with_ties=FALSE) %>%
    ungroup() %>% select(id, disability_year = interview_year)

  last_obs <- d %>%
    filter(!is.na(interview_year)) %>%
    group_by(id) %>% summarise(last_year = max(interview_year, na.rm=TRUE), .groups="drop")

  surv <- eligible %>%
    left_join(disability, by="id") %>%
    left_join(last_obs, by="id") %>%
    mutate(
      event_disability = as.integer(!is.na(disability_year)),
      event_death      = as.integer(!is.na(death_year)),
      event_dfs        = as.integer(event_disability == 1 | event_death == 1),
      # End time: min of (disability onset, death, last obs); all in year units
      t_end = pmin(disability_year, death_year, last_year, na.rm=TRUE),
      time_dfs = t_end - baseline_year
    ) %>%
    filter(!is.na(time_dfs) & time_dfs > 0)

  cat(sprintf("  Analytic survival N = %s, DFS events = %s (%.1f%%)\n",
              format(nrow(surv), big.mark=","),
              format(sum(surv$event_dfs), big.mark=","),
              100 * mean(surv$event_dfs)))
  cat(sprintf("  Median follow-up = %.1f years\n", median(surv$time_dfs)))

  saveRDS(surv, file.path(DATA_DIR, "klosa_analytic.rds"))
  surv
}

# ============================================================================
# STEP 4: Cox M4 interaction on KLoSA
# ============================================================================

step4_cox <- function(surv) {
  cat("\n", strrep("=", 70), "\n"); cat("Step 4: Cox M4 interaction\n"); cat(strrep("=", 70), "\n\n")

  dat <- surv %>%
    filter(!is.na(ic_cfa) & !is.na(env_formative) & !is.na(ic_env_cont) &
             !is.na(age) & !is.na(female) & !is.na(chronic_n))

  cat(sprintf("  Complete-case N = %s, events = %s\n",
              format(nrow(dat), big.mark=","),
              format(sum(dat$event_dfs), big.mark=",")))

  if (nrow(dat) < 500 || sum(dat$event_dfs) < 100) {
    cat("  WARNING: sample/events low — replication will be underpowered.\n")
  }

  fit <- coxph(
    Surv(time_dfs, event_dfs) ~ ic_cfa + env_formative + ic_env_cont +
      age + female + chronic_n,
    data = dat
  )
  print(summary(fit))

  s <- summary(fit)$coefficients %>% as.data.frame() %>%
    tibble::rownames_to_column("variable")

  # Export all model terms (not just ic_env_cont) for transparency
  out <- s %>% transmute(
    cohort="KLoSA", model="M4_Interaction_reduced", variable,
    log_hr = coef, se_log_hr = `se(coef)`,
    hr = `exp(coef)`, p_value = `Pr(>|z|)`,
    hr_lower = exp(log_hr - 1.96 * se_log_hr),
    hr_upper = exp(log_hr + 1.96 * se_log_hr),
    n = nrow(dat), n_events = sum(dat$event_dfs)
  )
  write.csv(out, file.path(RESULTS_DIR, "klosa_replication_m4.csv"), row.names=FALSE)

  int_row <- out %>% filter(variable == "ic_env_cont")
  cat(sprintf("\n  *** KLoSA IC×Env interaction HR = %.3f (95%% CI %.3f-%.3f), p = %.4f ***\n",
              int_row$hr, int_row$hr_lower, int_row$hr_upper, int_row$p_value))
  cat("\n  (Qualitative replication — reduced 5-indicator IC, 3-indicator Env)\n")

  out
}

# ============================================================================
# MAIN
# ============================================================================

main <- function() {
  t0 <- Sys.time()
  long <- step1_load()
  d    <- step2_composites(long)
  surv <- step3_survival(d)
  out  <- step4_cox(surv)
  cat(sprintf("\nTotal runtime: %.1f min\n",
              as.numeric(Sys.time() - t0, units="mins")))
}

if (sys.nframe() == 0) main()
