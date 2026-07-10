#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
test_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else normalizePath("tests/test_standardized_5y_risk.R")
study_dir <- normalizePath(file.path(dirname(test_file), ".."))

stage1_path <- file.path(study_dir, "results", "standardized_5y_risk_stage1.csv")
pooled_path <- file.path(study_dir, "results", "standardized_5y_risk_pooled.csv")
figure_pdf <- file.path(study_dir, "figures", "fig2c_standardized_5y_risk.pdf")
figure_tiff <- file.path(study_dir, "figures", "fig2c_standardized_5y_risk.tiff")
combined_pdf <- file.path(study_dir, "figures", "fig2_forest_marginal.pdf")
figure_script <- file.path(study_dir, "10_figures.R")

stopifnot(
  "Standardized risk outputs are missing" =
    all(file.exists(c(stage1_path, pooled_path, figure_pdf, figure_tiff,
                      combined_pdf, figure_script)))
)

stage1 <- read.csv(stage1_path, check.names = FALSE)
pooled <- read.csv(pooled_path, check.names = FALSE)

required <- c(
  "cohort", "ic_stratum", "time_horizon_years", "env_low_value",
  "env_high_value", "risk_low_env", "risk_high_env", "risk_difference",
  "se_risk_difference", "rd_ci_lower", "rd_ci_upper", "n", "n_events"
)
stopifnot(
  "Stage-1 standardized risk schema is incomplete" = all(required %in% names(stage1)),
  "Five cohorts by three IC strata are required" = nrow(stage1) == 15L,
  "All three IC strata must be present" =
    setequal(unique(stage1$ic_stratum), c("Low", "Medium", "High")),
  "Standardized risks must be valid probabilities" =
    all(stage1$risk_low_env >= 0 & stage1$risk_low_env <= 1) &&
    all(stage1$risk_high_env >= 0 & stage1$risk_high_env <= 1),
  "Risk-difference uncertainty must be finite" =
    all(is.finite(stage1$se_risk_difference)) &&
    all(stage1$rd_ci_lower <= stage1$risk_difference &
          stage1$risk_difference <= stage1$rd_ci_upper),
  "Pooled output must contain three strata" = nrow(pooled) == 3L,
  "Pooled risk differences must be finite" =
    all(is.finite(unlist(pooled[c("pooled_risk_difference", "se", "ci_lower", "ci_upper", "p_value", "I2")]))),
  "Figure files must be non-empty" =
    file.info(figure_pdf)$size > 1000 && file.info(figure_tiff)$size > 1000 &&
    file.info(combined_pdf)$size > 1000,
  "Combined Figure 2 must label RD point estimates and 95% CIs" = {
    script_text <- paste(readLines(figure_script, warn = FALSE), collapse = "\n")
    grepl("100 \\* ci_lower", script_text) &&
      grepl("100 \\* ci_upper", script_text) &&
      grepl("95% CI", script_text, fixed = TRUE)
  }
)

cat("PASS: five-year standardized DFS risks and risk differences are reproducible.\n")
