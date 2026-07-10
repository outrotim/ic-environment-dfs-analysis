#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
test_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else normalizePath("tests/test_ipcw_joint_sensitivity.R")
study_dir <- normalizePath(file.path(dirname(test_file), ".."))
results_dir <- file.path(study_dir, "results")

paths <- c(
  diagnostics = file.path(results_dir, "ipcw_weights_diagnostics.csv"),
  pooled = file.path(results_dir, "ipcw_stage2_pooled.csv"),
  joint_stage1 = file.path(results_dir, "joint_confounding_stage1.csv"),
  joint_stage2 = file.path(results_dir, "joint_confounding_stage2.csv")
)
stopifnot("Required enhanced sensitivity outputs are missing" = all(file.exists(paths)))

diagnostics <- read.csv(paths[["diagnostics"]], check.names = FALSE)
pooled <- read.csv(paths[["pooled"]], check.names = FALSE)
joint_stage1 <- read.csv(paths[["joint_stage1"]], check.names = FALSE)
joint_stage2 <- read.csv(paths[["joint_stage2"]], check.names = FALSE)

stopifnot(
  "All five cohorts must contribute to enhanced IPCW" =
    nrow(diagnostics) == 5L && pooled$k[[1]] == 5L,
  "IPCW diagnostics must record the denominator formula" =
    "denominator_formula" %in% names(diagnostics),
  "IPCW must condition on baseline IC" =
    all(grepl("ic_cfa", diagnostics$denominator_formula, fixed = TRUE)),
  "IPCW must condition on baseline environment" =
    all(grepl("env_formative", diagnostics$denominator_formula, fixed = TRUE)),
  "IPCW must condition on marital status" =
    all(grepl("marital", diagnostics$denominator_formula, fixed = TRUE)),
  "Education must be used where harmonised" =
    sum(grepl("edu_years", diagnostics$denominator_formula, fixed = TRUE)) == 4L,
  "Truncated IPCW diagnostics must be finite" =
    all(is.finite(diagnostics$sw_tr_mean)) && all(is.finite(diagnostics$sw_tr_max)),
  "Enhanced IPCW estimates must be finite" =
    all(is.finite(unlist(pooled[c("hr", "hr_lower", "hr_upper", "p_value", "i2")]))),
  "Joint confounding comparison must use four education-complete cohorts" =
    length(unique(joint_stage1$cohort)) == 4L && all(joint_stage2$k == 4L),
  "Both same-sample models must be present" =
    setequal(unique(joint_stage2$model), c("base_same_sample", "joint_education_marital")),
  "Joint confounding estimates must be finite" =
    all(is.finite(unlist(joint_stage2[c("hr", "hr_lower", "hr_upper", "p_value", "I2")])))
)

cat("PASS: enhanced IPCW and same-sample joint confounding analyses are reproducible.\n")
