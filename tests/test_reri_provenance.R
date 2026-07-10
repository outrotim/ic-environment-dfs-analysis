#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
test_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else normalizePath("tests/test_reri_provenance.R")
study_dir <- normalizePath(file.path(dirname(test_file), ".."))

stage1_path <- file.path(study_dir, "results", "reri_stage1.csv")
stage2_path <- file.path(study_dir, "results", "reri_stage2_pooled.csv")

stopifnot(
  "RERI stage-1 output is missing" = file.exists(stage1_path),
  "RERI pooled output is missing" = file.exists(stage2_path)
)

stage1 <- read.csv(stage1_path, check.names = FALSE)
stage2 <- read.csv(stage2_path, check.names = FALSE)

required_stage1 <- c(
  "cohort", "model_source", "reference_group", "reri", "se", "ci_lower",
  "ci_upper", "p_value", "n", "n_events"
)
required_stage2 <- c(
  "model_source", "estimand", "method", "k", "pooled_reri", "se",
  "ci_lower", "ci_upper", "p_value", "tau2", "I2"
)

stopifnot(
  "RERI stage-1 schema is incomplete" = all(required_stage1 %in% names(stage1)),
  "RERI pooled schema is incomplete" = all(required_stage2 %in% names(stage2)),
  "RERI must be computed from cohort-specific M5 models" = all(stage1$model_source == "M5_6group"),
  "RERI pooled provenance must name M5" = identical(stage2$model_source[[1]], "M5_6group"),
  "All five cohorts must contribute" = nrow(stage1) == 5L && stage2$k[[1]] == 5L,
  "RERI uncertainty must be finite" = all(is.finite(stage1$se)) && all(is.finite(stage2$se)),
  "RERI confidence intervals must contain their estimates" =
    all(stage1$ci_lower <= stage1$reri & stage1$reri <= stage1$ci_upper) &&
    stage2$ci_lower[[1]] <= stage2$pooled_reri[[1]] &&
    stage2$pooled_reri[[1]] <= stage2$ci_upper[[1]],
  "RERI p-value must be valid" = stage2$p_value[[1]] >= 0 && stage2$p_value[[1]] <= 1
)

cat("PASS: RERI provenance and uncertainty are reproducible from cohort-specific M5 models.\n")
