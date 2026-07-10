#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
parse_bound <- function(prefix, default) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  as.integer(sub(paste0("^", prefix, "="), "", hit[[1]]))
}

from <- parse_bound("--from", 1L)
to <- parse_bound("--to", 19L)
if (is.na(from) || is.na(to) || from < 1L || to > 19L || from > to) {
  stop("Use --from=1..19 and --to=1..19 with from <= to")
}

scripts <- c(
  "code/01_data_load.R",
  "code/02_harmonization.R",
  "code/03_ic_cfa.R",
  "code/04_env_index.R",
  "code/05_descriptive.R",
  "code/06_ipd_ma.R",
  "code/07_trajectory.R",
  "code/08_msm.R",
  "code/09_sensitivity.R",
  "code/10_figures.R",
  "code/11_strobe_flowchart.R",
  "code/12_consistency_check.R",
  "code/13_table1_from_surv.R",
  "code/14_competing_risks.R",
  "code/15_objective_env.R",
  "code/16_marginal_effect_plot.R",
  "code/17_klosa_replication.R",
  "code/18_ipcw_sensitivity.R",
  "code/19_standardized_5y_risk.R"
)

Sys.setenv(STUDY17_PROJECT_DIR = normalizePath(getwd(), mustWork = TRUE))
for (index in seq.int(from, to)) {
  message(sprintf("[%02d/19] %s", index, scripts[[index]]))
  status <- system2("Rscript", scripts[[index]])
  if (!identical(status, 0L)) stop("Pipeline failed at ", scripts[[index]])
}
