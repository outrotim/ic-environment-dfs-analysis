#!/usr/bin/env Rscript
# ==============================================================================
# 19_standardized_5y_risk.R
# Study 17: Five-year standardized DFS event risks by IC stratum
#
# Estimand:
#   Within each cohort and IC tertile, standardize the 5-year probability of
#   disability or death after setting environment support to the cohort-specific
#   25th versus 75th percentile. Observed IC and covariate distributions are
#   retained within each stratum. Model-based uncertainty conditions on the
#   estimated baseline cumulative hazard and propagates coefficient covariance.
#
# Outputs:
#   results/standardized_5y_risk_stage1.csv
#   results/standardized_5y_risk_pooled.csv
#   figures/fig2c_standardized_5y_risk.pdf
#   figures/fig2c_standardized_5y_risk.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(metafor)
})

SCRIPT_DIR <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) normalizePath(Sys.getenv("STUDY17_PROJECT_DIR", unset = getwd()), mustWork = FALSE)
)
DATA_DIR <- file.path(SCRIPT_DIR, "data")
RESULTS_DIR <- file.path(SCRIPT_DIR, "results")
FIGURE_DIR <- file.path(SCRIPT_DIR, "figures")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

TIME_HORIZON <- 5

design_matrix <- function(fit, newdata) {
  matrix <- model.matrix(delete.response(terms(fit)), data = newdata)
  if ("(Intercept)" %in% colnames(matrix)) {
    matrix <- matrix[, colnames(matrix) != "(Intercept)", drop = FALSE]
  }
  matrix[, names(coef(fit)), drop = FALSE]
}

standardized_risk <- function(fit, newdata, baseline_hazard) {
  x <- design_matrix(fit, newdata)
  beta <- coef(fit)
  lp <- drop(x %*% beta)
  relative_hazard <- exp(lp)
  survival_probability <- exp(-baseline_hazard * relative_hazard)
  risk <- mean(1 - survival_probability)

  derivative_eta <- survival_probability * baseline_hazard * relative_hazard
  gradient <- colMeans(x * derivative_eta)
  variance <- drop(t(gradient) %*% vcov(fit) %*% gradient)

  list(
    risk = risk,
    se = sqrt(max(variance, 0)),
    gradient = gradient
  )
}

pool_logit_risk <- function(data, risk_column, se_column) {
  risk <- pmin(pmax(data[[risk_column]], 1e-6), 1 - 1e-6)
  yi <- qlogis(risk)
  sei <- data[[se_column]] / (risk * (1 - risk))
  fit <- rma(yi = yi, sei = sei, method = "REML")
  tibble(
    risk = plogis(as.numeric(fit$beta)),
    ci_lower = plogis(as.numeric(fit$ci.lb)),
    ci_upper = plogis(as.numeric(fit$ci.ub)),
    I2 = fit$I2
  )
}

run_stage1 <- function() {
  surv <- readRDS(file.path(DATA_DIR, "survival_cohort.rds")) %>%
    filter(
      !is.na(ic_cfa_c) & !is.na(env_form_c) & !is.na(ic_env_cont) &
        !is.na(ic_cat) & !is.na(age) & !is.na(female) & !is.na(chronic_n) &
        time_dfs > 0
    ) %>%
    mutate(ic_cat = factor(ic_cat, levels = c("Low", "Medium", "High")))

  rows <- list()
  for (coh in sort(unique(surv$cohort))) {
    dat <- surv %>% filter(cohort == coh) %>% droplevels()
    fit <- coxph(
      Surv(time_dfs, event_dfs) ~
        ic_cfa_c + env_form_c + ic_env_cont + age + female + chronic_n,
      data = dat,
      x = TRUE
    )

    baseline <- basehaz(fit, centered = FALSE) %>% filter(time <= TIME_HORIZON)
    if (nrow(baseline) == 0) stop("No baseline hazard before year 5 in ", coh)
    baseline_hazard <- tail(baseline$hazard, 1)
    env_values <- quantile(dat$env_form_c, probs = c(0.25, 0.75), na.rm = TRUE)

    for (stratum in c("Low", "Medium", "High")) {
      target <- dat %>% filter(ic_cat == stratum)
      low_env <- target %>%
        mutate(
          env_form_c = unname(env_values[[1]]),
          ic_env_cont = ic_cfa_c * env_form_c
        )
      high_env <- target %>%
        mutate(
          env_form_c = unname(env_values[[2]]),
          ic_env_cont = ic_cfa_c * env_form_c
        )

      low <- standardized_risk(fit, low_env, baseline_hazard)
      high <- standardized_risk(fit, high_env, baseline_hazard)
      gradient_difference <- low$gradient - high$gradient
      variance_difference <- drop(
        t(gradient_difference) %*% vcov(fit) %*% gradient_difference
      )
      se_difference <- sqrt(max(variance_difference, 0))
      risk_difference <- low$risk - high$risk

      rows[[paste(coh, stratum, sep = "_")]] <- tibble(
        cohort = coh,
        ic_stratum = stratum,
        outcome = "DFS event (first ADL disability or death)",
        time_horizon_years = TIME_HORIZON,
        standardization = "Observed covariates within IC stratum; environment set to cohort-specific P25 vs P75",
        env_low_value = unname(env_values[[1]]),
        env_high_value = unname(env_values[[2]]),
        risk_low_env = low$risk,
        se_risk_low_env = low$se,
        risk_low_env_ci_lower = pmax(0, low$risk - 1.96 * low$se),
        risk_low_env_ci_upper = pmin(1, low$risk + 1.96 * low$se),
        risk_high_env = high$risk,
        se_risk_high_env = high$se,
        risk_high_env_ci_lower = pmax(0, high$risk - 1.96 * high$se),
        risk_high_env_ci_upper = pmin(1, high$risk + 1.96 * high$se),
        risk_difference = risk_difference,
        se_risk_difference = se_difference,
        rd_ci_lower = risk_difference - 1.96 * se_difference,
        rd_ci_upper = risk_difference + 1.96 * se_difference,
        p_value = 2 * pnorm(abs(risk_difference / se_difference), lower.tail = FALSE),
        n = nrow(target),
        n_events = sum(target$event_dfs),
        baseline_hazard_5y = baseline_hazard
      )
    }
  }

  stage1 <- bind_rows(rows)
  write_csv(stage1, file.path(RESULTS_DIR, "standardized_5y_risk_stage1.csv"))
  stage1
}

run_stage2 <- function(stage1) {
  pooled_rows <- list()
  for (stratum in c("Low", "Medium", "High")) {
    data <- stage1 %>% filter(ic_stratum == stratum)
    rd_fit <- rma(
      yi = risk_difference,
      sei = se_risk_difference,
      data = data,
      method = "REML"
    )
    low <- pool_logit_risk(data, "risk_low_env", "se_risk_low_env")
    high <- pool_logit_risk(data, "risk_high_env", "se_risk_high_env")

    pooled_rows[[stratum]] <- tibble(
      ic_stratum = stratum,
      outcome = "DFS event (first ADL disability or death)",
      time_horizon_years = TIME_HORIZON,
      environment_contrast = "cohort-specific P25 vs P75",
      k = rd_fit$k,
      pooled_risk_low_env = low$risk,
      risk_low_env_ci_lower = low$ci_lower,
      risk_low_env_ci_upper = low$ci_upper,
      pooled_risk_high_env = high$risk,
      risk_high_env_ci_lower = high$ci_lower,
      risk_high_env_ci_upper = high$ci_upper,
      pooled_risk_difference = as.numeric(rd_fit$beta),
      se = as.numeric(rd_fit$se),
      ci_lower = as.numeric(rd_fit$ci.lb),
      ci_upper = as.numeric(rd_fit$ci.ub),
      p_value = as.numeric(rd_fit$pval),
      tau2 = rd_fit$tau2,
      I2 = rd_fit$I2
    )
  }

  pooled <- bind_rows(pooled_rows) %>%
    mutate(ic_stratum = factor(ic_stratum, levels = c("Low", "Medium", "High"))) %>%
    arrange(ic_stratum) %>%
    mutate(ic_stratum = as.character(ic_stratum))
  write_csv(pooled, file.path(RESULTS_DIR, "standardized_5y_risk_pooled.csv"))
  pooled
}

plot_results <- function(pooled) {
  plot_data <- bind_rows(
    pooled %>% transmute(
      ic_stratum,
      environment = "Lower support (P25)",
      risk = pooled_risk_low_env,
      ci_lower = risk_low_env_ci_lower,
      ci_upper = risk_low_env_ci_upper
    ),
    pooled %>% transmute(
      ic_stratum,
      environment = "Higher support (P75)",
      risk = pooled_risk_high_env,
      ci_lower = risk_high_env_ci_lower,
      ci_upper = risk_high_env_ci_upper
    )
  ) %>%
    mutate(
      ic_stratum = factor(ic_stratum, levels = c("Low", "Medium", "High")),
      environment = factor(
        environment,
        levels = c("Lower support (P25)", "Higher support (P75)")
      )
    )

  rd_labels <- pooled %>%
    mutate(
      ic_stratum = factor(ic_stratum, levels = c("Low", "Medium", "High")),
      label = sprintf("RD %.1f pp", 100 * pooled_risk_difference),
      y = pmax(pooled_risk_low_env, pooled_risk_high_env) + 0.025
    )

  figure <- ggplot(
    plot_data,
    aes(x = ic_stratum, y = risk, color = environment, group = environment)
  ) +
    geom_line(linewidth = 0.7, position = position_dodge(width = 0.15)) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.08,
      linewidth = 0.55,
      position = position_dodge(width = 0.15)
    ) +
    geom_point(size = 2.5, position = position_dodge(width = 0.15)) +
    geom_text(
      data = rd_labels,
      aes(x = ic_stratum, y = y, label = label),
      inherit.aes = FALSE,
      size = 2.8,
      color = "grey20"
    ) +
    scale_color_manual(
      values = c("Lower support (P25)" = "#D55E00", "Higher support (P75)" = "#0072B2"),
      name = NULL
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.14))
    ) +
    labs(
      x = "Intrinsic capacity tertile",
      y = "Standardized 5-year DFS event risk"
    ) +
    theme_minimal(base_size = 10, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(color = "grey20"),
      axis.title = element_text(color = "grey10"),
      legend.position = "bottom",
      plot.margin = margin(8, 12, 8, 8)
    )

  ggsave(
    file.path(FIGURE_DIR, "fig2c_standardized_5y_risk.pdf"),
    figure,
    width = 7.2,
    height = 4.6,
    device = cairo_pdf
  )
  ggsave(
    file.path(FIGURE_DIR, "fig2c_standardized_5y_risk.tiff"),
    figure,
    width = 7.2,
    height = 4.6,
    dpi = 300,
    compression = "lzw"
  )
}

main <- function() {
  stage1 <- run_stage1()
  pooled <- run_stage2(stage1)
  plot_results(pooled)
  print(
    pooled %>%
      transmute(
        ic_stratum,
        risk_low_env = pooled_risk_low_env,
        risk_high_env = pooled_risk_high_env,
        risk_difference = pooled_risk_difference,
        ci_lower,
        ci_upper,
        p_value,
        I2
      ) %>%
      mutate(across(where(is.numeric), ~ round(., 4)))
  )
}

if (sys.nframe() == 0) main()
