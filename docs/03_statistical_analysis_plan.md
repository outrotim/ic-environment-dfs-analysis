# Statistical Analysis Plan (SAP)

## Intrinsic Capacity–Environment Fit and Disability-Free Survival in Older Adults: A Harmonized Individual Participant Data Meta-Analysis

> **Version**: 1.1 (Skeleton) | **Date**: 2026-03-11 | **Status**: DRAFT

---

## 1. Study Design Overview

**Design**: Prospective cohort study using harmonized individual participant data (IPD) from five core international aging cohorts, with one additional sensitivity cohort.

**Core Cohorts (5)**: HRS (US), ELSA (UK), CHARLS (China), SHARE (Europe, 28 countries), MHAS (Mexico).

**Sensitivity Cohort**: KLoSA (South Korea) — lacks gait speed; locomotion domain uses grip strength only.

**Excluded**: LASI (India) — only Wave 1 cross-sectional data available; no mortality data.

**Data Source**: Local harmonized datasets from Gateway to Global Aging Data.
**Data Path**: `${STUDY17_RAW_DATA_DIR}/`

**Population**: Community-dwelling adults aged ≥60 years at baseline with no pre-existing severe ADL disability (≤1 ADL limitation at baseline).

**Primary Objective**: To determine whether environmental support modifies the association between intrinsic capacity (IC) and disability-free survival (DFS), and whether a capacity–environment mismatch index predicts functional outcomes beyond IC alone.

---

## 2. Sample Size Considerations

### 2.1 Expected Analytic Sample

| Cohort | Eligible (≥60, community) | Expected after exclusions | Follow-up waves | Tier |
|--------|--------------------------|--------------------------|-----------------|------|
| HRS | ~12,000 | ~9,000–10,000 | 8–9 waves (2006–2022) | Core |
| ELSA | ~8,000 | ~6,000–7,000 | 9–10 waves (2002–2020) | Core |
| CHARLS | ~10,000 | ~7,000–8,000 | 4–5 waves (2011–2020) | Core |
| SHARE | ~40,000 | ~25,000–30,000 | 6–8 waves (2004–2022) | Core |
| MHAS | ~8,000 | ~5,000–6,000 | 5–6 waves (2001–2018) | Core |
| **Core Total** | **~78,000** | **~52,000–61,000** | — | — |
| *KLoSA* | *~6,000* | *~4,000–5,000* | *7–8 waves (2006–2020)* | *Sensitivity* |

### 2.2 Power Considerations

With N ≈ 55,000 (5 core cohorts) and expected DFS event rate of 25–35% over follow-up:
- Power > 99% to detect main effects of IC and Env on DFS (HR ≥ 1.10)
- Power > 95% to detect multiplicative interaction (HR ratio ≥ 1.15) at α = 0.05
- Power > 80% for subgroup analyses (by sex, education tertile) with N ≥ 5,000/group
- Sufficient power for trajectory class comparisons with ≥3 classes of N ≥ 2,000 each

No formal sample size calculation is required as this is a secondary analysis of existing data; the above provides reassurance that interactions are detectable.

---

## 3. Variable Construction

### 3.1 Intrinsic Capacity (IC) Score

#### 3.1.1 Indicator Selection
Nine indicators across WHO 5 domains (all available 4/4 cohorts):

| Domain | Indicators | Scoring |
|--------|-----------|---------|
| Cognition | Immediate word recall (0–10), Delayed word recall (0–10), Serial 7s (0–5) | Higher = better |
| Psychological | Depression score (CES-D 8/10 or EURO-D) | Reverse-coded: Higher = better |
| Sensory | Self-rated vision (1–5), Self-rated hearing (1–5) | Reverse-coded: Higher = better |
| Locomotion | Grip strength (kg), Gait speed (m/s, where available) | Higher = better |
| Vitality | BMI (U-shaped transformation), Self-rated health (1–5) | Higher = better |

#### 3.1.2 Multi-Group Confirmatory Factor Analysis (MG-CFA)

**Model specification**:
- Single second-order factor (IC) loading on 5 domain-level first-order factors
- Alternative: Bifactor model with general IC factor + domain-specific factors

**Grouping variable**: Cohort (HRS, ELSA, CHARLS, SHARE, MHAS; KLoSA in sensitivity)

**Measurement invariance testing** (sequential):

| Level | Constraint | Test | Criterion |
|-------|-----------|------|-----------|
| Configural | Same factor structure | Fit indices | CFI > 0.90, RMSEA < 0.08 |
| Metric | Equal factor loadings | ΔCFI | |ΔCFI| < 0.010 |
| Scalar | Equal intercepts | ΔCFI | |ΔCFI| < 0.010 |
| Strict (optional) | Equal residual variances | ΔCFI | |ΔCFI| < 0.010 |

If full scalar invariance fails, partial invariance will be tested by freeing the most non-invariant intercepts (identified via modification indices) until ΔCFI criterion is met.

**Output**: Factor score for each individual at each wave (EAP or regression-based scoring).

**Software**: R `lavaan` package with `cfa()` and `measurementInvariance()`.

#### 3.1.3 Depression Score Cross-Calibration

Since SHARE uses EURO-D (12 items) while others use CES-D variants:

**Primary approach**: Extract 5–6 conceptually overlapping items (depressed mood, sleep disturbance, loss of interest, fatigue, loneliness, crying) from both scales → create a reduced harmonized depression indicator.

**Sensitivity**: IRT co-calibration using common items as linking anchors across CES-D and EURO-D (2PL or graded response model).

### 3.2 Environment Support Index

#### 3.2.1 Core Variables (4/4 cohorts)

| Variable | Operationalization |
|----------|-------------------|
| Social participation | Count of social activities (0–8+), standardized |
| Social support | Frequency of emotional/instrumental support from family/friends, standardized |
| Living arrangement | Coded as: alone (0), with spouse only (1), with children/others (2) |
| Financial strain | Difficulty paying bills / meeting needs (reverse-coded), standardized |

#### 3.2.2 Construction

**Primary**: Latent variable via CFA (single factor, 4 indicators) → factor score per individual per wave.

**Sensitivity**: Simple sum of z-scores (equally weighted).

**Sensitivity with loneliness**: Add UCLA-3 loneliness score in HRS + ELSA + CHARLS only (3-cohort analysis).

### 3.3 Capacity–Environment Mismatch Index

#### 3.3.1 Continuous Version (Primary Inference)
- Product interaction term: IC_score × Env_score in regression models
- Both variables mean-centered within cohort before computing interaction

#### 3.3.2 Categorical Version (Primary Visualization)
- IC: tertiles (Low / Medium / High) within each cohort
- Env: median split (Low / High) within each cohort
- Cross-classification: 6 groups (3 × 2)
- Key contrasts:
  - "Double jeopardy": Low IC + Low Env (reference)
  - "Compensated": Low IC + High Env
  - "At-risk": High IC + Low Env
  - "Resilient": High IC + High Env

---

## 4. Outcome Definitions

### 4.1 Primary Outcome: Disability-Free Survival (DFS)

**Event**: First occurrence of either:
1. **Sustained ADL/IADL disability**: ≥1 ADL limitation reported at ≥2 consecutive waves, OR
2. **All-cause mortality**

**Time origin**: Date of baseline assessment.
**Time scale**: Years from baseline.
**Censoring**: Last known alive and disability-free assessment date.

### 4.2 Secondary Outcomes

| Outcome | Definition | Type |
|---------|-----------|------|
| All-cause mortality | Death from any cause | Time-to-event |
| Care dependence | ADL disability + receiving help from others | Time-to-event |
| Frailty onset | Transition from non-frail/pre-frail to frail (modified Fried) | Time-to-event |
| Dementia-free survival (cognitive subgroup) | Dementia diagnosis or cognitive score below cohort-specific threshold, or death | Time-to-event |

---

## 5. Statistical Analysis

### 5.1 Descriptive Analysis

1. Baseline characteristics by cohort (Table 1): mean ± SD or median [IQR] for continuous, N (%) for categorical
2. IC score distribution by cohort, sex, age group (Figure: violin/ridgeline plots)
3. Env score distribution by cohort
4. Cross-tabulation of IC × Env categories by cohort
5. Kaplan-Meier curves for DFS by IC tertile, Env group, and 6-group cross-classification (Figure)

### 5.2 Layer 2: Two-Stage IPD Meta-Analysis (Primary Analysis)

#### Stage 1: Within-Cohort Models

For each cohort k (k = 1, ..., 5), fit:

**Model A (Continuous interaction)**:
```
h_k(t | X) = h_0k(t) × exp(β1·IC + β2·Env + β3·IC×Env + γ·Z)
```
where Z = age, sex, education, marital status, chronic disease count, cohort-specific covariates.

**Model B (Categorical groups)**:
```
h_k(t | X) = h_0k(t) × exp(Σ β_g·I(group=g) + γ·Z)
```
where g ∈ {Low-IC/High-Env, Med-IC/Low-Env, Med-IC/High-Env, High-IC/Low-Env, High-IC/High-Env}, reference = Low-IC/Low-Env.

**Model specification**:
- Cox proportional hazards (primary)
- Flexible parametric survival model (Royston-Parmar, sensitivity) for time-varying HRs
- Proportional hazards assumption tested via Schoenfeld residuals and log-log plots

#### Stage 2: Random-Effects Meta-Analysis

Pool cohort-specific estimates using restricted maximum likelihood (REML) random-effects meta-analysis:

- Pool β3 (interaction HR) from Model A across 5 core cohorts
- Pool group-specific HRs from Model B across 5 core cohorts
- Report: pooled HR [95% CI], I², τ², prediction interval
- Forest plots for each parameter

**Heterogeneity assessment**: I² statistic, Cochran's Q, τ² estimation.

**For SHARE**: Optionally treat each country (or country group) as a separate unit in Stage 1, yielding ~30 estimates for Stage 2 → more granular meta-analysis.

#### 5.2.1 Additive Interaction Measures

In addition to multiplicative interaction (HR for IC×Env), compute:
- **RERI** (Relative Excess Risk due to Interaction): RERI = HR_11 − HR_10 − HR_01 + 1
- **AP** (Attributable Proportion due to Interaction): AP = RERI / HR_11
- **SI** (Synergy Index): SI = (HR_11 − 1) / [(HR_10 − 1) + (HR_01 − 1)]

where subscripts denote low IC (1/0) and low Env (1/0).

Bootstrap 95% CIs for additive measures (1000 replicates).

### 5.3 Layer 3: Longitudinal Trajectory Analysis

#### 5.3.1 IC Trajectories

**Method**: Group-Based Trajectory Modeling (GBTM) or Latent Class Growth Analysis (LCGA)

**Specification**:
- Outcome: IC score across waves
- Polynomial order: linear, quadratic, or cubic (selected by BIC)
- Number of classes: 2–6, selected by BIC, entropy, average posterior probability (>0.70), and clinical interpretability
- Covariates: age at baseline, sex (as predictors of class membership)

**Expected trajectories**: "Stable high", "Gradual decline", "Rapid decline", "Late decline", "Improving" (pending data)

#### 5.3.2 Env Trajectories

Same approach as IC trajectories applied to Env score.

#### 5.3.3 Co-Trajectory Analysis

**Method A**: Dual trajectory modeling (joint GBTM for IC and Env simultaneously)
**Method B**: Cross-classify IC trajectory class × Env trajectory class → combined trajectory typology

**Association with DFS**: Cox model with trajectory class as exposure:
```
h(t) = h_0(t) × exp(Σ β_j·I(co-trajectory=j) + γ·Z)
```

**Incremental predictive value**: Compare C-statistic (Harrell's C) of:
1. Model with baseline IC + Env only
2. Model with IC trajectory class
3. Model with IC trajectory + Env trajectory
4. Model with co-trajectory class

Test differences using likelihood ratio test and ΔC-statistic (Pencina method).

### 5.4 Layer 4: Causal Inference — Marginal Structural Models (MSM)

#### 5.4.1 Rationale
IC and Env at each wave may be affected by prior DFS status (reverse causation) and by time-varying confounders (e.g., new chronic disease) that are themselves affected by prior IC/Env. Standard time-varying Cox models may introduce collider bias.

#### 5.4.2 Model Specification

**Step 1: Estimate stabilized inverse probability of treatment weights (IPTW)**

At each wave t, fit:
```
logit[P(Env_t = high | Env_{t-1}, IC_{t-1}, L_t, L_{t-1})]
```
and
```
logit[P(Env_t = high | Env_{t-1})]
```

Stabilized weight:
```
SW_t = Π_{s=1}^{t} [P(Env_s | Env_{s-1})] / [P(Env_s | Env_{s-1}, IC_{s-1}, L_s, L_{s-1})]
```

where L = time-varying confounders (chronic disease count, ADL status, hospitalization).

**Step 2: Weighted pooled logistic regression** (approximating continuous-time MSM):
```
logit[P(Y_{t+1} = 1 | Ȳ_t = 0)] = α_0 + α_1·t + β1·IC_t + β2·Env_t + β3·IC_t×Env_t + γ·V
```
where V = baseline covariates, Y = DFS event, Ȳ_t = event history.

**Weight diagnostics**:
- Distribution of weights (mean, SD, min, max, 1st/99th percentiles)
- Truncation at 1st and 99th percentiles if extreme weights observed
- Covariate balance check (standardized mean differences before/after weighting)

**Step 3: Comparison with naive time-varying Cox**
Report both MSM-weighted and unweighted estimates to quantify bias direction and magnitude.

#### 5.4.3 Alternative Causal Methods (Sensitivity)
- **Parametric g-formula**: Monte Carlo simulation-based estimation of counterfactual DFS under different Env intervention scenarios
- **g-estimation** of structural nested failure time models
- **Target trial emulation**: Define hypothetical intervention (e.g., "high vs low Env from baseline onward") and estimate per-protocol effect

---

## 6. Sensitivity Analyses

| # | Analysis | Rationale |
|---|----------|-----------|
| S1 | IC sum score (z-score average) vs latent score | Measurement robustness |
| S2 | Env sum score vs latent score | Measurement robustness |
| S3 | Exclude baseline ADL ≥ 1 (stricter exclusion) | Reduce reverse causation |
| S4 | Fine-Gray competing risk model (death as competing event for disability) | Competing risks |
| S5 | Stratified by sex | Effect modification |
| S6 | Stratified by education (tertiles) | Effect modification |
| S7 | Stratified by age group (60–69, 70–79, 80+) | Effect modification |
| S8 | Leave-one-cohort-out meta-analysis | Influence of single cohort |
| S8b | Original 4-cohort analysis (exclude MHAS) | Consistency with original design |
| S8c | Include KLoSA (grip-only locomotion domain) | Geographic coverage extension |
| S9 | 5-variable Env index with loneliness (HRS + ELSA + CHARLS only) | Extended environment |
| S10 | Different IC categorization thresholds (median, quartiles) | Dose-response |
| S11 | Multiple imputation (20 datasets, MICE by chained equations) | Missing data |
| S12 | Complete case analysis | Bias assessment |
| S13 | E-value for unmeasured confounding | Unmeasured confounding |
| S14 | One-stage IPD-MA (mixed-effects Cox) vs two-stage | Modeling approach |
| S15 | SHARE country-level analysis (selected large countries) | Within-SHARE heterogeneity |
| S16 | Restrict to ≥65 years | Age threshold sensitivity |
| S17 | Exclude first 2 years of follow-up | Reverse causation |

---

## 7. Subgroup Analyses (Pre-Specified)

| Subgroup | Justification |
|----------|--------------|
| Sex (male vs female) | IC components and social environments differ by sex |
| Age group (60–69, 70–79, ≥80) | Compensation effect may vary by age |
| Education (low vs medium vs high) | Cognitive reserve and resource access |
| Urban vs rural | Environmental support structures differ |
| Cohort/region (East Asia vs Europe vs Americas vs UK) | Cultural and system differences |
| Baseline frailty status (non-frail vs pre-frail) | Effect may be stronger in pre-frail |
| Multimorbidity (0–1 vs ≥2 chronic conditions) | Disease burden interaction |

---

## 8. Missing Data

### 8.1 Primary Strategy: Multiple Imputation

- Method: Multiple Imputation by Chained Equations (MICE)
- Number of imputations: M = 20
- Imputation model: includes all analysis variables + auxiliary variables (variables associated with missingness but not in analysis model)
- Performed separately within each cohort (preserving cohort-specific distributions)
- Predictive mean matching (PMM) for continuous, logistic for binary, multinomial logit for categorical
- Convergence assessed via trace plots

### 8.2 Secondary Strategies

- Complete case analysis (sensitivity S12)
- Inverse probability of censoring weighting (IPCW) for informative dropout
- Pattern mixture models (if missing data exceeds 30% for any key variable)

---

## 9. Reporting

### 9.1 Tables

| Table | Content |
|-------|---------|
| Table 1 | Baseline characteristics by cohort |
| Table 2 | IC domain scores and IC composite score by cohort |
| Table 3 | Environment support indicators by cohort |
| Table 4 | Two-stage IPD-MA: main effects and interaction (continuous) |
| Table 5 | Two-stage IPD-MA: 6-group categorical analysis |
| Table 6 | Additive interaction measures (RERI, AP, SI) |
| Table 7 | Trajectory class descriptions and prevalence |
| Table 8 | Co-trajectory classes and DFS association |
| Table 9 | MSM results vs naive estimates comparison |
| eTable 1–N | Sensitivity and subgroup analyses |

### 9.2 Figures

| Figure | Content |
|--------|---------|
| Figure 1 | Study flow diagram (STROBE) |
| Figure 2 | Conceptual framework: IC × Environment → Functional Ability (WHO model) |
| Figure 3 | KM curves for DFS by 4-group IC×Env classification |
| Figure 4 | Forest plot: pooled HRs from two-stage IPD-MA |
| Figure 5 | IC and Env trajectory classes (spaghetti + mean trajectory) |
| Figure 6 | Heat map: co-trajectory × DFS risk |
| Figure 7 | Additive interaction (RERI) visualization |
| eFigure 1–N | Measurement invariance, weight diagnostics, subgroup forests |

### 9.3 Reporting Guidelines

- **STROBE** for observational studies
- **PRISMA-IPD** for IPD meta-analysis components
- **GRIPS** for risk prediction components (if applicable)
- **RECORD** for routinely collected health data

---

## 10. Software

| Task | Software | Packages |
|------|----------|----------|
| Data management | R 4.3+ / Python 3.10+ | `tidyverse`, `data.table`, `pandas` |
| CFA / measurement invariance | R | `lavaan`, `semTools` |
| IRT calibration | R | `mirt`, `ltm` |
| Survival analysis | R | `survival`, `survminer`, `flexsurv` |
| IPD meta-analysis | R | `ipdmetan` (Stata), `metafor`, `meta` |
| Trajectory modeling | R / SAS | `lcmm`, `flexmix`, `crimCV` / PROC TRAJ |
| MSM / IPTW | R | `ipw`, `WeightIt`, `cobalt`, `survey` |
| G-formula | R | `gfoRmula`, `ltmle` |
| Multiple imputation | R | `mice`, `Amelia` |
| Competing risks | R | `cmprsk`, `tidycmprsk` |
| Visualization | R | `ggplot2`, `patchwork`, `forestplot` |

---

## 11. Timeline (Projected)

| Phase | Duration | Milestones |
|-------|----------|-----------|
| Protocol & registration | 2–3 weeks | PROSPERO/OSF registration (data already obtained locally) |
| Data loading & cleaning | 2–3 weeks | Harmonized analytic dataset from local files |
| Measurement alignment (Layer 1) | 3–4 weeks | IC/Env scores, invariance report |
| Primary analysis (Layer 2) | 3–4 weeks | IPD-MA results, main tables/figures |
| Trajectory analysis (Layer 3) | 3–4 weeks | Trajectory classes, co-trajectory results |
| MSM analysis (Layer 4) | 3–4 weeks | Causal estimates, weight diagnostics |
| Sensitivity & subgroup | 2–3 weeks | All sensitivity results |
| Manuscript drafting | 4–6 weeks | Full manuscript + supplement |
| Internal review & revision | 2–3 weeks | Final version |
| **Total** | **~6–8 months** | — |

---

## 12. Ethical Considerations

- All data are de-identified and publicly available for research use
- Data use agreements will be obtained for each cohort as required
- No individual-level data will be shared; only aggregate results reported
- Study protocol will be pre-registered on OSF or PROSPERO (as appropriate)
- IRB exemption or approval will be obtained from the lead institution as required by local regulations

---

## Appendix A: DAG (Directed Acyclic Graph)

```
                    ┌─────────────┐
                    │ Baseline     │
                    │ Confounders  │
                    │ (age, sex,   │
                    │  education,  │
                    │  SES, comor) │
                    └──────┬──────┘
                           │
              ┌────────────┼────────────┐
              ▼            ▼            ▼
        ┌──────────┐ ┌──────────┐ ┌──────────────────┐
        │   IC_t   │ │  Env_t   │ │ Time-varying     │
        │(Intrinsic│ │(Environ- │ │ confounders L_t  │
        │ Capacity)│ │  ment)   │ │ (new disease,    │
        └────┬─────┘ └────┬─────┘ │  hospitalization) │
             │             │       └────────┬─────────┘
             │             │                │
             ▼             ▼                │
        ┌────────────────────────┐          │
        │   IC_t × Env_t         │◄─────────┘
        │   (Interaction)        │
        └───────────┬────────────┘
                    │
                    ▼
        ┌────────────────────────┐
        │  Disability-Free       │
        │  Survival (DFS)        │
        └────────────────────────┘
```

**Key causal assumptions**:
1. No unmeasured confounding (conditional on Z) — assessed by E-value
2. Positivity: All IC × Env combinations observed in data
3. Consistency: Well-defined exposure (IC score, Env score)
4. No interference between individuals (SUTVA)
5. For MSM: correct specification of weight models

---

## Appendix B: Planned Code Repository Structure

```
study17_ic_environment_ipd/
├── 01_study_schema.md
├── 02_variable_harmonization_map.md
├── 03_statistical_analysis_plan.md
├── PROGRESS.md
├── data/
│   ├── raw/           # Original downloaded data (not tracked in git)
│   ├── harmonized/    # Harmonized analytic datasets
│   └── codebooks/     # Variable codebooks
├── scripts/
│   ├── 01_data_load.R           # Load from local files (no download needed)
│   ├── 02_harmonization.R
│   ├── 03_ic_cfa.R
│   ├── 04_env_index.R
│   ├── 05_descriptive.R
│   ├── 06_ipd_ma.R
│   ├── 07_trajectory.R
│   ├── 08_msm.R
│   ├── 09_sensitivity.R
│   └── 10_figures.R
├── results/
├── figures/
└── manuscript/
    └── manuscript_draft.md
```
