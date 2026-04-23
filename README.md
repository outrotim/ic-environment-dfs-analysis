# IC–Environment Interaction and Disability-Free Survival — Analysis Code

> **Companion code repository for the manuscript:**
> *Environmental support delays disability onset most when intrinsic capacity is preserved: a harmonised analysis of 102,006 older adults across five international cohorts with qualitative replication in a sixth.*

This repository contains the full analytic pipeline (18 R scripts) used to produce the variable harmonisation, primary two-stage IPD meta-analysis, causal marginal structural model (MSM/IPTW), competing-risk decomposition, KLoSA qualitative replication, IPCW loss-to-follow-up sensitivity, trajectory analysis, and all 20 pre-specified sensitivity analyses reported in the manuscript.

**Repository purpose**: methodological transparency and independent re-analysis by readers with access to the source cohorts. It is **not** a plug-and-play pipeline and does **not** include data.

---

## Repository contents

### Root
| File | Purpose |
|---|---|
| `README.md` | This file |
| `LICENSE` | MIT License (code) |
| `sessionInfo.txt` | Exact R and package versions used |
| `.gitignore` | Excludes data, results, figures, cache |

### Numbered R scripts (pipeline order)

| # | Script | Purpose |
|---|---|---|
| 01 | `01_data_load.R` | Load 5 cohort CSVs (HRS, ELSA, CHARLS, SHARE, MHAS) into a pooled long-format dataset |
| 02 | `02_harmonization.R` | Harmonise IC indicators (z-scores, reverse-coding) and environmental variables across cohorts |
| 03 | `03_ic_cfa.R` | Multi-group confirmatory factor analysis (MG-CFA) of IC; sequential measurement invariance testing |
| 04 | `04_env_index.R` | Formative environment composite (social participation + financial + living arrangement) |
| 05 | `05_descriptive.R` | Baseline Table 1 (by cohort; by IC × Env group) and variable availability matrix |
| 06 | `06_ipd_ma.R` | **Primary analysis**: two-stage IPD meta-analysis — cohort-specific Cox models (M1–M5) and REML random-effects pooling |
| 07 | `07_trajectory.R` | Latent-class growth analysis (LCGA) for IC and environment trajectories |
| 08 | `08_msm.R` | **Causal analysis**: marginal structural model with stabilised IPTW; time-varying confounding adjustment |
| 09 | `09_sensitivity.R` | 20 pre-specified sensitivity analyses |
| 10 | `10_figures.R` | Publication figures (Figures 1–3 + eFigures 1–2) |
| 11 | `11_strobe_flowchart.R` | STROBE-IPD participant flow diagram |
| 12 | `12_consistency_check.R` | Numerical consistency audit against manuscript claims |
| 13 | `13_table1_from_surv.R` | Regenerate Table 1 from the survival analytic dataset |
| 14 | `14_competing_risks.R` | Cause-specific hazard (CSH) decomposition: CSH-Disability vs CSH-Death |
| 15 | `15_objective_env.R` | Structural-only environment proxy (financial + living, excluding social participation) |
| 16 | `16_marginal_effect_plot.R` | Marginal effect of environment across IC continuum (delta method + REML pooling) |
| 17 | `17_klosa_replication.R` | Qualitative external replication in KLoSA (reduced 5-indicator IC, 2-factor CFA) |
| 18 | `18_ipcw_sensitivity.R` | IPCW sensitivity analysis for the 44,555 eligible baseline participants excluded for loss-to-follow-up |

---

## Data availability (layered statement)

This repository is **code-only**. Individual participant data are not redistributable under the respective data-use agreements (DUAs).

### Source cohort access

| Cohort | Access route | Access level |
|---|---|---|
| Health and Retirement Study (HRS) | https://hrs.isr.umich.edu | Registered researchers (no fee) |
| English Longitudinal Study of Ageing (ELSA) | https://www.elsa-project.ac.uk | Registered researchers |
| China Health and Retirement Longitudinal Study (CHARLS) | https://charls.pku.edu.cn | Registered researchers |
| Survey of Health, Ageing and Retirement in Europe (SHARE) | https://share-eric.eu | Research contract |
| Mexican Health and Aging Study (MHAS) | https://www.mhasweb.org | Registration required |
| Korean Longitudinal Study of Aging (KLoSA) | https://survey.keis.or.kr | Registered researchers |

All main analyses used **Harmonized** files distributed via the Gateway to Global Aging Data (https://g2aging.org).

### What is **not** in this repository

- Individual participant data in any format (prohibited by cohort DUAs)
- Intermediate analytic files (`data/*.rds`) derived from protected data
- Raw cohort CSV files

---

## How to reproduce the analyses

Prerequisites:
1. R ≥ 4.5 with packages listed in `sessionInfo.txt` (tidyverse, survival, metafor, lavaan, haven, lcmm).
2. Access to the six cohort Harmonized data files.
3. Familiarity with each cohort's variable naming conventions.

Pipeline outline:

```r
# Set working directory to the cloned repository root
# setwd("path/to/ic-environment-dfs-analysis")

# Place raw cohort CSVs at data/raw/{hrs,elsa,charls,share,mhas}.csv
# (scripts expect this layout)

# Step 1: Load cohort data and pool into long format
source("01_data_load.R")

# Step 2: Harmonise IC + environment indicators
source("02_harmonization.R")

# Step 3: MG-CFA with invariance testing
source("03_ic_cfa.R")

# Step 4: Environment formative composite
source("04_env_index.R")

# Step 5: Descriptive tables
source("05_descriptive.R")

# Step 6: Primary IPD meta-analysis (Cox M1-M5)
source("06_ipd_ma.R")

# Step 7: LCGA trajectory analysis
source("07_trajectory.R")

# Step 8: MSM with IPTW
source("08_msm.R")

# Steps 9-18: Sensitivity, figures, supplementary analyses
# (run individually as needed)
```

---

## Important caveats for reusers

1. **Measurement invariance**: The three-factor IC CFA achieved scalar (not strict) invariance across the five main cohorts. Users extending this IC composite to new populations should re-fit the MG-CFA and verify at least scalar invariance before cross-population comparisons.
2. **Proportional hazards**: The M4 global PH test was violated in SHARE (p<0·001) while the other four cohorts passed; individual key terms (IC, environment, interaction) passed in SHARE. The two-stage random-effects meta-analytic framework absorbs between-cohort PH heterogeneity.
3. **KLoSA reduced indicator set**: The Harmonized KLoSA E.2 release lacks cognitive indicators; `17_klosa_replication.R` uses a reduced 5-indicator IC with a 2-factor CFA. Direct translation of the 7-indicator composite to KLoSA is not supported.
4. **Environment formative composite**: The index is formative (not reflective); users should not apply internal-consistency reliability metrics (e.g., Cronbach's α) designed for reflective scales.
5. **Selection bias diagnostics**: 44,555 of 146,561 eligible baseline participants (30%) had no post-baseline follow-up and were excluded from the analytic cohort. The IPCW sensitivity (`18_ipcw_sensitivity.R`) shows the pooled interaction HR is not materially biased by this exclusion (HR 0·898 unweighted → 0·895 IPCW-adjusted).

---

## License

- **Code** (`R/*.R`, `*.txt`): [MIT License](LICENSE).
- **Documentation** (`README.md`): [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

## Citation

If you use or adapt this code, please cite the companion manuscript:

```
[Authors]. Environmental support delays disability onset most when intrinsic
capacity is preserved: a harmonised analysis of 102,006 older adults across
five international cohorts with qualitative replication in a sixth.
[Journal] [Year]; [Volume]: [Pages]. DOI: [DOI]
```

A `CITATION.cff` file will be added upon manuscript acceptance.

---

## Acknowledgments

We acknowledge the data providers for HRS, ELSA, CHARLS, SHARE, MHAS, KLoSA, and the Gateway to Global Aging Data. Full cohort-specific acknowledgments (including grant numbers) are provided in the companion manuscript.

---

## Contact

Contact details for the corresponding author appear in the published manuscript. Issues with the code may be opened via GitHub Issues.
