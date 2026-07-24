# Analysis of a Contaminated Phase II RCT for Antihypertensive Efficacy

**Inferencia 2026** — Department of Statistics, Presidency University, Kolkata

**Authors:** Pratik Majumdar, Sourav Dey, Sampad Nag

## Overview

This project evaluates whether a new antihypertensive drug (**Trx**) reduces
systolic blood pressure (SBP) compared to placebo, using a simulated Phase II
RCT dataset deliberately constructed to mimic real-world data problems:

- **Device switch** at Site 2 partway through the trial
- **Differential dropout** across sites
- **Protocol deviations** concentrated at Site 3
- **Missing values** in `age` and `bmi`

The goal was to work through these "contamination" issues the way an
analyst would with real trial data handed off with no documentation, and
determine whether the treatment effect estimate still holds up.

## Key Result

A linear regression of change in SBP (baseline − Week 12) on treatment,
adjusting for site, age, sex, BMI, and baseline SBP, found a statistically
significant treatment effect:

- **Estimate:** 8.80 mmHg greater reduction with Trx vs placebo
- **p < 0.001**, 95% CI [6.41, 11.18]
- Effect size (Cohen's d): 1.32

The result was **robust** to:
- excluding the site with protocol deviations (Site 3) — estimate shifted only ~2.6%
- the device switch at Site 2, which showed no significant measurement shift (t-test p = 0.47)
- the missing-data mechanism, which was consistent with Missing At Random (MAR),
  since dropout was strongly predicted by site (esp. Site 3) rather than
  unobserved factors

## Repository Structure

```
.
├── R/
│   └── analysis.R                  # Full analysis script (cleaning → EDA → inference → sensitivity)
├── data/
│   └── raw/
│       └── trialx_hackathon_data.csv   # Raw simulated trial data (180 obs, 13 vars)
├── output/
│   ├── figures/                    # Histograms, boxplots, density plots
│   └── tables/                     # Balance table, regression output, etc.
├── report/
│   ├── Inferencia_Presentation_2026.lyx   # Beamer slide deck (LyX source)
│   ├── Inferencia_Paper_Presentation.lyx  # Companion paper/report (LyX source)
│   └── Inferencia_Abstract_2026.pdf       # Conference abstract
└── README.md
```

## Data Dictionary

| Variable          | Description                                          |
|-------------------|-------------------------------------------------------|
| `patient_id`       | Unique patient identifier                             |
| `site`             | Clinical site (1, 2, or 3)                             |
| `treatment`        | `Trx` or `Placebo`                                     |
| `age`              | Age in years (some missing)                            |
| `sex`              | Male or Female                                          |
| `bmi`              | Body Mass Index (some missing)                          |
| `sbp_0`            | Systolic BP at baseline (Week 0)                        |
| `sbp_4`            | Systolic BP at Week 4                                    |
| `sbp_8`            | Systolic BP at Week 8                                     |
| `sbp_12`           | Systolic BP at Week 12                                     |
| `dropout_week`     | Week of dropout; `NA` = completed trial                    |
| `device_flag`      | `1` = original device, `2` = new device (Site 2 only)      |
| `deviation_flag`   | `1` = protocol deviation recorded                            |

## Methodology

1. **Data cleaning** — identified genuine missingness in `age`/`bmi` (vs. structural NAs elsewhere); imputed via truncated normal distribution after confirming normality (Shapiro-Wilk).
2. **Outlier detection** — Tukey's IQR method on continuous variables.
3. **Balance checks** — baseline covariates (age, sex, BMI) compared across arms; no significant imbalance found.
4. **Device-switch test** — two-sample t-test comparing pre/post-switch SBP at Site 2; no significant shift.
5. **Dropout mechanism** — logistic regression of dropout on covariates to assess MAR vs MNAR.
6. **Primary inference** — linear regression with robust (HC) standard errors on ΔSBP, adjusting for site, age, sex, BMI, baseline SBP.
7. **Multiplicity-corrected site-level tests** — Bonferroni-adjusted per-site treatment effect tests.
8. **Sensitivity analysis** — re-fit the primary model excluding Site 3 to check robustness to the flagged site.

## Reproducing the Analysis

Requires R (≥ 4.0) with the following packages:

```r
install.packages(c(
  "truncnorm", "tidyverse", "car", "lmtest",
  "naniar", "tableone", "sandwich", "multcomp",
  "ggplot2", "patchwork"
))
```

Then, from the repository root:

```r
source("R/analysis.R")
```

The script expects the raw data at `data/raw/trialx_hackathon_data.csv`.

## Limitations

- Dataset is simulated, designed for analytical practice rather than a real clinical trial.
- Heterogeneity and elevated dropout at Site 3 remain a caveat for interpretation, even though the primary conclusion is robust to its exclusion.

## Citation

Majumdar, P., Dey, S., & Nag, S. (2026). *Analysis of a Contaminated Phase II
RCT for Antihypertensive Efficacy*. Presented at Inferencia 2026, Department
of Statistics, Presidency University, Kolkata.
