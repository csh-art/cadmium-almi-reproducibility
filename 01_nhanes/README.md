# 01 — NHANES 2011–2018 analysis

## What is included

This module contains:
- a reconstructed raw-XPT-to-analysis pipeline;
- the exact raw-file hashes and variable dictionaries used to reconstruct the early stages;
- archived audit outputs from the original analysis;
- aggregate Table 1, Table 2, and sensitivity source data;
- verbatim copies of the originally archived final table/Figure 2 code.

## Important provenance note

The original project retained excellent audit outputs and frozen objects, but the earliest data-construction commands were not preserved as standalone `.R` scripts. Therefore scripts `01_build_core_dataset.R` through `03_finalize_analysis_dataset.R` are explicitly labeled **reconstructed**. They were reverse-checked against the archived analysis and the exact original XPT inputs.

## Independent reconstruction checks already completed

Using the exact 52 uploaded XPT files:
- all 52 SHA-256 hashes matched the archived manifest;
- strict analytic N reproduced exactly: **7,628**;
- cycle-specific strict N reproduced exactly: **2,604 / 1,523 / 1,337 / 2,164**;
- all covariate nonmissing counts reproduced exactly;
- alcohol categories reproduced exactly: **3,315 / 1,598 / 2,257 / 458 missing**;
- physical-activity categories reproduced exactly: **5,207 / 2,381 / 40 missing**;
- hypertension, diabetes, and hyperlipidemia category counts reproduced exactly;
- Model 1–4 and S1–S4 complete-case sample sizes reproduced exactly;
- 2017–2018 serum creatinine calibration reproduced the archived mean exactly using `1.051 × LBXSCR − 0.06945` before the 2021 race-free CKD-EPI equation.

See `audit/RECONSTRUCTION_VALIDATION_REPORT.md`.

## How to run in R

From the repository root:

```r
Sys.setenv(NHANES_RAW_XPT_DIR = "E:/path/to/01_raw_xpt")
source("01_nhanes/scripts/08_run_all.R", encoding = "UTF-8")
```

Required packages for the reconstructed analytical pipeline:
- `haven`
- `survey`

The archived final DOCX/figure generation scripts use additional packages such as `ggplot2`, `scales`, `officer`, and `flextable`.

## Data not committed to GitHub

Do not commit:
- raw XPT files;
- participant-level analytic RDS files;
- 45–50 MB fitted-model objects.

These are recreated locally. Aggregate source data and audit outputs are committed instead.
