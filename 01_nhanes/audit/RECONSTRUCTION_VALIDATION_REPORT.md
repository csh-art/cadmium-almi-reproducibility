# NHANES reconstruction validation report

Validation source: exact uploaded `01_raw_xpt.zip` plus archived project audit files.

## Raw-file integrity

- Expected XPT files: **52**
- Uploaded XPT files checked: **52**
- SHA-256 matches to archived copy manifest: **52/52**
- Hash status: **PASS**

## Core analytic sample

| Check | Reconstructed | Archived |
|---|---:|---:|
| 2011–2012 strict N | 2,604 | 2,604 |
| 2013–2014 strict N | 1,523 | 1,523 |
| 2015–2016 strict N | 1,337 | 1,337 |
| 2017–2018 strict N | 2,164 | 2,164 |
| Overall strict N | 7,628 | 7,628 |

Strict limb validity was reproduced using `DXAEXSTS == 1`, all four limb tissue-invalidity codes equal to 0, complete four-limb lean-mass measurements, valid measured height, cadmium, survey weight, strata, and PSU.

## Covariate nonmissing counts

All archived counts were reproduced exactly, including PIR 7,012; smoking 7,624; alcohol 7,170; physical activity 7,588; hypertension 7,621; diabetes 7,622; hyperlipidemia 7,594; eGFR 7,497; cotinine 7,539; and blood lead 7,628.

## Category-level checks

- Alcohol: none/rare 3,315; less than weekly 1,598; weekly or more 2,257; missing 458.
- Physical activity: meets guideline 5,207; below guideline 2,381; missing 40.
- Hypertension: no 5,815; yes 1,806; missing 7.
- Diabetes: no 6,921; borderline 140; yes 561; missing 6.
- Hyperlipidemia: no 5,691; yes 1,903; missing 34.

Hyperlipidemia was reproduced as self-reported high cholesterol (`BPQ080 == 1`) **or current prescribed cholesterol medication (`BPQ100D == 1`)**.

## Model sample-size checks

- Model 1: 7,628
- Model 2: 7,628
- Model 3: 6,575
- Model 4: 6,536
- Sensitivity 1 (eGFR): 6,423
- Sensitivity 2 (cotinine): 6,460
- Sensitivity 3 (lead): 6,536
- Sensitivity 4 (all three): 6,420

## Creatinine calibration check

For NHANES 2017–2018, converting Cobas 6000 creatinine to the earlier DxC 660i scale using the documented backward equation

`calibrated creatinine = 1.051 × LBXSCR − 0.06945`

reproduced the archived means exactly:
- raw creatinine mean: 0.8478699340 mg/dL;
- calibrated mean: 0.8216613007 mg/dL;
- calibrated 2021 CKD-EPI eGFR mean: 105.6730466 mL/min/1.73m².

## Remaining execution check

The container used to assemble this repository does not contain R. Therefore the final `survey` model coefficients were not independently re-executed here. Instead, the repository contains:
1. archived exact aggregate regression outputs;
2. reconstructed R scripts with locked expected-value assertions;
3. verbatim original final Figure 2/table-generation code.

Run `01_nhanes/scripts/08_run_all.R` in the archived R environment (R 4.6.0, `survey` 4.5) before the public GitHub/Zenodo release. The run should reproduce the archived Model 3 beta of −0.1596438517 and other locked results.
