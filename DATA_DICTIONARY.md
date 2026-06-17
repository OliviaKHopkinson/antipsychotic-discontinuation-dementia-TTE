# Data dictionary for main analysis scripts

This document summarises the expected structure of the input datasets used by the main analysis scripts.

## Dataset structure

The analysis datasets should be in long format, with one row per patient per 28-day follow-up interval.

## Core identifiers and time variables

| Variable | Description |
|---|---|
| `patid` | Patient identifier |
| `t` | Follow-up interval number; each interval represents 28 days |
| `t_squared` | Squared follow-up interval term |
| `T0` | Trial baseline date |
| `startdate` | Start of available follow-up |
| `enddate` | End of available follow-up |
| `months_since_diagnosis` | Time from dementia diagnosis to T0 |
| `age_T0` | Age at trial baseline |
| `gender` | Sex/gender variable used in models |

## Treatment and censoring variables

| Variable | Description |
|---|---|
| `treatment_startdate` | Start of qualifying antipsychotic treatment episode |
| `treatment_enddate` | End of qualifying antipsychotic treatment episode |
| `dose_reduction_date` | Date of first prescription following T0 with a dose reduction of at least 50% from baseline, same antipsychotic drug substance required |
| `reinitiation_date` | Date of antipsychotic reinitiation after discontinuation |

## Outcome variables

| Variable | Description |
|---|---|
| `stroke_outcome_date` | Stroke outcome date |
| `pneumonia_outcome_date` | Pneumonia outcome date |
| `fracture_outcome_date` | Fracture outcome date |
| `delirium_outcome_date` | Delirium outcome date |
| `negcontrol_outcome_date` | Negative control outcome date |
| `deathdate_adjusted` | Death date used for mortality analyses, utilises ONS data |

## Covariates

Baseline covariates are denoted with the suffix `_baseline`. Time-varying covariates use the same root variable names without the `_baseline` suffix. All of these are binary variables. Comorbidities are 1 if there is any recording of that comorbidity in primary or secondary care in all history. Codelists are available in this repository.

Medication use covariates are 1 if there was a treatment of that medication in the prior 28 days.

Examples include comorbidities such as `angina_baseline`, `anxiety_baseline`, `diabetes_baseline`, and `kidney_disease_baseline`, and medication variables such as `antidepressants_baseline`, `antiplatelets_baseline`, `dementia_drugs_baseline`, and `oral_anticoagulants_baseline`.

## Data sharing note

The original CPRD-linked patient-level datasets cannot be shared publicly. This dictionary is provided so that researchers with appropriate approvals can understand the required input structure.


## Model covariates

The broad candidate covariate set includes baseline and time-varying clinical conditions and medications. Final censoring model terms may differ between outcomes because model fit was optimised separately for each outcome and treatment strategy, with covariates retained or removed based on diagnostics including inspection of standard errors.


## Associated study

These variables relate to the main analysis code for:

Hopkinson OK, Chobanov J, Goordeen D, Man KKC, Wong ICK, Reeve E, Bell JS, Wei L, Howard R, Lau WCY. Abrupt discontinuation and tapering of antipsychotics in older adults with dementia: a target trial emulation study. The Lancet Healthy Longevity. 2026;7(5). doi:10.1016/j.lanhl.2026.100852
