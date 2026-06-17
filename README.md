# Code for main analysis for 'Abrupt discontinuation and tapering of antipsychotics in older adults with dementia: a target trial emulation study'

**Author:** Olivia Hopkinson  
**Contact:** olivia.hopkinson.22@ucl.ac.uk 
**GitHub:** @OliviaKHopkinson

## Overview

This repository contains the analysis code for the main analysis for the study of antipsychotic discontinuation strategies among people living with dementia, published in The Lancet Healthy Longevity..

The main analysis uses a clone-censor-weight target trial emulation framework to compare antipsychotic continuation with two discontinuation strategies:

1. Abrupt discontinuation versus continuation
2. Tapering discontinuation versus continuation

The analysis is repeated among people eligible after 12 weeks and 24 weeks of continuous antipsychotic use. Follow-up is structured in 28-day intervals. Inverse probability of censoring weights are used to account for artificial censoring induced by assignment to treatment strategies.

## Main analysis scripts

The commented analysis scripts are:

- `07_tapering_trial_analysis.R`  
  Main analysis code for the tapering discontinuation versus continuation comparison.

- `08_abrupt_discontinuation_trial_analysis.R`  
  Main analysis code for the abrupt discontinuation versus continuation comparison.

These scripts are intended to document the primary analysis workflow used for the study. They include analyses for stroke, pneumonia, fracture, delirium, death, and the negative control outcome where applicable.

## Important data access note

The individual-level CPRD-linked data used in the study cannot be shared publicly because of data governance restrictions. The code is therefore provided to document the analysis workflow and support transparency. Researchers with appropriate data access should adapt the file paths and input datasets to match their approved data environment.

No individual-level patient data are included in this repository.

## Expected input data

The scripts expect long-format analytic datasets with one row per patient per 28-day follow-up interval. The original analysis used two input files:

- `processed_data/combined/combined_12weeks_longformat.txt`
- `processed_data/combined/combined_24weeks_longformat.txt`

Key variables expected by the scripts include:

- `patid`: patient identifier
- `t`: follow-up interval number
- `t_squared`: squared follow-up interval term
- `T0`: trial baseline date
- `startdate` and `enddate`: available follow-up period
- `treatment_startdate` and `treatment_enddate`: qualifying antipsychotic treatment episode dates
- `dose_reduction_date`: date of first prescription with at least a 50% dose reduction compared with baseline, where relevant. Otherwise NA.
- `reinitiation_date`: date of antipsychotic reinitiation, where relevant. Otherwise NA.
- `stroke_outcome_date`, `pneumonia_outcome_date`, `fracture_outcome_date`, `delirium_outcome_date`, `negcontrol_outcome_date`, and `deathdate_adjusted`: first date of outcome after T0 where relevant, otherwise NA
- baseline covariates ending in `_baseline`
- time-varying covariates with the same root names but without the `_baseline` suffix

## Analysis summary

Each analysis script follows the same broad structure:

1. Load packages and input data
2. Convert dates and numeric variables to the required formats
3. Create outcome indicators for each 28-day interval
4. Create continuation and discontinuation/tapering clones
5. Apply treatment-strategy-specific artificial censoring rules
6. Estimate inverse probability of censoring weights
7. Truncate weights at the 99th percentile
8. Fit weighted outcome models
9. Estimate risks and risk differences
10. Run patient-level bootstrap analyses for uncertainty estimation


## Outcome-specific model specification

The scripts use the same overall clone-censor-weight workflow across outcomes, but censoring model fit was optimised separately for each outcome and treatment strategy. Candidate covariates were retained or removed based on model diagnostics, including inspection of standard errors, to avoid unstable estimates from sparse or poorly estimated covariate effects. For that reason, the covariate terms included in the IPCW models are not identical across all outcomes.

## Suggested citation

When using or adapting this code, please cite the associated study and acknowledge Olivia Hopkinson as the author of the analysis code.

Hopkinson OK, Chobanov J, Goordeen D, Man KKC, Wong ICK, Reeve E, Bell JS, Wei L, Howard R, Lau WCY. Abrupt discontinuation and tapering of antipsychotics in older adults with dementia: a target trial emulation study. The Lancet Healthy Longevity. 2026;7(5). doi:10.1016/j.lanhl.2026.100852

A RIS export of the reference is included at `references/hopkinson_2026_lancet_healthy_longevity.ris`. A GitHub citation file is also included as `CITATION.cff`.

## Licence

Add the intended licence here before making the repository public. For academic code, common options include MIT, GPL-3.0, or CC BY 4.0 for documentation. Check institutional or funder requirements before selecting a licence.
