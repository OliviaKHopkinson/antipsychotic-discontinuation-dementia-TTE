################################################################################
# 07_tapering_trial_analysis.R
#
# Main analysis code for study of antipsychotic discontinuation
# strategies among people living with dementia.
#
# Author: Olivia Hopkinson
# Contact: olivia.hopkinson.22@ucl.ac.uk
# GitHub/handle: @OliviaKHopkinson
#
# Analysis: Tapering discontinuation versus continuation
#
# Purpose:
#   This script forms part of the main analysis code for the associated study.
#   It implements a clone-censor-weight target trial emulation using 28-day
#   follow-up intervals, inverse probability of censoring weights, and weighted
#   outcome models.
#
# Data access note:
#   The individual-level CPRD-linked data used for the study cannot be shared
#   publicly because of data governance restrictions. This script is provided to
#   document the analysis workflow and to support transparency and reproducibility
#   for researchers with appropriate data access and equivalently structured data.
#
# Associated study:
#   Hopkinson OK, Chobanov J, Goordeen D, Man KKC, Wong ICK, Reeve E, Bell JS,
#   Wei L, Howard R, Lau WCY. Abrupt discontinuation and tapering of
#   antipsychotics in older adults with dementia: a target trial emulation
#   study. The Lancet Healthy Longevity. 2026;7(5).
#   doi:10.1016/j.lanhl.2026.100852
################################################################################


################################################################################
# 0. Setup ----------------------------------------------------------------------
#
# Load packages, define common covariate lists, and
# read the analysis datasets. The expected input data are in long format with one
# row per patient per 28-day follow-up interval.
################################################################################

library(tidyverse)
library(data.table)
library(doseminer)
library(stringr)
library(tidyr)
library(fuzzyjoin)
library(geepack)
library(cobalt)
library(survey)
library(tableone)
library(lmtest)
library(sandwich)
library(splitstackshape)
library(boot)
library(pbapply)

# helper function to get the model stats for each IPCW model
ipcw_model_counts <- function(df,
                              id = "patid",
                              time = "t",
                              event = "censored") {
  
  df %>%
    ungroup() %>%
    summarise(
      n_patients = n_distinct(.data[[id]]),
      
      # total person-time (person-months)
      person_months = n(),
      
      # number of censoring events
      n_events = sum(.data[[event]] == 1, na.rm = TRUE),
      
      median_fu_months = median(
        ave(.data[[time]], .data[[id]], FUN = length),
        na.rm = TRUE
      ),
      max_fu_months = max(.data[[time]], na.rm = TRUE)
    )
}


set.seed(7)

# turn off scientific notation
options(scipen = 999) 
baseline_covariates <-  c("patid", "angina_baseline", "anxiety_baseline", "cancer_baseline",
                          "copd_baseline", "atrial_fibrillation_baseline", "depression_baseline",
                          "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline",
                          "heart_hypertension_baseline", "kidney_disease_baseline", 
                          "osteoporosis_baseline", "rheumatoid_arthritis_baseline", "venous_thromboembolism_baseline",
                          "ventricular_arrhythmia_baseline", "ace_inhibitors_baseline", "antidepressants_baseline",
                          "antiepileptics_baseline", "antiplatelets_baseline", "anxiolytics_baseline",
                          "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline", "dementia_drugs_baseline",
                          "diabetes_drug_baseline", "diuretics_baseline", "lipid_regulating_drugs_baseline",
                          "oral_anticoagulants_baseline", "mi_baseline")    


################################################################################
# Expected input data structure -------------------------------------------------
#
# The two input files used below contain trial-specific long-format datasets for
# people eligible after 12 or 24 weeks of continuous antipsychotic use.
#
# Key variables assumed by the analysis include:
#   patid                  Patient identifier
#   t                      Follow-up interval number; each interval is 28 days
#   t_squared              Squared follow-up interval term for flexible time
#   T0                     Trial baseline date
#   startdate, enddate     Available follow-up period
#   treatment_startdate    Start of qualifying antipsychotic treatment episode
#   treatment_enddate      End of qualifying antipsychotic treatment episode
#   dose_reduction_date    Date of dose reduction of >= 50%, where relevant
#   reinitiation_date      Date of antipsychotic reinitiation, where relevant
#   *_outcome_date         Outcome-specific event dates
#
# Baseline covariates end in _baseline. Time-varying covariates use the same root
# variable names without the _baseline suffix and are updated over follow-up.
################################################################################

combined_12weeks <- fread("processed_data/combined/combined_12weeks_longformat.txt", colClasses="character")
combined_12weeks <- combined_12weeks %>%
  mutate(t = as.numeric(t),
         diagnosis_date = as.Date(diagnosis_date, format="%Y-%m-%d"),
         age_at_diagnosis = as.numeric(age_at_diagnosis),
         deathdate_adjusted = as.Date(deathdate_adjusted, format="%Y-%m-%d"),
         stroke_outcome_date = as.Date(stroke_outcome_date, format="%Y-%m-%d"),
         pneumonia_outcome_date = as.Date(pneumonia_outcome_date, format="%Y-%m-%d"),
         fracture_outcome_date = as.Date(fracture_outcome_date, format="%Y-%m-%d"),
         delirium_outcome_date = as.Date(delirium_outcome_date, format="%Y-%m-%d"),
         negcontrol_outcome_date = as.Date(negcontrol_outcome_date, format="%Y-%m-%d"),
         treatment_startdate = as.Date(treatment_startdate, format="%Y-%m-%d"),
         treatment_enddate = as.Date(treatment_enddate, format="%Y-%m-%d"),
         T0 = as.Date(T0, format="%Y-%m-%d"),
         startdate = as.Date(startdate, format="%Y-%m-%d"),
         enddate = as.Date(enddate, format="%Y-%m-%d"),
         t_squared = as.numeric(t_squared),
         reinitiation_date = as.Date(reinitiation_date, format="%Y-%m-%d"),
         dose_reduction_date = as.Date(dose_reduction_date, format="%Y-%m-%d"),
         t_squared = as.numeric(t_squared),
         months_since_diagnosis = as.numeric(months_since_diagnosis),
         age_T0 = ceiling(as.numeric(age_T0)))

combined_24weeks <- fread("processed_data/combined/combined_24weeks_longformat.txt", colClasses="character")
combined_24weeks <- combined_24weeks %>%
  mutate(t = as.numeric(t),
         diagnosis_date = as.Date(diagnosis_date, format="%Y-%m-%d"),
         age_at_diagnosis = as.numeric(age_at_diagnosis),
         deathdate_adjusted = as.Date(deathdate_adjusted, format="%Y-%m-%d"),
         stroke_outcome_date = as.Date(stroke_outcome_date, format="%Y-%m-%d"),
         pneumonia_outcome_date = as.Date(pneumonia_outcome_date, format="%Y-%m-%d"),
         fracture_outcome_date = as.Date(fracture_outcome_date, format="%Y-%m-%d"),
         delirium_outcome_date = as.Date(delirium_outcome_date, format="%Y-%m-%d"),
         treatment_startdate = as.Date(treatment_startdate, format="%Y-%m-%d"),
         treatment_enddate = as.Date(treatment_enddate, format="%Y-%m-%d"),
         T0 = as.Date(T0, format="%Y-%m-%d"),
         startdate = as.Date(startdate, format="%Y-%m-%d"),
         enddate = as.Date(enddate, format="%Y-%m-%d"),
         t_squared = as.numeric(t_squared),
         reinitiation_date = as.Date(reinitiation_date, format="%Y-%m-%d"),
         dose_reduction_date = as.Date(dose_reduction_date, format="%Y-%m-%d"),
         t_squared = as.numeric(t_squared),
         months_since_diagnosis = as.numeric(months_since_diagnosis),
         age_T0 = ceiling(as.numeric(age_T0)))


################################################################################
# Analysis workflow applied to each outcome ------------------------------------
#
# Each outcome analysis follows the same broad target trial emulation workflow:
#   1. Create the outcome indicator for the relevant 28-day interval.
#   2. Derive discontinuation, dose reduction, reinitiation, and follow-up dates.
#   3. Create continuation and tapering discontinuation clones.
#   4. Apply treatment-strategy-specific artificial censoring rules.
#   5. Estimate inverse probability of censoring weights.
#   6. Truncate weights at the 99th percentile.
#   7. Fit weighted outcome models.
#   8. Estimate risks, risk differences, and bootstrap confidence intervals.
#
# Model fitting note:
#   The broad candidate covariate set was consistent across analyses, but model
#   fit was optimised separately for each outcome and treatment strategy. Terms
#   were included or removed based on model diagnostics, including standard
#   errors, to avoid unstable model fits while retaining clinically relevant
#   adjustment where possible.
################################################################################

################################################################################
# STROKE OUTCOME
################################################################################

#### stroke outcome analysis
tapering_stroke_trials <- function(long_format, time_period, full_sample=FALSE) {
    
    long_format <- long_format %>%
      mutate(outcome_date = ifelse(stroke_outcome_date > enddate & !is.na(enddate), NA, stroke_outcome_date)) %>%
      mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
      mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
      mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
      mutate(outcome = ifelse(is.na(outcome), 0, outcome))
    
    # add discontinuation date
    long_format <- long_format %>%
      mutate(discontinuation_date = treatment_enddate + 28)
    
    long_format <- long_format %>%
      mutate(discontinuation_date = case_when(
        (!is.na(enddate) & discontinuation_date > enddate) |
          (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
        TRUE ~ discontinuation_date
      ))        
    
    
    # create a variable for maximum followup
    long_format <- long_format %>%
      group_by(patid) %>%
      mutate(maximum_followup = min(stroke_outcome_date, enddate, na.rm=TRUE))
    
    # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
    
    long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                            long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
    
    
    long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                             long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
    
    
    long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                          long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
    
    
    
    long_format <- long_format %>%
      mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
      mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
      mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
    
    # assuming all available followup so take from the overall long format
    long_format <- long_format %>%
      mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
      mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
    
    ################################################################################
    # Continuation arm
    ################################################################################
    # continuation clone are only censored due to discontinuation
    continuation_arm <- long_format %>%
      mutate(treatment = "continuation") %>%
      mutate(censoring_date = discontinuation_date) %>%
      mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
      mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
      mutate(censored = ifelse(is.na(censored), 0, censored))
    
    # end of followup is the earliest of censoring date or outcome date or enddate
    continuation_arm <- continuation_arm %>%
      mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
      mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
      mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
    
    # if outcome date is after the censoring date, set to 0
    continuation_arm <- continuation_arm %>%
      mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                                censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                              0, outcome))
    
    # in cases where censoring date is after the outcome date, set censor month to 0 for that month
    # this is because  followup ends at that point so ignore the censoring event
    continuation_arm <- continuation_arm %>%
      arrange(patid, t) %>%
      mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                                 censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                               0, censored))
    
    # keep rows with followup set to 0
    continuation_arm <- continuation_arm %>%
      mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
      filter(followup == 1)
    
    # create stabilized IPCW weights for the continuation arm
    ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                   # baseline conditions
                                   angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                   cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                   heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                   kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                   venous_thromboembolism_baseline + mi_baseline + 
                                   # time-varying conditions
                                   angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                   diabetes + heartdisease + heart_failure + heart_hypertension +
                                   kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                   # baseline medications
                                   ace_inhibitors_baseline  + antidepressants_baseline +
                                   antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                   arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                   dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                   lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                   # time-varying medications
                                   ace_inhibitors + antidepressants + antiepileptics + 
                                   antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                   dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                   oral_anticoagulants,
                                 data = continuation_arm, family=quasibinomial())
    
    # checked model summary, all have low standard errors
    continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
    
    # get the cumulative probability of remaining uncensored
    continuation_arm <- continuation_arm %>%
      group_by(patid) %>%
      arrange(patid, t) %>%
      mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
      ungroup()
    
    continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
    print("continuation arm unstabilized weights")
    print(summary(continuation_arm$unstabilized_weight))
    
    baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                            # baseline conditions
                                            angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                            cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                            heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                            kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                            venous_thromboembolism_baseline + mi_baseline + 
                                            # baseline medications
                                            ace_inhibitors_baseline + antidepressants_baseline +
                                            antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                            arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                            dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                            lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                          data = continuation_arm, family=quasibinomial())
    
    # checked model summary, all have low standard errors
    continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
    continuation_arm <- continuation_arm %>%
      arrange(patid, t) %>%
      group_by(patid) %>%
      mutate(cum_numerator = cumprod(numerator))
    
    continuation_arm <- continuation_arm %>%
      mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
    
    print("continuation arm stabilized weights")
    print(summary(continuation_arm$stabilized_weight))
    
    # truncate the weights at 99th percentile based on uncensored population
    continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
    upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
    continuation_arm <- continuation_arm %>%
      mutate(stabilized_weight_untruncated = stabilized_weight) %>%
      mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
    
    lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
    continuation_arm <- continuation_arm %>%
      mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
      mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
    
    ################################################################################
    # Tapering arm
    ################################################################################
    # censor at earliest of:
    # 1. reinitiation
    # 2. end of grace period if failure to discontinue
    # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
    tapering_arm <- long_format %>%
      mutate(treatment = "tapering") %>%
      group_by(patid) %>%
      mutate(discontinued = cummax(discontinue)) %>%
      dplyr::select(-discontinue) %>%
      ungroup()
    
    
    # filter out discontinuation dates and dose reduction dates that are not within grace period
    tapering_arm <- tapering_arm %>%
      mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
             dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
    
    # create variables for 'no discontinuation date' and no dose reduction date
    tapering_arm <- tapering_arm %>%
      group_by(patid) %>%
      mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
      fill(no_discontinuation_date, .direction="down") %>%
      fill(no_discontinuation_date, .direction="up") %>%
      ungroup() %>%
      group_by(patid) %>%
      mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                       is.na(dose_reduction_date), discontinuation_date, NA)))
    
    tapering_arm <- tapering_arm %>%
      mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
    
    # followup is earliest of outcome date, enddate, censoring date
    tapering_arm <- tapering_arm %>%
      mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
      mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
      mutate(followup = ifelse(t <= followup_end, 1, 0))
    
    # if outcome date is after the censoring date, set to 0
    tapering_arm <- tapering_arm %>%
      mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
      mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                              & !is.na(outcome_month), 0, outcome))
    
    # create binary variables for each of the reasons for censoring
    tapering_arm <- tapering_arm %>%
      mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
      mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
      mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
      mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
      mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
      mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
      mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
             no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
             reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
    
    # if reinitiation censor is in the same month as another censor, set to 0 as the person
    # would be censored by the other two reasons first
    tapering_arm <- tapering_arm %>%
      mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                            reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
    
    
    # in cases where the outcome occurs in the same month as the censor and the outcome
    # is BEFORE censoring, set the censor variable to 0
    tapering_arm <- tapering_arm %>%
      group_by(patid) %>%
      mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                  !is.na(outcome_month) & !is.na(censoring_month) &
                                                  outcome_date < censoring_date & t == outcome_month,
             0, no_discontinuation_censor)) %>%
      mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                                  !is.na(outcome_month) & !is.na(censoring_month) &
                                                  outcome_date < censoring_date & t == outcome_month,
             0, no_dose_reduction_censor)) %>%
      mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                                  !is.na(outcome_month) & !is.na(censoring_month) &
                                                  outcome_date < censoring_date & t == outcome_month,
             0, reinitiation_censor))
      
    # keep data where followup == 1
    tapering_arm <- tapering_arm %>%
      filter(followup == 1)
    
    # first set of weights, failure to discontinuation, which only applies to t == 7
    grace_period_data <- tapering_arm %>%
      filter(t == 7) %>%
      mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
    
    ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                      # baseline conditions
                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                      venous_thromboembolism_baseline + mi_baseline + 
                                      # time-varying conditions
                                      anxiety + atrial_fibrillation + cancer + depression +
                                      diabetes + heartdisease + heart_failure + 
                                      kidney_disease + osteoporosis + venous_thromboembolism + mi +
                                      # baseline medications
                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                      # time-varying medications
                                      ace_inhibitors + antidepressants + antiepileptics + 
                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                      oral_anticoagulants,
                                    data = grace_period_data, family=quasibinomial())
    
    # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
    grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
    
    # as there is only one time point, no need to do cumultive probabilities
    grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
    print("discontinuation: no discontinuation unstabilized weights")
    print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
    
    baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                                # baseline conditions
                                                angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                venous_thromboembolism_baseline + mi_baseline + 
                                                # baseline medications
                                                ace_inhibitors_baseline + antidepressants_baseline +
                                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                              data = grace_period_data, family=quasibinomial())
    
    grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
    grace_period_data <- grace_period_data %>%
      mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
    print("discontinuation: no discontinuation stabilized weights")
    print(summary(grace_period_data$stabilized_weight_no_discontinuation))
    
    grace_period_data <- grace_period_data %>% 
      dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
    tapering_arm <- tapering_arm %>%
      left_join(grace_period_data, by=c("patid", "t"))
    
    ## no dose reduction weights which should only be applied to those who have not discontinued
    # those who have discontinued without being censored cannot be censored due to no dose reduction
    # must also only occur during grace period
    not_discontinued_data <- tapering_arm %>%
      filter(t < 7) %>%
      filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
    
    not_discontinued_data <- not_discontinued_data %>%
      mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
    
    
    ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                        # baseline conditions
                                                        angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                        cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                        heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                        kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                        venous_thromboembolism_baseline + mi_baseline + 
                                                        # time-varying conditions
                                                        angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                        diabetes + heartdisease + heart_failure + heart_hypertension +
                                                        kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                        # baseline medications
                                                        ace_inhibitors_baseline  + antidepressants_baseline +
                                                        antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                        arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                        dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                        lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                        # time-varying medications
                                                        ace_inhibitors + antidepressants + antiepileptics + 
                                                        antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                        dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                        oral_anticoagulants,
                                 data = not_discontinued_data, family=quasibinomial())
    
    # all covariates have reasonable standard errors apart from rheumatoid arthritis
    not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
    
    # need to calculate the cumulative probabilities of the numerator and denominator
    # get the cumulative probability of remaining uncensored
    not_discontinued_data <- not_discontinued_data %>%
      group_by(patid) %>%
      arrange(patid, t) %>%
      mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
      ungroup()
    
    not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
    print("discontinuation no dose reduction unstabilized weights")
    print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
    
    baseline_ps_model_discontinuation_no_dose_reduction <-  
      glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
            # baseline conditions
            angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
            cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
            heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
            kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
            venous_thromboembolism_baseline + mi_baseline + 
            # baseline medications
            ace_inhibitors_baseline + antidepressants_baseline +
            antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
            arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
            dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
            lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
          data = not_discontinued_data, family=quasibinomial())
     
    
    not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
    
    # get the cumulative probabilities for the numerator
    not_discontinued_data <- not_discontinued_data %>%
      group_by(patid) %>%
      arrange(patid, t) %>%
      mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
      ungroup()
    
    
    not_discontinued_data <- not_discontinued_data %>%
      mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
    print("discontinuation: no dose reduction stabilized weights")
    print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
    
    not_discontinued_data <- not_discontinued_data %>%
      dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
    
    tapering_arm <- tapering_arm %>%
      left_join(not_discontinued_data, by=c("patid", "t"))
    
    # reinitiation weights, which are only applied to those who have discontinued  
      reinitiation_weight_data <- tapering_arm %>%
      filter(discontinued == 1)
    
    reinitiation_weight_data <- reinitiation_weight_data %>%
      mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
    
    reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                      # baseline conditions
                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                      venous_thromboembolism_baseline + mi_baseline + 
                                      # time-varying conditions
                                      anxiety + atrial_fibrillation + depression +
                                      diabetes + kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                      # baseline medications
                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                      # time-varying medications
                                      ace_inhibitors + antidepressants + antiepileptics + 
                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                      oral_anticoagulants,
                                    data = reinitiation_weight_data, family=quasibinomial())
    
    # very high standard errors for angina, cancer, copd, heart disease, heart failure, hypertension
    reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
    
    # get the cumulative probabilities
    reinitiation_weight_data <- reinitiation_weight_data %>%
      group_by(patid) %>%
      arrange(patid, t) %>%
      mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
      ungroup()
    
    reinitiation_weight_data <- reinitiation_weight_data %>%
      mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
    print("reinitiation weights unstabilized")
    print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
    
    reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # baseline medications
                                    ace_inhibitors_baseline + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                  data = reinitiation_weight_data, family=quasibinomial())
    
    reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
    # get the cumulative probabilities
    reinitiation_weight_data <- reinitiation_weight_data %>%
      group_by(patid) %>%
      arrange(patid, t) %>%
      mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
      ungroup()
    
    reinitiation_weight_data <- reinitiation_weight_data %>%
      mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
    print("reinitiation weights stabilized")
    print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
    
    reinitiation_weight_data <- reinitiation_weight_data %>%
      dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
    
    tapering_arm <- tapering_arm %>%
      left_join(reinitiation_weight_data, by=c("patid", "t"))
    
    # fill each of the weights downwards
    tapering_arm <- tapering_arm %>%
      group_by(patid) %>%
      arrange(patid, t) %>%
      fill(
        stabilized_weight_no_discontinuation,
        unstabilized_weight_no_discontinuation,
        stabilized_weight_no_dose_reduction,
        unstabilized_weight_no_dose_reduction,
        stabilized_weight_reinitiation,
        unstabilized_weight_reinitiation,
        .direction = "down"
      ) %>%
      ungroup()
    
    # fill in the rest in with 1s
    tapering_arm <- tapering_arm %>%
      mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                           stabilized_weight_no_discontinuation),
             unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                             unstabilized_weight_no_discontinuation),
             stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                          stabilized_weight_no_dose_reduction),
             unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                            unstabilized_weight_no_dose_reduction),
             stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                     stabilized_weight_reinitiation),
             unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                       unstabilized_weight_reinitiation))
    
    # multiple all of the weights together
    tapering_arm <- tapering_arm %>%
      mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
               stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
             unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
               unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
    
    # truncate the stabilized weights at the 99th percentile
    quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
    tapering_arm <- tapering_arm %>%
      mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
    print("summary of truncated stabilized weights for tapering arm:")
    print(summary(tapering_arm$stabilized_weight))
    
    # truncate the unstabilized weights at the 99th percentile
    upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
    
    tapering_arm <- tapering_arm %>%
      mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
    print("summary of truncated unstabilized weights for tapering arm:")
    print(summary(tapering_arm$unstabilized_weight))
    
    tapering_arm <- tapering_arm %>%
      mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                                 1, 0)) %>%
      mutate(censored = ifelse(uncensored == 0, 1, 0))   
    
    vars <- c("stabilized_weight", "stabilized_weight_untruncated")
    
    # function to compute the summary stats
    get_summary_stats <- function(x) {
      c(
        mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        p1 = quantile(x, 0.01, na.rm = TRUE),
        p50 = quantile(x, 0.50, na.rm = TRUE),
        p99 = quantile(x, 0.99, na.rm = TRUE),
        max = max(x, na.rm = TRUE)
      )
    }
    
    summary_tapering <- do.call(rbind, lapply(vars, function(v) {
      stats <- get_summary_stats(tapering_arm[[v]])
      data.frame(variable = v, t(stats))
    }))
    summary_tapering <- summary_tapering %>%
      mutate(arm = "tapering")
    
    summary_continuation <- do.call(rbind, lapply(vars, function(v) {
      stats <- get_summary_stats(continuation_arm[[v]])
      data.frame(variable = v, t(stats))
    }))
    summary_continuation <- summary_continuation %>%
      mutate(arm = "continuation")
    
    summary_weights <- rbind(summary_continuation, summary_tapering)
    
    fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/stroke.csv"))
    ################################################################################
    # check covariate balance after grace period
    ################################################################################
    covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                    "cancer_baseline", "copd_baseline", "depression_baseline",
                    "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                    "heart_hypertension_baseline", "kidney_disease_baseline",
                    "ace_inhibitors_baseline", "antidepressants_baseline", 
                    "antiplatelets_baseline", "anxiolytics_baseline",
                    "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                    "dementia_drugs_baseline", "diabetes_drug_baseline",
                    "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                    "oral_anticoagulants_baseline", "osteoporosis_baseline",
                    "venous_thromboembolism_baseline", 
                    "rheumatoid_arthritis_baseline", "mi_baseline",
                    "months_since_diagnosis", "age_T0", "gender")
    
    
    continuation_t7 <- continuation_arm %>%
      filter(t == 7 & uncensored == 1) %>%
      dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
    
    tapering_t7 <- tapering_arm %>%
      filter(t == 7 & uncensored == 1) %>%
      dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
    
    
    combined_t7 <- rbind(continuation_t7, tapering_t7)
    
    check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
    
    # Calculate the balance table for the subset
    balance_table <- svyCreateTableOne(vars = covariates, 
                                       strata = "treatment", 
                                       data = check_weighting, 
                                       test = FALSE, 
                                       smd = TRUE)
    
    table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
    
    table_one_df$variable <- rownames(table_one_df)
    
    # balance table before weighting
    pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
    
    # Calculate the balance table for the subset
    pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                  strata = "treatment", 
                                                  data = pre_weighting_smd, 
                                                  test = FALSE, 
                                                  smd = TRUE)
    
    
    pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
    pre_weighting_table$variable <- rownames(pre_weighting_table)
    
    
    if(full_sample == TRUE){
      fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/stroke.csv"))
      fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_stroke.csv"))
    }
    
    ################################################################################
    # Outcome model
    ################################################################################
    
    continuation_df <- continuation_arm %>%
      dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                    age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
    
    tapering_df <- tapering_arm %>%
      dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                    age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
    
    combined_data <- rbind(continuation_df, tapering_df) %>%
      filter(censored == 0)
    
    
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/stroke.csv" ))
  
    
    event_by_t <- combined_data %>%
      group_by(t) %>%
      summarise(
        n_at_risk = n_distinct(patid),
        n_events  = sum(outcome == 1, na.rm = TRUE),
        .groups = "drop"
      )
    
    fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/stroke.csv" ))
    
    
    # patient years by arm
    py_by_arm <- combined_data %>%
      group_by(treatment) %>%
      summarise(
        patient_years = ceiling(sum(stabilized_weight * (1/12))),
        .groups = "drop"
      )
    print(py_by_arm)
    
    # number of outcomes
    weighted_events <- with(
      combined_data,
      tapply(stabilized_weight * outcome, treatment, sum)
    )
    print(weighted_events)
    
    crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
    # to create a survival curve, create a model with an interaction between t and treatment
    # and t_squared and treatment
    survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                  + t*treatment + t_squared*treatment + 
                                  + gender + age_T0 + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                  atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline +
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
   
   continuation_df <- long_format %>%
     filter(followup_curve == 1) %>%
     mutate(treatment = "continuation")
   
   tapering_df <- long_format %>%
     filter(followup_curve == 1) %>%
     mutate(treatment = "tapering")
   
   # create hazards under each strategy
   tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
   continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
   
   # calculate survival from the cumulative of 1-hazard
   tapering_df <- tapering_df %>%
     group_by(patid) %>%
     mutate(survival1 = cumprod(1-hazard1)) %>%
     ungroup()
   
   continuation_df <- continuation_df %>%
     group_by(patid) %>%
     mutate(survival0 = cumprod(1-hazard0)) %>%
     ungroup()
   
   # estimate risks as 1-survival
   tapering_df <- tapering_df %>%
     mutate(risk1 = 1-survival1)
   
   continuation_df <- continuation_df %>%
     mutate(risk0 = 1-survival0)
   
   # calculate mean for each time point
   tapering_risk <- tapering_df %>%
     group_by(t, treatment) %>%
     summarise(risk = mean(risk1),
               survival = mean(survival1),
               hazard = mean(hazard1))
   
   continuation_risk <- continuation_df %>%
     group_by(t, treatment) %>%
     summarise(risk = mean(risk0),
               survival = mean(survival0),
               hazard = mean(hazard0))
   
   risk_curves <- full_join(tapering_risk, continuation_risk)
   
   risk_curves_wider <- risk_curves %>%
     pivot_wider(
       id_cols = t,
       names_from = treatment,
       values_from = c(risk, hazard, survival),
       names_sep = "_"
     )
   
   risk_curves_wider <- risk_curves_wider %>%
     mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
   
   if(full_sample == TRUE){
     fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/stroke.csv"))
   }
   
   # hazard ratio is average of hazard ratios over time
   hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
   
   tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
   continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
   
   tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
   continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
   
   # calculate risk ratios and differences
   risk_ratio_12 <- tapering_risk_12/continuation_risk_12
   risk_diff_12 <- tapering_risk_12-continuation_risk_12
   
   risk_ratio_24 <- tapering_risk_24/continuation_risk_24
   risk_diff_24 <- tapering_risk_24-continuation_risk_24

   risk_df <- data.frame(
     estimate_hr = hazard_ratio_estimate,
     risk_ratio_12months = risk_ratio_12,
     risk_ratio_24months = risk_ratio_24,
     risk_diff_12months = risk_diff_12,
     risk_diff_24months = risk_diff_24,
     abs_risk_continuation_12 = continuation_risk_12,
     abs_risk_continuation_24 = continuation_risk_24,
     abs_risk_tapering_12 = tapering_risk_12, 
     abs_risk_tapering_24 = tapering_risk_24
   )
   rownames(risk_df) <- NULL
   
   if(full_sample == TRUE){
     fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/stroke.csv"))
   }
   
   risk_curves <- risk_curves %>%
     mutate(risk=risk*100)
   
   risk_curves <- rbind(
     risk_curves,
     data.frame(t = 0, risk = 0, treatment = "continuation"),
     data.frame(t = 0, risk = 0, treatment = "tapering")
   )
   
   # Plot the survival curves
   survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
     geom_line(linewidth=1)+
     labs(
       title = "Stroke Cumulative Incidence Curves by Treatment",
       x = "Time (months)",
       y = "Cumulative Incidence (%)",
       color = "Treatment"
     ) +
     theme_minimal(base_size = 14) +
     theme(
       axis.text = element_text(size = 14),     # tick label size
       axis.title = element_text(size = 16),    # axis title size
       axis.ticks = element_line(linewidth = 0.8),
       axis.ticks.length = unit(0.25, "cm"),
       legend.title = element_text(size = 14),
       legend.text = element_text(size = 13),
       plot.title = element_text(size = 16, face = "bold")
     ) +
     # vertical line at end of grace period
     geom_vline(xintercept = 6, 
                linetype = "dashed", 
                linewidth = 1,
                color = "black") +
     ylim(0, 32)
   
   if(full_sample == TRUE){
     ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/stroke.png"), plot = survival_curve_plot,
                   dpi=300, width = 8, height = 6)
   }
   
   return(risk_df)
   
}


################################################################################
# Main full-sample analyses -----------------------------------------------------
#
# The calls below run each outcome analysis using the full analytic dataset. These
# are the primary analyses reported for the study. Bootstrap analyses are run in
# the following section to estimate uncertainty around risk differences.
################################################################################

tapering_stroke_trials(combined_12weeks, '12weeks', full_sample = TRUE)
tapering_stroke_trials(combined_24weeks, '24weeks', full_sample = TRUE)


################################################################################
# PNEUMONIA OUTCOME
################################################################################

#### pneumonia outcome analysis


tapering_pneumonia_trials <- function(long_format, time_period, full_sample=FALSE) {
  
  long_format <- long_format %>%
    mutate(outcome_date = ifelse(pneumonia_outcome_date > enddate & !is.na(enddate), NA, pneumonia_outcome_date)) %>%
    mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
    mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
    mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
    mutate(outcome = ifelse(is.na(outcome), 0, outcome))
  
  # add discontinuation date
  long_format <- long_format %>%
    mutate(discontinuation_date = treatment_enddate + 28)
  
  long_format <- long_format %>%
    mutate(discontinuation_date = case_when(
      (!is.na(enddate) & discontinuation_date > enddate) |
        (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
      TRUE ~ discontinuation_date
    ))
  
  
  # create a variable for maximum followup
  long_format <- long_format %>%
    group_by(patid) %>%
    mutate(maximum_followup = min(pneumonia_outcome_date, enddate, na.rm=TRUE))
  
  # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
  
  long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                    long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                     long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                  long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
  
  
  
  long_format <- long_format %>%
    mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
    mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
    mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
  
  # assuming all available followup so take from the overall long format
  long_format <- long_format %>%
    mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
    mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
  
  ################################################################################
  # Continuation arm
  ################################################################################
  # continuation clone are only censored due to discontinuation
  continuation_arm <- long_format %>%
    mutate(treatment = "continuation") %>%
    mutate(censoring_date = discontinuation_date) %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
    mutate(censored = ifelse(is.na(censored), 0, censored))
  
  # end of followup is the earliest of censoring date or outcome date or enddate
  continuation_arm <- continuation_arm %>%
    mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
    mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
    mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  continuation_arm <- continuation_arm %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                              censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                            0, outcome))
  
  # in cases where censoring date is after the outcome date, set censor month to 0 for that month
  # this is because  followup ends at that point so ignore the censoring event
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                               censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                             0, censored))
  
  # keep rows with followup set to 0
  continuation_arm <- continuation_arm %>%
    mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
    filter(followup == 1)
  
  # create stabilized IPCW weights for the continuation arm
  ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                 # baseline conditions
                                 angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                 cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                 heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                 kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                 venous_thromboembolism_baseline + mi_baseline + 
                                 # time-varying conditions
                                 angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                 diabetes + heartdisease + heart_failure + heart_hypertension +
                                 kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                 # baseline medications
                                 ace_inhibitors_baseline  + antidepressants_baseline +
                                 antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                 arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                 dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                 lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                 # time-varying medications
                                 ace_inhibitors + antidepressants + antiepileptics + 
                                 antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                 dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                 oral_anticoagulants,
                               data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  
  continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
  
  # get the cumulative probability of remaining uncensored
  continuation_arm <- continuation_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
  print("continuation arm unstabilized weights")
  print(summary(continuation_arm$unstabilized_weight))
  
  baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                          # baseline conditions
                                          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                          venous_thromboembolism_baseline + mi_baseline + 
                                          # baseline medications
                                          ace_inhibitors_baseline + antidepressants_baseline +
                                          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                        data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    group_by(patid) %>%
    mutate(cum_numerator = cumprod(numerator))
  
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
  
  print("continuation arm stabilized weights")
  print(summary(continuation_arm$stabilized_weight))
  
  # truncate the weights at 99th percentile based on uncensored population
  continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
  upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight) %>%
    mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
  
  lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
  
  ################################################################################
  # Tapering arm
  ################################################################################
  # censor at earliest of:
  # 1. reinitiation
  # 2. end of grace period if failure to discontinue
  # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
  tapering_arm <- long_format %>%
    mutate(treatment = "tapering") %>%
    group_by(patid) %>%
    mutate(discontinued = cummax(discontinue)) %>%
    dplyr::select(-discontinue) %>%
    ungroup()
  
  
  # filter out discontinuation dates and dose reduction dates that are not within grace period
  tapering_arm <- tapering_arm %>%
    mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
           dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
  
  # create variables for 'no discontinuation date' and no dose reduction date
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
    fill(no_discontinuation_date, .direction="down") %>%
    fill(no_discontinuation_date, .direction="up") %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                     is.na(dose_reduction_date), discontinuation_date, NA)))
  
  tapering_arm <- tapering_arm %>%
    mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
  
  # followup is earliest of outcome date, enddate, censoring date
  tapering_arm <- tapering_arm %>%
    mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
    mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
    mutate(followup = ifelse(t <= followup_end, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  tapering_arm <- tapering_arm %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                            & !is.na(outcome_month), 0, outcome))
  
  # create binary variables for each of the reasons for censoring
  tapering_arm <- tapering_arm %>%
    mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
    mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
    mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
    mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
    mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
    mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
    mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
           no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
           reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
  
  # if reinitiation censor is in the same month as another censor, set to 0 as the person
  # would be censored by the other two reasons first
  tapering_arm <- tapering_arm %>%
    mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                          reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
  
  
  # in cases where the outcome occurs in the same month as the censor and the outcome
  # is BEFORE censoring, set the censor variable to 0
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                !is.na(outcome_month) & !is.na(censoring_month) &
                                                outcome_date < censoring_date & t == outcome_month,
                                              0, no_discontinuation_censor)) %>%
    mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                               !is.na(outcome_month) & !is.na(censoring_month) &
                                               outcome_date < censoring_date & t == outcome_month,
                                             0, no_dose_reduction_censor)) %>%
    mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                          !is.na(outcome_month) & !is.na(censoring_month) &
                                          outcome_date < censoring_date & t == outcome_month,
                                        0, reinitiation_censor))
  
  # keep data where followup == 1
  tapering_arm <- tapering_arm %>%
    filter(followup == 1)
  
  # first set of weights, failure to discontinuation, which only applies to t == 7
  grace_period_data <- tapering_arm %>%
    filter(t == 7) %>%
    mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
  
  ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + cancer + depression +
                                    diabetes + heartdisease + heart_failure + 
                                    kidney_disease + osteoporosis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = grace_period_data, family=quasibinomial())
  
  # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
  
  grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
  
  # as there is only one time point, no need to do cumultive probabilities
  grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
  print("discontinuation: no discontinuation unstabilized weights")
  print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
  
  baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                              # baseline conditions
                                              angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                              cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                              heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                              kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                              venous_thromboembolism_baseline + mi_baseline + 
                                              # baseline medications
                                              ace_inhibitors_baseline + antidepressants_baseline +
                                              antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                              arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                              dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                              lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                            data = grace_period_data, family=quasibinomial())
  
  grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
  grace_period_data <- grace_period_data %>%
    mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
  print("discontinuation: no discontinuation stabilized weights")
  print(summary(grace_period_data$stabilized_weight_no_discontinuation))
  
  grace_period_data <- grace_period_data %>% 
    dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
  tapering_arm <- tapering_arm %>%
    left_join(grace_period_data, by=c("patid", "t"))
  
  ## no dose reduction weights which should only be applied to those who have not discontinued
  # those who have discontinued without being censored cannot be censored due to no dose reduction
  # must also only occur during grace period
  not_discontinued_data <- tapering_arm %>%
    filter(t < 7) %>%
    filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
  
  
  ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                      # baseline conditions
                                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                      venous_thromboembolism_baseline + mi_baseline + 
                                                      # time-varying conditions
                                                      angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                      diabetes + heartdisease + heart_failure + heart_hypertension +
                                                      kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                      # baseline medications
                                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                      # time-varying medications
                                                      ace_inhibitors + antidepressants + antiepileptics + 
                                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                      oral_anticoagulants,
                                                    data = not_discontinued_data, family=quasibinomial())
  # all covariates have reasonable standard errors apart from rheumatoid arthritis
  
  not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
  
  # need to calculate the cumulative probabilities of the numerator and denominator
  # get the cumulative probability of remaining uncensored
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
  print("discontinuation no dose reduction unstabilized weights")
  print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
  
  baseline_ps_model_discontinuation_no_dose_reduction <-  
    glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
          # baseline conditions
          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
          venous_thromboembolism_baseline + mi_baseline + 
          # baseline medications
          ace_inhibitors_baseline + antidepressants_baseline +
          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
        data = not_discontinued_data, family=quasibinomial())
  
  
  not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
  
  # get the cumulative probabilities for the numerator
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
    ungroup()
  
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
  print("discontinuation: no dose reduction stabilized weights")
  print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
  
  not_discontinued_data <- not_discontinued_data %>%
    dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
  
  tapering_arm <- tapering_arm %>%
    left_join(not_discontinued_data, by=c("patid", "t"))
  
  # reinitiation weights, which are only applied to those who have discontinued  
  reinitiation_weight_data <- tapering_arm %>%
    filter(discontinued == 1)
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
  
  reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + depression +
                                    diabetes + kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = reinitiation_weight_data, family=quasibinomial())
  
  # very high standard errors for angina, cancer, copd, heart disease, heart failure, hypertension
  reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
  
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
  print("reinitiation weights unstabilized")
  print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
  
  reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                  cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                  kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline + 
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
  print("reinitiation weights stabilized")
  print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
  
  tapering_arm <- tapering_arm %>%
    left_join(reinitiation_weight_data, by=c("patid", "t"))
  
  # fill each of the weights downwards
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    fill(
      stabilized_weight_no_discontinuation,
      unstabilized_weight_no_discontinuation,
      stabilized_weight_no_dose_reduction,
      unstabilized_weight_no_dose_reduction,
      stabilized_weight_reinitiation,
      unstabilized_weight_reinitiation,
      .direction = "down"
    ) %>%
    ungroup()
  
  # fill in the rest in with 1s
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                         stabilized_weight_no_discontinuation),
           unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                           unstabilized_weight_no_discontinuation),
           stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                        stabilized_weight_no_dose_reduction),
           unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                          unstabilized_weight_no_dose_reduction),
           stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                   stabilized_weight_reinitiation),
           unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                     unstabilized_weight_reinitiation))
  
  # multiple all of the weights together
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
             stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
           unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
             unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
  
  # truncate the stabilized weights at the 99th percentile
  quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
  print("summary of truncated stabilized weights for tapering arm:")
  print(summary(tapering_arm$stabilized_weight))
  
  # truncate the unstabilized weights at the 99th percentile
  upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
  
  tapering_arm <- tapering_arm %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
  print("summary of truncated unstabilized weights for tapering arm:")
  print(summary(tapering_arm$unstabilized_weight))
  
  tapering_arm <- tapering_arm %>%
    mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                               1, 0)) %>%
    mutate(censored = ifelse(uncensored == 0, 1, 0))   
  
  vars <- c("stabilized_weight", "stabilized_weight_untruncated")
  
  # Function to compute the summary stats
  get_summary_stats <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p1 = quantile(x, 0.01, na.rm = TRUE),
      p50 = quantile(x, 0.50, na.rm = TRUE),
      p99 = quantile(x, 0.99, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  summary_tapering <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(tapering_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_tapering <- summary_tapering %>%
    mutate(arm = "tapering")
  
  summary_continuation <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(continuation_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_continuation <- summary_continuation %>%
    mutate(arm = "continuation")
  
  summary_weights <- rbind(summary_continuation, summary_tapering)
  
  fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/pneumonia.csv"))
  ################################################################################
  # check covariate balance after grace period
  ################################################################################
  covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                  "cancer_baseline", "copd_baseline", "depression_baseline",
                  "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                  "heart_hypertension_baseline", "kidney_disease_baseline",
                  "ace_inhibitors_baseline", "antidepressants_baseline", 
                  "antiplatelets_baseline", "anxiolytics_baseline",
                  "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                  "dementia_drugs_baseline", "diabetes_drug_baseline",
                  "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                  "oral_anticoagulants_baseline", "osteoporosis_baseline",
                  "venous_thromboembolism_baseline", 
                  "rheumatoid_arthritis_baseline", "mi_baseline",
                  "months_since_diagnosis", "age_T0", "gender")
  
  
  continuation_t7 <- continuation_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  tapering_t7 <- tapering_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  
  combined_t7 <- rbind(continuation_t7, tapering_t7)
  
  check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
  
  # Calculate the balance table for the subset
  balance_table <- svyCreateTableOne(vars = covariates, 
                                     strata = "treatment", 
                                     data = check_weighting, 
                                     test = FALSE, 
                                     smd = TRUE)
  
  table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  
  table_one_df$variable <- rownames(table_one_df)
  
  # balance table before weighting
  pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
  
  # Calculate the balance table for the subset
  pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                strata = "treatment", 
                                                data = pre_weighting_smd, 
                                                test = FALSE, 
                                                smd = TRUE)
  
  
  pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  pre_weighting_table$variable <- rownames(pre_weighting_table)
  
  
  if(full_sample == TRUE){
    fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pneumonia.csv"))
    fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_pneumonia.csv"))
  }
  
  ################################################################################
  # Outcome model
  ################################################################################
  
  continuation_df <- continuation_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  tapering_df <- tapering_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  combined_data <- rbind(continuation_df, tapering_df) %>%
    filter(censored == 0)
  
  
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/pneumonia.csv" ))
  
  
  # patient years by arm
  py_by_arm <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      # each row = 1 month, so 1/12 year
      patient_years = ceiling(sum(stabilized_weight * (1/12))),
      .groups = "drop"
    )
  print(py_by_arm)
  
  # number of outcomes
  weighted_events <- with(
    combined_data,
    tapply(stabilized_weight * outcome, treatment, sum)
  )
  print(weighted_events)
  
  
  crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
  # to create a survival curve, create a model with an interaction between t and treatment
  # and t_squared and treatment
  survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                + t*treatment + t_squared*treatment + 
                                + gender + age_T0 + months_since_diagnosis +
                                # baseline conditions
                                angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                venous_thromboembolism_baseline + mi_baseline +
                                # baseline medications
                                ace_inhibitors_baseline + antidepressants_baseline +
                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                              data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
  
  continuation_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "continuation")
  
  tapering_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "tapering")
  
  # create hazards under each strategy
  tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
  continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
  
  # calculate survival from the cumulative of 1-hazard
  tapering_df <- tapering_df %>%
    group_by(patid) %>%
    mutate(survival1 = cumprod(1-hazard1)) %>%
    ungroup()
  
  continuation_df <- continuation_df %>%
    group_by(patid) %>%
    mutate(survival0 = cumprod(1-hazard0)) %>%
    ungroup()
  
  # estimate risks as 1-survival
  tapering_df <- tapering_df %>%
    mutate(risk1 = 1-survival1)
  
  continuation_df <- continuation_df %>%
    mutate(risk0 = 1-survival0)
  
  # calculate mean for each time point
  tapering_risk <- tapering_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk1),
              survival = mean(survival1),
              hazard = mean(hazard1))
  
  continuation_risk <- continuation_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk0),
              survival = mean(survival0),
              hazard = mean(hazard0))
  
  risk_curves <- full_join(tapering_risk, continuation_risk)
  
  risk_curves_wider <- risk_curves %>%
    pivot_wider(
      id_cols = t,
      names_from = treatment,
      values_from = c(risk, hazard, survival),
      names_sep = "_"
    )
  
  risk_curves_wider <- risk_curves_wider %>%
    mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
  
  if(full_sample == TRUE){
    fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/pneumonia.csv"))
  }
  
  # hazard ratio is average of hazard ratios over time
  hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
  
  tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
  continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
  
  tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
  continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
  
  # calculate risk ratios and differences
  risk_ratio_12 <- tapering_risk_12/continuation_risk_12
  risk_diff_12 <- tapering_risk_12-continuation_risk_12
  
  risk_ratio_24 <- tapering_risk_24/continuation_risk_24
  risk_diff_24 <- tapering_risk_24-continuation_risk_24
  
  
  risk_df <- data.frame(
    estimate_hr = hazard_ratio_estimate,
    risk_ratio_12months = risk_ratio_12,
    risk_ratio_24months = risk_ratio_24,
    risk_diff_12months = risk_diff_12,
    risk_diff_24months = risk_diff_24,
    abs_risk_continuation_12 = continuation_risk_12,
    abs_risk_continuation_24 = continuation_risk_24,
    abs_risk_tapering_12 = tapering_risk_12, 
    abs_risk_tapering_24 = tapering_risk_24
  )
  
  rownames(risk_df) <- NULL
  
  if(full_sample == TRUE){
    fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/pneumonia.csv"))
  }
  
  risk_curves <- risk_curves %>%
    mutate(risk=risk*100)
  
  risk_curves <- rbind(
    risk_curves,
    data.frame(t = 0, risk = 0, treatment = "continuation"),
    data.frame(t = 0, risk = 0, treatment = "tapering")
  )
  
  # Plot the survival curves
  survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
    geom_line(linewidth=1)+
    labs(
      title = "Pneumonia Cumulative Incidence Curves by Treatment",
      x = "Time (months)",
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),     # tick label size
      axis.title = element_text(size = 16),    # axis title size
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    # vertical line at end of grace period
    geom_vline(xintercept = 6, 
               linetype = "dashed", 
               linewidth = 1,
               color = "black") 
  
  if(full_sample == TRUE){
    ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/pneumonia.png"), plot = survival_curve_plot,
           dpi=300, width = 8, height = 6)
  }
  
  return(risk_df)
  
}


tapering_pneumonia_trials_24weeks <- function(long_format, time_period, full_sample=FALSE) {
  
  long_format <- long_format %>%
    mutate(outcome_date = ifelse(pneumonia_outcome_date > enddate & !is.na(enddate), NA, pneumonia_outcome_date)) %>%
    mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
    mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
    mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
    mutate(outcome = ifelse(is.na(outcome), 0, outcome))
  
  # add discontinuation date
  long_format <- long_format %>%
    mutate(discontinuation_date = treatment_enddate + 28)
  
  long_format <- long_format %>%
    mutate(discontinuation_date = case_when(
      (!is.na(enddate) & discontinuation_date > enddate) |
        (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
      TRUE ~ discontinuation_date
    ))
  
  
  # create a variable for maximum followup
  long_format <- long_format %>%
    group_by(patid) %>%
    mutate(maximum_followup = min(pneumonia_outcome_date, enddate, na.rm=TRUE))
  
  # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
  
  long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                    long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                     long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                  long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
  
  
  
  long_format <- long_format %>%
    mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
    mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
    mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
  
  # assuming all available followup so take from the overall long format
  long_format <- long_format %>%
    mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
    mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
  
  ################################################################################
  # Continuation arm
  ################################################################################
  # continuation clone are only censored due to discontinuation
  continuation_arm <- long_format %>%
    mutate(treatment = "continuation") %>%
    mutate(censoring_date = discontinuation_date) %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
    mutate(censored = ifelse(is.na(censored), 0, censored))
  
  # end of followup is the earliest of censoring date or outcome date or enddate
  continuation_arm <- continuation_arm %>%
    mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
    mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
    mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  continuation_arm <- continuation_arm %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                              censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                            0, outcome))
  
  # in cases where censoring date is after the outcome date, set censor month to 0 for that month
  # this is because  followup ends at that point so ignore the censoring event
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                               censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                             0, censored))
  
  # keep rows with followup set to 0
  continuation_arm <- continuation_arm %>%
    mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
    filter(followup == 1)
  
  # create stabilized IPCW weights for the continuation arm
  ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                 # baseline conditions
                                 angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                 cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                 heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                 kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                 venous_thromboembolism_baseline + mi_baseline + 
                                 # time-varying conditions
                                 angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                 diabetes + heartdisease + heart_failure + heart_hypertension +
                                 kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                 # baseline medications
                                 ace_inhibitors_baseline  + antidepressants_baseline +
                                 antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                 arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                 dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                 lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                 # time-varying medications
                                 ace_inhibitors + antidepressants + antiepileptics + 
                                 antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                 dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                 oral_anticoagulants,
                               data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  
  continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
  
  # get the cumulative probability of remaining uncensored
  continuation_arm <- continuation_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
  print("continuation arm unstabilized weights")
  print(summary(continuation_arm$unstabilized_weight))
  
  baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                          # baseline conditions
                                          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                          venous_thromboembolism_baseline + mi_baseline + 
                                          # baseline medications
                                          ace_inhibitors_baseline + antidepressants_baseline +
                                          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                        data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    group_by(patid) %>%
    mutate(cum_numerator = cumprod(numerator))
  
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
  
  print("continuation arm stabilized weights")
  print(summary(continuation_arm$stabilized_weight))
  
  # truncate the weights at 99th percentile based on uncensored population
  continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
  upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight) %>%
    mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
  
  lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
  
  ################################################################################
  # Tapering arm
  ################################################################################
  # censor at earliest of:
  # 1. reinitiation
  # 2. end of grace period if failure to discontinue
  # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
  tapering_arm <- long_format %>%
    mutate(treatment = "tapering") %>%
    group_by(patid) %>%
    mutate(discontinued = cummax(discontinue)) %>%
    dplyr::select(-discontinue) %>%
    ungroup()
  
  
  # filter out discontinuation dates and dose reduction dates that are not within grace period
  tapering_arm <- tapering_arm %>%
    mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
           dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
  
  # create variables for 'no discontinuation date' and no dose reduction date
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
    fill(no_discontinuation_date, .direction="down") %>%
    fill(no_discontinuation_date, .direction="up") %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                     is.na(dose_reduction_date), discontinuation_date, NA)))
  
  tapering_arm <- tapering_arm %>%
    mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
  
  # followup is earliest of outcome date, enddate, censoring date
  tapering_arm <- tapering_arm %>%
    mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
    mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
    mutate(followup = ifelse(t <= followup_end, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  tapering_arm <- tapering_arm %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                            & !is.na(outcome_month), 0, outcome))
  
  # create binary variables for each of the reasons for censoring
  tapering_arm <- tapering_arm %>%
    mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
    mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
    mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
    mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
    mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
    mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
    mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
           no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
           reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
  
  # if reinitiation censor is in the same month as another censor, set to 0 as the person
  # would be censored by the other two reasons first
  tapering_arm <- tapering_arm %>%
    mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                          reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
  
  
  # in cases where the outcome occurs in the same month as the censor and the outcome
  # is BEFORE censoring, set the censor variable to 0
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                !is.na(outcome_month) & !is.na(censoring_month) &
                                                outcome_date < censoring_date & t == outcome_month,
                                              0, no_discontinuation_censor)) %>%
    mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                               !is.na(outcome_month) & !is.na(censoring_month) &
                                               outcome_date < censoring_date & t == outcome_month,
                                             0, no_dose_reduction_censor)) %>%
    mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                          !is.na(outcome_month) & !is.na(censoring_month) &
                                          outcome_date < censoring_date & t == outcome_month,
                                        0, reinitiation_censor))
  
  # keep data where followup == 1
  tapering_arm <- tapering_arm %>%
    filter(followup == 1)
  
  # first set of weights, failure to discontinuation, which only applies to t == 7
  grace_period_data <- tapering_arm %>%
    filter(t == 7) %>%
    mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
  
  ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    cancer + heart_failure + 
                                    kidney_disease + 
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = grace_period_data, family=quasibinomial())
  
  # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
  
  grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
  
  # as there is only one time point, no need to do cumultive probabilities
  grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
  print("discontinuation: no discontinuation unstabilized weights")
  print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
  
  baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                              # baseline conditions
                                              angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                              cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                              heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                              kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                              venous_thromboembolism_baseline + mi_baseline + 
                                              # baseline medications
                                              ace_inhibitors_baseline + antidepressants_baseline +
                                              antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                              arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                              dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                              lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                            data = grace_period_data, family=quasibinomial())
  
  grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
  grace_period_data <- grace_period_data %>%
    mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
  print("discontinuation: no discontinuation stabilized weights")
  print(summary(grace_period_data$stabilized_weight_no_discontinuation))
  
  grace_period_data <- grace_period_data %>% 
    dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
  tapering_arm <- tapering_arm %>%
    left_join(grace_period_data, by=c("patid", "t"))
  
  ## no dose reduction weights which should only be applied to those who have not discontinued
  # those who have discontinued without being censored cannot be censored due to no dose reduction
  # must also only occur during grace period
  not_discontinued_data <- tapering_arm %>%
    filter(t < 7) %>%
    filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
  
  
  ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                      # baseline conditions
                                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                      venous_thromboembolism_baseline + mi_baseline + 
                                                      # time-varying conditions
                                                      angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                      diabetes + heartdisease + heart_failure + heart_hypertension +
                                                      kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                      # baseline medications
                                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                      # time-varying medications
                                                      ace_inhibitors + antidepressants + antiepileptics + 
                                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                      oral_anticoagulants,
                                                    data = not_discontinued_data, family=quasibinomial())
  # all covariates have reasonable standard errors apart from rheumatoid arthritis
  
  not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
  
  # need to calculate the cumulative probabilities of the numerator and denominator
  # get the cumulative probability of remaining uncensored
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
  print("discontinuation no dose reduction unstabilized weights")
  print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
  
  baseline_ps_model_discontinuation_no_dose_reduction <-  
    glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
          # baseline conditions
          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
          venous_thromboembolism_baseline + mi_baseline + 
          # baseline medications
          ace_inhibitors_baseline + antidepressants_baseline +
          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
        data = not_discontinued_data, family=quasibinomial())
  
  
  not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
  
  # get the cumulative probabilities for the numerator
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
    ungroup()
  
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
  print("discontinuation: no dose reduction stabilized weights")
  print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
  
  not_discontinued_data <- not_discontinued_data %>%
    dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
  
  tapering_arm <- tapering_arm %>%
    left_join(not_discontinued_data, by=c("patid", "t"))
  
  # reinitiation weights, which are only applied to those who have discontinued  
  reinitiation_weight_data <- tapering_arm %>%
    filter(discontinued == 1)
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
  
  reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + depression +
                                    diabetes + kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = reinitiation_weight_data, family=quasibinomial())
  
  # very high standard errors for angina, cancer, copd, heart disease, heart failure, hypertension
  reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
  
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
  print("reinitiation weights unstabilized")
  print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
  
  reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                  cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                  kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline + 
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
  print("reinitiation weights stabilized")
  print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
  
  tapering_arm <- tapering_arm %>%
    left_join(reinitiation_weight_data, by=c("patid", "t"))
  
  # fill each of the weights downwards
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    fill(
      stabilized_weight_no_discontinuation,
      unstabilized_weight_no_discontinuation,
      stabilized_weight_no_dose_reduction,
      unstabilized_weight_no_dose_reduction,
      stabilized_weight_reinitiation,
      unstabilized_weight_reinitiation,
      .direction = "down"
    ) %>%
    ungroup()
  
  # fill in the rest in with 1s
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                         stabilized_weight_no_discontinuation),
           unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                           unstabilized_weight_no_discontinuation),
           stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                        stabilized_weight_no_dose_reduction),
           unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                          unstabilized_weight_no_dose_reduction),
           stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                   stabilized_weight_reinitiation),
           unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                     unstabilized_weight_reinitiation))
  
  # multiple all of the weights together
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
             stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
           unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
             unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
  
  # truncate the stabilized weights at the 99th percentile
  quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
  print("summary of truncated stabilized weights for tapering arm:")
  print(summary(tapering_arm$stabilized_weight))
  
  # truncate the unstabilized weights at the 99th percentile
  upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
  
  tapering_arm <- tapering_arm %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
  print("summary of truncated unstabilized weights for tapering arm:")
  print(summary(tapering_arm$unstabilized_weight))
  
  tapering_arm <- tapering_arm %>%
    mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                               1, 0)) %>%
    mutate(censored = ifelse(uncensored == 0, 1, 0))   
  
  vars <- c("stabilized_weight", "stabilized_weight_untruncated")
  
  # Function to compute the summary stats
  get_summary_stats <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p1 = quantile(x, 0.01, na.rm = TRUE),
      p50 = quantile(x, 0.50, na.rm = TRUE),
      p99 = quantile(x, 0.99, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  summary_tapering <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(tapering_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_tapering <- summary_tapering %>%
    mutate(arm = "tapering")
  
  summary_continuation <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(continuation_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_continuation <- summary_continuation %>%
    mutate(arm = "continuation")
  
  summary_weights <- rbind(summary_continuation, summary_tapering)
  
  fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/pneumonia.csv"))
  
  ################################################################################
  # check covariate balance after grace period
  ################################################################################
  covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                  "cancer_baseline", "copd_baseline", "depression_baseline",
                  "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                  "heart_hypertension_baseline", "kidney_disease_baseline",
                  "ace_inhibitors_baseline", "antidepressants_baseline", 
                  "antiplatelets_baseline", "anxiolytics_baseline",
                  "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                  "dementia_drugs_baseline", "diabetes_drug_baseline",
                  "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                  "oral_anticoagulants_baseline", "osteoporosis_baseline",
                  "venous_thromboembolism_baseline", 
                  "rheumatoid_arthritis_baseline", "mi_baseline",
                  "months_since_diagnosis", "age_T0", "gender")
  
  
  continuation_t7 <- continuation_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  tapering_t7 <- tapering_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  
  combined_t7 <- rbind(continuation_t7, tapering_t7)
  
  check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
  
  # Calculate the balance table for the subset
  balance_table <- svyCreateTableOne(vars = covariates, 
                                     strata = "treatment", 
                                     data = check_weighting, 
                                     test = FALSE, 
                                     smd = TRUE)
  
  table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  
  table_one_df$variable <- rownames(table_one_df)
  
  # balance table before weighting
  pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
  
  # Calculate the balance table for the subset
  pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                strata = "treatment", 
                                                data = pre_weighting_smd, 
                                                test = FALSE, 
                                                smd = TRUE)
  
  
  pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  pre_weighting_table$variable <- rownames(pre_weighting_table)
  
  
  if(full_sample == TRUE){
    fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pneumonia.csv"))
    fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_pneumonia.csv"))
  }
  
  ################################################################################
  # Outcome model
  ################################################################################
  
  continuation_df <- continuation_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  tapering_df <- tapering_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  combined_data <- rbind(continuation_df, tapering_df) %>%
    filter(censored == 0)
  
  
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/pneumonia.csv" ))
  
  
  # patient years by arm
  py_by_arm <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      # each row = 1 month, so 1/12 year
      patient_years = ceiling(sum(stabilized_weight * (1/12))),
      .groups = "drop"
    )
  print(py_by_arm)
  
  # number of outcomes
  weighted_events <- with(
    combined_data,
    tapply(stabilized_weight * outcome, treatment, sum)
  )
  print(weighted_events)
  
  
  crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
  # to create a survival curve, create a model with an interaction between t and treatment
  # and t_squared and treatment
  survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                + t*treatment + t_squared*treatment + 
                                + gender + age_T0 + months_since_diagnosis +
                                # baseline conditions
                                angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                venous_thromboembolism_baseline + mi_baseline +
                                # baseline medications
                                ace_inhibitors_baseline + antidepressants_baseline +
                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                              data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
  
  continuation_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "continuation")
  
  tapering_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "tapering")
  
  # create hazards under each strategy
  tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
  continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
  
  # calculate survival from the cumulative of 1-hazard
  tapering_df <- tapering_df %>%
    group_by(patid) %>%
    mutate(survival1 = cumprod(1-hazard1)) %>%
    ungroup()
  
  continuation_df <- continuation_df %>%
    group_by(patid) %>%
    mutate(survival0 = cumprod(1-hazard0)) %>%
    ungroup()
  
  # estimate risks as 1-survival
  tapering_df <- tapering_df %>%
    mutate(risk1 = 1-survival1)
  
  continuation_df <- continuation_df %>%
    mutate(risk0 = 1-survival0)
  
  # calculate mean for each time point
  tapering_risk <- tapering_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk1),
              survival = mean(survival1),
              hazard = mean(hazard1))
  
  continuation_risk <- continuation_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk0),
              survival = mean(survival0),
              hazard = mean(hazard0))
  
  risk_curves <- full_join(tapering_risk, continuation_risk)
  
  risk_curves_wider <- risk_curves %>%
    pivot_wider(
      id_cols = t,
      names_from = treatment,
      values_from = c(risk, hazard, survival),
      names_sep = "_"
    )
  
  risk_curves_wider <- risk_curves_wider %>%
    mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
  
  if(full_sample == TRUE){
    fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/pneumonia.csv"))
  }
  
  # hazard ratio is average of hazard ratios over time
  hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
  
  tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
  continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
  
  tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
  continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
  
  # calculate risk ratios and differences
  risk_ratio_12 <- tapering_risk_12/continuation_risk_12
  risk_diff_12 <- tapering_risk_12-continuation_risk_12
  
  risk_ratio_24 <- tapering_risk_24/continuation_risk_24
  risk_diff_24 <- tapering_risk_24-continuation_risk_24
  
  
  risk_df <- data.frame(
    estimate_hr = hazard_ratio_estimate,
    risk_ratio_12months = risk_ratio_12,
    risk_ratio_24months = risk_ratio_24,
    risk_diff_12months = risk_diff_12,
    risk_diff_24months = risk_diff_24,
    abs_risk_continuation_12 = continuation_risk_12,
    abs_risk_continuation_24 = continuation_risk_24,
    abs_risk_tapering_12 = tapering_risk_12, 
    abs_risk_tapering_24 = tapering_risk_24
  )
  
  rownames(risk_df) <- NULL
  
  if(full_sample == TRUE){
    fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/pneumonia.csv"))
  }
  
  risk_curves <- risk_curves %>%
    mutate(risk=risk*100)
  
  risk_curves <- rbind(
    risk_curves,
    data.frame(t = 0, risk = 0, treatment = "continuation"),
    data.frame(t = 0, risk = 0, treatment = "tapering")
  )
  
  # Plot the survival curves
  survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
    geom_line(linewidth=1)+
    labs(
      title = "Pneumonia Cumulative Incidence Curves by Treatment",
      x = "Time (months)",
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),     # tick label size
      axis.title = element_text(size = 16),    # axis title size
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    # vertical line at end of grace period
    geom_vline(xintercept = 6, 
               linetype = "dashed", 
               linewidth = 1,
               color = "black") 
  if(full_sample == TRUE){
    ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/pneumonia.png"), plot = survival_curve_plot,
           dpi=300, width = 8, height = 6)
  }
  
  return(risk_df)
  
}

tapering_pneumonia_trials(combined_12weeks, '12weeks', full_sample = TRUE)
tapering_pneumonia_trials_24weeks(combined_24weeks, '24weeks', full_sample = TRUE)


################################################################################
# FRACTURE OUTCOME
################################################################################

#### fracture outcome analysis


tapering_fracture_trials <- function(long_format, time_period, full_sample=FALSE) {
  
  long_format <- long_format %>%
    mutate(outcome_date = ifelse(fracture_outcome_date > enddate & !is.na(enddate), NA, fracture_outcome_date)) %>%
    mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
    mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
    mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
    mutate(outcome = ifelse(is.na(outcome), 0, outcome))
  
  # add discontinuation date
  long_format <- long_format %>%
    mutate(discontinuation_date = treatment_enddate + 28)
  
  long_format <- long_format %>%
    mutate(discontinuation_date = case_when(
      (!is.na(enddate) & discontinuation_date > enddate) |
        (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
      TRUE ~ discontinuation_date
    ))
  
  
  # create a variable for maximum followup
  long_format <- long_format %>%
    group_by(patid) %>%
    mutate(maximum_followup = min(fracture_outcome_date, enddate, na.rm=TRUE))
  
  # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
  
  long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                    long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                     long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                  long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
  
  
  
  long_format <- long_format %>%
    mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
    mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
    mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
  
  # assuming all available followup so take from the overall long format
  long_format <- long_format %>%
    mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
    mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
  
  ################################################################################
  # Continuation arm
  ################################################################################
  # continuation clone are only censored due to discontinuation
  continuation_arm <- long_format %>%
    mutate(treatment = "continuation") %>%
    mutate(censoring_date = discontinuation_date) %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
    mutate(censored = ifelse(is.na(censored), 0, censored))
  
  # end of followup is the earliest of censoring date or outcome date or enddate
  continuation_arm <- continuation_arm %>%
    mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
    mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
    mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  continuation_arm <- continuation_arm %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                              censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                            0, outcome))
  
  # in cases where censoring date is after the outcome date, set censor month to 0 for that month
  # this is because  followup ends at that point so ignore the censoring event
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                               censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                             0, censored))
  
  # keep rows with followup set to 0
  continuation_arm <- continuation_arm %>%
    mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
    filter(followup == 1)
  
  # create stabilized IPCW weights for the continuation arm
  ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                 # baseline conditions
                                 angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                 cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                 heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                 kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                 venous_thromboembolism_baseline + mi_baseline + 
                                 # time-varying conditions
                                 angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                 diabetes + heartdisease + heart_failure + heart_hypertension +
                                 kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                 # baseline medications
                                 ace_inhibitors_baseline  + antidepressants_baseline +
                                 antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                 arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                 dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                 lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                 # time-varying medications
                                 ace_inhibitors + antidepressants + antiepileptics + 
                                 antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                 dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                 oral_anticoagulants,
                               data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  
  continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
  
  # get the cumulative probability of remaining uncensored
  continuation_arm <- continuation_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
  print("continuation arm unstabilized weights")
  print(summary(continuation_arm$unstabilized_weight))
  
  baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                          # baseline conditions
                                          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                          venous_thromboembolism_baseline + mi_baseline + 
                                          # baseline medications
                                          ace_inhibitors_baseline + antidepressants_baseline +
                                          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                        data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    group_by(patid) %>%
    mutate(cum_numerator = cumprod(numerator))
  
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
  
  print("continuation arm stabilized weights")
  print(summary(continuation_arm$stabilized_weight))
  
  # truncate the weights at 99th percentile based on uncensored population
  continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
  upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight) %>%
    mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
  
  lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
  
  ################################################################################
  # Tapering arm
  ################################################################################
  # censor at earliest of:
  # 1. reinitiation
  # 2. end of grace period if failure to discontinue
  # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
  tapering_arm <- long_format %>%
    mutate(treatment = "tapering") %>%
    group_by(patid) %>%
    mutate(discontinued = cummax(discontinue)) %>%
    dplyr::select(-discontinue) %>%
    ungroup()
  
  
  # filter out discontinuation dates and dose reduction dates that are not within grace period
  tapering_arm <- tapering_arm %>%
    mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
           dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
  
  # create variables for 'no discontinuation date' and no dose reduction date
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
    fill(no_discontinuation_date, .direction="down") %>%
    fill(no_discontinuation_date, .direction="up") %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                     is.na(dose_reduction_date), discontinuation_date, NA)))
  
  tapering_arm <- tapering_arm %>%
    mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
  
  # followup is earliest of outcome date, enddate, censoring date
  tapering_arm <- tapering_arm %>%
    mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
    mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
    mutate(followup = ifelse(t <= followup_end, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  tapering_arm <- tapering_arm %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                            & !is.na(outcome_month), 0, outcome))
  
  # create binary variables for each of the reasons for censoring
  tapering_arm <- tapering_arm %>%
    mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
    mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
    mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
    mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
    mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
    mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
    mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
           no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
           reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
  
  # if reinitiation censor is in the same month as another censor, set to 0 as the person
  # would be censored by the other two reasons first
  tapering_arm <- tapering_arm %>%
    mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                          reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
  
  
  # in cases where the outcome occurs in the same month as the censor and the outcome
  # is BEFORE censoring, set the censor variable to 0
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                !is.na(outcome_month) & !is.na(censoring_month) &
                                                outcome_date < censoring_date & t == outcome_month,
                                              0, no_discontinuation_censor)) %>%
    mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                               !is.na(outcome_month) & !is.na(censoring_month) &
                                               outcome_date < censoring_date & t == outcome_month,
                                             0, no_dose_reduction_censor)) %>%
    mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                          !is.na(outcome_month) & !is.na(censoring_month) &
                                          outcome_date < censoring_date & t == outcome_month,
                                        0, reinitiation_censor))
  
  # keep data where followup == 1
  tapering_arm <- tapering_arm %>%
    filter(followup == 1)
  
  # first set of weights, failure to discontinuation, which only applies to t == 7
  grace_period_data <- tapering_arm %>%
    filter(t == 7) %>%
    mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
  
  ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + cancer + depression +
                                    diabetes + heartdisease + heart_failure + 
                                    kidney_disease + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = grace_period_data, family=quasibinomial())
  
  # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
  
  grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
  
  # as there is only one time point, no need to do cumultive probabilities
  grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
  print("discontinuation: no discontinuation unstabilized weights")
  print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
  
  baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                              # baseline conditions
                                              angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                              cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                              heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                              kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                              venous_thromboembolism_baseline + mi_baseline + 
                                              # baseline medications
                                              ace_inhibitors_baseline + antidepressants_baseline +
                                              antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                              arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                              dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                              lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                            data = grace_period_data, family=quasibinomial())
  
  grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
  grace_period_data <- grace_period_data %>%
    mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
  print("discontinuation: no discontinuation stabilized weights")
  print(summary(grace_period_data$stabilized_weight_no_discontinuation))
  
  grace_period_data <- grace_period_data %>% 
    dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
  tapering_arm <- tapering_arm %>%
    left_join(grace_period_data, by=c("patid", "t"))
  
  ## no dose reduction weights which should only be applied to those who have not discontinued
  # those who have discontinued without being censored cannot be censored due to no dose reduction
  # must also only occur during grace period
  not_discontinued_data <- tapering_arm %>%
    filter(t < 7) %>%
    filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
  
  
  ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                      # baseline conditions
                                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                      venous_thromboembolism_baseline + mi_baseline + 
                                                      # time-varying conditions
                                                      angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                      diabetes + heartdisease + heart_failure + heart_hypertension +
                                                      kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                      # baseline medications
                                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                      # time-varying medications
                                                      ace_inhibitors + antidepressants + antiepileptics + 
                                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                      oral_anticoagulants,
                                                    data = not_discontinued_data, family=quasibinomial())
  # all covariates have reasonable standard errors apart from rheumatoid arthritis
  
  not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
  
  # need to calculate the cumulative probabilities of the numerator and denominator
  # get the cumulative probability of remaining uncensored
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
  print("discontinuation no dose reduction unstabilized weights")
  print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
  
  baseline_ps_model_discontinuation_no_dose_reduction <-  
    glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
          # baseline conditions
          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
          venous_thromboembolism_baseline + mi_baseline + 
          # baseline medications
          ace_inhibitors_baseline + antidepressants_baseline +
          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
        data = not_discontinued_data, family=quasibinomial())
  
  
  not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
  
  # get the cumulative probabilities for the numerator
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
    ungroup()
  
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
  print("discontinuation: no dose reduction stabilized weights")
  print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
  
  not_discontinued_data <- not_discontinued_data %>%
    dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
  
  tapering_arm <- tapering_arm %>%
    left_join(not_discontinued_data, by=c("patid", "t"))
  
  # reinitiation weights, which are only applied to those who have discontinued  
  reinitiation_weight_data <- tapering_arm %>%
    filter(discontinued == 1)
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
  
  reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + depression +
                                    diabetes + kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = reinitiation_weight_data, family=quasibinomial())
  
  # very high standard errors for angina, cancer, copd, heart disease, heart failure, hypertension
  reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
  
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
  print("reinitiation weights unstabilized")
  print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
  
  reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                  cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                  kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline + 
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
  print("reinitiation weights stabilized")
  print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
  
  tapering_arm <- tapering_arm %>%
    left_join(reinitiation_weight_data, by=c("patid", "t"))
  
  # fill each of the weights downwards
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    fill(
      stabilized_weight_no_discontinuation,
      unstabilized_weight_no_discontinuation,
      stabilized_weight_no_dose_reduction,
      unstabilized_weight_no_dose_reduction,
      stabilized_weight_reinitiation,
      unstabilized_weight_reinitiation,
      .direction = "down"
    ) %>%
    ungroup()
  
  # fill in the rest in with 1s
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                         stabilized_weight_no_discontinuation),
           unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                           unstabilized_weight_no_discontinuation),
           stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                        stabilized_weight_no_dose_reduction),
           unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                          unstabilized_weight_no_dose_reduction),
           stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                   stabilized_weight_reinitiation),
           unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                     unstabilized_weight_reinitiation))
  
  # multiple all of the weights together
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
             stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
           unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
             unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
  
  # truncate the stabilized weights at the 99th percentile
  quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
  print("summary of truncated stabilized weights for tapering arm:")
  print(summary(tapering_arm$stabilized_weight))
  
  # truncate the unstabilized weights at the 99th percentile
  upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
  
  tapering_arm <- tapering_arm %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
  print("summary of truncated unstabilized weights for tapering arm:")
  print(summary(tapering_arm$unstabilized_weight))
  
  tapering_arm <- tapering_arm %>%
    mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                               1, 0)) %>%
    mutate(censored = ifelse(uncensored == 0, 1, 0))   
  
  vars <- c("stabilized_weight", "stabilized_weight_untruncated")
  
  # Function to compute the summary stats
  get_summary_stats <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p1 = quantile(x, 0.01, na.rm = TRUE),
      p50 = quantile(x, 0.50, na.rm = TRUE),
      p99 = quantile(x, 0.99, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  summary_tapering <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(tapering_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_tapering <- summary_tapering %>%
    mutate(arm = "tapering")
  
  summary_continuation <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(continuation_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_continuation <- summary_continuation %>%
    mutate(arm = "continuation")
  
  summary_weights <- rbind(summary_continuation, summary_tapering)
  
  fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/fracture.csv"))
  
  ################################################################################
  # check covariate balance after grace period
  ################################################################################
  covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                  "cancer_baseline", "copd_baseline", "depression_baseline",
                  "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                  "heart_hypertension_baseline", "kidney_disease_baseline",
                  "ace_inhibitors_baseline", "antidepressants_baseline", 
                  "antiplatelets_baseline", "anxiolytics_baseline",
                  "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                  "dementia_drugs_baseline", "diabetes_drug_baseline",
                  "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                  "oral_anticoagulants_baseline", "osteoporosis_baseline",
                  "venous_thromboembolism_baseline", 
                  "rheumatoid_arthritis_baseline", "mi_baseline",
                  "months_since_diagnosis", "age_T0", "gender")
  
  
  continuation_t7 <- continuation_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  tapering_t7 <- tapering_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  
  combined_t7 <- rbind(continuation_t7, tapering_t7)
  
  check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
  
  # Calculate the balance table for the subset
  balance_table <- svyCreateTableOne(vars = covariates, 
                                     strata = "treatment", 
                                     data = check_weighting, 
                                     test = FALSE, 
                                     smd = TRUE)
  
  table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  
  table_one_df$variable <- rownames(table_one_df)
  
  # balance table before weighting
  pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
  
  # Calculate the balance table for the subset
  pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                strata = "treatment", 
                                                data = pre_weighting_smd, 
                                                test = FALSE, 
                                                smd = TRUE)
  
  
  pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  pre_weighting_table$variable <- rownames(pre_weighting_table)
  
  
  if(full_sample == TRUE){
    fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/fracture.csv"))
    fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_fracture.csv"))
  }
  
  ################################################################################
  # Outcome model
  ################################################################################
  
  continuation_df <- continuation_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  tapering_df <- tapering_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  combined_data <- rbind(continuation_df, tapering_df) %>%
    filter(censored == 0)
  
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/fracture.csv" ))
  
  # patient years by arm
  py_by_arm <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      # each row = 1 month, so 1/12 year
      patient_years = ceiling(sum(stabilized_weight * (1/12))),
      .groups = "drop"
    )
  print(py_by_arm)
  
  # number of outcomes
  weighted_events <- with(
    combined_data,
    tapply(stabilized_weight * outcome, treatment, sum)
  )
  print(weighted_events)
  
  crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
  # to create a survival curve, create a model with an interaction between t and treatment
  # and t_squared and treatment
  survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                + t*treatment + t_squared*treatment + 
                                + gender + age_T0 + months_since_diagnosis +
                                # baseline conditions
                                angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                venous_thromboembolism_baseline + mi_baseline +
                                # baseline medications
                                ace_inhibitors_baseline + antidepressants_baseline +
                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                              data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
  
  continuation_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "continuation")
  
  tapering_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "tapering")
  
  # create hazards under each strategy
  tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
  continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
  
  # calculate survival from the cumulative of 1-hazard
  tapering_df <- tapering_df %>%
    group_by(patid) %>%
    mutate(survival1 = cumprod(1-hazard1)) %>%
    ungroup()
  
  continuation_df <- continuation_df %>%
    group_by(patid) %>%
    mutate(survival0 = cumprod(1-hazard0)) %>%
    ungroup()
  
  # estimate risks as 1-survival
  tapering_df <- tapering_df %>%
    mutate(risk1 = 1-survival1)
  
  continuation_df <- continuation_df %>%
    mutate(risk0 = 1-survival0)
  
  # calculate mean for each time point
  tapering_risk <- tapering_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk1),
              survival = mean(survival1),
              hazard = mean(hazard1))
  
  continuation_risk <- continuation_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk0),
              survival = mean(survival0),
              hazard = mean(hazard0))
  
  risk_curves <- full_join(tapering_risk, continuation_risk)
  
  risk_curves_wider <- risk_curves %>%
    pivot_wider(
      id_cols = t,
      names_from = treatment,
      values_from = c(risk, hazard, survival),
      names_sep = "_"
    )
  
  risk_curves_wider <- risk_curves_wider %>%
    mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
  
  if(full_sample == TRUE){
    fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/fracture.csv"))
  }
  
  # hazard ratio is average of hazard ratios over time
  hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
  
  tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
  continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
  
  tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
  continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
  
  # calculate risk ratios and differences
  risk_ratio_12 <- tapering_risk_12/continuation_risk_12
  risk_diff_12 <- tapering_risk_12-continuation_risk_12
  
  risk_ratio_24 <- tapering_risk_24/continuation_risk_24
  risk_diff_24 <- tapering_risk_24-continuation_risk_24
  
  
  risk_df <- data.frame(
    estimate_hr = hazard_ratio_estimate,
    risk_ratio_12months = risk_ratio_12,
    risk_ratio_24months = risk_ratio_24,
    risk_diff_12months = risk_diff_12,
    risk_diff_24months = risk_diff_24,
    abs_risk_continuation_12 = continuation_risk_12,
    abs_risk_continuation_24 = continuation_risk_24,
    abs_risk_tapering_12 = tapering_risk_12, 
    abs_risk_tapering_24 = tapering_risk_24
  )
  rownames(risk_df) <- NULL
  
  if(full_sample == TRUE){
    fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/fracture.csv"))
  }
  
  risk_curves <- risk_curves %>%
    mutate(risk=risk*100)
  
  risk_curves <- rbind(
    risk_curves,
    data.frame(t = 0, risk = 0, treatment = "continuation"),
    data.frame(t = 0, risk = 0, treatment = "tapering")
  )
  
  # Plot the survival curves
  survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
    geom_line(linewidth=1)+
    labs(
      title = "Fracture Cumulative Incidence Curves by Treatment",
      x = "Time (months)",
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),     # tick label size
      axis.title = element_text(size = 16),    # axis title size
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    # vertical line at end of grace period
    geom_vline(xintercept = 6, 
               linetype = "dashed", 
               linewidth = 1,
               color = "black") 
  
  if(full_sample == TRUE){
    ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/fracture.png"), plot = survival_curve_plot,
           dpi=300, width = 8, height = 6)
  }
  
  return(risk_df)
  
}

tapering_fracture_trials(combined_12weeks, '12weeks', full_sample = TRUE)
tapering_fracture_trials(combined_24weeks, '24weeks', full_sample = TRUE)

################################################################################
# DELIRIUM OUTCOME
################################################################################

#### delirium outcome analysis


tapering_delirium_trials <- function(long_format, time_period, full_sample=FALSE) {
  
  long_format <- long_format %>%
    mutate(outcome_date = ifelse(delirium_outcome_date > enddate & !is.na(enddate), NA, delirium_outcome_date)) %>%
    mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
    mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
    mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
    mutate(outcome = ifelse(is.na(outcome), 0, outcome))
  
  # add discontinuation date
  long_format <- long_format %>%
    mutate(discontinuation_date = treatment_enddate + 28)
  
  long_format <- long_format %>%
    mutate(discontinuation_date = case_when(
      (!is.na(enddate) & discontinuation_date > enddate) |
        (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
      TRUE ~ discontinuation_date
    ))
  
  
  # create a variable for maximum followup
  long_format <- long_format %>%
    group_by(patid) %>%
    mutate(maximum_followup = min(delirium_outcome_date, enddate, na.rm=TRUE))
  
  # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
  
  long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                    long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                     long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                  long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
  
  
  
  long_format <- long_format %>%
    mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
    mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
    mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
  
  # assuming all available followup so take from the overall long format
  long_format <- long_format %>%
    mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
    mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
  
  ################################################################################
  # Continuation arm
  ################################################################################
  # continuation clone are only censored due to discontinuation
  continuation_arm <- long_format %>%
    mutate(treatment = "continuation") %>%
    mutate(censoring_date = discontinuation_date) %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
    mutate(censored = ifelse(is.na(censored), 0, censored))
  
  # end of followup is the earliest of censoring date or outcome date or enddate
  continuation_arm <- continuation_arm %>%
    mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
    mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
    mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  continuation_arm <- continuation_arm %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                              censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                            0, outcome))
  
  # in cases where censoring date is after the outcome date, set censor month to 0 for that month
  # this is because  followup ends at that point so ignore the censoring event
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                               censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                             0, censored))
  
  # keep rows with followup set to 0
  continuation_arm <- continuation_arm %>%
    mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
    filter(followup == 1)
  
  # create stabilized IPCW weights for the continuation arm
  ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                 # baseline conditions
                                 angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                 cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                 heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                 kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                 venous_thromboembolism_baseline + mi_baseline + 
                                 # time-varying conditions
                                 angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                 diabetes + heartdisease + heart_failure + heart_hypertension +
                                 kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                 # baseline medications
                                 ace_inhibitors_baseline  + antidepressants_baseline +
                                 antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                 arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                 dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                 lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                 # time-varying medications
                                 ace_inhibitors + antidepressants + antiepileptics + 
                                 antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                 dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                 oral_anticoagulants,
                               data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  
  continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
  
  # get the cumulative probability of remaining uncensored
  continuation_arm <- continuation_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
  print("continuation arm unstabilized weights")
  print(summary(continuation_arm$unstabilized_weight))
  
  baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                          # baseline conditions
                                          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                          venous_thromboembolism_baseline + mi_baseline + 
                                          # baseline medications
                                          ace_inhibitors_baseline + antidepressants_baseline +
                                          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                        data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    group_by(patid) %>%
    mutate(cum_numerator = cumprod(numerator))
  
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
  
  print("continuation arm stabilized weights")
  print(summary(continuation_arm$stabilized_weight))
  
  # truncate the weights at 99th percentile based on uncensored population
  continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
  upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight) %>%
    mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
  
  lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
  
  ################################################################################
  # Tapering arm
  ################################################################################
  # censor at earliest of:
  # 1. reinitiation
  # 2. end of grace period if failure to discontinue
  # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
  tapering_arm <- long_format %>%
    mutate(treatment = "tapering") %>%
    group_by(patid) %>%
    mutate(discontinued = cummax(discontinue)) %>%
    dplyr::select(-discontinue) %>%
    ungroup()
  
  
  # filter out discontinuation dates and dose reduction dates that are not within grace period
  tapering_arm <- tapering_arm %>%
    mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
           dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
  
  # create variables for 'no discontinuation date' and no dose reduction date
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
    fill(no_discontinuation_date, .direction="down") %>%
    fill(no_discontinuation_date, .direction="up") %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                     is.na(dose_reduction_date), discontinuation_date, NA)))
  
  tapering_arm <- tapering_arm %>%
    mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
  
  # followup is earliest of outcome date, enddate, censoring date
  tapering_arm <- tapering_arm %>%
    mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
    mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
    mutate(followup = ifelse(t <= followup_end, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  tapering_arm <- tapering_arm %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                            & !is.na(outcome_month), 0, outcome))
  
  # create binary variables for each of the reasons for censoring
  tapering_arm <- tapering_arm %>%
    mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
    mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
    mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
    mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
    mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
    mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
    mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
           no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
           reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
  
  # if reinitiation censor is in the same month as another censor, set to 0 as the person
  # would be censored by the other two reasons first
  tapering_arm <- tapering_arm %>%
    mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                          reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
  
  
  # in cases where the outcome occurs in the same month as the censor and the outcome
  # is BEFORE censoring, set the censor variable to 0
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                !is.na(outcome_month) & !is.na(censoring_month) &
                                                outcome_date < censoring_date & t == outcome_month,
                                              0, no_discontinuation_censor)) %>%
    mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                               !is.na(outcome_month) & !is.na(censoring_month) &
                                               outcome_date < censoring_date & t == outcome_month,
                                             0, no_dose_reduction_censor)) %>%
    mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                          !is.na(outcome_month) & !is.na(censoring_month) &
                                          outcome_date < censoring_date & t == outcome_month,
                                        0, reinitiation_censor))
  
  # keep data where followup == 1
  tapering_arm <- tapering_arm %>%
    filter(followup == 1)
  
  # first set of weights, failure to discontinuation, which only applies to t == 7
  grace_period_data <- tapering_arm %>%
    filter(t == 7) %>%
    mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
  
  ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + cancer + depression +
                                    diabetes + heartdisease + heart_failure + 
                                    kidney_disease + osteoporosis + venous_thromboembolism +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = grace_period_data, family=quasibinomial())
  
  # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
  
  grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
  
  # as there is only one time point, no need to do cumulative probabilities
  grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
  print("discontinuation: no discontinuation unstabilized weights")
  print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
  
  baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                              # baseline conditions
                                              angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                              cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                              heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                              kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                              venous_thromboembolism_baseline + mi_baseline + 
                                              # baseline medications
                                              ace_inhibitors_baseline + antidepressants_baseline +
                                              antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                              arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                              dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                              lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                            data = grace_period_data, family=quasibinomial())
  
  grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
  grace_period_data <- grace_period_data %>%
    mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
  print("discontinuation: no discontinuation stabilized weights")
  print(summary(grace_period_data$stabilized_weight_no_discontinuation))
  
  grace_period_data <- grace_period_data %>% 
    dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
  tapering_arm <- tapering_arm %>%
    left_join(grace_period_data, by=c("patid", "t"))
  
  ## no dose reduction weights which should only be applied to those who have not discontinued
  # those who have discontinued without being censored cannot be censored due to no dose reduction
  # must also only occur during grace period
  not_discontinued_data <- tapering_arm %>%
    filter(t < 7) %>%
    filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
  
  
  ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                      # baseline conditions
                                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                      venous_thromboembolism_baseline + mi_baseline + 
                                                      # time-varying conditions
                                                      angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                      diabetes + heartdisease + heart_failure + heart_hypertension +
                                                      kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                      # baseline medications
                                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                      # time-varying medications
                                                      ace_inhibitors + antidepressants + antiepileptics + 
                                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                      oral_anticoagulants,
                                                    data = not_discontinued_data, family=quasibinomial())
  # all covariates have reasonable standard errors apart from rheumatoid arthritis
  
  not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
  
  # need to calculate the cumulative probabilities of the numerator and denominator
  # get the cumulative probability of remaining uncensored
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
  print("discontinuation no dose reduction unstabilized weights")
  print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
  
  baseline_ps_model_discontinuation_no_dose_reduction <-  
    glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
          # baseline conditions
          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
          venous_thromboembolism_baseline + mi_baseline + 
          # baseline medications
          ace_inhibitors_baseline + antidepressants_baseline +
          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
        data = not_discontinued_data, family=quasibinomial())
  
  
  not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
  
  # get the cumulative probabilities for the numerator
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
    ungroup()
  
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
  print("discontinuation: no dose reduction stabilized weights")
  print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
  
  not_discontinued_data <- not_discontinued_data %>%
    dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
  
  tapering_arm <- tapering_arm %>%
    left_join(not_discontinued_data, by=c("patid", "t"))
  
  # reinitiation weights, which are only applied to those who have discontinued  
  reinitiation_weight_data <- tapering_arm %>%
    filter(discontinued == 1)
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
  
  reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + depression +
                                    diabetes + kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = reinitiation_weight_data, family=quasibinomial())
  
  # very high standard errors for angina, cancer, copd, heart disease, heart failure, hypertension
  reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
  
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
  print("reinitiation weights unstabilized")
  print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
  
  reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                  cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                  kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline + 
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
  print("reinitiation weights stabilized")
  print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
  
  tapering_arm <- tapering_arm %>%
    left_join(reinitiation_weight_data, by=c("patid", "t"))
  
  # fill each of the weights downwards
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    fill(
      stabilized_weight_no_discontinuation,
      unstabilized_weight_no_discontinuation,
      stabilized_weight_no_dose_reduction,
      unstabilized_weight_no_dose_reduction,
      stabilized_weight_reinitiation,
      unstabilized_weight_reinitiation,
      .direction = "down"
    ) %>%
    ungroup()
  
  # fill in the rest in with 1s
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                         stabilized_weight_no_discontinuation),
           unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                           unstabilized_weight_no_discontinuation),
           stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                        stabilized_weight_no_dose_reduction),
           unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                          unstabilized_weight_no_dose_reduction),
           stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                   stabilized_weight_reinitiation),
           unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                     unstabilized_weight_reinitiation))
  
  # multiple all of the weights together
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
             stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
           unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
             unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
  
  # truncate the stabilized weights at the 99th percentile
  quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
  print("summary of truncated stabilized weights for tapering arm:")
  print(summary(tapering_arm$stabilized_weight))
  
  # truncate the unstabilized weights at the 99th percentile
  upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
  
  tapering_arm <- tapering_arm %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
  print("summary of truncated unstabilized weights for tapering arm:")
  print(summary(tapering_arm$unstabilized_weight))
  
  tapering_arm <- tapering_arm %>%
    mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                               1, 0)) %>%
    mutate(censored = ifelse(uncensored == 0, 1, 0))   
  
  vars <- c("stabilized_weight", "stabilized_weight_untruncated")
  
  # Function to compute the summary stats
  get_summary_stats <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p1 = quantile(x, 0.01, na.rm = TRUE),
      p50 = quantile(x, 0.50, na.rm = TRUE),
      p99 = quantile(x, 0.99, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  summary_tapering <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(tapering_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_tapering <- summary_tapering %>%
    mutate(arm = "tapering")
  
  summary_continuation <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(continuation_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_continuation <- summary_continuation %>%
    mutate(arm = "continuation")
  
  summary_weights <- rbind(summary_continuation, summary_tapering)
  
  fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/delirium.csv"))
  
  ################################################################################
  # check covariate balance after grace period
  ################################################################################
  covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                  "cancer_baseline", "copd_baseline", "depression_baseline",
                  "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                  "heart_hypertension_baseline", "kidney_disease_baseline",
                  "ace_inhibitors_baseline", "antidepressants_baseline", 
                  "antiplatelets_baseline", "anxiolytics_baseline",
                  "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                  "dementia_drugs_baseline", "diabetes_drug_baseline",
                  "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                  "oral_anticoagulants_baseline", "osteoporosis_baseline",
                  "venous_thromboembolism_baseline", 
                  "rheumatoid_arthritis_baseline", "mi_baseline",
                  "months_since_diagnosis", "age_T0", "gender")
  
  
  continuation_t7 <- continuation_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  tapering_t7 <- tapering_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  
  combined_t7 <- rbind(continuation_t7, tapering_t7)
  
  check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
  
  # Calculate the balance table for the subset
  balance_table <- svyCreateTableOne(vars = covariates, 
                                     strata = "treatment", 
                                     data = check_weighting, 
                                     test = FALSE, 
                                     smd = TRUE)
  
  table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  
  table_one_df$variable <- rownames(table_one_df)
  
  # balance table before weighting
  pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
  
  # Calculate the balance table for the subset
  pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                strata = "treatment", 
                                                data = pre_weighting_smd, 
                                                test = FALSE, 
                                                smd = TRUE)
  
  
  pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  pre_weighting_table$variable <- rownames(pre_weighting_table)
  
  
  if(full_sample == TRUE){
    fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/delirium.csv"))
    fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_delirium.csv"))
  }
  
  ################################################################################
  # Outcome model
  ################################################################################
  
  continuation_df <- continuation_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  tapering_df <- tapering_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  combined_data <- rbind(continuation_df, tapering_df) %>%
    filter(censored == 0)
  
  
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/delirium.csv" ))
  
  
  # patient years by arm
  py_by_arm <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      # each row = 1 month, so 1/12 year
      patient_years = ceiling(sum(stabilized_weight * (1/12))),
      .groups = "drop"
    )
  print(py_by_arm)
  
  # number of outcomes
  weighted_events <- with(
    combined_data,
    tapply(stabilized_weight * outcome, treatment, sum)
  )
  print(weighted_events)
  
  crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
  # to create a survival curve, create a model with an interaction between t and treatment
  # and t_squared and treatment
  survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                + t*treatment + t_squared*treatment + 
                                + gender + age_T0 + months_since_diagnosis +
                                # baseline conditions
                                angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                venous_thromboembolism_baseline + mi_baseline +
                                # baseline medications
                                ace_inhibitors_baseline + antidepressants_baseline +
                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                              data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
  
  continuation_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "continuation")
  
  tapering_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "tapering")
  
  # create hazards under each strategy
  tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
  continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
  
  # calculate survival from the cumulative of 1-hazard
  tapering_df <- tapering_df %>%
    group_by(patid) %>%
    mutate(survival1 = cumprod(1-hazard1)) %>%
    ungroup()
  
  continuation_df <- continuation_df %>%
    group_by(patid) %>%
    mutate(survival0 = cumprod(1-hazard0)) %>%
    ungroup()
  
  # estimate risks as 1-survival
  tapering_df <- tapering_df %>%
    mutate(risk1 = 1-survival1)
  
  continuation_df <- continuation_df %>%
    mutate(risk0 = 1-survival0)
  
  # calculate mean for each time point
  tapering_risk <- tapering_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk1),
              survival = mean(survival1),
              hazard = mean(hazard1))
  
  continuation_risk <- continuation_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk0),
              survival = mean(survival0),
              hazard = mean(hazard0))
  
  risk_curves <- full_join(tapering_risk, continuation_risk)
  
  risk_curves_wider <- risk_curves %>%
    pivot_wider(
      id_cols = t,
      names_from = treatment,
      values_from = c(risk, hazard, survival),
      names_sep = "_"
    )
  
  risk_curves_wider <- risk_curves_wider %>%
    mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
  
  if(full_sample == TRUE){
    fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/delirium.csv"))
  }
  
  # hazard ratio is average of hazard ratios over time
  hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
  
  tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
  continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
  
  tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
  continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
  
  # calculate risk ratios and differences
  risk_ratio_12 <- tapering_risk_12/continuation_risk_12
  risk_diff_12 <- tapering_risk_12-continuation_risk_12
  
  risk_ratio_24 <- tapering_risk_24/continuation_risk_24
  risk_diff_24 <- tapering_risk_24-continuation_risk_24
  
  
  risk_df <- data.frame(
    estimate_hr = hazard_ratio_estimate,
    risk_ratio_12months = risk_ratio_12,
    risk_ratio_24months = risk_ratio_24,
    risk_diff_12months = risk_diff_12,
    risk_diff_24months = risk_diff_24,
    abs_risk_continuation_12 = continuation_risk_12,
    abs_risk_continuation_24 = continuation_risk_24,
    abs_risk_tapering_12 = tapering_risk_12, 
    abs_risk_tapering_24 = tapering_risk_24
  )
  rownames(risk_df) <- NULL
  
  
  if(full_sample == TRUE){
    fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/delirium.csv"))
  }
  
  risk_curves <- risk_curves %>%
    mutate(risk=risk*100)
  
  risk_curves <- rbind(
    risk_curves,
    data.frame(t = 0, risk = 0, treatment = "continuation"),
    data.frame(t = 0, risk = 0, treatment = "tapering")
  )
  
  # Plot the survival curves
  survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
    geom_line(linewidth=1)+
    labs(
      title = "Delirium Cumulative Incidence Curves by Treatment",
      x = "Time (months)",
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),     # tick label size
      axis.title = element_text(size = 16),    # axis title size
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    # vertical line at end of grace period
    geom_vline(xintercept = 6, 
               linetype = "dashed", 
               linewidth = 1,
               color = "black") 
  if(full_sample == TRUE){
    ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/delirium.png"), plot = survival_curve_plot,
           dpi=300, width = 8, height = 6)
  }
  
  return(risk_df)
  
}


tapering_delirium_trials_24weeks <- function(long_format, time_period, full_sample=FALSE) {
  
  long_format <- long_format %>%
    mutate(outcome_date = ifelse(delirium_outcome_date > enddate & !is.na(enddate), NA, delirium_outcome_date)) %>%
    mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
    mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
    mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
    mutate(outcome = ifelse(is.na(outcome), 0, outcome))
  
  # add discontinuation date
  long_format <- long_format %>%
    mutate(discontinuation_date = treatment_enddate + 28)
  
  long_format <- long_format %>%
    mutate(discontinuation_date = case_when(
      (!is.na(enddate) & discontinuation_date > enddate) |
        (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
      TRUE ~ discontinuation_date
    ))
  
  
  # create a variable for maximum followup
  long_format <- long_format %>%
    group_by(patid) %>%
    mutate(maximum_followup = min(delirium_outcome_date, enddate, na.rm=TRUE))
  
  # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
  
  long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                    long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                     long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                  long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
  
  
  
  long_format <- long_format %>%
    mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
    mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
    mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
  
  # assuming all available followup so take from the overall long format
  long_format <- long_format %>%
    mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
    mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
  
  ################################################################################
  # Continuation arm
  ################################################################################
  # continuation clone are only censored due to discontinuation
  continuation_arm <- long_format %>%
    mutate(treatment = "continuation") %>%
    mutate(censoring_date = discontinuation_date) %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
    mutate(censored = ifelse(is.na(censored), 0, censored))
  
  # end of followup is the earliest of censoring date or outcome date or enddate
  continuation_arm <- continuation_arm %>%
    mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
    mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
    mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  continuation_arm <- continuation_arm %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                              censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                            0, outcome))
  
  # in cases where censoring date is after the outcome date, set censor month to 0 for that month
  # this is because  followup ends at that point so ignore the censoring event
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                               censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                             0, censored))
  
  # keep rows with followup set to 0
  continuation_arm <- continuation_arm %>%
    mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
    filter(followup == 1)
  
  # create stabilized IPCW weights for the continuation arm
  ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                 # baseline conditions
                                 angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                 cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                 heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                 kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                 venous_thromboembolism_baseline + mi_baseline + 
                                 # time-varying conditions
                                 angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                 diabetes + heartdisease + heart_failure + heart_hypertension +
                                 kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                 # baseline medications
                                 ace_inhibitors_baseline  + antidepressants_baseline +
                                 antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                 arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                 dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                 lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                 # time-varying medications
                                 ace_inhibitors + antidepressants + antiepileptics + 
                                 antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                 dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                 oral_anticoagulants,
                               data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  
  continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
  
  # get the cumulative probability of remaining uncensored
  continuation_arm <- continuation_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
  print("continuation arm unstabilized weights")
  print(summary(continuation_arm$unstabilized_weight))
  
  baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                          # baseline conditions
                                          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                          venous_thromboembolism_baseline + mi_baseline + 
                                          # baseline medications
                                          ace_inhibitors_baseline + antidepressants_baseline +
                                          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                        data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    group_by(patid) %>%
    mutate(cum_numerator = cumprod(numerator))
  
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
  
  print("continuation arm stabilized weights")
  print(summary(continuation_arm$stabilized_weight))
  
  # truncate the weights at 99th percentile based on uncensored population
  continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
  upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight) %>%
    mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
  
  lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
  
  ################################################################################
  # Tapering arm
  ################################################################################
  # censor at earliest of:
  # 1. reinitiation
  # 2. end of grace period if failure to discontinue
  # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
  tapering_arm <- long_format %>%
    mutate(treatment = "tapering") %>%
    group_by(patid) %>%
    mutate(discontinued = cummax(discontinue)) %>%
    dplyr::select(-discontinue) %>%
    ungroup()
  
  
  # filter out discontinuation dates and dose reduction dates that are not within grace period
  tapering_arm <- tapering_arm %>%
    mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
           dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
  
  # create variables for 'no discontinuation date' and no dose reduction date
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
    fill(no_discontinuation_date, .direction="down") %>%
    fill(no_discontinuation_date, .direction="up") %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                     is.na(dose_reduction_date), discontinuation_date, NA)))
  
  tapering_arm <- tapering_arm %>%
    mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
  
  # followup is earliest of outcome date, enddate, censoring date
  tapering_arm <- tapering_arm %>%
    mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
    mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
    mutate(followup = ifelse(t <= followup_end, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  tapering_arm <- tapering_arm %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                            & !is.na(outcome_month), 0, outcome))
  
  # create binary variables for each of the reasons for censoring
  tapering_arm <- tapering_arm %>%
    mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
    mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
    mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
    mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
    mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
    mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
    mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
           no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
           reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
  
  # if reinitiation censor is in the same month as another censor, set to 0 as the person
  # would be censored by the other two reasons first
  tapering_arm <- tapering_arm %>%
    mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                          reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
  
  
  # in cases where the outcome occurs in the same month as the censor and the outcome
  # is BEFORE censoring, set the censor variable to 0
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                !is.na(outcome_month) & !is.na(censoring_month) &
                                                outcome_date < censoring_date & t == outcome_month,
                                              0, no_discontinuation_censor)) %>%
    mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                               !is.na(outcome_month) & !is.na(censoring_month) &
                                               outcome_date < censoring_date & t == outcome_month,
                                             0, no_dose_reduction_censor)) %>%
    mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                          !is.na(outcome_month) & !is.na(censoring_month) &
                                          outcome_date < censoring_date & t == outcome_month,
                                        0, reinitiation_censor))
  
  # keep data where followup == 1
  tapering_arm <- tapering_arm %>%
    filter(followup == 1)
  
  # first set of weights, failure to discontinuation, which only applies to t == 7
  grace_period_data <- tapering_arm %>%
    filter(t == 7) %>%
    mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
  
  ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    cancer  + heart_failure + 
                                    kidney_disease + heart_hypertension + 
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = grace_period_data, family=quasibinomial())
  
  # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
  
  grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
  
  # as there is only one time point, no need to do cumulative probabilities
  grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
  print("discontinuation: no discontinuation unstabilized weights")
  print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
  
  baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                              # baseline conditions
                                              angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                              cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                              heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                              kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                              venous_thromboembolism_baseline + mi_baseline + 
                                              # baseline medications
                                              ace_inhibitors_baseline + antidepressants_baseline +
                                              antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                              arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                              dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                              lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                            data = grace_period_data, family=quasibinomial())
  
  grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
  grace_period_data <- grace_period_data %>%
    mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
  print("discontinuation: no discontinuation stabilized weights")
  print(summary(grace_period_data$stabilized_weight_no_discontinuation))
  
  grace_period_data <- grace_period_data %>% 
    dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
  tapering_arm <- tapering_arm %>%
    left_join(grace_period_data, by=c("patid", "t"))
  
  ## no dose reduction weights which should only be applied to those who have not discontinued
  # those who have discontinued without being censored cannot be censored due to no dose reduction
  # must also only occur during grace period
  not_discontinued_data <- tapering_arm %>%
    filter(t < 7) %>%
    filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
  
  
  ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                      # baseline conditions
                                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                      venous_thromboembolism_baseline + mi_baseline + 
                                                      # time-varying conditions
                                                      angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                      diabetes + heartdisease + heart_failure + heart_hypertension +
                                                      kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                      # baseline medications
                                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                      # time-varying medications
                                                      ace_inhibitors + antidepressants + antiepileptics + 
                                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                      oral_anticoagulants,
                                                    data = not_discontinued_data, family=quasibinomial())
  # all covariates have reasonable standard errors apart from rheumatoid arthritis
  
  not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
  
  # need to calculate the cumulative probabilities of the numerator and denominator
  # get the cumulative probability of remaining uncensored
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
  print("discontinuation no dose reduction unstabilized weights")
  print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
  
  baseline_ps_model_discontinuation_no_dose_reduction <-  
    glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
          # baseline conditions
          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
          venous_thromboembolism_baseline + mi_baseline + 
          # baseline medications
          ace_inhibitors_baseline + antidepressants_baseline +
          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
        data = not_discontinued_data, family=quasibinomial())
  
  
  not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
  
  # get the cumulative probabilities for the numerator
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
    ungroup()
  
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
  print("discontinuation: no dose reduction stabilized weights")
  print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
  
  not_discontinued_data <- not_discontinued_data %>%
    dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
  
  tapering_arm <- tapering_arm %>%
    left_join(not_discontinued_data, by=c("patid", "t"))
  
  # reinitiation weights, which are only applied to those who have discontinued  
  reinitiation_weight_data <- tapering_arm %>%
    filter(discontinued == 1)
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
  
  reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + depression +
                                    diabetes + kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = reinitiation_weight_data, family=quasibinomial())
  
  # very high standard errors for angina, cancer, copd, heart disease, heart failure, hypertension
  reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
  
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
  print("reinitiation weights unstabilized")
  print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
  
  reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                  cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                  kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline + 
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
  print("reinitiation weights stabilized")
  print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
  
  tapering_arm <- tapering_arm %>%
    left_join(reinitiation_weight_data, by=c("patid", "t"))
  
  # fill each of the weights downwards
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    fill(
      stabilized_weight_no_discontinuation,
      unstabilized_weight_no_discontinuation,
      stabilized_weight_no_dose_reduction,
      unstabilized_weight_no_dose_reduction,
      stabilized_weight_reinitiation,
      unstabilized_weight_reinitiation,
      .direction = "down"
    ) %>%
    ungroup()
  
  # fill in the rest in with 1s
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                         stabilized_weight_no_discontinuation),
           unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                           unstabilized_weight_no_discontinuation),
           stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                        stabilized_weight_no_dose_reduction),
           unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                          unstabilized_weight_no_dose_reduction),
           stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                   stabilized_weight_reinitiation),
           unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                     unstabilized_weight_reinitiation))
  
  # multiple all of the weights together
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
             stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
           unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
             unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
  
  # truncate the stabilized weights at the 99th percentile
  quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
  print("summary of truncated stabilized weights for tapering arm:")
  print(summary(tapering_arm$stabilized_weight))
  
  # truncate the unstabilized weights at the 99th percentile
  upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
  
  tapering_arm <- tapering_arm %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
  print("summary of truncated unstabilized weights for tapering arm:")
  print(summary(tapering_arm$unstabilized_weight))
  
  tapering_arm <- tapering_arm %>%
    mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                               1, 0)) %>%
    mutate(censored = ifelse(uncensored == 0, 1, 0))   
  
  vars <- c("stabilized_weight", "stabilized_weight_untruncated")
  
  # Function to compute the summary stats
  get_summary_stats <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p1 = quantile(x, 0.01, na.rm = TRUE),
      p50 = quantile(x, 0.50, na.rm = TRUE),
      p99 = quantile(x, 0.99, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  summary_tapering <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(tapering_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_tapering <- summary_tapering %>%
    mutate(arm = "tapering")
  
  summary_continuation <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(continuation_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_continuation <- summary_continuation %>%
    mutate(arm = "continuation")
  
  summary_weights <- rbind(summary_continuation, summary_tapering)
  
  fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/delirium.csv"))
  
  ################################################################################
  # check covariate balance after grace period
  ################################################################################
  covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                  "cancer_baseline", "copd_baseline", "depression_baseline",
                  "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                  "heart_hypertension_baseline", "kidney_disease_baseline",
                  "ace_inhibitors_baseline", "antidepressants_baseline", 
                  "antiplatelets_baseline", "anxiolytics_baseline",
                  "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                  "dementia_drugs_baseline", "diabetes_drug_baseline",
                  "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                  "oral_anticoagulants_baseline", "osteoporosis_baseline",
                  "venous_thromboembolism_baseline", 
                  "rheumatoid_arthritis_baseline", "mi_baseline",
                  "months_since_diagnosis", "age_T0", "gender")
  
  
  continuation_t7 <- continuation_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  tapering_t7 <- tapering_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  
  combined_t7 <- rbind(continuation_t7, tapering_t7)
  
  check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
  
  # Calculate the balance table for the subset
  balance_table <- svyCreateTableOne(vars = covariates, 
                                     strata = "treatment", 
                                     data = check_weighting, 
                                     test = FALSE, 
                                     smd = TRUE)
  
  table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  
  table_one_df$variable <- rownames(table_one_df)
  
  # balance table before weighting
  pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
  
  # Calculate the balance table for the subset
  pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                strata = "treatment", 
                                                data = pre_weighting_smd, 
                                                test = FALSE, 
                                                smd = TRUE)
  
  
  pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  pre_weighting_table$variable <- rownames(pre_weighting_table)
  
  
  if(full_sample == TRUE){
    fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/delirium.csv"))
    fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_delirium.csv"))
  }
  
  ################################################################################
  # Outcome model
  ################################################################################
  
  continuation_df <- continuation_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  tapering_df <- tapering_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  combined_data <- rbind(continuation_df, tapering_df) %>%
    filter(censored == 0)
  
  
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/delirium.csv" ))
  
  # patient years by arm
  py_by_arm <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      # each row = 1 month, so 1/12 year
      patient_years = ceiling(sum(stabilized_weight * (1/12))),
      .groups = "drop"
    )
  print(py_by_arm)
  
  # number of outcomes
  weighted_events <- with(
    combined_data,
    tapply(stabilized_weight * outcome, treatment, sum)
  )
  print(weighted_events)
  
  crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
  # to create a survival curve, create a model with an interaction between t and treatment
  # and t_squared and treatment
  survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                + t*treatment + t_squared*treatment + 
                                + gender + age_T0 + months_since_diagnosis +
                                # baseline conditions
                                angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                venous_thromboembolism_baseline + mi_baseline +
                                # baseline medications
                                ace_inhibitors_baseline + antidepressants_baseline +
                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                              data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
  
  continuation_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "continuation")
  
  tapering_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "tapering")
  
  # create hazards under each strategy
  tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
  continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
  
  # calculate survival from the cumulative of 1-hazard
  tapering_df <- tapering_df %>%
    group_by(patid) %>%
    mutate(survival1 = cumprod(1-hazard1)) %>%
    ungroup()
  
  continuation_df <- continuation_df %>%
    group_by(patid) %>%
    mutate(survival0 = cumprod(1-hazard0)) %>%
    ungroup()
  
  # estimate risks as 1-survival
  tapering_df <- tapering_df %>%
    mutate(risk1 = 1-survival1)
  
  continuation_df <- continuation_df %>%
    mutate(risk0 = 1-survival0)
  
  # calculate mean for each time point
  tapering_risk <- tapering_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk1),
              survival = mean(survival1),
              hazard = mean(hazard1))
  
  continuation_risk <- continuation_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk0),
              survival = mean(survival0),
              hazard = mean(hazard0))
  
  risk_curves <- full_join(tapering_risk, continuation_risk)
  
  risk_curves_wider <- risk_curves %>%
    pivot_wider(
      id_cols = t,
      names_from = treatment,
      values_from = c(risk, hazard, survival),
      names_sep = "_"
    )
  
  risk_curves_wider <- risk_curves_wider %>%
    mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
  
  if(full_sample == TRUE){
    fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/delirium.csv"))
  }
  
  # hazard ratio is average of hazard ratios over time
  hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
  
  tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
  continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
  
  tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
  continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
  
  # calculate risk ratios and differences
  risk_ratio_12 <- tapering_risk_12/continuation_risk_12
  risk_diff_12 <- tapering_risk_12-continuation_risk_12
  
  risk_ratio_24 <- tapering_risk_24/continuation_risk_24
  risk_diff_24 <- tapering_risk_24-continuation_risk_24
  
  
  risk_df <- data.frame(
    estimate_hr = hazard_ratio_estimate,
    risk_ratio_12months = risk_ratio_12,
    risk_ratio_24months = risk_ratio_24,
    risk_diff_12months = risk_diff_12,
    risk_diff_24months = risk_diff_24,
    abs_risk_continuation_12 = continuation_risk_12,
    abs_risk_continuation_24 = continuation_risk_24,
    abs_risk_tapering_12 = tapering_risk_12, 
    abs_risk_tapering_24 = tapering_risk_24
  )
  rownames(risk_df) <- NULL
  
  if(full_sample == TRUE){
    fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/delirium.csv"))
  }
  
  risk_curves <- risk_curves %>%
    mutate(risk=risk*100)
  
  risk_curves <- rbind(
    risk_curves,
    data.frame(t = 0, risk = 0, treatment = "continuation"),
    data.frame(t = 0, risk = 0, treatment = "tapering")
  )
  
  # Plot the survival curves
  survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
    geom_line(linewidth=1)+
    labs(
      title = "Delirium Cumulative Incidence Curves by Treatment",
      x = "Time (months)",
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),     # tick label size
      axis.title = element_text(size = 16),    # axis title size
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    # vertical line at end of grace period
    geom_vline(xintercept = 6, 
               linetype = "dashed", 
               linewidth = 1,
               color = "black") 
  
  if(full_sample == TRUE){
    ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/delirium.png"), plot = survival_curve_plot,
           dpi=300, width = 8, height = 6)
  }
  
  return(risk_df)
  
}

tapering_delirium_trials(combined_12weeks, '12weeks', full_sample = TRUE)

####### delirium, issues with covariate balance
tapering_delirium_trials_24weeks(combined_24weeks, '24weeks', full_sample = TRUE)


################################################################################
# DEATH OUTCOME
################################################################################

#### death outcome analysis


tapering_death_trials_12weeks <- function(long_format, time_period, full_sample=FALSE) {
  
  long_format <- long_format %>%
    mutate(outcome_date = ifelse(deathdate_adjusted > enddate & !is.na(enddate), NA, deathdate_adjusted)) %>%
    mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
    mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
    mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
    mutate(outcome = ifelse(is.na(outcome), 0, outcome))
  
  # add discontinuation date
  long_format <- long_format %>%
    mutate(discontinuation_date = treatment_enddate + 28)
  
  long_format <- long_format %>%
    mutate(discontinuation_date = case_when(
      (!is.na(enddate) & discontinuation_date > enddate) |
        (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
      TRUE ~ discontinuation_date
    ))
  
  
  # create a variable for maximum followup
  long_format <- long_format %>%
    group_by(patid) %>%
    mutate(maximum_followup = min(deathdate_adjusted, enddate, na.rm=TRUE))
  
  # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
  
  long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                    long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                     long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                  long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
  
  
  
  long_format <- long_format %>%
    mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
    mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
    mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
  
  # assuming all available followup so take from the overall long format
  long_format <- long_format %>%
    mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
    mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
  
  ################################################################################
  # Continuation arm
  ################################################################################
  # continuation clone are only censored due to discontinuation
  continuation_arm <- long_format %>%
    mutate(treatment = "continuation") %>%
    mutate(censoring_date = discontinuation_date) %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
    mutate(censored = ifelse(is.na(censored), 0, censored))
  
  # end of followup is the earliest of censoring date or outcome date or enddate
  continuation_arm <- continuation_arm %>%
    mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
    mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
    mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  continuation_arm <- continuation_arm %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                              censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                            0, outcome))
  
  # in cases where censoring date is after the outcome date, set censor month to 0 for that month
  # this is because  followup ends at that point so ignore the censoring event
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                               censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                             0, censored))
  
  # keep rows with followup set to 0
  continuation_arm <- continuation_arm %>%
    mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
    filter(followup == 1)
  
  # create stabilized IPCW weights for the continuation arm
  ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                 # baseline conditions
                                 angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                 cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                 heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                 kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                 venous_thromboembolism_baseline + mi_baseline + 
                                 # time-varying conditions
                                 angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                 diabetes + heartdisease + heart_failure + heart_hypertension +
                                 kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                 # baseline medications
                                 ace_inhibitors_baseline  + antidepressants_baseline +
                                 antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                 arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                 dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                 lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                 # time-varying medications
                                 ace_inhibitors + antidepressants + antiepileptics + 
                                 antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                 dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                 oral_anticoagulants,
                               data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  
  continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
  
  # get the cumulative probability of remaining uncensored
  continuation_arm <- continuation_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
  print("continuation arm unstabilized weights")
  print(summary(continuation_arm$unstabilized_weight))
  
  baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                          # baseline conditions
                                          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                          venous_thromboembolism_baseline + mi_baseline + 
                                          # baseline medications
                                          ace_inhibitors_baseline + antidepressants_baseline +
                                          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                        data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    group_by(patid) %>%
    mutate(cum_numerator = cumprod(numerator))
  
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
  
  print("continuation arm stabilized weights")
  print(summary(continuation_arm$stabilized_weight))
  
  # truncate the weights at 99th percentile based on uncensored population
  continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
  upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight) %>%
    mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
  
  lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
  
  ################################################################################
  # Tapering arm
  ################################################################################
  # censor at earliest of:
  # 1. reinitiation
  # 2. end of grace period if failure to discontinue
  # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
  tapering_arm <- long_format %>%
    mutate(treatment = "tapering") %>%
    group_by(patid) %>%
    mutate(discontinued = cummax(discontinue)) %>%
    dplyr::select(-discontinue) %>%
    ungroup()
  
  
  # filter out discontinuation dates and dose reduction dates that are not within grace period
  tapering_arm <- tapering_arm %>%
    mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
           dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
  
  # create variables for 'no discontinuation date' and no dose reduction date
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
    fill(no_discontinuation_date, .direction="down") %>%
    fill(no_discontinuation_date, .direction="up") %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                     is.na(dose_reduction_date), discontinuation_date, NA)))
  
  tapering_arm <- tapering_arm %>%
    mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
  
  # followup is earliest of outcome date, enddate, censoring date
  tapering_arm <- tapering_arm %>%
    mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
    mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
    mutate(followup = ifelse(t <= followup_end, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  tapering_arm <- tapering_arm %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                            & !is.na(outcome_month), 0, outcome))
  
  # create binary variables for each of the reasons for censoring
  tapering_arm <- tapering_arm %>%
    mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
    mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
    mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
    mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
    mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
    mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
    mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
           no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
           reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
  
  # if reinitiation censor is in the same month as another censor, set to 0 as the person
  # would be censored by the other two reasons first
  tapering_arm <- tapering_arm %>%
    mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                          reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
  
  
  # in cases where the outcome occurs in the same month as the censor and the outcome
  # is BEFORE censoring, set the censor variable to 0
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                !is.na(outcome_month) & !is.na(censoring_month) &
                                                outcome_date < censoring_date & t == outcome_month,
                                              0, no_discontinuation_censor)) %>%
    mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                               !is.na(outcome_month) & !is.na(censoring_month) &
                                               outcome_date < censoring_date & t == outcome_month,
                                             0, no_dose_reduction_censor)) %>%
    mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                          !is.na(outcome_month) & !is.na(censoring_month) &
                                          outcome_date < censoring_date & t == outcome_month,
                                        0, reinitiation_censor))
  
  # keep data where followup == 1
  tapering_arm <- tapering_arm %>%
    filter(followup == 1)
  
  # first set of weights, failure to discontinuation, which only applies to t == 7
  grace_period_data <- tapering_arm %>%
    filter(t == 7) %>%
    mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
  
  ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + cancer + depression +
                                    diabetes + heartdisease + heart_failure + 
                                    kidney_disease + osteoporosis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = grace_period_data, family=quasibinomial())
  
  # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
  
  grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
  
  # as there is only one time point, no need to do cumultive probabilities
  grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
  print("discontinuation: no discontinuation unstabilized weights")
  print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
  
  baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                              # baseline conditions
                                              angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                              cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                              heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                              kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                              venous_thromboembolism_baseline + mi_baseline + 
                                              # baseline medications
                                              ace_inhibitors_baseline + antidepressants_baseline +
                                              antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                              arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                              dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                              lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                            data = grace_period_data, family=quasibinomial())
  
  grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
  grace_period_data <- grace_period_data %>%
    mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
  print("discontinuation: no discontinuation stabilized weights")
  print(summary(grace_period_data$stabilized_weight_no_discontinuation))
  
  grace_period_data <- grace_period_data %>% 
    dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
  tapering_arm <- tapering_arm %>%
    left_join(grace_period_data, by=c("patid", "t"))
  
  ## no dose reduction weights which should only be applied to those who have not discontinued
  # those who have discontinued without being censored cannot be censored due to no dose reduction
  # must also only occur during grace period
  not_discontinued_data <- tapering_arm %>%
    filter(t < 7) %>%
    filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
  
  
  ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                      # baseline conditions
                                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                      venous_thromboembolism_baseline + mi_baseline + 
                                                      # time-varying conditions
                                                      angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                      diabetes + heartdisease + heart_failure + heart_hypertension +
                                                      kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                      # baseline medications
                                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                      # time-varying medications
                                                      ace_inhibitors + antidepressants + antiepileptics + 
                                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                      oral_anticoagulants,
                                                    data = not_discontinued_data, family=quasibinomial())
  # all covariates have reasonable standard errors apart from rheumatoid arthritis
  
  not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
  
  # need to calculate the cumulative probabilities of the numerator and denominator
  # get the cumulative probability of remaining uncensored
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
  print("discontinuation no dose reduction unstabilized weights")
  print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
  
  baseline_ps_model_discontinuation_no_dose_reduction <-  
    glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
          # baseline conditions
          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
          venous_thromboembolism_baseline + mi_baseline + 
          # baseline medications
          ace_inhibitors_baseline + antidepressants_baseline +
          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
        data = not_discontinued_data, family=quasibinomial())
  
  
  not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
  
  # get the cumulative probabilities for the numerator
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
    ungroup()
  
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
  print("discontinuation: no dose reduction stabilized weights")
  print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
  
  not_discontinued_data <- not_discontinued_data %>%
    dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
  
  tapering_arm <- tapering_arm %>%
    left_join(not_discontinued_data, by=c("patid", "t"))
  
  # reinitiation weights, which are only applied to those who have discontinued  
  reinitiation_weight_data <- tapering_arm %>%
    filter(discontinued == 1)
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
  
  reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + depression +
                                    diabetes + kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = reinitiation_weight_data, family=quasibinomial())
  
  # very high standard errors for angina, cancer, copd, heart disease, heart failure, hypertension
  reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
  
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
  print("reinitiation weights unstabilized")
  print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
  
  reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                  cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                  kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline + 
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
  print("reinitiation weights stabilized")
  print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
  
  tapering_arm <- tapering_arm %>%
    left_join(reinitiation_weight_data, by=c("patid", "t"))
  
  # fill each of the weights downwards
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    fill(
      stabilized_weight_no_discontinuation,
      unstabilized_weight_no_discontinuation,
      stabilized_weight_no_dose_reduction,
      unstabilized_weight_no_dose_reduction,
      stabilized_weight_reinitiation,
      unstabilized_weight_reinitiation,
      .direction = "down"
    ) %>%
    ungroup()
  
  # fill in the rest in with 1s
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                         stabilized_weight_no_discontinuation),
           unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                           unstabilized_weight_no_discontinuation),
           stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                        stabilized_weight_no_dose_reduction),
           unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                          unstabilized_weight_no_dose_reduction),
           stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                   stabilized_weight_reinitiation),
           unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                     unstabilized_weight_reinitiation))
  
  # multiple all of the weights together
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
             stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
           unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
             unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
  
  # truncate the stabilized weights at the 99th percentile
  quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
  print("summary of truncated stabilized weights for tapering arm:")
  print(summary(tapering_arm$stabilized_weight))
  
  # truncate the unstabilized weights at the 99th percentile
  upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
  
  tapering_arm <- tapering_arm %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
  print("summary of truncated unstabilized weights for tapering arm:")
  print(summary(tapering_arm$unstabilized_weight))
  
  tapering_arm <- tapering_arm %>%
    mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                               1, 0)) %>%
    mutate(censored = ifelse(uncensored == 0, 1, 0))   
  
  vars <- c("stabilized_weight", "stabilized_weight_untruncated")
  
  # Function to compute the summary stats
  get_summary_stats <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p1 = quantile(x, 0.01, na.rm = TRUE),
      p50 = quantile(x, 0.50, na.rm = TRUE),
      p99 = quantile(x, 0.99, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  summary_tapering <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(tapering_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_tapering <- summary_tapering %>%
    mutate(arm = "tapering")
  
  summary_continuation <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(continuation_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_continuation <- summary_continuation %>%
    mutate(arm = "continuation")
  
  summary_weights <- rbind(summary_continuation, summary_tapering)
  
  fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/death.csv"))
  
  ################################################################################
  # check covariate balance after grace period
  ################################################################################
  covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                  "cancer_baseline", "copd_baseline", "depression_baseline",
                  "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                  "heart_hypertension_baseline", "kidney_disease_baseline",
                  "ace_inhibitors_baseline", "antidepressants_baseline", 
                  "antiplatelets_baseline", "anxiolytics_baseline",
                  "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                  "dementia_drugs_baseline", "diabetes_drug_baseline",
                  "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                  "oral_anticoagulants_baseline", "osteoporosis_baseline",
                  "venous_thromboembolism_baseline", 
                  "rheumatoid_arthritis_baseline", "mi_baseline",
                  "months_since_diagnosis", "age_T0", "gender")
  
  
  continuation_t7 <- continuation_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  tapering_t7 <- tapering_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  
  combined_t7 <- rbind(continuation_t7, tapering_t7)
  
  check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
  
  # Calculate the balance table for the subset
  balance_table <- svyCreateTableOne(vars = covariates, 
                                     strata = "treatment", 
                                     data = check_weighting, 
                                     test = FALSE, 
                                     smd = TRUE)
  
  table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  
  table_one_df$variable <- rownames(table_one_df)
  
  # balance table before weighting
  pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
  
  # Calculate the balance table for the subset
  pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                strata = "treatment", 
                                                data = pre_weighting_smd, 
                                                test = FALSE, 
                                                smd = TRUE)
  
  
  pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  pre_weighting_table$variable <- rownames(pre_weighting_table)
  
  
  if(full_sample == TRUE){
    fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/death.csv"))
    fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_death.csv"))
  }
  
  ################################################################################
  # Outcome model
  ################################################################################
  
  continuation_df <- continuation_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  tapering_df <- tapering_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  combined_data <- rbind(continuation_df, tapering_df) %>%
    filter(censored == 0)
  
  
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/death.csv" ))
  
  
  # patient years by arm
  py_by_arm <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      # each row = 1 month, so 1/12 year
      patient_years = ceiling(sum(stabilized_weight * (1/12))),
      .groups = "drop"
    )
  print(py_by_arm)
  
  # number of outcomes
  weighted_events <- with(
    combined_data,
    tapply(stabilized_weight * outcome, treatment, sum)
  )
  print(weighted_events)
  
  crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
  # to create a survival curve, create a model with an interaction between t and treatment
  # and t_squared and treatment
  survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                + t*treatment + t_squared*treatment + 
                                + gender + age_T0 + months_since_diagnosis +
                                # baseline conditions
                                angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                venous_thromboembolism_baseline + mi_baseline +
                                # baseline medications
                                ace_inhibitors_baseline + antidepressants_baseline +
                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                              data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
  
  continuation_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "continuation")
  
  tapering_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "tapering")
  
  # create hazards under each strategy
  tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
  continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
  
  # calculate survival from the cumulative of 1-hazard
  tapering_df <- tapering_df %>%
    group_by(patid) %>%
    mutate(survival1 = cumprod(1-hazard1)) %>%
    ungroup()
  
  continuation_df <- continuation_df %>%
    group_by(patid) %>%
    mutate(survival0 = cumprod(1-hazard0)) %>%
    ungroup()
  
  # estimate risks as 1-survival
  tapering_df <- tapering_df %>%
    mutate(risk1 = 1-survival1)
  
  continuation_df <- continuation_df %>%
    mutate(risk0 = 1-survival0)
  
  # calculate mean for each time point
  tapering_risk <- tapering_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk1),
              survival = mean(survival1),
              hazard = mean(hazard1))
  
  continuation_risk <- continuation_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk0),
              survival = mean(survival0),
              hazard = mean(hazard0))
  
  risk_curves <- full_join(tapering_risk, continuation_risk)
  
  risk_curves_wider <- risk_curves %>%
    pivot_wider(
      id_cols = t,
      names_from = treatment,
      values_from = c(risk, hazard, survival),
      names_sep = "_"
    )
  
  risk_curves_wider <- risk_curves_wider %>%
    mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
  
  if(full_sample == TRUE){
    fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/death.csv"))
  }
  
  # hazard ratio is average of hazard ratios over time
  hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
  
  tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
  continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
  
  tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
  continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
  
  # calculate risk ratios and differences
  risk_ratio_12 <- tapering_risk_12/continuation_risk_12
  risk_diff_12 <- tapering_risk_12-continuation_risk_12
  
  risk_ratio_24 <- tapering_risk_24/continuation_risk_24
  risk_diff_24 <- tapering_risk_24-continuation_risk_24
  
  
  risk_df <- data.frame(
    estimate_hr = hazard_ratio_estimate,
    risk_ratio_12months = risk_ratio_12,
    risk_ratio_24months = risk_ratio_24,
    risk_diff_12months = risk_diff_12,
    risk_diff_24months = risk_diff_24,
    abs_risk_continuation_12 = continuation_risk_12,
    abs_risk_continuation_24 = continuation_risk_24,
    abs_risk_tapering_12 = tapering_risk_12, 
    abs_risk_tapering_24 = tapering_risk_24
  )
  rownames(risk_df) <- NULL
  
  
  if(full_sample == TRUE){
    fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/death.csv"))
  }
  
  risk_curves <- risk_curves %>%
    mutate(risk=risk*100)
  
  risk_curves <- rbind(
    risk_curves,
    data.frame(t = 0, risk = 0, treatment = "continuation"),
    data.frame(t = 0, risk = 0, treatment = "tapering")
  )
  
  # Plot the survival curves
  survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
    geom_line(linewidth=1)+
    labs(
      title = "Death Cumulative Incidence Curves by Treatment",
      x = "Time (months)",
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),     # tick label size
      axis.title = element_text(size = 16),    # axis title size
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    # vertical line at end of grace period
    geom_vline(xintercept = 6, 
               linetype = "dashed", 
               linewidth = 1,
               color = "black") 
  
  if(full_sample == TRUE){
    ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/death.png"), plot = survival_curve_plot,
           dpi=300, width = 8, height = 6)
  }
  
  return(risk_df)
  
}


tapering_death_trials_24weeks <- function(long_format, time_period, full_sample=FALSE) {
  
  long_format <- long_format %>%
    mutate(outcome_date = ifelse(deathdate_adjusted > enddate & !is.na(enddate), NA, deathdate_adjusted)) %>%
    mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
    mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
    mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
    mutate(outcome = ifelse(is.na(outcome), 0, outcome))
  
  # add discontinuation date
  long_format <- long_format %>%
    mutate(discontinuation_date = treatment_enddate + 28)
  
  long_format <- long_format %>%
    mutate(discontinuation_date = case_when(
      (!is.na(enddate) & discontinuation_date > enddate) |
        (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
      TRUE ~ discontinuation_date
    ))
  
  
  # create a variable for maximum followup
  long_format <- long_format %>%
    group_by(patid) %>%
    mutate(maximum_followup = min(deathdate_adjusted, enddate, na.rm=TRUE))
  
  # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
  
  long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                    long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                     long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                  long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
  
  
  
  long_format <- long_format %>%
    mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
    mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
    mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
  
  # assuming all available followup so take from the overall long format
  long_format <- long_format %>%
    mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
    mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
  
  ################################################################################
  # Continuation arm
  ################################################################################
  # continuation clone are only censored due to discontinuation
  continuation_arm <- long_format %>%
    mutate(treatment = "continuation") %>%
    mutate(censoring_date = discontinuation_date) %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
    mutate(censored = ifelse(is.na(censored), 0, censored))
  
  # end of followup is the earliest of censoring date or outcome date or enddate
  continuation_arm <- continuation_arm %>%
    mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
    mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
    mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  continuation_arm <- continuation_arm %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                              censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                            0, outcome))
  
  # in cases where censoring date is after the outcome date, set censor month to 0 for that month
  # this is because  followup ends at that point so ignore the censoring event
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                               censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                             0, censored))
  
  # keep rows with followup set to 0
  continuation_arm <- continuation_arm %>%
    mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
    filter(followup == 1)
  
  # create stabilized IPCW weights for the continuation arm
  ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                 # baseline conditions
                                 angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                 cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                 heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                 kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                 venous_thromboembolism_baseline + mi_baseline + 
                                 # time-varying conditions
                                 angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                 diabetes + heartdisease + heart_failure + heart_hypertension +
                                 kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                 # baseline medications
                                 ace_inhibitors_baseline  + antidepressants_baseline +
                                 antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                 arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                 dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                 lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                 # time-varying medications
                                 ace_inhibitors + antidepressants + antiepileptics + 
                                 antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                 dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                 oral_anticoagulants,
                               data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  
  continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
  
  # get the cumulative probability of remaining uncensored
  continuation_arm <- continuation_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
  print("continuation arm unstabilized weights")
  print(summary(continuation_arm$unstabilized_weight))
  
  baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                          # baseline conditions
                                          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                          venous_thromboembolism_baseline + mi_baseline + 
                                          # baseline medications
                                          ace_inhibitors_baseline + antidepressants_baseline +
                                          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                        data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    group_by(patid) %>%
    mutate(cum_numerator = cumprod(numerator))
  
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
  
  print("continuation arm stabilized weights")
  print(summary(continuation_arm$stabilized_weight))
  
  # truncate the weights at 99th percentile based on uncensored population
  continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
  upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight) %>%
    mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
  
  lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
  
  ################################################################################
  # Tapering arm
  ################################################################################
  # censor at earliest of:
  # 1. reinitiation
  # 2. end of grace period if failure to discontinue
  # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
  tapering_arm <- long_format %>%
    mutate(treatment = "tapering") %>%
    group_by(patid) %>%
    mutate(discontinued = cummax(discontinue)) %>%
    dplyr::select(-discontinue) %>%
    ungroup()
  
  
  # filter out discontinuation dates and dose reduction dates that are not within grace period
  tapering_arm <- tapering_arm %>%
    mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
           dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
  
  # create variables for 'no discontinuation date' and no dose reduction date
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
    fill(no_discontinuation_date, .direction="down") %>%
    fill(no_discontinuation_date, .direction="up") %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                     is.na(dose_reduction_date), discontinuation_date, NA)))
  
  tapering_arm <- tapering_arm %>%
    mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
  
  # followup is earliest of outcome date, enddate, censoring date
  tapering_arm <- tapering_arm %>%
    mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
    mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
    mutate(followup = ifelse(t <= followup_end, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  tapering_arm <- tapering_arm %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                            & !is.na(outcome_month), 0, outcome))
  
  # create binary variables for each of the reasons for censoring
  tapering_arm <- tapering_arm %>%
    mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
    mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
    mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
    mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
    mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
    mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
    mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
           no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
           reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
  
  # if reinitiation censor is in the same month as another censor, set to 0 as the person
  # would be censored by the other two reasons first
  tapering_arm <- tapering_arm %>%
    mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                          reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
  
  
  # in cases where the outcome occurs in the same month as the censor and the outcome
  # is BEFORE censoring, set the censor variable to 0
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                !is.na(outcome_month) & !is.na(censoring_month) &
                                                outcome_date < censoring_date & t == outcome_month,
                                              0, no_discontinuation_censor)) %>%
    mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                               !is.na(outcome_month) & !is.na(censoring_month) &
                                               outcome_date < censoring_date & t == outcome_month,
                                             0, no_dose_reduction_censor)) %>%
    mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                          !is.na(outcome_month) & !is.na(censoring_month) &
                                          outcome_date < censoring_date & t == outcome_month,
                                        0, reinitiation_censor))
  
  # keep data where followup == 1
  tapering_arm <- tapering_arm %>%
    filter(followup == 1)
  
  # first set of weights, failure to discontinuation, which only applies to t == 7
  grace_period_data <- tapering_arm %>%
    filter(t == 7) %>%
    mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
  
  ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    cancer + heart_failure + 
                                    kidney_disease +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = grace_period_data, family=quasibinomial())
  
  # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
  
  grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
  
  # as there is only one time point, no need to do cumultive probabilities
  grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
  print("discontinuation: no discontinuation unstabilized weights")
  print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
  
  baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                              # baseline conditions
                                              angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                              cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                              heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                              kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                              venous_thromboembolism_baseline + mi_baseline + 
                                              # baseline medications
                                              ace_inhibitors_baseline + antidepressants_baseline +
                                              antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                              arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                              dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                              lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                            data = grace_period_data, family=quasibinomial())
  
  grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
  grace_period_data <- grace_period_data %>%
    mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
  print("discontinuation: no discontinuation stabilized weights")
  print(summary(grace_period_data$stabilized_weight_no_discontinuation))
  
  grace_period_data <- grace_period_data %>% 
    dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
  tapering_arm <- tapering_arm %>%
    left_join(grace_period_data, by=c("patid", "t"))
  
  ## no dose reduction weights which should only be applied to those who have not discontinued
  # those who have discontinued without being censored cannot be censored due to no dose reduction
  # must also only occur during grace period
  not_discontinued_data <- tapering_arm %>%
    filter(t < 7) %>%
    filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
  
  
  ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                      # baseline conditions
                                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                      venous_thromboembolism_baseline + mi_baseline + 
                                                      # time-varying conditions
                                                      angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                      diabetes + heartdisease + heart_failure + heart_hypertension +
                                                      kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                      # baseline medications
                                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                      # time-varying medications
                                                      ace_inhibitors + antidepressants + antiepileptics + 
                                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                      oral_anticoagulants,
                                                    data = not_discontinued_data, family=quasibinomial())
  # all covariates have reasonable standard errors apart from rheumatoid arthritis
  
  not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
  
  # need to calculate the cumulative probabilities of the numerator and denominator
  # get the cumulative probability of remaining uncensored
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
  print("discontinuation no dose reduction unstabilized weights")
  print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
  
  baseline_ps_model_discontinuation_no_dose_reduction <-  
    glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
          # baseline conditions
          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
          venous_thromboembolism_baseline + mi_baseline + 
          # baseline medications
          ace_inhibitors_baseline + antidepressants_baseline +
          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
        data = not_discontinued_data, family=quasibinomial())
  
  
  not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
  
  # get the cumulative probabilities for the numerator
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
    ungroup()
  
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
  print("discontinuation: no dose reduction stabilized weights")
  print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
  
  not_discontinued_data <- not_discontinued_data %>%
    dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
  
  tapering_arm <- tapering_arm %>%
    left_join(not_discontinued_data, by=c("patid", "t"))
  
  # reinitiation weights, which are only applied to those who have discontinued  
  reinitiation_weight_data <- tapering_arm %>%
    filter(discontinued == 1)
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
  
  reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + depression +
                                    diabetes + kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = reinitiation_weight_data, family=quasibinomial())
  
  # very high standard errors for angina, cancer, copd, heart disease, heart failure, hypertension
  reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
  
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
  print("reinitiation weights unstabilized")
  print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
  
  reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                  cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                  kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline + 
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
  print("reinitiation weights stabilized")
  print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
  
  tapering_arm <- tapering_arm %>%
    left_join(reinitiation_weight_data, by=c("patid", "t"))
  
  # fill each of the weights downwards
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    fill(
      stabilized_weight_no_discontinuation,
      unstabilized_weight_no_discontinuation,
      stabilized_weight_no_dose_reduction,
      unstabilized_weight_no_dose_reduction,
      stabilized_weight_reinitiation,
      unstabilized_weight_reinitiation,
      .direction = "down"
    ) %>%
    ungroup()
  
  # fill in the rest in with 1s
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                         stabilized_weight_no_discontinuation),
           unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                           unstabilized_weight_no_discontinuation),
           stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                        stabilized_weight_no_dose_reduction),
           unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                          unstabilized_weight_no_dose_reduction),
           stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                   stabilized_weight_reinitiation),
           unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                     unstabilized_weight_reinitiation))
  
  # multiple all of the weights together
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
             stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
           unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
             unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
  
  # truncate the stabilized weights at the 99th percentile
  quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
  print("summary of truncated stabilized weights for tapering arm:")
  print(summary(tapering_arm$stabilized_weight))
  
  # truncate the unstabilized weights at the 99th percentile
  upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
  
  tapering_arm <- tapering_arm %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
  print("summary of truncated unstabilized weights for tapering arm:")
  print(summary(tapering_arm$unstabilized_weight))
  
  tapering_arm <- tapering_arm %>%
    mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                               1, 0)) %>%
    mutate(censored = ifelse(uncensored == 0, 1, 0))   
  
  vars <- c("stabilized_weight", "stabilized_weight_untruncated")
  
  # Function to compute the summary stats
  get_summary_stats <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p1 = quantile(x, 0.01, na.rm = TRUE),
      p50 = quantile(x, 0.50, na.rm = TRUE),
      p99 = quantile(x, 0.99, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  summary_tapering <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(tapering_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_tapering <- summary_tapering %>%
    mutate(arm = "tapering")
  
  summary_continuation <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(continuation_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_continuation <- summary_continuation %>%
    mutate(arm = "continuation")
  
  summary_weights <- rbind(summary_continuation, summary_tapering)
  
  fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/death.csv"))
  
  ################################################################################
  # check covariate balance after grace period
  ################################################################################
  covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                  "cancer_baseline", "copd_baseline", "depression_baseline",
                  "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                  "heart_hypertension_baseline", "kidney_disease_baseline",
                  "ace_inhibitors_baseline", "antidepressants_baseline", 
                  "antiplatelets_baseline", "anxiolytics_baseline",
                  "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                  "dementia_drugs_baseline", "diabetes_drug_baseline",
                  "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                  "oral_anticoagulants_baseline", "osteoporosis_baseline",
                  "venous_thromboembolism_baseline", 
                  "rheumatoid_arthritis_baseline", "mi_baseline",
                  "months_since_diagnosis", "age_T0", "gender")
  
  
  continuation_t7 <- continuation_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  tapering_t7 <- tapering_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  
  combined_t7 <- rbind(continuation_t7, tapering_t7)
  
  check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
  
  # Calculate the balance table for the subset
  balance_table <- svyCreateTableOne(vars = covariates, 
                                     strata = "treatment", 
                                     data = check_weighting, 
                                     test = FALSE, 
                                     smd = TRUE)
  
  table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  
  table_one_df$variable <- rownames(table_one_df)
  
  # balance table before weighting
  pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
  
  # Calculate the balance table for the subset
  pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                strata = "treatment", 
                                                data = pre_weighting_smd, 
                                                test = FALSE, 
                                                smd = TRUE)
  
  
  pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  pre_weighting_table$variable <- rownames(pre_weighting_table)
  
  
  if(full_sample == TRUE){
    fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/death.csv"))
    fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_death.csv"))
  }
  
  ################################################################################
  # Outcome model
  ################################################################################
  
  continuation_df <- continuation_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  tapering_df <- tapering_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  combined_data <- rbind(continuation_df, tapering_df) %>%
    filter(censored == 0)
  
  
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/death.csv" ))
  
  
  # patient years by arm
  py_by_arm <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      # each row = 1 month, so 1/12 year
      patient_years = ceiling(sum(stabilized_weight * (1/12))),
      .groups = "drop"
    )
  print(py_by_arm)
  
  # number of outcomes
  weighted_events <- with(
    combined_data,
    tapply(stabilized_weight * outcome, treatment, sum)
  )
  print(weighted_events)
  
  crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
  # to create a survival curve, create a model with an interaction between t and treatment
  # and t_squared and treatment
  survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                + t*treatment + t_squared*treatment + 
                                + gender + age_T0 + months_since_diagnosis +
                                # baseline conditions
                                angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                venous_thromboembolism_baseline + mi_baseline +
                                # baseline medications
                                ace_inhibitors_baseline + antidepressants_baseline +
                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                              data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
  
  continuation_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "continuation")
  
  tapering_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "tapering")
  
  # create hazards under each strategy
  tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
  continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
  
  # calculate survival from the cumulative of 1-hazard
  tapering_df <- tapering_df %>%
    group_by(patid) %>%
    mutate(survival1 = cumprod(1-hazard1)) %>%
    ungroup()
  
  continuation_df <- continuation_df %>%
    group_by(patid) %>%
    mutate(survival0 = cumprod(1-hazard0)) %>%
    ungroup()
  
  # estimate risks as 1-survival
  tapering_df <- tapering_df %>%
    mutate(risk1 = 1-survival1)
  
  continuation_df <- continuation_df %>%
    mutate(risk0 = 1-survival0)
  
  # calculate mean for each time point
  tapering_risk <- tapering_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk1),
              survival = mean(survival1),
              hazard = mean(hazard1))
  
  continuation_risk <- continuation_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk0),
              survival = mean(survival0),
              hazard = mean(hazard0))
  
  risk_curves <- full_join(tapering_risk, continuation_risk)
  
  risk_curves_wider <- risk_curves %>%
    pivot_wider(
      id_cols = t,
      names_from = treatment,
      values_from = c(risk, hazard, survival),
      names_sep = "_"
    )
  
  risk_curves_wider <- risk_curves_wider %>%
    mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
  
  if(full_sample == TRUE){
    fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/death.csv"))
  }
  
  # hazard ratio is average of hazard ratios over time
  hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
  
  tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
  continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
  
  tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
  continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
  
  # calculate risk ratios and differences
  risk_ratio_12 <- tapering_risk_12/continuation_risk_12
  risk_diff_12 <- tapering_risk_12-continuation_risk_12
  
  risk_ratio_24 <- tapering_risk_24/continuation_risk_24
  risk_diff_24 <- tapering_risk_24-continuation_risk_24
  
  
  risk_df <- data.frame(
    estimate_hr = hazard_ratio_estimate,
    risk_ratio_12months = risk_ratio_12,
    risk_ratio_24months = risk_ratio_24,
    risk_diff_12months = risk_diff_12,
    risk_diff_24months = risk_diff_24,
    abs_risk_continuation_12 = continuation_risk_12,
    abs_risk_continuation_24 = continuation_risk_24,
    abs_risk_tapering_12 = tapering_risk_12, 
    abs_risk_tapering_24 = tapering_risk_24
  )
  rownames(risk_df) <- NULL
  
  
  if(full_sample == TRUE){
    fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/death.csv"))
  }
  
  risk_curves <- risk_curves %>%
    mutate(risk=risk*100)
  
  risk_curves <- rbind(
    risk_curves,
    data.frame(t = 0, risk = 0, treatment = "continuation"),
    data.frame(t = 0, risk = 0, treatment = "tapering")
  )
  # Plot the survival curves
  survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
    geom_line(linewidth=1)+
    labs(
      title = "Death Cumulative Incidence Curves by Treatment",
      x = "Time (months)",
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),     # tick label size
      axis.title = element_text(size = 16),    # axis title size
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    # vertical line at end of grace period
    geom_vline(xintercept = 6, 
               linetype = "dashed", 
               linewidth = 1,
               color = "black") 
  if(full_sample == TRUE){
    ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/death.png"), plot = survival_curve_plot,
           dpi=300, width = 8, height = 6)
  }
  
  return(risk_df)
  
}

tapering_death_trials_12weeks(combined_12weeks, '12weeks', full_sample = TRUE)
tapering_death_trials_24weeks(combined_24weeks, '24weeks', full_sample = TRUE)

################################################################################
# NEGATIVE CONTROL OUTCOME
################################################################################

#### negcontrol outcome analysis
tapering_negcontrol_trials <- function(long_format, time_period, full_sample=FALSE) {
  
  long_format <- long_format %>%
    mutate(outcome_date = ifelse(negcontrol_outcome_date > enddate & !is.na(enddate), NA, negcontrol_outcome_date)) %>%
    mutate(outcome_date = as.Date(outcome_date, origin="1970-01-01")) %>%
    mutate(outcome_month = ceiling((as.numeric(outcome_date - T0) + 1)/28)) %>%
    mutate(outcome = ifelse(t == outcome_month, 1, 0)) %>%
    mutate(outcome = ifelse(is.na(outcome), 0, outcome))
  
  # add discontinuation date
  long_format <- long_format %>%
    mutate(discontinuation_date = treatment_enddate + 28)
  
  long_format <- long_format %>%
    mutate(discontinuation_date = case_when(
      (!is.na(enddate) & discontinuation_date > enddate) |
        (!is.na(outcome_date) & discontinuation_date > outcome_date) ~ as.Date(NA),
      TRUE ~ discontinuation_date
    ))
  
  
  # create a variable for maximum followup
  long_format <- long_format %>%
    group_by(patid) %>%
    mutate(maximum_followup = min(negcontrol_outcome_date, enddate, na.rm=TRUE))
  
  # remove dose_reduction_dates, discontinuation_dates, and reinitiation_dates that happen after maximum_followup
  
  long_format$dose_reduction_date[!is.na(long_format$dose_reduction_date) & 
                                    long_format$maximum_followup <= long_format$dose_reduction_date] <- NA # Removes dose reduction dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$discontinuation_date[!is.na(long_format$discontinuation_date) & 
                                     long_format$maximum_followup <= long_format$discontinuation_date] <- NA # Removes discontinuation dates on or after maximum follow-up (should not be included in analyses)
  
  
  long_format$reinitiation_date[!is.na(long_format$reinitiation_date) & 
                                  long_format$maximum_followup <= long_format$reinitiation_date] <- NA # Removes reinitiation dates on or after maximum follow-up (should not be included in analyses)
  
  
  
  long_format <- long_format %>%
    mutate(discontinuation_month = ceiling((as.numeric(discontinuation_date - T0) + 1)/28)) %>%
    mutate(discontinue = ifelse(t == discontinuation_month, 1, 0)) %>%
    mutate(discontinue = ifelse(is.na(discontinue), 0, discontinue))
  
  # assuming all available followup so take from the overall long format
  long_format <- long_format %>%
    mutate(followup_curve = ceiling((as.numeric(maximum_followup - T0) + 1)/28)) %>%
    mutate(followup_curve = ifelse(t <= followup_curve, 1, 0))
  
  ################################################################################
  # Continuation arm
  ################################################################################
  # continuation clone are only censored due to discontinuation
  continuation_arm <- long_format %>%
    mutate(treatment = "continuation") %>%
    mutate(censoring_date = discontinuation_date) %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(censored = ifelse(t == censoring_month, 1, 0)) %>%
    mutate(censored = ifelse(is.na(censored), 0, censored))
  
  # end of followup is the earliest of censoring date or outcome date or enddate
  continuation_arm <- continuation_arm %>%
    mutate(end_of_followup_date = pmin(outcome_date, censoring_date, enddate, na.rm=TRUE)) %>%
    mutate(end_of_followup_month = ceiling((as.numeric(end_of_followup_date - T0) + 1)/28)) %>%
    mutate(followup = ifelse(t <= end_of_followup_month, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  continuation_arm <- continuation_arm %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >=
                              censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                            0, outcome))
  
  # in cases where censoring date is after the outcome date, set censor month to 0 for that month
  # this is because  followup ends at that point so ignore the censoring event
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    mutate(censored = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date <
                               censoring_date & !is.na(censoring_month) & !is.na(outcome_month),
                             0, censored))
  
  # keep rows with followup set to 0
  continuation_arm <- continuation_arm %>%
    mutate(uncensored = ifelse(censored == 1, 0, 1)) %>%
    filter(followup == 1)
  
  # create stabilized IPCW weights for the continuation arm
  ps_model_continuation <- glm(uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                 # baseline conditions
                                 angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                 cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                 heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                 kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                 venous_thromboembolism_baseline + mi_baseline + 
                                 # time-varying conditions
                                 angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                 diabetes + heartdisease + heart_failure + heart_hypertension +
                                 kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                 # baseline medications
                                 ace_inhibitors_baseline  + antidepressants_baseline +
                                 antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                 arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                 dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                 lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                 # time-varying medications
                                 ace_inhibitors + antidepressants + antiepileptics + 
                                 antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                 dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                 oral_anticoagulants,
                               data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  
  continuation_arm$prob_uncensored <- predict(ps_model_continuation, type="response")
  
  # get the cumulative probability of remaining uncensored
  continuation_arm <- continuation_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  continuation_arm$unstabilized_weight <- 1/continuation_arm$cum_prob_uncensored
  print("continuation arm unstabilized weights")
  print(summary(continuation_arm$unstabilized_weight))
  
  baseline_ps_model_continuation <- glm(uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
                                          # baseline conditions
                                          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                          venous_thromboembolism_baseline + mi_baseline + 
                                          # baseline medications
                                          ace_inhibitors_baseline + antidepressants_baseline +
                                          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                        data = continuation_arm, family=quasibinomial())
  
  # checked model summary, all have low standard errors
  continuation_arm$numerator <- predict(baseline_ps_model_continuation, type="response")
  continuation_arm <- continuation_arm %>%
    arrange(patid, t) %>%
    group_by(patid) %>%
    mutate(cum_numerator = cumprod(numerator))
  
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight = cum_numerator/cum_prob_uncensored)
  
  print("continuation arm stabilized weights")
  print(summary(continuation_arm$stabilized_weight))
  
  # truncate the weights at 99th percentile based on uncensored population
  continuation_uncensored <- continuation_arm %>% filter(censored !=0 )
  upper_bound_stab <- quantile(continuation_uncensored$stabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight) %>%
    mutate(stabilized_weight = pmin(stabilized_weight, upper_bound_stab))
  
  lower_bound_unstab <- quantile(continuation_uncensored$unstabilized_weight, probs=0.99, na.rm=TRUE)
  continuation_arm <- continuation_arm %>%
    mutate(unstabilized_weight_untruncated = unstabilized_weight) %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight, lower_bound_unstab))
  
  ################################################################################
  # Tapering arm
  ################################################################################
  # censor at earliest of:
  # 1. reinitiation
  # 2. end of grace period if failure to discontinue
  # 3. discontinuation date is discontinuation date is during grace period and no dose reduction first
  tapering_arm <- long_format %>%
    mutate(treatment = "tapering") %>%
    group_by(patid) %>%
    mutate(discontinued = cummax(discontinue)) %>%
    dplyr::select(-discontinue) %>%
    ungroup()
  
  
  # filter out discontinuation dates and dose reduction dates that are not within grace period
  tapering_arm <- tapering_arm %>%
    mutate(discontinuation_date = as.Date(ifelse(discontinuation_date >= T0 + 168, NA, discontinuation_date)),
           dose_reduction_date = as.Date(ifelse(dose_reduction_date >= T0 + 168, NA, dose_reduction_date)))
  
  # create variables for 'no discontinuation date' and no dose reduction date
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_date = as.Date(ifelse(lag(discontinued) == 0 & t == 7, T0 + 168, as.Date(NA)))) %>%
    fill(no_discontinuation_date, .direction="down") %>%
    fill(no_discontinuation_date, .direction="up") %>%
    ungroup() %>%
    group_by(patid) %>%
    mutate(no_dose_reduction_date = as.Date(ifelse(discontinuation_date < dose_reduction_date |
                                                     is.na(dose_reduction_date), discontinuation_date, NA)))
  
  tapering_arm <- tapering_arm %>%
    mutate(censoring_date = pmin(no_discontinuation_date, no_dose_reduction_date, reinitiation_date, na.rm=TRUE))
  
  # followup is earliest of outcome date, enddate, censoring date
  tapering_arm <- tapering_arm %>%
    mutate(followup_end_date = pmin(outcome_date, enddate, censoring_date, na.rm=TRUE)) %>%
    mutate(followup_end = ceiling((as.numeric(followup_end_date - T0) + 1) / 28)) %>%
    mutate(followup = ifelse(t <= followup_end, 1, 0))
  
  # if outcome date is after the censoring date, set to 0
  tapering_arm <- tapering_arm %>%
    mutate(censoring_month = ceiling((as.numeric(censoring_date - T0) + 1) / 28)) %>%
    mutate(outcome = ifelse(t == censoring_month & outcome_month == censoring_month & outcome_date >= censoring_date & !is.na(censoring_month)
                            & !is.na(outcome_month), 0, outcome))
  
  # create binary variables for each of the reasons for censoring
  tapering_arm <- tapering_arm %>%
    mutate(no_discontinuation_censor = ceiling((as.numeric(no_discontinuation_date - T0) + 1)/ 28)) %>%
    mutate(no_discontinuation_censor = ifelse(t < no_discontinuation_censor, 0, 1)) %>%
    mutate(no_dose_reduction_censor = ceiling((as.numeric(no_dose_reduction_date - T0) + 1)/ 28)) %>%
    mutate(no_dose_reduction_censor = ifelse(t < no_dose_reduction_censor, 0, 1)) %>%
    mutate(reinitiation_censor = ceiling((as.numeric(reinitiation_date - T0) + 1)/ 28)) %>%
    mutate(reinitiation_censor = ifelse(t < reinitiation_censor, 0, 1)) %>%
    mutate(no_discontinuation_censor = ifelse(is.na(no_discontinuation_censor), 0, no_discontinuation_censor),
           no_dose_reduction_censor = ifelse(is.na(no_dose_reduction_censor), 0, no_dose_reduction_censor),
           reinitiation_censor = ifelse(is.na(reinitiation_censor), 0, reinitiation_censor))
  
  # if reinitiation censor is in the same month as another censor, set to 0 as the person
  # would be censored by the other two reasons first
  tapering_arm <- tapering_arm %>%
    mutate(reinitiation_censor = ifelse(reinitiation_censor == 1 & no_discontinuation_censor == 1 |
                                          reinitiation_censor == 1 & no_dose_reduction_censor == 1, 0, reinitiation_censor))
  
  
  # in cases where the outcome occurs in the same month as the censor and the outcome
  # is BEFORE censoring, set the censor variable to 0
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    mutate(no_discontinuation_censor = ifelse(outcome_month == censoring_month & 
                                                !is.na(outcome_month) & !is.na(censoring_month) &
                                                outcome_date < censoring_date & t == outcome_month,
                                              0, no_discontinuation_censor)) %>%
    mutate(no_dose_reduction_censor = ifelse(outcome_month == censoring_month & 
                                               !is.na(outcome_month) & !is.na(censoring_month) &
                                               outcome_date < censoring_date & t == outcome_month,
                                             0, no_dose_reduction_censor)) %>%
    mutate(reinitiation_censor = ifelse(outcome_month == censoring_month & 
                                          !is.na(outcome_month) & !is.na(censoring_month) &
                                          outcome_date < censoring_date & t == outcome_month,
                                        0, reinitiation_censor))
  
  # keep data where followup == 1
  tapering_arm <- tapering_arm %>%
    filter(followup == 1)
  
  # first set of weights, failure to discontinuation, which only applies to t == 7
  grace_period_data <- tapering_arm %>%
    filter(t == 7) %>%
    mutate(remain_uncensored = ifelse(no_discontinuation_censor == 1, 0, 1))
  
  ps_model_discontinuation <- glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + cancer + depression +
                                    diabetes + heartdisease + heart_failure + 
                                    kidney_disease + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = grace_period_data, family=quasibinomial())
  
  # created models with all variables but high standard errors for angina, copd, rheumatoid arthritis, hypertension
  
  grace_period_data$prob_uncensored <- predict(ps_model_discontinuation, type="response")
  
  # as there is only one time point, no need to do cumultive probabilities
  grace_period_data$unstabilized_weight_no_discontinuation <- 1/grace_period_data$prob_uncensored
  print("discontinuation: no discontinuation unstabilized weights")
  print(summary(grace_period_data$unstabilized_weight_no_discontinuation))
  
  baseline_ps_model_discontinuation <-  glm(remain_uncensored ~ age_T0 + gender + months_since_diagnosis +
                                              # baseline conditions
                                              angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                              cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                              heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                              kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                              venous_thromboembolism_baseline + mi_baseline + 
                                              # baseline medications
                                              ace_inhibitors_baseline + antidepressants_baseline +
                                              antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                              arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                              dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                              lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                            data = grace_period_data, family=quasibinomial())
  
  grace_period_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation, type="response")
  grace_period_data <- grace_period_data %>%
    mutate(stabilized_weight_no_discontinuation = treatment_allocation_numerator/prob_uncensored)
  print("discontinuation: no discontinuation stabilized weights")
  print(summary(grace_period_data$stabilized_weight_no_discontinuation))
  
  grace_period_data <- grace_period_data %>% 
    dplyr::select(patid, t, stabilized_weight_no_discontinuation, unstabilized_weight_no_discontinuation)
  
  tapering_arm <- tapering_arm %>%
    left_join(grace_period_data, by=c("patid", "t"))
  
  ## no dose reduction weights which should only be applied to those who have not discontinued
  # those who have discontinued without being censored cannot be censored due to no dose reduction
  # must also only occur during grace period
  not_discontinued_data <- tapering_arm %>%
    filter(t < 7) %>%
    filter(lag(discontinued) == 0  | is.na(lag(discontinued)))
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(remain_uncensored = ifelse(no_dose_reduction_censor == 1, 0, 1))
  
  
  ps_model_discontinuation_no_dose_reduction <- glm(remain_uncensored ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                                      # baseline conditions
                                                      angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                                      cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                                      heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                                      kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                                      venous_thromboembolism_baseline + mi_baseline + 
                                                      # time-varying conditions
                                                      angina + anxiety + atrial_fibrillation + cancer + copd + depression +
                                                      diabetes + heartdisease + heart_failure + heart_hypertension +
                                                      kidney_disease + osteoporosis  + venous_thromboembolism + mi +
                                                      # baseline medications
                                                      ace_inhibitors_baseline  + antidepressants_baseline +
                                                      antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                                      arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                                      dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                                      lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                                      # time-varying medications
                                                      ace_inhibitors + antidepressants + antiepileptics + 
                                                      antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                                      dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                                      oral_anticoagulants,
                                                    data = not_discontinued_data, family=quasibinomial())
  # all covariates have reasonable standard errors apart from rheumatoid arthritis
  
  not_discontinued_data$prob_uncensored <- predict(ps_model_discontinuation_no_dose_reduction, type="response")
  
  # need to calculate the cumulative probabilities of the numerator and denominator
  # get the cumulative probability of remaining uncensored
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(prob_uncensored)) %>%
    ungroup()
  
  not_discontinued_data$unstabilized_weight_no_dose_reduction <- 1/not_discontinued_data$cum_prob_uncensored
  print("discontinuation no dose reduction unstabilized weights")
  print(summary(not_discontinued_data$unstabilized_weight_no_dose_reduction))
  
  baseline_ps_model_discontinuation_no_dose_reduction <-  
    glm(remain_uncensored ~  t + t_squared + age_T0 + gender + months_since_diagnosis +
          # baseline conditions
          angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
          cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
          heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
          kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
          venous_thromboembolism_baseline + mi_baseline + 
          # baseline medications
          ace_inhibitors_baseline + antidepressants_baseline +
          antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
          arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
          dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
          lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
        data = not_discontinued_data, family=quasibinomial())
  
  
  not_discontinued_data$treatment_allocation_numerator <- predict(baseline_ps_model_discontinuation_no_dose_reduction, type="response")
  
  # get the cumulative probabilities for the numerator
  not_discontinued_data <- not_discontinued_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_numerator = cumprod(treatment_allocation_numerator)) %>%
    ungroup()
  
  
  not_discontinued_data <- not_discontinued_data %>%
    mutate(stabilized_weight_no_dose_reduction = cum_prob_numerator/cum_prob_uncensored)
  print("discontinuation: no dose reduction stabilized weights")
  print(summary(not_discontinued_data$stabilized_weight_no_dose_reduction)) 
  
  not_discontinued_data <- not_discontinued_data %>%
    dplyr::select(patid, t, stabilized_weight_no_dose_reduction, unstabilized_weight_no_dose_reduction)
  
  tapering_arm <- tapering_arm %>%
    left_join(not_discontinued_data, by=c("patid", "t"))
  
  # reinitiation weights, which are only applied to those who have discontinued  
  reinitiation_weight_data <- tapering_arm %>%
    filter(discontinued == 1)
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(remain_discontinued = ifelse(reinitiation_censor == 1, 0, 1))
  
  reinitiation_denominator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                    # baseline conditions
                                    angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                    cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                    heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                    kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                    venous_thromboembolism_baseline + mi_baseline + 
                                    # time-varying conditions
                                    anxiety + atrial_fibrillation + depression +
                                   kidney_disease + osteoporosis + rheumatoid_arthritis + venous_thromboembolism + mi +
                                    # baseline medications
                                    ace_inhibitors_baseline  + antidepressants_baseline +
                                    antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                    arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                    dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                    lipid_regulating_drugs_baseline + oral_anticoagulants_baseline +
                                    # time-varying medications
                                    ace_inhibitors + antidepressants + antiepileptics + 
                                    antiplatelets + anxiolytics + arbs + beta_blockers + calcium_blockers +
                                    dementia_drugs + diabetes_drugs + diuretics + lipid_regulating_drugs +
                                    oral_anticoagulants,
                                  data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_denom <- predict(reinitiation_denominator, type="response")
  
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_prob_uncensored = cumprod(reinitiation_denom)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(unstabilized_weight_reinitiation = 1/cum_prob_uncensored)
  print("reinitiation weights unstabilized")
  print(summary(reinitiation_weight_data$unstabilized_weight_reinitiation))
  
  reinitiation_numerator <- glm(remain_discontinued ~ t + t_squared + age_T0 + gender + months_since_diagnosis +
                                  # baseline conditions
                                  angina_baseline + anxiety_baseline + atrial_fibrillation_baseline +
                                  cancer_baseline + copd_baseline + depression_baseline + diabetes_baseline +
                                  heartdisease_baseline + heart_failure_baseline + heart_hypertension_baseline +
                                  kidney_disease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                  venous_thromboembolism_baseline + mi_baseline + 
                                  # baseline medications
                                  ace_inhibitors_baseline + antidepressants_baseline +
                                  antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                  arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                  dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                  lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                                data = reinitiation_weight_data, family=quasibinomial())
  
  reinitiation_weight_data$reinitiation_numerator <- predict(reinitiation_numerator, type="response")
  # get the cumulative probabilities
  reinitiation_weight_data <- reinitiation_weight_data %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    mutate(cum_reinitiation_numerator = cumprod(reinitiation_numerator)) %>%
    ungroup()
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    mutate(stabilized_weight_reinitiation = cum_reinitiation_numerator/cum_prob_uncensored)
  print("reinitiation weights stabilized")
  print(summary(reinitiation_weight_data$stabilized_weight_reinitiation))
  
  reinitiation_weight_data <- reinitiation_weight_data %>%
    dplyr::select(patid, t, stabilized_weight_reinitiation, unstabilized_weight_reinitiation)
  
  tapering_arm <- tapering_arm %>%
    left_join(reinitiation_weight_data, by=c("patid", "t"))
  
  # fill each of the weights downwards
  tapering_arm <- tapering_arm %>%
    group_by(patid) %>%
    arrange(patid, t) %>%
    fill(
      stabilized_weight_no_discontinuation,
      unstabilized_weight_no_discontinuation,
      stabilized_weight_no_dose_reduction,
      unstabilized_weight_no_dose_reduction,
      stabilized_weight_reinitiation,
      unstabilized_weight_reinitiation,
      .direction = "down"
    ) %>%
    ungroup()
  
  # fill in the rest in with 1s
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_no_discontinuation = ifelse(is.na(stabilized_weight_no_discontinuation), 1, 
                                                         stabilized_weight_no_discontinuation),
           unstabilized_weight_no_discontinuation = ifelse(is.na(unstabilized_weight_no_discontinuation), 1,
                                                           unstabilized_weight_no_discontinuation),
           stabilized_weight_no_dose_reduction = ifelse(is.na(stabilized_weight_no_dose_reduction), 1,
                                                        stabilized_weight_no_dose_reduction),
           unstabilized_weight_no_dose_reduction = ifelse(is.na(unstabilized_weight_no_dose_reduction), 1,
                                                          unstabilized_weight_no_dose_reduction),
           stabilized_weight_reinitiation = ifelse(is.na(stabilized_weight_reinitiation), 1,
                                                   stabilized_weight_reinitiation),
           unstabilized_weight_reinitiation = ifelse(is.na(unstabilized_weight_reinitiation), 1,
                                                     unstabilized_weight_reinitiation))
  
  # multiple all of the weights together
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight_untruncated = stabilized_weight_no_discontinuation *
             stabilized_weight_no_dose_reduction * stabilized_weight_reinitiation,
           unstabilized_weight_untruncated = unstabilized_weight_no_discontinuation *
             unstabilized_weight_no_dose_reduction * unstabilized_weight_reinitiation)
  
  # truncate the stabilized weights at the 99th percentile
  quantile_cutoff <- quantile(tapering_arm$stabilized_weight_untruncated, 0.99, na.rm=TRUE)
  tapering_arm <- tapering_arm %>%
    mutate(stabilized_weight = pmin(stabilized_weight_untruncated, quantile_cutoff))
  print("summary of truncated stabilized weights for tapering arm:")
  print(summary(tapering_arm$stabilized_weight))
  
  # truncate the unstabilized weights at the 99th percentile
  upper_bound_unstab <- quantile(tapering_arm$unstabilized_weight_untruncated, probs=0.99)
  
  tapering_arm <- tapering_arm %>%
    mutate(unstabilized_weight = pmin(unstabilized_weight_untruncated, upper_bound_unstab))
  print("summary of truncated unstabilized weights for tapering arm:")
  print(summary(tapering_arm$unstabilized_weight))
  
  tapering_arm <- tapering_arm %>%
    mutate(uncensored = ifelse(no_discontinuation_censor == 0 & no_dose_reduction_censor == 0 & reinitiation_censor == 0,
                               1, 0)) %>%
    mutate(censored = ifelse(uncensored == 0, 1, 0))   
  
  vars <- c("stabilized_weight", "stabilized_weight_untruncated")
  
  # Function to compute the summary stats
  get_summary_stats <- function(x) {
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p1 = quantile(x, 0.01, na.rm = TRUE),
      p50 = quantile(x, 0.50, na.rm = TRUE),
      p99 = quantile(x, 0.99, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  summary_tapering <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(tapering_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_tapering <- summary_tapering %>%
    mutate(arm = "tapering")
  
  summary_continuation <- do.call(rbind, lapply(vars, function(v) {
    stats <- get_summary_stats(continuation_arm[[v]])
    data.frame(variable = v, t(stats))
  }))
  summary_continuation <- summary_continuation %>%
    mutate(arm = "continuation")
  
  summary_weights <- rbind(summary_continuation, summary_tapering)
  
  fwrite(summary_weights, paste0("analysis_results/Tapering Trials/", time_period, "/weights_distribution/negcontrol.csv"))
  
  ################################################################################
  # check covariate balance after grace period
  ################################################################################
  covariates <- c("angina_baseline", "anxiety_baseline", "atrial_fibrillation_baseline", 
                  "cancer_baseline", "copd_baseline", "depression_baseline",
                  "diabetes_baseline", "heartdisease_baseline", "heart_failure_baseline", 
                  "heart_hypertension_baseline", "kidney_disease_baseline",
                  "ace_inhibitors_baseline", "antidepressants_baseline", 
                  "antiplatelets_baseline", "anxiolytics_baseline",
                  "arbs_baseline", "beta_blockers_baseline", "calcium_blockers_baseline",
                  "dementia_drugs_baseline", "diabetes_drug_baseline",
                  "diuretics_baseline", "lipid_regulating_drugs_baseline", 
                  "oral_anticoagulants_baseline", "osteoporosis_baseline",
                  "venous_thromboembolism_baseline", 
                  "rheumatoid_arthritis_baseline", "mi_baseline",
                  "months_since_diagnosis", "age_T0", "gender")
  
  
  continuation_t7 <- continuation_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  tapering_t7 <- tapering_arm %>%
    filter(t == 7 & uncensored == 1) %>%
    dplyr::select(patid, all_of(covariates), unstabilized_weight, treatment)
  
  
  combined_t7 <- rbind(continuation_t7, tapering_t7)
  
  check_weighting <- svydesign(ids = ~1, data = combined_t7, weights = ~unstabilized_weight)
  
  # Calculate the balance table for the subset
  balance_table <- svyCreateTableOne(vars = covariates, 
                                     strata = "treatment", 
                                     data = check_weighting, 
                                     test = FALSE, 
                                     smd = TRUE)
  
  table_one_df <- as.data.frame(print(balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  
  table_one_df$variable <- rownames(table_one_df)
  
  # balance table before weighting
  pre_weighting_smd <- svydesign(ids = ~1, data = combined_t7)
  
  # Calculate the balance table for the subset
  pre_weight_balance_table <- svyCreateTableOne(vars = covariates, 
                                                strata = "treatment", 
                                                data = pre_weighting_smd, 
                                                test = FALSE, 
                                                smd = TRUE)
  
  
  pre_weighting_table <- as.data.frame(print(pre_weight_balance_table, quote = FALSE, noSpaces = TRUE, smd = TRUE))
  pre_weighting_table$variable <- rownames(pre_weighting_table)
  
  
  if(full_sample == TRUE){
    fwrite(table_one_df, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/negcontrol.csv"))
    fwrite(pre_weighting_table, paste0("analysis_results/Tapering Trials/", time_period, "/balance_checks/pre_weighting_negcontrol.csv"))
  }
  
  ################################################################################
  # Outcome model
  ################################################################################
  
  continuation_df <- continuation_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  tapering_df <- tapering_arm %>%
    dplyr::select(all_of(baseline_covariates), t, t_squared, treatment, outcome_date, outcome_month, outcome,
                  age_T0, gender, months_since_diagnosis, stabilized_weight, unstabilized_weight, censored)
  
  combined_data <- rbind(continuation_df, tapering_df) %>%
    filter(censored == 0)
  
  
  event_by_t <- combined_data %>%
    group_by(t) %>%
    summarise(
      n_at_risk = n_distinct(patid),
      n_events  = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(event_by_t,paste0("analysis_results/Tapering Trials/", time_period, "/number_at_risk/neg_control.csv" ))
  
  
  # patient years by arm
  py_by_arm <- combined_data %>%
    group_by(treatment) %>%
    summarise(
      # each row = 1 month, so 1/12 year
      patient_years = ceiling(sum(stabilized_weight * (1/12))),
      .groups = "drop"
    )
  print(py_by_arm)
  
  # number of outcomes
  weighted_events <- with(
    combined_data,
    tapply(stabilized_weight * outcome, treatment, sum)
  )
  print(weighted_events)
  
  crude_model <- glm(outcome ~ t + t_squared + treatment, data = combined_data, family=quasibinomial())
  
  # to create a survival curve, create a model with an interaction between t and treatment
  # and t_squared and treatment
  survival_curve_model <- glm(outcome ~ t + t_squared + treatment +
                                + t*treatment + t_squared*treatment + 
                                + gender + age_T0 + months_since_diagnosis +
                                # baseline conditions
                                angina_baseline + anxiety_baseline + cancer_baseline + copd_baseline + 
                                atrial_fibrillation_baseline + depression_baseline + diabetes_baseline +
                                heartdisease_baseline + osteoporosis_baseline + rheumatoid_arthritis_baseline +
                                venous_thromboembolism_baseline + mi_baseline +
                                # baseline medications
                                ace_inhibitors_baseline + antidepressants_baseline +
                                antiepileptics_baseline + antiplatelets_baseline + anxiolytics_baseline +
                                arbs_baseline + beta_blockers_baseline + calcium_blockers_baseline +
                                dementia_drugs_baseline + diabetes_drug_baseline + diuretics_baseline + 
                                lipid_regulating_drugs_baseline + oral_anticoagulants_baseline,
                              data = combined_data, weights = stabilized_weight, family=quasibinomial())
  
  
  continuation_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "continuation")
  
  tapering_df <- long_format %>%
    filter(followup_curve == 1) %>%
    mutate(treatment = "tapering")
  
  # create hazards under each strategy
  tapering_df$hazard1 <- predict(survival_curve_model, newdata=tapering_df, type="response")
  continuation_df$hazard0 <- predict(survival_curve_model, newdata=continuation_df, type="response")
  
  # calculate survival from the cumulative of 1-hazard
  tapering_df <- tapering_df %>%
    group_by(patid) %>%
    mutate(survival1 = cumprod(1-hazard1)) %>%
    ungroup()
  
  continuation_df <- continuation_df %>%
    group_by(patid) %>%
    mutate(survival0 = cumprod(1-hazard0)) %>%
    ungroup()
  
  # estimate risks as 1-survival
  tapering_df <- tapering_df %>%
    mutate(risk1 = 1-survival1)
  
  continuation_df <- continuation_df %>%
    mutate(risk0 = 1-survival0)
  
  # calculate mean for each time point
  tapering_risk <- tapering_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk1),
              survival = mean(survival1),
              hazard = mean(hazard1))
  
  continuation_risk <- continuation_df %>%
    group_by(t, treatment) %>%
    summarise(risk = mean(risk0),
              survival = mean(survival0),
              hazard = mean(hazard0))
  
  risk_curves <- full_join(tapering_risk, continuation_risk)
  
  risk_curves_wider <- risk_curves %>%
    pivot_wider(
      id_cols = t,
      names_from = treatment,
      values_from = c(risk, hazard, survival),
      names_sep = "_"
    )
  
  risk_curves_wider <- risk_curves_wider %>%
    mutate(hazard_ratio = log(survival_tapering)/log(survival_continuation))
  
  if(full_sample == TRUE){
    fwrite(risk_curves_wider, paste0("analysis_results/Tapering Trials/", time_period, "/outcome_model/negcontrol.csv"))
  }
  
  # hazard ratio is average of hazard ratios over time
  hazard_ratio_estimate <- mean(risk_curves_wider$hazard_ratio)
  
  tapering_risk_12 <- tapering_risk[tapering_risk$t == 12,]$risk * 100
  continuation_risk_12 <- continuation_risk[continuation_risk$t == 12,]$risk * 100
  
  tapering_risk_24 <- tapering_risk[tapering_risk$t == 24,]$risk * 100
  continuation_risk_24 <- continuation_risk[continuation_risk$t == 24,]$risk * 100
  
  # calculate risk ratios and differences
  risk_ratio_12 <- tapering_risk_12/continuation_risk_12
  risk_diff_12 <- tapering_risk_12-continuation_risk_12
  
  risk_ratio_24 <- tapering_risk_24/continuation_risk_24
  risk_diff_24 <- tapering_risk_24-continuation_risk_24
  
  
  risk_df <- data.frame(
    estimate_hr = hazard_ratio_estimate,
    risk_ratio_12months = risk_ratio_12,
    risk_ratio_24months = risk_ratio_24,
    risk_diff_12months = risk_diff_12,
    risk_diff_24months = risk_diff_24,
    abs_risk_continuation_12 = continuation_risk_12,
    abs_risk_continuation_24 = continuation_risk_24,
    abs_risk_tapering_12 = tapering_risk_12, 
    abs_risk_tapering_24 = tapering_risk_24
  )
  rownames(risk_df) <- NULL
  
  if(full_sample == TRUE){
    fwrite(risk_df, paste0("analysis_results/Tapering Trials/", time_period, "/risk_estimates/negcontrol.csv"))
  }
  
  risk_curves <- risk_curves %>%
    mutate(risk=risk*100)
  
  risk_curves <- rbind(
    risk_curves,
    data.frame(t = 0, risk = 0, treatment = "continuation"),
    data.frame(t = 0, risk = 0, treatment = "tapering")
  )
  
  # Plot the survival curves
  survival_curve_plot <- ggplot(risk_curves, aes(x = t, y = risk, color = treatment)) +
    geom_line(linewidth=1)+
    labs(
      title = "Skin Conditions (Negative Control Outcome) Cumulative Incidence Curves by Treatment",
      x = "Time (months)",
      y = "Cumulative Incidence (%)",
      color = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),     # tick label size
      axis.title = element_text(size = 16),    # axis title size
      axis.ticks = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    # vertical line at end of grace period
    geom_vline(xintercept = 6, 
               linetype = "dashed", 
               linewidth = 1,
               color = "black") 
  
  if(full_sample == TRUE){
    ggsave(paste0("analysis_results/Tapering Trials/", time_period, "/survival_curves/negcontrol.png"), plot = survival_curve_plot,
           dpi=300, width = 8, height = 6)
  }
  
  return(risk_df)
  
}

tapering_negcontrol_trials(combined_12weeks, '12weeks', full_sample = TRUE)


################################################################################
# BOOTSTRAP CONFIDENCE INTERVALS
################################################################################

################################################################################
# Bootstrap estimation ----------------------------------------------------------
#
# Non-parametric bootstrap resampling is used to estimate uncertainty around the
# absolute risk estimates and risk differences. Resampling is performed at the
# patient level so that all follow-up intervals for a sampled patient remain
# together. New patient identifiers are created within each bootstrap sample to
# avoid duplicate IDs after sampling with replacement.
#
# The code below first runs the 12-week eligibility analyses and then the 24-week
# eligibility analyses for the tapering discontinuation versus continuation comparison.
################################################################################


n_boot <- 500

stroke_bootstrapped_estimates <- data.frame()
pneumonia_bootstrapped_estimates <- data.frame()
fracture_bootstrapped_estimates <- data.frame()
death_bootstrapped_estimates_12weeks <- data.frame()
death_bootstrapped_estimates_24weeks <- data.frame()
delirium_bootstrapped_estimates <- data.frame()
negcontrol_bootstrapped_estimates <- data.frame()

# Pre-split data by patient_id (MUCH faster for repeated lookups)
df_split <- split(combined_12weeks, combined_12weeks$patid)
patids <- names(df_split)
n_patients <- length(patids)

for (i in 1:n_boot){
  print("bootstrap sample #")
  print(i)
  
  sampled_ids <- sample(patids, size = n_patients, replace = TRUE)
  
  sample_data <- do.call(rbind, lapply(seq_along(sampled_ids), function(j) {
    patient_data <- df_split[[sampled_ids[j]]]
    patient_data$bootstrap_rep <- j
    patient_data$patid<- paste0(sampled_ids[j], "_", j)  # Make unique if needed
    patient_data
  }))
    
    
  stroke_estimate <- tapering_stroke_trials(sample_data, '12weeks', full_sample=FALSE)
  stroke_bootstrapped_estimates <- rbind(stroke_estimate, stroke_bootstrapped_estimates)
  
  pneumonia_estimate <- tapering_pneumonia_trials(sample_data, '12weeks', full_sample=FALSE)
  pneumonia_bootstrapped_estimates <- rbind(pneumonia_estimate, pneumonia_bootstrapped_estimates)
  
  fracture_estimate <- tapering_fracture_trials(sample_data, '12weeks', full_sample=FALSE)
  fracture_bootstrapped_estimates <- rbind(fracture_estimate, fracture_bootstrapped_estimates)
  
  death_estimate <- tapering_death_trials_12weeks(sample_data, '12weeks', full_sample=FALSE)
  death_bootstrapped_estimates_12weeks <- rbind(death_estimate, death_bootstrapped_estimates_12weeks)
  
  delirium_estimate <- tapering_delirium_trials(sample_data, '12weeks', full_sample=FALSE)
  delirium_bootstrapped_estimates <- rbind(delirium_estimate, delirium_bootstrapped_estimates)
  
  fwrite(stroke_bootstrapped_estimates, 
         "analysis_results/Tapering Trials/12weeks/bootstrapped_stroke.csv", row.names=FALSE)
  fwrite(pneumonia_bootstrapped_estimates, 
         "analysis_results/Tapering Trials/12weeks/bootstrapped_pneumonia.csv", row.names=FALSE)
  fwrite(fracture_bootstrapped_estimates, 
         "analysis_results/Tapering Trials/12weeks/bootstrapped_fracture.csv", row.names=FALSE)
  fwrite(death_bootstrapped_estimates_12weeks, 
        "analysis_results/Tapering Trials/12weeks/bootstrapped_death.csv", row.names=FALSE)
  fwrite(delirium_bootstrapped_estimates,
        "analysis_results/Tapering Trials/12weeks/bootstrapped_delirium.csv", row.names=FALSE)
  
  gc()
}

##### 24 weeks
n_boot <- 500

stroke_bootstrapped_estimates <- data.frame()
pneumonia_bootstrapped_estimates <- data.frame()
fracture_bootstrapped_estimates <- data.frame()
death_bootstrapped_estimates_24weeks <- data.frame()
delirium_bootstrapped_estimates <- data.frame()

# Pre-split data by patient_id (MUCH faster for repeated lookups)
df_split <- split(combined_24weeks, combined_24weeks$patid)
patids <- names(df_split)
n_patients <- length(patids)


for (i in 1:n_boot){
  print("bootstrap sample #")
  print(i)
  
  sampled_ids <- sample(patids, size = n_patients, replace = TRUE)
  
  sample_data <- do.call(rbind, lapply(seq_along(sampled_ids), function(j) {
    patient_data <- df_split[[sampled_ids[j]]]
    patient_data$bootstrap_rep <- j
    patient_data$patid<- paste0(sampled_ids[j], "_", j)  # Make unique if needed
    patient_data
  }))
  
  
  stroke_estimate <- tapering_stroke_trials(sample_data, '24weeks', full_sample=FALSE)
  stroke_bootstrapped_estimates <- rbind(stroke_estimate, stroke_bootstrapped_estimates)
  
  pneumonia_estimate <- tapering_pneumonia_trials_24weeks(sample_data, '24weeks', full_sample=FALSE)
  pneumonia_bootstrapped_estimates <- rbind(pneumonia_estimate, pneumonia_bootstrapped_estimates)
  
  fracture_estimate <- tapering_fracture_trials(sample_data, '24weeks', full_sample=FALSE)
  fracture_bootstrapped_estimates <- rbind(fracture_estimate, fracture_bootstrapped_estimates)
  
  death_estimate <- tapering_death_trials_24weeks(sample_data, '24weeks', full_sample=FALSE)
  death_bootstrapped_estimates_24weeks <- rbind(death_estimate, death_bootstrapped_estimates_24weeks)
  
  delirium_estimate <- tapering_delirium_trials_24weeks(sample_data, '24weeks', full_sample=FALSE)
  delirium_bootstrapped_estimates <- rbind(delirium_estimate, delirium_bootstrapped_estimates)
  
  fwrite(stroke_bootstrapped_estimates, 
         "analysis_results/Tapering Trials/24weeks/bootstrapped_stroke.csv", row.names=FALSE)
  fwrite(pneumonia_bootstrapped_estimates, 
         "analysis_results/Tapering Trials/24weeks/bootstrapped_pneumonia.csv", row.names=FALSE)
  fwrite(fracture_bootstrapped_estimates, 
         "analysis_results/Tapering Trials/24weeks/bootstrapped_fracture.csv", row.names=FALSE)
  fwrite(death_bootstrapped_estimates_24weeks, 
         "analysis_results/Tapering Trials/24weeks/bootstrapped_death.csv", row.names=FALSE)
   fwrite(delirium_bootstrapped_estimates,
         "analysis_results/Tapering Trials/24weeks/bootstrapped_delirium.csv", row.names=FALSE)
  
  gc()
}




