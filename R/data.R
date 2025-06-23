#' UKPDS Model Coefficients for Individual Traits
#'
#' @description
#' A 3D array containing a set of coefficients for the UKPDS risk equations.
#' This array specifically holds the coefficients related to individual-level
#' patient characteristics and biomarkers. It is a subset of a larger
#' coefficient set, separated for use in the vectorized C++ simulation functions.
#'
#' @format A 3D numeric array with dimensions 62 x 25 x 1.
#' \describe{
#'   \item{Dimension 1 (Rows)}{62 model predictors/parameters (e.g., `"age"`, `"a1c"`, `"female"`).}
#'   \item{Dimension 2 (Columns)}{25 health outcome equations (e.g., `"hba1c"`, `"mi1_male"`, `"death_nhne"`).}
#'   \item{Dimension 3 (Slices)}{1 bootstrap replicate (`"boot_rep_1"`).}
#' }
#'
#' @source
#' These coefficients are based on the UKPDS Outcomes Model 2 risk equations.
#'
"a_coef_ukpds_ind_traits"

#' UKPDS Model Ancillary Coefficients (lambda, rho, death)
#'
#' @description
#' A 3D array containing a set of ancillary coefficients for the UKPDS risk
#' equations. This array specifically holds the `lambda` (location/intercept),
#' `rho` (shape), and `death` parameters, which are used alongside the main
#' coefficients in the simulation functions.
#'
#' @format A 3D numeric array with dimensions 3 x 25 x 1.
#' \describe{
#'   \item{Dimension 1 (Rows)}{3 ancillary parameters: `"lambda"`, `"rho"`, `"death"`.}
#'   \item{Dimension 2 (Columns)}{25 health outcome equations (e.g., `"hba1c"`, `"mi1_male"`, `"death_nhne"`).}
#'   \item{Dimension 3 (Slices)}{1 bootstrap replicate (`"boot_rep_1"`).}
#' }
#'
#' @source
#' These coefficients are based on the UKPDS Outcomes Model 2 risk equations.
#'
"a_coef_ukpds_other_ind_traits"

#' UKPDS Model Bootstrap Replicates for PSA
#'
#' @description
#' A data frame containing 4,001 sets of bootstrapped coefficients for the UKPDS
#' risk equations. Each row represents a complete parameter set for one iteration
#' of a Probabilistic Sensitivity Analysis (PSA).
#'
#' @format A tibble (data frame) with 4001 rows and 48 columns:
#' \describe{
#'   \item{Rows}{Each row is a single bootstrap replicate.}
#'   \item{Columns}{Each column represents a specific model parameter, such as
#'     `re_hb_fem` (the female coefficient for the HbA1c risk equation) or
#'     `re_sbp_cons` (the intercept for the SBP risk equation).}
#' }
#'
#' @source
#' Generated from a bootstrapping procedure of the original UKPDS Outcomes
#' Model 2 risk equations.
#'
"df_ukpds_bootstraps"

#' Base UKPDS Model Coefficients
#'
#' @description
#' A data frame containing the base coefficients for the UKPDS Outcomes Model 2
#' risk equations. This data frame serves as a template that can be updated with
#' bootstrapped parameter sets for use in a Probabilistic Sensitivity Analysis (PSA).
#'
#' @format A tibble (data frame) with 65 rows and 26 columns:
#' \describe{
#'   \item{Parameter}{A character column listing the model predictors or parameters
#'     (e.g., `"age"`, `"female"`, `"lambda"`).}
#'   \item{hba1c}{Coefficients for the HbA1c risk equation.}
#'   \item{sbp}{Coefficients for the Systolic Blood Pressure risk equation.}
#'   \item{ldl}{Coefficients for the LDL Cholesterol risk equation.}
#'   \item{hdl}{Coefficients for the HDL Cholesterol risk equation.}
#'   \item{bmi}{Coefficients for the Body Mass Index risk equation.}
#'   \item{heart_rate}{Coefficients for the heart rate risk equation.}
#'   \item{wbc}{Coefficients for the White Blood Cell count risk equation.}
#'   \item{haem}{Coefficients for the Haemoglobin risk equation.}
#'   \item{chf}{Coefficients for the Congestive Heart Failure risk equation.}
#'   \item{ihd}{Coefficients for the Ischemic Heart Disease risk equation.}
#'   \item{mi1_male}{Coefficients for the first Myocardial Infarction (male) risk equation.}
#'   \item{mi1_female}{Coefficients for the first Myocardial Infarction (female) risk equation.}
#'   \item{mi2}{Coefficients for a subsequent Myocardial Infarction risk equation.}
#'   \item{stroke_1}{Coefficients for the first stroke risk equation.}
#'   \item{stroke_2}{Coefficients for a subsequent stroke risk equation.}
#'   \item{blindness}{Coefficients for the blindness risk equation.}
#'   \item{ulcer}{Coefficients for the foot ulcer risk equation.}
#'   \item{amp1_no_ulcer}{Coefficients for the first amputation (without ulcer history) risk equation.}
#'   \item{amp1_yes_ulcer}{Coefficients for the first amputation (with ulcer history) risk equation.}
#'   \item{amp2}{Coefficients for a subsequent amputation risk equation.}
#'   \item{esrd}{Coefficients for the End-Stage Renal Disease risk equation.}
#'   \item{death_nhne}{Coefficients for the death (no history, no event) risk equation.}
#'   \item{death_1st_event}{Coefficients for the death (first event) risk equation.}
#'   \item{death_yhne}{Coefficients for the death (yes history, no event) risk equation.}
#'   \item{death_yhye}{Coefficients for the death (yes history, yes event) risk equation.}
#' }
#'
#' @source
#' The base coefficients are derived from the UKPDS Outcomes Model 2.
#'
"df_ukpds_coef"

#' Baseline Patient Population for UKPDS Simulation
#'
#' @description
#' A matrix containing the baseline characteristics for a cohort of simulated
#' patients. Each row represents a single patient at the start of the simulation
#' (cycle 0), and each column represents a specific demographic, biomarker, or
#' health history variable.
#'
#' @format A matrix with 250,000 rows and 66 columns.
#' \describe{
#'   \item{id}{A unique identifier for each patient.}
#'   \item{age}{Patient's age at the start of the simulation.}
#'   \item{age_diag}{Patient's age at diagnosis.}
#'   \item{black}{Indicator for Black ethnicity.}
#'   \item{indian}{Indicator for South Asian ethnicity.}
#'   \item{female}{Indicator for female sex.}
#'   \item{diab_dur}{Duration of diabetes in years.}
#'   \item{diab_dur_log}{Natural log of diabetes duration.}
#'   \item{smoke}{Indicator for current smoking status.}
#'   \item{a1c}{Baseline HbA1c value.}
#'   \item{a1c_lag}{Lagged (previous cycle) HbA1c value.}
#'   \item{a1c_first}{HbA1c value at the start of the simulation (cycle 0).}
#'   \item{bmi}{Baseline Body Mass Index.}
#'   \item{bmi_lt18_5}{Indicator for BMI less than 18.5.}
#'   \item{bmi_gte25}{Indicator for BMI greater than or equal to 25.}
#'   \item{bmi_lag}{Lagged (previous cycle) BMI value.}
#'   \item{bmi_first}{BMI at the start of the simulation (cycle 0).}
#'   \item{egfr}{Estimated Glomerular Filtration Rate.}
#'   \item{egfr_lt60}{Indicator for eGFR less than 60.}
#'   \item{egfr_gte60}{Indicator for eGFR greater than or equal to 60.}
#'   \item{egfr_real}{eGFR value used in some internal calculations.}
#'   \item{hdl}{High-Density Lipoprotein cholesterol.}
#'   \item{hdl_lag}{Lagged (previous cycle) HDL value.}
#'   \item{hdl_first}{HDL at the start of the simulation (cycle 0).}
#'   \item{hdl_real}{HDL value used in some internal calculations.}
#'   \item{heart_rate}{Heart rate in beats per minute.}
#'   \item{heart_rate_lag}{Lagged (previous cycle) heart rate value.}
#'   \item{heart_rate_first}{Heart rate at the start of the simulation (cycle 0).}
#'   \item{heart_rate_real}{Heart rate value used in some internal calculations.}
#'   \item{ldl}{Low-Density Lipoprotein cholesterol.}
#'   \item{ldl_gt35}{Indicator for LDL greater than 3.5 mmol/L (spline term).}
#'   \item{ldl_lag}{Lagged (previous cycle) LDL value.}
#'   \item{ldl_first}{LDL at the start of the simulation (cycle 0).}
#'   \item{ldl_real}{LDL value used in some internal calculations.}
#'   \item{albumin_mm}{Urinary albumin concentration.}
#'   \item{sbp}{Systolic Blood Pressure.}
#'   \item{sbp_lag}{Lagged (previous cycle) SBP value.}
#'   \item{sbp_first}{SBP at the start of the simulation (cycle 0).}
#'   \item{sbp_real}{SBP value used in some internal calculations.}
#'   \item{wbc}{White Blood Cell count.}
#'   \item{wbc_lag}{Lagged (previous cycle) WBC value.}
#'   \item{wbc_first}{WBC at the start of the simulation (cycle 0).}
#'   \item{heamo}{Haemoglobin level.}
#'   \item{heamo_first}{Haemoglobin at the start of the simulation (cycle 0).}
#'   \item{amp_event}{Indicator for a new amputation event in the current cycle.}
#'   \item{amp_event2}{Indicator for a second amputation event.}
#'   \item{amp_hist}{Indicator for a history of amputation.}
#'   \item{atria_fib}{Indicator for atrial fibrillation.}
#'   \item{blindness_event}{Indicator for a new blindness event.}
#'   \item{blindness_hist}{Indicator for a history of blindness.}
#'   \item{chf_event}{Indicator for a new Congestive Heart Failure event.}
#'   \item{chf_hist}{Indicator for a history of Congestive Heart Failure.}
#'   \item{esrd_event}{Indicator for a new End-Stage Renal Disease event.}
#'   \item{esrd_hist}{Indicator for a history of End-Stage Renal Disease.}
#'   \item{ihd_event}{Indicator for a new Ischemic Heart Disease event.}
#'   \item{ihd_hist}{Indicator for a history of Ischemic Heart Disease.}
#'   \item{mi_event}{Indicator for a new Myocardial Infarction event.}
#'   \item{mi_hist}{Indicator for a history of Myocardial Infarction.}
#'   \item{pvd_event}{Indicator for a Peripheral Vascular Disease event.}
#'   \item{stroke_event}{Indicator for a new stroke event.}
#'   \item{stroke_hist}{Indicator for a history of stroke.}
#'   \item{ulcer_event}{Indicator for a new foot ulcer event.}
#'   \item{ulcer_hist}{Indicator for a history of foot ulcer.}
#'   \item{lambda}{Placeholder for the lambda parameter.}
#'   \item{rho}{Placeholder for the rho parameter.}
#'   \item{death}{Indicator for death.}
#' }
#'
#' @source
#' A synthetic dataset created to be representative of a population with
#' Type 2 Diabetes for use with the UKPDS Outcomes Model 2.
#'
"m_ukpds_pop"

#' UKPDS Model Predictor and Parameter Names
#'
#' @description
#' A character vector containing the names of all 65 predictors, health state
#' variables, and ancillary parameters used across the UKPDS Outcomes Model 2
#' risk equations.
#'
#' This vector defines the row names for the main coefficient matrix
#' (`df_ukpds_coef`) and the column names for the patient characteristics
#' matrix (`m_ukpds_pop`), ensuring consistency between them.
#'
#' @format A character vector with 65 elements:
#' \describe{
#'   \item{Demographics & Baseline}{`age`, `age_diag`, `black`, `indian`, `female`, `diab_dur`, `diab_dur_log`, `smoke`}
#'   \item{Biomarkers & Clinical Measurements}{
#'     `a1c`, `a1c_lag`, `a1c_first`,
#'     `bmi`, `bmi_lt18_5`, `bmi_gte25`, `bmi_lag`, `bmi_first`,
#'     `egfr`, `egfr_lt60`, `egfr_gte60`, `egfr_real`,
#'     `hdl`, `hdl_lag`, `hdl_first`, `hdl_real`,
#'     `heart_rate`, `heart_rate_lag`, `heart_rate_first`, `heart_rate_real`,
#'     `ldl`, `ldl_gt35`, `ldl_lag`, `ldl_first`, `ldl_real`,
#'     `albumin_mm`,
#'     `sbp`, `sbp_lag`, `sbp_first`, `sbp_real`,
#'     `wbc`, `wbc_lag`, `wbc_first`,
#'     `heamo`, `heamo_first`
#'   }
#'   \item{Health Event & History Flags}{
#'     `amp_event`, `amp_event2`, `amp_hist`,
#'     `atria_fib`,
#'     `blindness_event`, `blindness_hist`,
#'     `chf_event`, `chf_hist`,
#'     `esrd_event`, `esrd_hist`,
#'     `ihd_event`, `ihd_hist`,
#'     `mi_event`, `mi_hist`,
#'     `pvd_event`,
#'     `stroke_event`, `stroke_hist`,
#'     `ulcer_event`, `ulcer_hist`
#'   }
#'   \item{Ancillary Model Parameters}{`lambda`, `rho`, `death`}
#' }
#'
#' @source
#' The variable names are derived from the UKPDS Outcomes Model 2.
#'
"v_coef_names"

#' UKPDS Model Outcome and Risk Equation Names
#'
#' @description
#' A character vector containing the names of all 25 health outcomes and
#' biomarker risk equations from the UKPDS Outcomes Model 2.
#'
#' This vector defines the column names for the main coefficient matrix
#' (`df_ukpds_coef`), specifying which set of coefficients corresponds to each
#' health event or biomarker progression model.
#'
#' @format A character vector with 25 elements, categorized as follows:
#' \describe{
#'   \item{Biomarker Progression Models}{
#'     `hba1c`, `sbp`, `ldl`, `hdl`, `bmi`, `heart_rate`, `wbc`, `haem`
#'   }
#'   \item{Health Event Models (Non-fatal)}{
#'     `chf`, `ihd`,
#'     `mi1_male`, `mi1_female`, `mi2`,
#'     `stroke_1`, `stroke_2`,
#'     `blindness`,
#'     `ulcer`,
#'     `amp1_no_ulcer`, `amp1_yes_ulcer`, `amp2`,
#'     `esrd`
#'   }
#'   \item{Mortality Models}{
#'     `death_nhne`, `death_1st_event`, `death_yhne`, `death_yhye`
#'   }
#' }
#'
#' @source
#' The names are derived from the UKPDS Outcomes Model 2 risk equations.
#'
"v_factors_names"
