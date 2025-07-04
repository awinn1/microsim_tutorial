# This is Wael's version of the model

rm(list = ls())

library(pacman)
pacman::p_load(haven,
               Rcpp,
               RcppArmadillo,
               profvis,
               readr)

# Step 1: Import the matrix of coefficients ####
# Load necessary libraries
# Read the coefficient matrix from a CSV or RData file
df_UKPDS_coef <- readr::read_csv("data/ukpds_coef.csv") # Load coefficient matrix from CSV

# Replace NAs with 0s to avoid missing values in calculations
df_UKPDS_coef[is.na(df_UKPDS_coef)] <- 0  

# Extract parameter names (used as row names)
v_coef_names <- df_UKPDS_coef$Parameter # Get row names from the 'Parameter' column

# Determine the number of parameters (rows)
n_coef <-  length(v_coef_names)  # Count the number of parameters

# Extract factor names (used as column names), excluding the first column
v_factors_names <- colnames(df_UKPDS_coef[-1])  # Get column names excluding 'Parameter'

# Determine the number of factors (columns)
n_equa <- length(v_factors_names)  # Count the number of factors

# Allow for bootstrapped coefficients 
n_boot <- 1

# Create an array that holds onto everything!
a_coef_ukpds <- array(
  data = NA,
  dim = c(n_coef, n_equa, n_boot),
  dimnames = list(v_coef_names, v_factors_names, paste0("boot_rep_", 1:n_boot))
)

dimnames(a_coef_ukpds)

# Fill in the array with coefficients from the dataset
a_coef_ukpds[, v_factors_names, 1] <- as.matrix(
  df_UKPDS_coef[, v_factors_names, drop = FALSE]
)

# Extract individual characteristics at initial boot strap slice (Aaron, why?)
a_coef_ukpds_ind_traits<- a_coef_ukpds[1:62, , "boot_rep_1", drop = FALSE]
a_coef_ukpds_other_ind_traits<- a_coef_ukpds[63:65, , ,drop = FALSE]

# Step 2: Load the patient dataset and save as matrix
m_ukpds_pop <- read_csv("data/population.csv") |> 
  as.matrix()

# show the names of the variables and rows
print(dimnames(m_ukpds_pop)) 

seed    <- 1234                     # random number generator state
num_i <- 250000                          # number of simulated individuals
# Define the number of time points
num_cycles <- 20                    # maximum length of a simulation 
set.seed(seed)    # set the seed to ensure reproducible samples below
v_ids <- paste("id",   1:num_i,    sep ="_")
v_cycle_nms <- paste("cycle", 0:num_cycles, sep ="_")

# Create AN ARRAY with columns for each variable, row for each person, and a 
# slice for each period
a_all_ind_traits <- array(   
  data = NA, 
  dim = c(num_i, n_coef, num_cycles + 1),
  dimnames = list(v_ids, v_coef_names , v_cycle_nms )  
)
# need this to be the same number of columns as the coefficient table is long/rows
print(dim(a_all_ind_traits)) # to verify the dimensions
print(dimnames(a_all_ind_traits)) # to verify the dimension names

# Fill patient trace 
## First confirm all trait names in m_ukpds_pop and a_all_ind_traits match
stopifnot(
  all(dimnames(a_all_ind_traits)[[2]] %in% colnames(m_ukpds_pop))
)
## Identify m_ukpds_pop columns with time-invariant data
v_t_invar_col_nms <- grep(
  pattern = "^age_diag$|^black$|^indian$|^female$|^smoke$|_first$",
  x = colnames(m_ukpds_pop), 
  ignore.case = FALSE,
  value = TRUE
)
## Identify m_ukpds_pop columns with time-varying data
v_t_var_col_nms <- setdiff(
  dimnames(a_all_ind_traits)[[2]], v_t_invar_col_nms
)
a_all_ind_traits[, v_t_var_col_nms, 1] <- m_ukpds_pop[
  , v_t_var_col_nms, drop = FALSE
]
a_all_ind_traits[, v_t_invar_col_nms, ] <- m_ukpds_pop[
  , v_t_invar_col_nms, drop = FALSE
]

# break array up into stuff individual traits
a_ind_traits <- a_all_ind_traits[,1:62,] 
# remaining array that captures lambda, rho and death
a_other_ind_traits <- a_all_ind_traits[,63:65,] 

# need this to be the same number of columns as the coefficient table is long/rows
print(dim(a_ind_traits)) # to verify the dimensions
print(dimnames(a_ind_traits)) # to verify the dimension names

# Step 3: Define functions for risk factor progression ####
# Function for linear progression of risk factors

#' @title Calculate Biomarkers 
#'
#' @details
#' This function calculates patient-specific factors to predict the time path of
#' a biomarker. 
#'
#' @param m_ind_traits A matrix containing patient characteristics at a specific
#' time step.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for 
#' calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of other coefficients used
#' for calculating risk.
#' @param biomarker_eq A character string specifying the health outcome equation
#' (e.g., "ihd").
#' 
#' @return Matrix containing the updated biomarker is stored.
#' @export
#' @examples
#' /dontrun{
#' # Calculate the updated value for the 'hba1c' biomarker
#' m_updated_hba1c <- biomarker(
#'   m_ind_traits = a_ind_traits[ , ,1],
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
#'   biomarker_eq = "hba1c"
#' )
#'
#' # The 'm_updated_hba1c' one-column matrix now contains the updated biomarker
#' values. You can inspect the result:
#' head(updated_hba1c)
#' }
biomarker_w <- function(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  m_updated_biomarker <- m_ind_traits %*%
    a_coef_ukpds_ind_traits[, biomarker_eq, 1] +
    a_coef_ukpds_other_ind_traits["lambda", biomarker_eq, 1]
  
  return(m_updated_biomarker)
}

# Step 4: Create a function to apply all risk factor models ####
# Combine risk factor functions into a single pipeline
# Update patient data over time

#' @title Update Multiple Biomarkers in a Transition Matrix
#'
#' @description
#' This function updates multiple biomarker values in the transition matrix for
#' a given time step.
#'
#' @param m_ind_traits_init A matrix containing patient characteristics at the
#' initial time step (obtained from `a_ind_traits[,,1]`).
#' @param m_ind_traits_nStep The patient matrix to update with new biomarker
#' values (obtained from `a_ind_traits[,,next_row]`).
#' @param m_ind_traits_step The patient matrix containing patient data at the
#' current time step (obtained from `a_ind_traits[,,time_step]`).
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for
#' calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of other coefficients used
#' for calculating risk (e.g., lambda).
#'
#' @return The updated transition matrix `m_ind_traits_nStep` with new biomarker
#' values.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have your a_ind_traits array
#' time_step <- 1
#' next_row <- 2
#'
#' m_ind_traits_init <- a_ind_traits[,,1]
#' m_ind_traits_nStep <- a_ind_traits[,,next_row]
#' m_ind_traits_step <- a_ind_traits[,,time_step]
#'
#' # Assuming a_coef_ukpds_ind_traits and a_coef_ukpds_other_ind_traits are defined
#'
#' updated_m_ind_traits_nStep <- update_all_biomarkers(
#'   m_ind_traits_init = m_ind_traits_init,
#'   m_ind_traits_nStep = m_ind_traits_nStep,
#'   m_ind_traits_step = m_ind_traits_step,
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
#' )
#' }
update_all_biomarkers_w <- function(
    m_ind_traits_init,
    m_ind_traits_nStep,
    m_ind_traits_step,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits) {
  # predict the next period (and perform transformations as needed)
  # the biomarkers use real values of variables, but the event equations use transformed variables
  m_ind_traits_nStep[, "a1c"] <- biomarker_w(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "hba1c")
  m_ind_traits_nStep[, "sbp_real"] <- biomarker_w(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "sbp")
  m_ind_traits_nStep[, "sbp"] <- m_ind_traits_nStep[, "sbp_real"] / 10
  m_ind_traits_nStep[, "ldl_real"] <- biomarker_w(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "ldl")
  m_ind_traits_nStep[, "ldl"] <- m_ind_traits_nStep[, "ldl_real"] * 10
  m_ind_traits_nStep[, "hdl_real"] <- biomarker_w(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "hdl")
  m_ind_traits_nStep[, "hdl"] <- m_ind_traits_nStep[, "hdl_real"] * 10
  m_ind_traits_nStep[, "bmi"] <- biomarker_w(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "bmi")
  m_ind_traits_nStep[, "heart_rate_real"] <- biomarker_w(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "heart_rate")
  m_ind_traits_nStep[, "heart_rate"] <- m_ind_traits_nStep[, "heart_rate_real"] / 10
  m_ind_traits_nStep[, "wbc"] <- biomarker_w(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "wbc")
  m_ind_traits_nStep[, "heamo"] <- biomarker_w(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "haem")
  
  # Update lag and first occurrence columns
  m_ind_traits_nStep[, "a1c_lag"] <- m_ind_traits_step[, "a1c"]
  m_ind_traits_nStep[, "bmi_lag"] <- m_ind_traits_step[, "bmi"]
  m_ind_traits_nStep[, "bmi_lt18_5"] <- as.integer(m_ind_traits_nStep[, "bmi"] < 18.5)
  m_ind_traits_nStep[, "bmi_gte25"] <- as.integer(m_ind_traits_nStep[, "bmi"] >= 25)
  
  m_ind_traits_nStep[, "hdl_lag"] <- m_ind_traits_step[, "hdl_real"]
  m_ind_traits_nStep[, "heart_rate_lag"] <- m_ind_traits_step[, "heart_rate_real"]
  
  # Check if this is functioning as a spline
  m_ind_traits_nStep[, "ldl_gt35"] <- as.integer(m_ind_traits_nStep[, "ldl_real"] > 35) / 10
  m_ind_traits_nStep[, "ldl_lag"] <- m_ind_traits_step[, "ldl_real"]
  m_ind_traits_nStep[, "sbp_lag"] <- m_ind_traits_step[, "sbp_real"]
  m_ind_traits_nStep[, "wbc_lag"] <- m_ind_traits_step[, "wbc"]
  
  # Update additional values
  m_ind_traits_nStep[, "egfr"] <- m_ind_traits_init[, "egfr"]
  m_ind_traits_nStep[, "egfr_real"] <- m_ind_traits_init[, "egfr_real"]
  m_ind_traits_nStep[, "egfr_lt60"] <- m_ind_traits_init[, "egfr_lt60"]
  m_ind_traits_nStep[, "egfr_gte60"] <- m_ind_traits_init[, "egfr_gte60"]
  m_ind_traits_nStep[, "albumin_mm"] <- m_ind_traits_init[, "albumin_mm"]
  
  # Return updated matrix
  return(m_ind_traits_nStep)
}

# Step 5: Define event functions (Weibull/Exponential and Logistic) ####
# Weibull distribution function for event occurrence
# Logistic regression for binary event prediction

#' @title Calculate Transition Probability Based on a Weibull Model and Update Patient
#' State 
#'
#' @details
#' This function calculates patient-specific factors, cumulative hazards, 
#' and the transition probability for a given health outcome (e.g., "ihd"). 
#' The function updates the provided `m_ind_traits` matrix with the event
#' occurrence at the specified time step.
#' Note: An exponential model is a special case of the Weibull model where the
#' shape parameter (ρ) is set to 1, meaning the hazard function remains constant
#' over time, resulting in a constant rate of event occurrence rather than a
#' time-dependent rate.
#'
#' @param m_ind_traits_step A matrix containing patient characteristics at a
#' specific time step (obtained from `a_ind_traits[,,time_step]`).
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for
#' calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of other coefficients used
#' for calculating risk (e.g., lambda).
#' @param health_outcome A character string specifying the health outcome
#' equation (e.g., "ihd").
#' 
#' @return One-column matrix specifying whether the event occurred.
#' @export
#' 
#' @examples
#' \dontrun{
#' time_step <- 1
#' health_outcome <- "ihd"
#'
#' # Extract the matrix for the specific time step
#' m_ind_traits_step <- a_ind_traits[,,time_step]
#'
#' # Calculate the probability of a health event
#' event_occurred <- weibull_event(
#'   m_ind_traits_step = m_ind_traits_step,
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
#'   health_outcome = health_outcome
#' )
#'
#' # Print whether the event occurred for each individual at this time step
#' print(head(event_occurred))
#' }
weibull_event_w <- function(
    m_ind_traits_step,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- m_ind_traits_step %*%
    a_coef_ukpds_ind_traits[, health_outcome, 1] +
    as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1])
  
  # Compute cumulative hazard at the current time step
  cum_hazard_t <- exp(patient_factors) *
    (m_ind_traits_step[, "diab_dur"]^
       (a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]))
  
  # Compute cumulative hazard at the next time step (by adding 1 year to
  # diabetes duration)
  cum_hazard_t1 <- exp(patient_factors) *
    ((m_ind_traits_step[, "diab_dur"] + 1)^
       (a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]))
  
  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(nrow(m_ind_traits_step))
  
  # Return the updated matrix
  return(event)
}

#' @title Calculate Transition Probability Based on a Logistic Regression and Update Patient State 
#'
#' @details
#' This function calculates patient-specific factors, cumulative hazards, 
#' and the transition probability for a given health outcome (e.g., "ihd"). 
#' The function updates the provided `m_ind_traits` matrix with the event occurrence 
#' at the specified time step.
#'
#' @param m_ind_traits_step A matrix containing patient characteristics at a
#' specific time step (obtained from `a_ind_traits[,,time_step]`).
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for
#' calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of other coefficients used
#' for calculating risk (e.g., lambda).
#' @param health_outcome A character string specifying the health outcome
#' equation (e.g., "ihd").
#'  
#' @return One-column matrix specifying whether the event occurred.
#' @export
#' 
#' @examples
#' \dontrun{
#' time_step <- 1
#' health_outcome <- "ulcer"
#'
#' # Extract the matrix for the specific time step
#' m_ind_traits_step <- a_ind_traits[,,time_step]
#'
#' # Calculate the probability of a health event
#' event_occurred <- logistic_event(
#'   m_ind_traits_step = m_ind_traits_step,
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
#'   health_outcome = health_outcome
#' )
#'
#' # Print whether the event occurred for each individual at this time step
#' print(head(event_occurred))
#' }
logistic_event_w <- function(
    m_ind_traits_step,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- m_ind_traits_step %*%
    a_coef_ukpds_ind_traits[, health_outcome, 1] +
    as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1])
  
  # Calculate transition probability
  trans_prob <- 1-(exp(-patient_factors)/(1+exp(-patient_factors)))^1
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(nrow(m_ind_traits_step))
  
  # Return the value
  return(event)
}

# Step 6: Initialize event and history variables ####
#' @title Update Health Events Over Time Steps
#'
#' @description This function updates health events in a patient data matrix
#' (`m_ind_traits_step`) by applying Weibull and logistic event functions in a
#' randomized order.
#'
#' @param m_ind_traits_init A matrix containing patient-level data at the
#' initial time step (obtained from `a_ind_traits[,,1]`).
#' @param m_ind_traits_step A matrix containing patient-level data at the
#' current time step to be updated (obtained from `a_ind_traits[,,time_step]`).
#' @param m_ind_traits_pStep A matrix containing patient-level data at the
#' previous time step (obtained from `a_ind_traits[,,max(1, time_step - 1)]`).
#' @param a_coef_ukpds_ind_traits A coefficient matrix used in Weibull and
#' logistic event calculations.
#' @param a_coef_ukpds_other_ind_traits A coefficient matrix containing other
#' parameters for Weibull and logistic event calculations (e.g., lambda, rho).
#'
#' @return Updated `m_ind_traits` matrix with event and history values updated
#' for the given time step.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming a_ind_traits and a_coef_ukpds_ind_traits are defined
#' time_step <- 2
#' m_ind_traits_init <- a_ind_traits[,,1]
#' m_ind_traits_step <- a_ind_traits[,,time_step]
#' m_ind_traits_pStep <- a_ind_traits[,,time_step - 1]
#'
#' updated_m_ind_traits_step <- update_health_events(
#'   m_ind_traits_init = m_ind_traits_init,
#'   m_ind_traits_step = m_ind_traits_step,
#'   m_ind_traits_pStep = m_ind_traits_pStep,
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
#' )
#' print(head(updated_m_ind_traits_step))
#' }
update_health_events_w <- function(
    m_ind_traits_init,
    m_ind_traits_step,
    m_ind_traits_pStep,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits) {
  # Ensure m_ind_traits_step remains a matrix
  if (!is.matrix(m_ind_traits_step)) {
    stop("m_ind_traits_step must be a matrix.")
  }
  
  # Initialize event variables and update history
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")
  
  # create event and history column names once and save each group in a vector
  v_event_cols      <- paste0(events, "_event")
  v_history_cols    <- paste0(events, "_hist")
  
  # Update history columns in one vectorized call
  m_ind_traits_step[, v_event_cols] <- 0
  m_ind_traits_step[, v_history_cols] <- pmax(
    m_ind_traits_pStep[, v_history_cols],
    m_ind_traits_pStep[, v_event_cols]
  )
  
  m_ind_traits_step[, "atria_fib"] <- m_ind_traits_init[, "atria_fib"]
  m_ind_traits_step[, "pvd_event"] <- 0
  m_ind_traits_step[, "amp_event2"] <- 0
  
  # Randomize event order
  randomized_events <- sample(events)
  
  for (event in randomized_events) {
    if (event == "amp") {
      amp1_no_ulcer <- weibull_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "amp1_no_ulcer"
      )
      amp1_yes_ulcer <- weibull_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "amp1_yes_ulcer"
      )
      
      m_ind_traits_step[, "amp_event"] <-
        (amp1_no_ulcer * (m_ind_traits_step[, "ulcer_hist"] == 0)) +
        (amp1_yes_ulcer * (m_ind_traits_step[, "ulcer_hist"] == 1))
      
      # Ensure that this is a new event
      m_ind_traits_step[, "amp_event"] <-
        m_ind_traits_step[, "amp_event"] *
        (m_ind_traits_step[, "amp_hist"] == 0)
      
      # Calculate amp2 event
      amp2 <- weibull_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "amp2"
      )
      
      m_ind_traits_step[, "amp_event2"] <- 0
      
      m_ind_traits_step[, "amp_event2"] <-
        amp2 * (m_ind_traits_step[, "amp_hist"] == 1)
      
    } else if (event == "mi") {
      mi1_male <- weibull_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "mi1_male"
      )
      mi1_female <- weibull_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "mi1_female"
      )
      
      m_ind_traits_step[, "mi_event"] <-
        (mi1_male * (m_ind_traits_step[, "female"] == 0)) +
        (mi1_female * (m_ind_traits_step[, "female"] == 1))
      
      m_ind_traits_step[, "mi_event"] <-
        m_ind_traits_step[, "mi_event"] *
        (m_ind_traits_step[, "mi_hist"] == 0)
      
      mi2 <- weibull_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "mi2"
      )
      
      m_ind_traits_step[, "mi_event"] <-
        (m_ind_traits_step[, "mi_hist"] == 0) *
        m_ind_traits_step[, "mi_event"] +
        (m_ind_traits_step[, "mi_hist"] == 1) * mi2
      
    } else if (event == "stroke") {
      stroke1 <- weibull_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "stroke_1"
      )
      stroke2 <- weibull_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "stroke_2"
      )
      
      m_ind_traits_step[, "stroke_event"] <-
        (stroke1 * (m_ind_traits_step[, "stroke_hist"] == 0)) +
        (stroke2 * (m_ind_traits_step[, "stroke_hist"] == 1))
      
      m_ind_traits_step[, "stroke_event"] <-
        m_ind_traits_step[, "stroke_event"] *
        (m_ind_traits_step[, "stroke_hist"] == 0)
      
    } else if (event == "ulcer") {
      m_ind_traits_step[, "ulcer_event"] <- logistic_event_w(
        m_ind_traits_step = m_ind_traits_step,
        a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
        health_outcome = "ulcer"
      )
      
      m_ind_traits_step[, "ulcer_event"] <- m_ind_traits_step[, "ulcer_event"] *
        (m_ind_traits_step[, "ulcer_hist"] == 0)
      
    } else {
      m_ind_traits_step[, paste0(event, "_event")] <-
        weibull_event_w(
          m_ind_traits_step = m_ind_traits_step,
          a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
          a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
          health_outcome = event
        )
      
      m_ind_traits_step[, paste0(event, "_event")] <-
        m_ind_traits_step[, paste0(event, "_event")] *
        (m_ind_traits_step[, paste0(event, "_hist")] == 0)
    }
  }
  
  return(m_ind_traits_step)
}

# Step 7: Define a mortality function ####
# Combine relevant event functions affecting mortality
# Estimate survival probability over time
#' @title Calculate Transition Probability Based on a Gompertz Model and Update
#' Patient State
#'
#' @details
#' This function calculates patient-specific factors, cumulative hazards,
#' and the transition probability for mortality using a Gompertz model.
#' The function determines whether the mortality event occurs at the specified
#' time step.
#'
#' @param m_ind_traits_step A matrix containing patient characteristics at a
#' specific time step (obtained from `a_ind_traits[,,time_step]`).
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for
#' calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of other coefficients used
#' for calculating risk (e.g., rho for Gompertz).
#' @param health_outcome A character string specifying the health outcome
#' equation (e.g., "death_yhne").
#'
#' @return A logical vector indicating whether the mortality event occurred for
#' each individual.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming a_ind_traits, a_coef_ukpds_ind_traits, and
#' # a_coef_ukpds_other_ind_traits are already defined.
#' time_step <- 1
#' health_outcome <- "death_yhne"
#'
#' # Extract the matrix for the specific time step
#' m_ind_traits_step <- a_ind_traits[,,time_step]
#'
#' # Calculate the probability of mortality using the Gompertz model
#' mortality_event <- gompertz_event(
#'   m_ind_traits_step = m_ind_traits_step,
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
#'   health_outcome = health_outcome
#' )
#'
#' # Print whether the mortality event occurred for each individual
#' print(head(mortality_event))
#' }
gompertz_event_w <- function(
    m_ind_traits_step,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- m_ind_traits_step %*%
    a_coef_ukpds_ind_traits[, health_outcome, 1] +
    as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1])
  
  # Compute cumulative hazard at the current time step
  cum_hazard_t <- (1/a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]) *
    exp(patient_factors) * (
      exp(m_ind_traits_step[, "age"] *
            (a_coef_ukpds_other_ind_traits["rho", health_outcome, 1])) -
        1)
  
  # Compute cumulative hazard at the next time step (by adding 1 year to
  # diabetes duration)
  cum_hazard_t1 <- (1/a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]) *
    exp(patient_factors) * (
      exp((m_ind_traits_step[, "age"] + 1) *
            (a_coef_ukpds_other_ind_traits["rho", health_outcome, 1])) -
        1)
  
  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(nrow(m_ind_traits_step))
  
  # Return the updated matrix
  return(event)
}

#' @title Calculate Mortality Events for a Given Time Step
#'
#' @description This function calculates mortality events for a given time step
#' based on new health events and medical history using Gompertz and logistic
#' models.
#'
#' @param m_ind_traits_step A matrix containing patient-level data at the
#' current time step (obtained from `a_ind_traits[,,time_step]`).
#' @param m_other_ind_traits_step A matrix containing mortality status at the
#' current time step (obtained from `a_other_ind_traits[,,time_step]`).
#' @param m_other_ind_traits_pStep A matrix containing mortality status at
#' the previous time step (obtained from
#' `a_other_ind_traits[,,max(1, time_step - 1)]`).
#' @param m_ind_traits_pStep A matrix containing patient-level data at the
#' previous time step (obtained from `a_ind_traits[,,max(1, time_step - 1)]`).
#' @param a_coef_ukpds_ind_traits A coefficient matrix used in Gompertz and
#' logistic event calculations.
#' @param a_coef_ukpds_other_ind_traits A coefficient matrix containing other
#' parameters for Gompertz and logistic event calculations (e.g., lambda, rho).
#'
#' @return Numeric vector representing the mortality and should be assigned
#' (saved) to column "death" in the matrix `m_other_ind_traits_step` or the
#' array.
#' @export
#' 
#' @examples
#' \dontrun{
#' time_step <- 1
#' m_ind_traits_step <- a_ind_traits[,,time_step]
#' m_other_ind_traits_step <- a_other_ind_traits[,,time_step]
#' m_other_ind_traits_pStep <- a_other_ind_traits[,,max(1, time_step - 1)]
#' m_ind_traits_pStep <- a_ind_traits[,,max(1, time_step - 1)]
#'
#' v_death_step <- mortality(
#'   m_ind_traits_step = m_ind_traits_step,
#'   m_other_ind_traits_step = m_other_ind_traits_step,
#'   m_other_ind_traits_pStep = m_other_ind_traits_pStep,
#'   m_ind_traits_pStep = m_ind_traits_pStep,
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
#' )
#' print(head(v_death_step))
#' }
mortality_w <- function(
    m_ind_traits_step, 
    m_other_ind_traits_step, 
    m_other_ind_traits_pStep, 
    m_ind_traits_pStep, 
    a_coef_ukpds_ind_traits, 
    a_coef_ukpds_other_ind_traits) {
  
  # Calculate new health event occurrence and prior history
  
  # Define events of interest
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")
  
  # Create vectors containing events and event-history names:
  v_event_cols <- paste0(events, "_event")
  v_hist_cols  <- paste0(events, "_hist")
  
  # Get the maximum across those columns, for the given time_step
  # Calculate any new health event
  new_event  <- max(m_ind_traits_step[, v_event_cols])
  # Calculate any prior history of health events
  any_history <- max(m_ind_traits_step[, v_hist_cols])
  
  
  # Determine event-history combinations
  nhne <- new_event == 0 & any_history == 0  # No history, no event
  yhne <- new_event == 0 & any_history == 1  # Yes history, no event
  nhye <- new_event == 1 & any_history == 0  # No history, new event
  yhye <- new_event == 1 & any_history == 1  # Yes history, new event
  
  # Mortality calculations using Gompertz and logistic models
  death_nhne <- gompertz_event_w(
    m_ind_traits_step = m_ind_traits_step, 
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits, 
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits, 
    health_outcome = "death_nhne"
  )
  
  death_yhne <- gompertz_event_w(
    m_ind_traits_step = m_ind_traits_step, 
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits, 
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
    health_outcome = "death_yhne"
  )
  
  death_nhye <- logistic_event_w(
    m_ind_traits_step = m_ind_traits_step,
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits, 
    health_outcome = "death_1st_event"
  )
  
  death_yhye <- logistic_event_w(
    m_ind_traits_step = m_ind_traits_step, 
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits, 
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits, 
    health_outcome = "death_yhye"
  )
  
  # Calculate new mortality status
  new_death <- nhne * death_nhne + 
    yhne * death_yhne + 
    nhye * death_nhye + 
    yhye * death_yhye
  
  # Update the mortality status in the matrix for the given time step
  m_deaths <- as.integer(new_death) + m_other_ind_traits_pStep[, "death"]
  
  return(m_deaths)
}

# Step 8: Simulate disease progression and mortality for all patients ####
# Initialize patient data
# Loop through time points to update risk factors and events
# Store results
# discount rate 


# Step 9: Simulate disease progression and mortality for 999 additional patients  ####
# Loop over all patients
# Store and summarize population-level results

discount_rate <-0.03
# qalys
q_baseline <- 0.785
q_blindness <- -0.074
q_amp<- -0.280
q_chf <- -0.108
q_esrd <- -0.204
q_ihd <- -0.090
q_mi <- -0.055
q_stroke <- -0.164
q_ulcer <- -0.170
# costs   
c_baseline <- 1990
c_blindness_e <- 4247
c_blindness_c <- 2206
c_amp_e<- 15153      
c_amp_c<- 5328
c_chf_e <- 5650
c_chf_c <- 4277
c_esrd_e <- 43359
c_esrd_c <- 43359
c_ihd_e <- 14001
c_ihd_c <- 3550
c_mi_e <- 9518
c_mi_c <- 3424
c_stroke_e <- 10755
c_stroke_c <- 3534
c_ulcer_e <- 7076
c_ulcer_c <- 1072

column_names <- list("cost", "qalys", "disc_costs", "disc_qalys")
patient_summary_file <- matrix(   
  data = NA, 
  nrow = num_i, 
  ncol = 4,
  dimnames = list(NULL,column_names)  # Removing row names saves 80% mem (16 MB)
)

ptm <- proc.time()

# Extract first slice, matrix of initial values
m_ind_traits_init <- a_ind_traits[, , 1]

# predict 2nd period biomarkers
a_ind_traits[, , 2] <-  update_all_biomarkers_w(
  m_ind_traits_init = m_ind_traits_init,
  m_ind_traits_nStep  = a_ind_traits[,, 2],
  m_ind_traits_step = a_ind_traits[,, 1],
  a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
  a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
)

# create a patient population 
for (time_step in 2:num_cycles) {
  a_ind_traits[, "age", time_step] <- a_ind_traits[
    ,"age", max(1, time_step - 1)
  ] + 1
  a_ind_traits[, "diab_dur", time_step] <- a_ind_traits[
    ,"diab_dur" , max(1, time_step - 1)
  ] + 1    
  a_ind_traits[, "diab_dur_log", time_step] <- log(
    a_ind_traits[,"diab_dur", time_step]
  )
  
  # ready to simulate 
  
  # event prediction at t
  a_ind_traits[,, time_step] <- update_health_events_w(
    m_ind_traits_init = m_ind_traits_init,
    m_ind_traits_step = a_ind_traits[,, time_step],
    m_ind_traits_pStep = a_ind_traits[,, time_step - 1],
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
  )
  
  # mortality prediction at t  time_step <- 1
  a_other_ind_traits[, "death" , time_step] <- mortality_w(
    m_ind_traits_step = a_ind_traits[,, time_step],
    m_ind_traits_pStep = a_ind_traits[,, max(1, time_step - 1)],
    m_other_ind_traits_step = a_other_ind_traits[,, time_step],
    m_other_ind_traits_pStep = a_other_ind_traits[,, max(1, time_step - 1)],
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
  )
  
  # predict the risk factors for the next cycle (t+1) 
  a_ind_traits[, , time_step + 1] <-  update_all_biomarkers_w(
    m_ind_traits_init = m_ind_traits_init,
    m_ind_traits_nStep  = a_ind_traits[,, time_step + 1],
    m_ind_traits_step = a_ind_traits[,, time_step],
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
  )
  
  # # provide some feedback to the console
  # if(time_step/2 == round(time_step/2,0)) {
  #   cat('\r', paste(time_step/num_cycles * 100, "% done", sep = " "))
  # }
}

# did it work?!?
# a_ind_traits[,,7]

(proc.time() - ptm)/60
# > (proc.time() - ptm)/60 [8 cycles]
#    user  system elapsed 
# 0.10830 0.03440 0.15965 
# > (proc.time() - ptm)/60 [20 cycles]
#      user    system   elapsed 
# 0.2674667 0.0596500 0.3432000 
# > (proc.time() - ptm)/60
#       user     system    elapsed [20 cycles]
# 0.26766667 0.05363333 0.33563333

# Summary matrix
m_summary <- matrix(   
  data = NA, 
  nrow = length(v_cycle_nms), 
  ncol = 4,   
  dimnames = list(v_cycle_nms, column_names)  
) 
m_summary <- m_summary[-nrow(m_summary), ]

for (time_step in 1:num_cycles) {
  m_summary[time_step, "cost"] <- 
    as.double(m_other_ind_traits[time_step,"death"]==0)*c_baseline + 
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step,"blindness_event"] * c_blindness_e +
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step,"blindness_hist"] * c_blindness_c +
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step, "amp_event"] * c_amp_e + 
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step, "amp_event2"] * c_amp_e + 
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step, "amp_hist"] * c_amp_c + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "chf_event"] * c_chf_e + 
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step, "chf_hist"] * c_chf_c + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "esrd_event"] * c_esrd_e + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "esrd_hist"] * c_esrd_c + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "ihd_event"] * c_ihd_e + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "ihd_hist"] * c_ihd_c + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "mi_event"] * c_mi_e + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "mi_hist"] * c_mi_c + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "stroke_event"] * c_stroke_e + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "stroke_hist"] * c_stroke_c + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "ulcer_event"] * c_ulcer_e + 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "ulcer_hist"] * c_ulcer_c
  
  m_summary[time_step, "qalys"]<- as.double(m_other_ind_traits[time_step,"death"]==0)*q_baseline + min(
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step,"blindness_event"] * q_blindness ,
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step,"blindness_hist"] * q_blindness ,
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step, "amp_event"] * q_amp , 
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step, "amp_event2"] * q_amp , 
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step, "amp_hist"] * q_amp ,
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "chf_event"] * q_chf , 
    as.double(m_other_ind_traits[time_step,"death"]==0)*m_ind_traits[time_step, "chf_hist"] * q_chf ,
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "esrd_event"] * q_esrd , 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "esrd_hist"] * q_esrd , 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "ihd_event"] * q_ihd , 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "ihd_hist"] * q_ihd , 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "mi_event"] * q_mi , 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "mi_hist"] * q_mi ,
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "stroke_event"] * q_stroke , 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "stroke_hist"] * q_stroke , 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "ulcer_event"] * q_ulcer , 
    as.double(m_other_ind_traits[time_step,"death"]==0)* m_ind_traits[time_step, "ulcer_hist"] * q_ulcer )
  m_summary[time_step, "disc_costs"] <- m_summary[time_step, "cost"] / (1 + discount_rate)^time_step
  m_summary[time_step, "disc_qalys"] <- m_summary[time_step, "qalys"] / (1 + discount_rate)^time_step
}

patient_summary_file[patient,"cost"]<-sum(m_summary[,"cost"])
patient_summary_file[patient,"disc_costs"]<-sum(m_summary[,"disc_costs"])
patient_summary_file[patient,"qalys"]<-sum(m_summary[,"qalys"])
patient_summary_file[patient,"disc_qalys"]<-sum(m_summary[,"disc_qalys"])

# Step 10: Summarize and visualize results
# Calculate summary statistics
# Generate plots to visualize disease progression trends