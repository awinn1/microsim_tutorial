mem.maxVSize()

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
  dim = c(num_i, n_coef, num_cycles+1),
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

#' Calculate Biomarkers 
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
biomarker <- function(
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

#' Update Multiple Biomarkers in a Transition Matrix
#'
#' @details
#' This function updates multiple biomarker values in the transition matrix for
#' a given time step.
#'
#' @param m_ind_traits_init A matrix containing patient characteristics at the
#' initial time step (obtained from `a_ind_traits[,,1]`).
#' @param m_ind_traits_row The patient matrix to update with new biomarker
#' values (obtained from `a_ind_traits[,,next_row]`).
#' @param m_ind_traits_step The patient matrix containing patient data at the
#' current time step (obtained from `a_ind_traits[,,time_step]`).
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for
#' calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of other coefficients used
#' for calculating risk (e.g., lambda).
#'
#' @return The updated transition matrix `m_ind_traits_row` with new biomarker
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
#' m_ind_traits_row <- a_ind_traits[,,next_row]
#' m_ind_traits_step <- a_ind_traits[,,time_step]
#'
#' # Assuming a_coef_ukpds_ind_traits and a_coef_ukpds_other_ind_traits are defined
#'
#' updated_m_ind_traits_row <- update_all_biomarkers(
#'   m_ind_traits_init = m_ind_traits_init,
#'   m_ind_traits_row = m_ind_traits_row,
#'   m_ind_traits_step = m_ind_traits_step,
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
#' )
#' }
update_all_biomarkers <- function(
    m_ind_traits_init,
    m_ind_traits_row,
    m_ind_traits_step,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits) {
  # predict the next period (and perform transformations as needed)
  # the biomarkers use real values of variables, but the event equations use transformed variables
  m_ind_traits_row[, "a1c"] <- biomarker(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "hba1c")
  m_ind_traits_row[, "sbp_real"] <- biomarker(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "sbp")
  m_ind_traits_row[, "sbp"] <- m_ind_traits_row[, "sbp_real"] / 10
  m_ind_traits_row[, "ldl_real"] <- biomarker(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "ldl")
  m_ind_traits_row[, "ldl"] <- m_ind_traits_row[, "ldl_real"] * 10
  m_ind_traits_row[, "hdl_real"] <- biomarker(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "hdl")
  m_ind_traits_row[, "hdl"] <- m_ind_traits_row[, "hdl_real"] * 10
  m_ind_traits_row[, "bmi"] <- biomarker(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "bmi")
  m_ind_traits_row[, "heart_rate_real"] <- biomarker(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "heart_rate")
  m_ind_traits_row[, "heart_rate"] <- m_ind_traits_row[, "heart_rate_real"] / 10
  m_ind_traits_row[, "wbc"] <- biomarker(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "wbc")
  m_ind_traits_row[, "heamo"] <- biomarker(m_ind_traits_step, a_coef_ukpds_ind_traits, a_coef_ukpds_other_ind_traits, biomarker_eq = "haem")
  
  # Update lag and first occurrence columns
  m_ind_traits_row[, "a1c_lag"] <- m_ind_traits_step[, "a1c"]
  m_ind_traits_row[, "bmi_lag"] <- m_ind_traits_step[, "bmi"]
  m_ind_traits_row[, "bmi_lt18_5"] <- as.integer(m_ind_traits_row[, "bmi"] < 18.5)
  m_ind_traits_row[, "bmi_gte25"] <- as.integer(m_ind_traits_row[, "bmi"] >= 25)
  
  m_ind_traits_row[, "hdl_lag"] <- m_ind_traits_step[, "hdl_real"]
  m_ind_traits_row[, "heart_rate_lag"] <- m_ind_traits_step[, "heart_rate_real"]
  
  # Check if this is functioning as a spline
  m_ind_traits_row[, "ldl_gt35"] <- as.integer(m_ind_traits_row[, "ldl_real"] > 35) / 10
  m_ind_traits_row[, "ldl_lag"] <- m_ind_traits_step[, "ldl_real"]
  m_ind_traits_row[, "sbp_lag"] <- m_ind_traits_step[, "sbp_real"]
  m_ind_traits_row[, "wbc_lag"] <- m_ind_traits_step[, "wbc"]
  
  # Update additional values
  m_ind_traits_row[, "egfr"] <- m_ind_traits_init[, "egfr"]
  m_ind_traits_row[, "egfr_real"] <- m_ind_traits_init[, "egfr_real"]
  m_ind_traits_row[, "egfr_lt60"] <- m_ind_traits_init[, "egfr_lt60"]
  m_ind_traits_row[, "egfr_gte60"] <- m_ind_traits_init[, "egfr_gte60"]
  m_ind_traits_row[, "albumin_mm"] <- m_ind_traits_init[, "albumin_mm"]
  
  # Return updated matrix
  return(m_ind_traits_row)
}

# Step 5: Define event functions (Weibull/Exponential and Logistic) ####
# Weibull distribution function for event occurrence
# Logistic regression for binary event prediction

#' Calculate Transition Probability Based on a Weibull Model and Update Patient
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
weibull_event <- function(
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
  event <- trans_prob > runif(1)
  
  # Return the updated matrix
  return(event)
}

#' Calculate Transition Probability Based on a Logistic Regression and Update Patient State 
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
#' health_outcome <- "ihd"
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
logistic_event <- function(
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
  event <- trans_prob > runif(1)
  
  # Return the value
  return(event)
}

#
# u <- logistic_event(m_ind_traits, a_coef_ukpds_ind_traits, health_outcome = "ulcer", health_event = "ulcer_event", time_step = 1)
#m_ind_traits <- weibull_event(m_ind_traits, a_coef_ukpds_ind_traits, health_outcome = "ihd", health_event = "ihd_event", time_step = 1)

# Step 6: Initialize event and history variables ####
#' @title Update Health Events Over Time Steps
#' @description This function updates health events in a patient data matrix (`m_ind_traits`) by applying Weibull 
#' and logistic event functions in a randomized order across multiple time steps.
#'
#' @param a_ind_traits An array containing patient-level data, including health event history.
#' @param a_coef_ukpds_ind_traits A coefficient matrix used in Weibull and logistic event calculations.
#' @param time_step An integer indicating the current time step to update events.
#'
#' @return Updated `m_ind_traits` matrix with event and history values updated for the given time step.
#'
#' @examples
#' \dontrun{
#' a_ind_traits <- update_health_events(m_ind_traits, a_coef_ukpds_ind_traits, time_step = 1)
#' }
#' 
#' @export
update_health_events <- function(a_ind_traits, a_coef_ukpds_ind_traits, time_step) {
  # Ensure m_ind_traits remains a matrix
  if (!is.array(a_ind_traits)) {
    stop("m_ind_traits must be an array.")
  }

  # Initialize event variables and update history
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")
  
  # create event and history column names once and save each group in a vector
  v_event_cols      <- paste0(events, "_event")
  v_history_cols    <- paste0(events, "_hist")
  
  # Update history columns in one vectorized call
  a_ind_traits[, v_event_cols , time_step] <- 0
  a_ind_traits[, v_history_cols, time_step] <- pmax(
    a_ind_traits[, v_history_cols, max(1, time_step - 1)],
    a_ind_traits[, v_event_cols, max(1, time_step - 1)] 
  )
  
  
  a_ind_traits[, "atria_fib", time_step] <- a_ind_traits[, "atria_fib", 1]
  a_ind_traits[, "pvd_event", time_step] <- 0
  a_ind_traits[, "amp_event2", time_step] <- 0
  # Randomize event order
  randomized_events <- sample(events)
  
  for (events in randomized_events) {
    if (events == "amp") {
      amp1_no_ulcer <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                                     health_outcome = "amp1_no_ulcer", health_event = "amp_event", time_step = time_step)
      amp1_yes_ulcer <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                                      health_outcome = "amp1_yes_ulcer", health_event = "amp_event", time_step = time_step)
      
      a_ind_traits[, "amp_event", time_step] <- 
        (amp1_no_ulcer * (a_ind_traits[, "ulcer_hist", time_step] == 0)) + 
        (amp1_yes_ulcer * (a_ind_traits[, "ulcer_hist", time_step] == 1))
      
      
      # Ensure that this is a new event
      a_ind_traits[, "amp_event", time_step] <- 
        a_ind_traits[, "amp_event", time_step] * (a_ind_traits[, "amp_hist", time_step] == 0)
      
      # Calculate amp2 event
      amp2 <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                            health_outcome = "amp2", health_event = "amp_event2", time_step = time_step)
      
      a_ind_traits[, "amp_event2", time_step] <- 0 
      
      a_ind_traits[, "amp_event2", time_step] <- 
        amp2 * (a_ind_traits[, "amp_hist", time_step] == 1)
      
      
      
      
    } else if (events == "mi") {
      mi1_male <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                                health_outcome = "mi1_male", health_event = "mi_event", time_step = time_step)
      mi1_female <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                                  health_outcome = "mi1_female", health_event = "mi_event", time_step = time_step)
      
      a_ind_traits[, "mi_event", time_step] <- 
        (mi1_male * (a_ind_traits[, "female", time_step] == 0)) + 
        (mi1_female * (a_ind_traits[, "female", time_step] == 1))
      
      a_ind_traits[, "mi_event", time_step] <- 
        a_ind_traits[, "mi_event", time_step] * (a_ind_traits[, "mi_hist", time_step] == 0)
      
      mi2 <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                           health_outcome = "mi2", health_event = "mi_event", time_step = time_step)
      
      a_ind_traits[, "mi_event", time_step] <- 
        (a_ind_traits[, "mi_hist", time_step] == 0) * a_ind_traits[, "mi_event", time_step] + 
        (a_ind_traits[, "mi_hist", time_step] == 1) * mi2
      
    } else if (events == "stroke") {
      stroke1 <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                               health_outcome = "stroke_1", health_event = "stroke_event", time_step = time_step)
      stroke2 <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                               health_outcome = "stroke_2", health_event = "stroke_event", time_step = time_step)
      
      a_ind_traits[, "stroke_event", time_step] <- 
        (stroke1 * (a_ind_traits[, "stroke_hist", time_step] == 0)) + 
        (stroke2 * (a_ind_traits[, "stroke_hist", time_step] == 1))
      
      a_ind_traits[, "stroke_event", time_step] <- 
        a_ind_traits[, "stroke_event", time_step] * (a_ind_traits[, "stroke_hist", time_step] == 0)
      
    } else if (events == "ulcer") {
      a_ind_traits[, "ulcer_event", time_step] <- 
        logistic_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                       health_outcome = "ulcer", health_event = "ulcer_event", time_step = time_step)
      
      a_ind_traits[, "ulcer_event", time_step] <- 
        a_ind_traits[, "ulcer_event", time_step] * (a_ind_traits[, "ulcer_hist", time_step] == 0)
      
    } else {
      a_ind_traits[, paste0(events, "_event"), time_step] <- 
        weibull_event(a_ind_traits, a_coef_ukpds_ind_traits, 
                      health_outcome = events, health_event = paste0(events, "_event"), time_step = time_step)
      
      a_ind_traits[, paste0(events, "_event"), time_step] <- 
        a_ind_traits[, paste0(events, "_event"), time_step] * 
        (a_ind_traits[, paste0(events, "_hist"), time_step] == 0)
    }
  }
  
  return(a_ind_traits)
}


# Step 7: Define a mortality function ####
# Combine relevant event functions affecting mortality
# Estimate survival probability over time
#' Calculate Transition Probability Based on a Gompertz Model and Update Patient State
#'
#'
#' This function calculates patient-specific factors, cumulative hazards, 
#' and the transition probability for mortality. 
#' The function updates the provided `m_ind_traits` matrix with the event occurrence 
#' at the specified time step.
#'
#' @param a_ind_traits An array containing patient characteristics over time.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#' 
#' @return The event occurrence stored.
#' @export
gompertz_event <- function(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome, health_event, time_step) {
  
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (a_ind_traits[, ,time_step] %*%  a_coef_ukpds_ind_traits[, health_outcome, 1] + 
                        as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1]) )
  
  # Compute cumulative hazard at the current time step
  cum_hazard_t <- (1/a_coef_ukpds_other_ind_traits["rho", health_outcome, 1])* exp(patient_factors) * (exp(a_ind_traits[, "age" , time_step]*(a_coef_ukpds_other_ind_traits["rho", health_outcome, 1])) -1 )
  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- (1/a_coef_ukpds_other_ind_traits["rho", health_outcome, 1])* exp(patient_factors) * (exp((a_ind_traits[ , "age" , time_step]+1)*(a_coef_ukpds_other_ind_traits["rho", health_outcome, 1])) -1 )
  
  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(1)
  
  # Return the updated matrix
  return(event)
}


#' @title Calculate Mortality Events for a Given Time Step
#' 
#' @description This function calculates mortality events for a given time step 
#' based on new health events and medical history using Gompertz and logistic models.
#'
#' @param a_ind_traits An array containing patient-level data, including health event history.
#' @param a_other_ind_traits An array containing lambda, rhos and death. 
#' @param a_coef_ukpds_ind_traits A coefficient matrix used in Gompertz and logistic event calculations.
#' @param time_step An integer specifying the time step at which mortality should be calculated.
#'
#' @return The updated `m_ind_traits` matrix with the mortality status recorded for the specified time step.
#'
#' @examples
#' \dontrun{
#' m_ind_traits <- mortality(m_ind_traits, a_coef_ukpds_ind_traits, time_step = 5)
#' }
#'
#' @export
mortality <- function(a_ind_traits, a_other_ind_traits, a_coef_ukpds_ind_traits, time_step) {
  
  # Calculate new health event occurrence and prior history
  
  # Define events of interest
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")
  
  # Create vectors containing events and event-history names:
  v_event_cols <- paste0(events, "_event")
  v_hist_cols  <- paste0(events, "_hist")
  
  # Get the maximum across those columns, for the given time_step
  # Calculate any new health event
  new_event  <- max(a_ind_traits[, v_event_cols,time_step])
  # Calculate any prior history of health events
  any_history <- max(a_ind_traits[, v_hist_cols, time_step])  

  
  # Determine event-history combinations
  nhne <- new_event == 0 & any_history == 0  # No history, no event
  yhne <- new_event == 0 & any_history == 1  # Yes history, no event
  nhye <- new_event == 1 & any_history == 0  # No history, new event
  yhye <- new_event == 1 & any_history == 1  # Yes history, new event
  
  # Mortality calculations using Gompertz and logistic models
  death_nhne <- gompertz_event(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome = "death_nhne", 
                               health_event = "death_nhne", time_step = time_step)
  
  death_yhne <- gompertz_event(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome = "death_yhne", 
                               health_event = "death_yhne", time_step = time_step)
  
  death_nhye <- logistic_event(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome = "death_1st_event", 
                               health_event = "death_nhye", time_step = time_step)
  
  death_yhye <- logistic_event(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome = "death_yhye", 
                               health_event = "death_yhye", time_step = time_step)
  
  # Calculate new mortality status
  new_death <- nhne * death_nhne + yhne * death_yhne + nhye * death_nhye + yhye * death_yhye
  
  # Update the mortality status in the matrix for the given time step
  a_other_ind_traits[, "death" , time_step] <- new_death + a_other_ind_traits[, "death", max(time_step - 1, 1)]
  
  return(a_other_ind_traits)
}


# Step 8: Simulate disease progression and mortality for all patients ####
# Initialize patient data
# Loop through time points to update risk factors and events
# Store results
# discount rate 


# Step 9: Simulate disease progression and mortality for 999 additional patients
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
  dimnames = list(c(1:250000),column_names)  
)

ptm <- proc.time()

# predict 2nd period biomarkers
a_ind_traits<- update_all_biomarkers(a_ind_traits, a_coef_ukpds_ind_traits, 
                                     time_step = 1, next_row = 2) 

#print(patient)
#create a patient population 
for (time_step in 4:num_cycles) {
  a_ind_traits[ ,"age",time_step]<-a_ind_traits[,"age", max(1,time_step-1)] +1
  a_ind_traits[ ,"diab_dur",time_step]<-a_ind_traits[,"diab_dur" , max(1,time_step-1)]+1    
  a_ind_traits[,"diab_dur_log",time_step]<- (log(a_ind_traits[,"diab_dur", time_step]))
  
  # a_ind_traits<- update_all_biomarkers(a_ind_traits, a_coef_ukpds_ind_traits, 
  #                                      time_step = time_step, next_row = time_step+1) 
  
  # ready to simulate 
  # event prediction at t
  a_ind_traits <- update_health_events(a_ind_traits, a_coef_ukpds_ind_traits, 
                                       time_step = time_step)
  # mortality prediction at t

  a_other_ind_traits <- mortality(a_ind_traits, a_other_ind_traits, a_coef_ukpds_ind_traits, time_step = time_step)
  #predict the risk factors for the next cycle (t+1) 

  a_ind_traits<- update_all_biomarkers(a_ind_traits, a_coef_ukpds_ind_traits, time_step = time_step, next_row = time_step+1) 

}

# did it work?!?
a_ind_traits[,,7]

(proc.time() - ptm)/60

 # m_ind_traits_new <- m_ind_traits[-nrow(m_ind_traits), ]
  
 

  m_summary <- matrix(   
    data = NA, 
    nrow = length(v_cycle_nms), 
    ncol = 4,   
    dimnames = list(v_cycle_nms,column_names)  
  ) 
  m_summary <- m_summary[-nrow(m_summary), ]

  for (time_step in 1:num_cycles) {
    m_summary[time_step, "cost"]<- as.double(m_other_ind_traits[time_step,"death"]==0)*c_baseline + 
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
}




# Step 10: Summarize and visualize results
# Calculate summary statistics
# Generate plots to visualize disease progression trends
