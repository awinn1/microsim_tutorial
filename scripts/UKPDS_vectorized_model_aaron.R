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
UKPDS_coef <- readr::read_csv("data/ukpds_coef.csv")  # Load coefficient matrix from CSV

# Replace NAs with 0s to avoid missing values in calculations
UKPDS_coef[is.na(UKPDS_coef)] <- 0  

# Extract parameter names (used as row names)
v_coef_names <- UKPDS_coef$Parameter # Get row names from the 'Parameter' column

# Determine the number of parameters (rows)
n_coef_names <-  length(v_coef_names)  # Count the number of parameters

# Extract factor names (used as column names), excluding the first column
v_factors_names <- as.vector(colnames(UKPDS_coef[-1]))  # Get column names excluding 'Parameter'

# Determine the number of factors (columns)
n_equa_names <- length(v_factors_names)  # Count the number of factors

# allow for bootstrapped coefficients 
boot <- 1
rep_names <- paste0("boot_rep_", 1:boot)

#create an array that holds onto everything!
a_coef_ukpds <- array(
  data = NA,
  dim = c(n_coef_names, n_equa_names, boot),
  dimnames = list(v_coef_names, v_factors_names, rep_names)
)

dimnames(a_coef_ukpds)
#fill in the array with coefficents from the dataset
a_coef_ukpds[,1,1]<-UKPDS_coef$hba1c 
a_coef_ukpds[,2,1]<-UKPDS_coef$sbp 
a_coef_ukpds[,3,1]<-UKPDS_coef$ldl 
a_coef_ukpds[,4,1]<-UKPDS_coef$hdl 
a_coef_ukpds[,5,1]<-UKPDS_coef$bmi 
a_coef_ukpds[,6,1]<-UKPDS_coef$heart_rate 
a_coef_ukpds[,7,1]<-UKPDS_coef$wbc 
a_coef_ukpds[,8,1]<-UKPDS_coef$haem 
a_coef_ukpds[,9,1]<-UKPDS_coef$chf 
a_coef_ukpds[,10,1]<-UKPDS_coef$ihd 
a_coef_ukpds[,11,1]<-UKPDS_coef$mi1_male 
a_coef_ukpds[,12,1]<-UKPDS_coef$mi1_female 
a_coef_ukpds[,13,1]<-UKPDS_coef$mi2 
a_coef_ukpds[,14,1]<-UKPDS_coef$stroke_1 
a_coef_ukpds[,15,1]<-UKPDS_coef$stroke_2 
a_coef_ukpds[,16,1]<-UKPDS_coef$blindness 
a_coef_ukpds[,17,1]<-UKPDS_coef$ulcer 
a_coef_ukpds[,18,1]<-UKPDS_coef$amp1_no_ulcer 
a_coef_ukpds[,19,1]<-UKPDS_coef$amp1_yes_ulcer 
a_coef_ukpds[,20,1]<-UKPDS_coef$amp2 
a_coef_ukpds[,21,1]<-UKPDS_coef$esrd 

a_coef_ukpds[,22,1]<-UKPDS_coef$death_nhne 
a_coef_ukpds[,23,1]<-UKPDS_coef$death_1st_event 
a_coef_ukpds[,24,1]<-UKPDS_coef$death_yhne 
a_coef_ukpds[,25,1]<-UKPDS_coef$death_yhye 

#print(dim(a_coef_ukpds_ind_traits)) # to verify the dimensions, 65 rows, 25 columns, 1 slice
#print(dimnames(a_coef_ukpds_ind_traits)) # to verify the dimension names

a_coef_ukpds_ind_traits<- a_coef_ukpds[1:62, , "boot_rep_1", drop = FALSE]
a_coef_ukpds_other_ind_traits<- a_coef_ukpds[63:65,, ,drop = FALSE]


# Step 2: Create the patient dataset
ukpds_pop <- read_csv("data/population.csv")
# show the names of the variables and rows
print(dimnames(ukpds_pop)) 

seed    <- 1234                     # random number generator state
num_i <- 250000                          # number of simulated individuals
# Define the number of time points
num_cycles <- 20                    # maximum length of a simulation 
set.seed(seed)    # set the seed to ensure reproducible samples below
ids <- paste("id",   1:num_i,    sep ="_")
cycles <- paste("cycle", 0:num_cycles, sep ="_")

nbr_coef_names <- n_coef_names
# Create AN ARRAY with columns for each variable, row for each person, and a 
# slice for each period
a_all_ind_traits <- array(   
  data = NA, 
  dim = c(num_i, nbr_coef_names, num_cycles+1),
  dimnames = list(ids, v_coef_names , cycles )  
)
# need this to be the same number of columns as the coefficient table is long/rows
print(dim(a_all_ind_traits)) # to verify the dimensions
print(dimnames(a_all_ind_traits)) # to verify the dimension names

#Fill patient trace 
a_all_ind_traits[,"age", 1] <- ukpds_pop$age
a_all_ind_traits[,"age_diag", ] <- ukpds_pop$age_diag
a_all_ind_traits[,"black", ] <- ukpds_pop$black
a_all_ind_traits[,"indian", ] <- ukpds_pop$indian
a_all_ind_traits[,"female", ] <- ukpds_pop$female
a_all_ind_traits[,"diab_dur", 1] <- ukpds_pop$diab_dur
a_all_ind_traits[,"diab_dur_log", 1] <- ukpds_pop$diab_dur_log
a_all_ind_traits[,"smoke", ] <- ukpds_pop$smoke
a_all_ind_traits[,"a1c", 1] <- ukpds_pop$a1c
a_all_ind_traits[,"a1c_lag", 1] <- ukpds_pop$a1c_lag
a_all_ind_traits[,"a1c_first", ] <- ukpds_pop$a1c_first
a_all_ind_traits[,"bmi", 1] <- ukpds_pop$bmi
a_all_ind_traits[,"bmi_lt18_5", 1] <- ukpds_pop$bmi_lt18_5
a_all_ind_traits[,"bmi_gte25", 1] <- ukpds_pop$bmi_gte25
a_all_ind_traits[,"bmi_lag", 1] <- ukpds_pop$bmi_lag
a_all_ind_traits[,"bmi_first", ] <- ukpds_pop$bmi_first
a_all_ind_traits[,"egfr", 1] <- ukpds_pop$egfr
a_all_ind_traits[,"egfr_lt60", 1] <- ukpds_pop$egfr_lt60
a_all_ind_traits[,"egfr_gte60", 1] <- ukpds_pop$egfr_gte60
a_all_ind_traits[,"egfr_real", 1] <- ukpds_pop$egfr_real
a_all_ind_traits[,"hdl", 1] <- ukpds_pop$hdl
a_all_ind_traits[,"hdl_lag", 1] <- ukpds_pop$hdl_lag
a_all_ind_traits[,"hdl_first", ] <- ukpds_pop$hdl_first
a_all_ind_traits[,"hdl_real", 1] <- ukpds_pop$hdl_real
a_all_ind_traits[,"heart_rate", 1] <- ukpds_pop$heart_rate
a_all_ind_traits[,"heart_rate_lag", 1] <- ukpds_pop$heart_rate_lag
a_all_ind_traits[,"heart_rate_first", ] <- ukpds_pop$heart_rate_first
a_all_ind_traits[,"heart_rate_real", 1] <- ukpds_pop$heart_rate_real
a_all_ind_traits[,"ldl", 1] <- ukpds_pop$ldl
a_all_ind_traits[,"ldl_gt35", 1] <- ukpds_pop$ldl_gt35
a_all_ind_traits[,"ldl_lag", 1] <- ukpds_pop$ldl_lag
a_all_ind_traits[,"ldl_first", ] <- ukpds_pop$ldl_first
a_all_ind_traits[,"ldl_real", 1] <- ukpds_pop$ldl_real
a_all_ind_traits[,"albumin_mm", 1] <- ukpds_pop$albumin_mm
a_all_ind_traits[,"sbp", 1] <- ukpds_pop$sbp
a_all_ind_traits[,"sbp_lag", 1] <- ukpds_pop$sbp_lag
a_all_ind_traits[,"sbp_first", ] <- ukpds_pop$sbp_first
a_all_ind_traits[,"sbp_real", 1] <- ukpds_pop$sbp_real
a_all_ind_traits[,"wbc", 1] <- ukpds_pop$wbc
a_all_ind_traits[,"wbc_lag", 1] <- ukpds_pop$wbc_lag
a_all_ind_traits[,"wbc_first", ] <- ukpds_pop$wbc_first
a_all_ind_traits[,"heamo", 1] <- ukpds_pop$heamo
a_all_ind_traits[,"heamo_first", ] <- ukpds_pop$heamo_first
a_all_ind_traits[,"amp_event", 1] <- ukpds_pop$amp_event
a_all_ind_traits[,"amp_event2", 1] <- ukpds_pop$amp_event2
a_all_ind_traits[,"amp_hist", 1] <- ukpds_pop$amp_hist
a_all_ind_traits[,"atria_fib", 1] <- ukpds_pop$atria_fib
a_all_ind_traits[,"blindness_event", 1] <- ukpds_pop$blindness_event
a_all_ind_traits[,"blindness_hist", 1] <- ukpds_pop$blindness_hist
a_all_ind_traits[,"chf_event", 1] <- ukpds_pop$chf_event
a_all_ind_traits[,"chf_hist", 1] <- ukpds_pop$chf_hist
a_all_ind_traits[,"esrd_event", 1] <- ukpds_pop$esrd_event
a_all_ind_traits[,"esrd_hist", 1] <- ukpds_pop$esrd_hist
a_all_ind_traits[,"ihd_event", 1] <- ukpds_pop$ihd_event
a_all_ind_traits[,"ihd_hist", 1] <- ukpds_pop$ihd_hist
a_all_ind_traits[,"mi_event", 1] <- ukpds_pop$mi_event
a_all_ind_traits[,"mi_hist", 1] <- ukpds_pop$mi_hist
a_all_ind_traits[,"pvd_event", 1] <- ukpds_pop$pvd_event
a_all_ind_traits[,"stroke_event", 1] <- ukpds_pop$stroke_event
a_all_ind_traits[,"stroke_hist", 1] <- ukpds_pop$stroke_hist
a_all_ind_traits[,"ulcer_event", 1] <- ukpds_pop$ulcer_event
a_all_ind_traits[,"ulcer_hist", 1] <- ukpds_pop$ulcer_hist
a_all_ind_traits[,"lambda", 1] <- ukpds_pop$lambda
a_all_ind_traits[,"rho", 1] <- ukpds_pop$rho
a_all_ind_traits[,"death", 1] <- ukpds_pop$death

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
#' This function calculates patient-specific factors to predict the time path of a biomarker. 
#'
#' @param a_ind_traits An array containing patient characteristics over time.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param biomarker_eq A character string specifying the health outcome equation (e.g., "ihd").
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#' 
#' @return The updated biomarker is stored.
#' @export
biomarker <- function(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq,  time_step) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  updated_biomarker <- (a_ind_traits[,,time_step] %*%  a_coef_ukpds_ind_traits[,  biomarker_eq, 1] + 
                          a_coef_ukpds_other_ind_traits["lambda",  biomarker_eq, 1] )
  
  return(updated_biomarker)
}



# Step 4: Create a function to apply all risk factor models ####
# Combine risk factor functions into a single pipeline
# Update patient data over time

#' Update Multiple Biomarkers in a Transition Matrix
#'
#' This function updates multiple biomarker values in the transition matrix for a given time step.
#'
#' @param a_ind_traits The patient trace, an array containing patient data with biomarker and event columns.
#' @param a_coef_ukpds_ind_traits A coefficient matrix containing biomarker and event equations.
#' @param time_step An integer representing the current time step.
#' @param next_row An integer indicating the row in `m_ind_traits` to update with new biomarker values.
#'
#' @return The updated transition matrix `m_ind_traits` with new biomarker values in the specified row.
#' 
#' @examples
#' # Example usage
#' m_ind_traits <- update_all_biomarkers(m_ind_traits, a_coef_ukpds_ind_traits, time_step = 1, next_row = 2)
#'
#' @export

update_all_biomarkers <- function(a_ind_traits, a_coef_ukpds_ind_traits, time_step, next_row) {
  # predict the next period (and perform transformations as needed)
  # the biomarkers use real values of variables, but the event equations use transformed variables
  a_ind_traits[, "a1c", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "hba1c", time_step = time_step)   
  a_ind_traits[, "sbp_real", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "sbp", time_step = time_step)   
  a_ind_traits[, "sbp", next_row] <- a_ind_traits[, "sbp_real", next_row] / 10   
  a_ind_traits[, "ldl_real", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "ldl", time_step = time_step)   
  a_ind_traits[, "ldl", next_row] <- a_ind_traits[, "ldl_real", next_row] * 10   
  a_ind_traits[, "hdl_real", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "hdl", time_step = time_step)   
  a_ind_traits[, "hdl", next_row] <- a_ind_traits[, "hdl_real", next_row] * 10    
  a_ind_traits[, "bmi", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "bmi", time_step = time_step)   
  a_ind_traits[, "heart_rate_real", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "heart_rate", time_step = time_step)   
  a_ind_traits[, "heart_rate", next_row] <- a_ind_traits[, "heart_rate_real", next_row] / 10   
  a_ind_traits[, "wbc", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "wbc", time_step = time_step)   
  a_ind_traits[, "heamo", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "haem", time_step = time_step)         
  
  # Update lag and first occurrence columns   
  a_ind_traits[, "a1c_lag", next_row] <- a_ind_traits[, "a1c", time_step]   
  a_ind_traits[, "bmi_lag", next_row] <- a_ind_traits[, "bmi", time_step]   
  a_ind_traits[, "bmi_lt18_5", next_row] <- as.integer(a_ind_traits[, "bmi", next_row] < 18.5)   
  a_ind_traits[, "bmi_gte25", next_row] <- as.integer(a_ind_traits[, "bmi", next_row] >= 25)         
  
  a_ind_traits[, "hdl_lag", next_row] <- a_ind_traits[, "hdl_real", time_step]      
  a_ind_traits[, "heart_rate_lag", next_row] <- a_ind_traits[, "heart_rate_real", time_step]   
  
  # Check if this is functioning as a spline   
  a_ind_traits[, "ldl_gt35", next_row] <- as.integer(a_ind_traits[, "ldl_real", next_row] > 35) / 10   
  a_ind_traits[, "ldl_lag", next_row] <- a_ind_traits[, "ldl_real", time_step]   
  a_ind_traits[, "sbp_lag", next_row] <- a_ind_traits[, "sbp_real", time_step]   
  a_ind_traits[, "wbc_lag", next_row] <- a_ind_traits[, "wbc", time_step]      
  
  # Update additional values   
  a_ind_traits[, "egfr", next_row] <- a_ind_traits[, "egfr", 1]   
  a_ind_traits[, "egfr_real", next_row] <- a_ind_traits[, "egfr_real", 1]   
  a_ind_traits[, "egfr_lt60", next_row] <- a_ind_traits[, "egfr_lt60", 1]   
  a_ind_traits[, "egfr_gte60", next_row] <- a_ind_traits[, "egfr_gte60", 1]   
  a_ind_traits[, "albumin_mm", next_row] <- a_ind_traits[, "albumin_mm", 1]  
  
  
  # Return updated matrix
  return(a_ind_traits)
}


# Step 5: Define event functions (Weibull/Exponential and Logistic) ####
# Weibull distribution function for event occurrence
# Logistic regression for binary event prediction


#' Calculate Transition Probability Based on a Weibull Model and Update Patient State
#' Note: An exponential model is a special case of the Weibull model where the shape 
#' parameter (ρ) is set to 1, meaning the hazard function remains constant over time, 
#' resulting in a constant rate of event occurrence rather than a time-dependent rate.
#'
#' This function calculates patient-specific factors, cumulative hazards, 
#' and the transition probability for a given health outcome (e.g., "ihd"). 
#' The function updates the provided `m_ind_traits` matrix with the event occurrence 
#' at the specified time step.
#'
#' @param m_ind_traits A matrix containing patient characteristics over time.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#' 
#' @return Whether the event occurred.
#' @export
weibull_event <- function(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome, health_event, time_step) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (a_ind_traits[,,time_step] %*%  a_coef_ukpds_ind_traits[, health_outcome, 1] + 
                        as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1]) )
  
  # Compute cumulative hazard at the current time step
  cum_hazard_t <- exp(patient_factors) * (a_ind_traits[, "diab_dur",time_step]^(a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]) )
  
  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- exp(patient_factors) * ((a_ind_traits[, "diab_dur",time_step] + 1)^(a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]) )
  
  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(nrow(a_ind_traits))
  
  # Return the updated matrix
  return(event)
}




#' Calculate Transition Probability Based on a Logistic Regression and Update Patient State 
#'
#' This function calculates patient-specific factors, cumulative hazards, 
#' and the transition probability for a given health outcome (e.g., "ihd"). 
#' The function updates the provided `m_ind_traits` matrix with the event occurrence 
#' at the specified time step.
#'
#' @param a_ind_traits A 3d array containing patient characteristics over time (slices).
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#' 
#' @return Whether the event occurred.
#' @export
logistic_event <- function(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome, health_event, time_step) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (a_ind_traits[,,time_step] %*%  a_coef_ukpds_ind_traits[, health_outcome, 1] + 
                        as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1]) )
  
  # Calculate transition probability
  trans_prob=1-(exp(-patient_factors)/(1+exp(-patient_factors)))^1
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(nrow(a_ind_traits))
  
  
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
  event <- trans_prob > runif(nrow(a_ind_traits))
  
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
# for (time_step in 2:num_cycles)
for (time_step in 2:5) {
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
# a_ind_traits[,,7]

(proc.time() - ptm)/60

# m_ind_traits_new <- m_ind_traits[-nrow(m_ind_traits), ]



m_summary <- matrix(   
  data = NA, 
  nrow = length(cycles), 
  ncol = 4,   
  dimnames = list(cycles,column_names)  
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