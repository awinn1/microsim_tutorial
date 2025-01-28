
library(pacman)
pacman::p_load(haven,
               Rcpp,
               RcppArmadillo,
               readr,
               dplyr,
               profvis)


# Step 1: Import the matrix of coefficients ####
# Load necessary libraries
# Read the coefficient matrix from a CSV or RData file
UKPDS_coef <- read_csv("data/ukpds_coef.csv")  # Load coefficient matrix from CSV

# Replace NAs with 0s to avoid missing values in calculations
UKPDS_coef[is.na(UKPDS_coef)] <- 0  

# Extract parameter names (used as row names)
v_coef_names <- as.vector(UKPDS_coef$Parameter)  # Get row names from the 'Parameter' column

# Determine the number of parameters (rows)
nbr_coef_names <- as.numeric(length(v_coef_names))  # Count the number of parameters

# Extract factor names (used as column names), excluding the first column
v_factors_names <- as.vector(colnames(UKPDS_coef[-1]))  # Get column names excluding 'Parameter'

# Determine the number of factors (columns)
nbr_equa_names <- as.numeric(length(v_factors_names))  # Count the number of factors

# allow for bootstrapped coefficients 
boot <- 1
rep_names <- paste0("boot_rep_", 1:boot)

#create an array that holds onto everything!
a_coef_C <- array(
  data = NA,
  dim = c(nbr_coef_names, nbr_equa_names, boot),
  dimnames = list(v_coef_names, v_factors_names, rep_names)
)
# what is the faster way to do this? 
#fill in the array with coefficents from the dataset
a_coef_C[,1,1]<-UKPDS_coef$hba1c
a_coef_C[,2,1]<-UKPDS_coef$sbp
a_coef_C[,3,1]<-UKPDS_coef$ldl
a_coef_C[,4,1]<-UKPDS_coef$hdl
a_coef_C[,5,1]<-UKPDS_coef$bmi
a_coef_C[,6,1]<-UKPDS_coef$heart_rate
a_coef_C[,7,1]<-UKPDS_coef$wbc
a_coef_C[,8,1]<-UKPDS_coef$haem
a_coef_C[,9,1]<-UKPDS_coef$chf
a_coef_C[,10,1]<-UKPDS_coef$ihd
a_coef_C[,11,1]<-UKPDS_coef$mi1_male
a_coef_C[,12,1]<-UKPDS_coef$mi1_female
a_coef_C[,13,1]<-UKPDS_coef$mi2
a_coef_C[,14,1]<-UKPDS_coef$stroke_1
a_coef_C[,15,1]<-UKPDS_coef$stroke_2
a_coef_C[,16,1]<-UKPDS_coef$blindness
a_coef_C[,17,1]<-UKPDS_coef$ulcer
a_coef_C[,18,1]<-UKPDS_coef$amp1_no_ulcer
a_coef_C[,19,1]<-UKPDS_coef$amp1_yes_ulcer
a_coef_C[,20,1]<-UKPDS_coef$amp2
a_coef_C[,21,1]<-UKPDS_coef$esrd
	
a_coef_C[,22,1]<-UKPDS_coef$death_nhne
a_coef_C[,23,1]<-UKPDS_coef$death_1st_event
a_coef_C[,24,1]<-UKPDS_coef$death_yhne
a_coef_C[,25,1]<-UKPDS_coef$death_yhye

print(dim(a_coef_C)) # to verify the dimensions, 65 rows, 25 columns, 1 slice
print(dimnames(a_coef_C)) # to verify the dimension names

# Step 2: Create the patient dataset
ukpds_pop <- read_csv("data/population.csv")
# show the names of the variables and rows
print(dimnames(ukpds_pop)) 

seed    <- 1234                     # random number generator state
num_i <- 4                          # number of simulated individuals
# Define the number of time points
num_cycles <- 20                    # maximum length of a simulation 
set.seed(seed)    # set the seed to ensure reproducible samples below
ids <- paste("id",   1:num_i,    sep ="_")
cycles <- paste("cycle", 0:num_cycles, sep ="_")

# Create a matrix with columns for each variable
m_TR <- matrix(   
  data = NA, 
  nrow = length(cycles), 
  ncol = nbr_coef_names,   
  dimnames = list(cycles,v_coef_names)  
)

# need this to be the same number of columns as the coefficient table is long/rows
print(dim(m_TR)) # to verify the dimensions
print(dimnames(m_TR)) # to verify the dimension names

#which patient to simulate
#' Initialize baseline values for multiple patients
#'
#' @param num_patients The total number of patients to process.
#' @param ukpds_pop A data frame containing patient characteristics.
#' @param m_TR A matrix to store patient data.
#' @return The updated matrix with initialized patient data.
#' @export
initialize_patients <- function(num_patients, ukpds_pop, m_TR) {
#  patient<- 1
  patient<- num_patients 
    m_TR[patient, "age"] <- as.matrix(ukpds_pop[patient, "age"])
    m_TR[patient, "age_diag"] <- as.matrix(ukpds_pop[patient, "age_diag"])
    m_TR[patient, "black"] <- as.matrix(ukpds_pop[patient, "black"])
    m_TR[patient, "indian"] <- as.matrix(ukpds_pop[patient, "indian"])
    m_TR[patient, "female"] <- as.matrix(ukpds_pop[patient, "female"])
    m_TR[patient, "diab_dur"] <- as.matrix(ukpds_pop[patient, "diab_dur"])
    m_TR[patient, "diab_dur_log"] <- as.matrix(ukpds_pop[patient, "diab_dur_log"])
    m_TR[patient, "smoke"] <- as.matrix(ukpds_pop[patient, "smoke"])
    m_TR[patient, "a1c"] <- as.matrix(ukpds_pop[patient, "a1c"])
    m_TR[patient, "a1c_lag"] <- as.matrix(ukpds_pop[patient, "a1c_lag"])
    m_TR[patient, "a1c_first"] <- as.matrix(ukpds_pop[patient, "a1c_first"])
    m_TR[patient, "bmi"] <- as.matrix(ukpds_pop[patient, "bmi"])
    m_TR[patient, "bmi_lt18_5"] <- as.matrix(ukpds_pop[patient, "bmi_lt18_5"])
    m_TR[patient, "bmi_gte25"] <- as.matrix(ukpds_pop[patient, "bmi_gte25"])
    m_TR[patient, "bmi_lag"] <- as.matrix(ukpds_pop[patient, "bmi_lag"])
    m_TR[patient, "bmi_first"] <- as.matrix(ukpds_pop[patient, "bmi_first"])
    m_TR[patient, "egfr"] <- as.matrix(ukpds_pop[patient, "egfr"])
    m_TR[patient, "egfr_lt60"] <- as.matrix(ukpds_pop[patient, "egfr_lt60"])
    m_TR[patient, "egfr_gte60"] <- as.matrix(ukpds_pop[patient, "egfr_gte60"])
    m_TR[patient, "egfr_real"] <- as.matrix(ukpds_pop[patient, "egfr_real"])
    m_TR[patient, "hdl"] <- as.matrix(ukpds_pop[patient, "hdl"])
    m_TR[patient, "hdl_lag"] <- as.matrix(ukpds_pop[patient, "hdl_lag"])
    m_TR[patient, "hdl_first"] <- as.matrix(ukpds_pop[patient, "hdl_first"])
    m_TR[patient, "hdl_real"] <- as.matrix(ukpds_pop[patient, "hdl_real"])
    m_TR[patient, "heart_rate"] <- as.matrix(ukpds_pop[patient, "heart_rate"])
    m_TR[patient, "heart_rate_lag"] <- as.matrix(ukpds_pop[patient, "heart_rate_lag"])
    m_TR[patient, "heart_rate_first"] <- as.matrix(ukpds_pop[patient, "heart_rate_first"])
    m_TR[patient, "heart_rate_real"] <- as.matrix(ukpds_pop[patient, "heart_rate_real"])
    m_TR[patient, "ldl"] <- as.matrix(ukpds_pop[patient, "ldl"])
    m_TR[patient, "ldl_gt35"] <- as.matrix(ukpds_pop[patient, "ldl_gt35"])
    m_TR[patient, "ldl_lag"] <- as.matrix(ukpds_pop[patient, "ldl_lag"])
    m_TR[patient, "ldl_first"] <- as.matrix(ukpds_pop[patient, "ldl_first"])
    m_TR[patient, "ldl_real"] <- as.matrix(ukpds_pop[patient, "ldl_real"])
    m_TR[patient, "albumin_mm"] <- as.matrix(ukpds_pop[patient, "albumin_mm"])
    m_TR[patient, "sbp"] <- as.matrix(ukpds_pop[patient, "sbp"])
    m_TR[patient, "sbp_lag"] <- as.matrix(ukpds_pop[patient, "sbp_lag"])
    m_TR[patient, "sbp_first"] <- as.matrix(ukpds_pop[patient, "sbp_first"])
    m_TR[patient, "sbp_real"] <- as.matrix(ukpds_pop[patient, "sbp_real"])
    m_TR[patient, "wbc"] <- as.matrix(ukpds_pop[patient, "wbc"])
    m_TR[patient, "wbc_lag"] <- as.matrix(ukpds_pop[patient, "wbc_lag"])
    m_TR[patient, "wbc_first"] <- as.matrix(ukpds_pop[patient, "wbc_first"])
    
    # Hardcoded hemoglobin values
    m_TR[patient, "heamo"] <- 15
    m_TR[patient, "heamo_first"] <- 15
    
    # Event history tracking
    event_vars <- c("amp_event", "amp_event2", "amp_hist", "atria_fib",
                    "blindness_event", "blindness_hist", "chf_event", "chf_hist",
                    "esrd_event", "esrd_hist", "ihd_event", "ihd_hist",
                    "mi_event", "mi_hist", "pvd_event", "stroke_event",
                    "stroke_hist", "ulcer_event", "ulcer_hist")
    
    for (var in event_vars) {
      m_TR[patient, var] <- as.matrix(ukpds_pop[patient, var])
    }
    
    # Set default values for lambda, rho, and death
    m_TR[patient, "lambda"] <- 0
    m_TR[patient, "rho"] <- 1
    m_TR[patient, "death"] <- 0
  
  return(m_TR)
}

# dimnames(m_TR)


# Step 3: Define functions for risk factor progression ####
# Function for linear progression of risk factors
#' Calculate Biomarkers 
#'
#' This function calculates patient-specific factors to predict the time path of a biomarker. 
#'
#' @param m_TR A matrix containing patient characteristics over time.
#' @param a_coef_C A 3D array of coefficients used for calculating risk.
#' @param biomarker_eq A character string specifying the health outcome equation (e.g., "ihd").
#' @param time_step An integer indicating the row in `m_TR` to use for calculations.
#' 
#' @return The updated `m_TR` matrix with the event occurrence stored.
#' @export
biomarker <- function(m_TR, a_coef_C, biomarker_eq,  time_step) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  updated_biomarker <- (m_TR[max(1,time_step-1), 1:62] %*%  a_coef_C[1:62,  biomarker_eq, 1] + 
                          as.vector(a_coef_C["lambda",  biomarker_eq, 1]) )
  
  return(updated_biomarker)
}



# Step 4: Create a function to apply all risk factor models ####
# Combine risk factor functions into a single pipeline
# Update patient data over time

#' Update Multiple Biomarkers in a Transition Matrix
#'
#' This function updates multiple biomarker values in the transition matrix for a given time step.
#'
#' @param m_TR The patient trace, a matrix containing patient data with biomarker and event columns.
#' @param a_coef_C A coefficient matrix containing biomarker and event equations.
#' @param time_step An integer representing the current time step.
#' @param next_row An integer indicating the row in `m_TR` to update with new biomarker values.
#'
#' @return The updated transition matrix `m_TR` with new biomarker values in the specified row.
#' 
#' @examples
#' # Example usage
#' m_TR <- all_biomarkers(m_TR, a_coef_C, time_step = 1, next_row = 2)
#'
#' @export

all_biomarkers <- function(m_TR, a_coef_C, time_step, next_row) {
  
# predict the next period (and perform transformations as needed)
  # the biomarkers use real values of variables, but the event equations use transformed variables
  m_TR[next_row, "a1c"] <- biomarker(m_TR, a_coef_C, biomarker_eq = "hba1c", time_step = time_step)
  m_TR[next_row, "sbp_real"] <- biomarker(m_TR, a_coef_C, biomarker_eq = "sbp", time_step = time_step)
    m_TR[next_row, "sbp"] <- m_TR[next_row, "sbp_real"] /10
  m_TR[next_row, "ldl_real"] <- biomarker(m_TR, a_coef_C, biomarker_eq = "ldl", time_step = time_step)
    m_TR[next_row, "ldl"] <- m_TR[next_row, "ldl_real"] * 10
  m_TR[next_row, "hdl_real"] <- biomarker(m_TR, a_coef_C, biomarker_eq = "hdl", time_step = time_step)
    m_TR[next_row, "hdl"] <- m_TR[next_row, "hdl_real"] * 10 
  m_TR[next_row, "bmi"] <- biomarker(m_TR, a_coef_C, biomarker_eq = "bmi", time_step = time_step)
  m_TR[next_row, "heart_rate_real"] <- biomarker(m_TR, a_coef_C, biomarker_eq = "heart_rate", time_step = time_step)
    m_TR[next_row, "heart_rate"] <- m_TR[next_row, "heart_rate_real"] /10
  m_TR[next_row, "wbc"] <- biomarker(m_TR, a_coef_C, biomarker_eq = "wbc", time_step = time_step)
  m_TR[next_row, "heamo"] <- biomarker(m_TR, a_coef_C, biomarker_eq = "haem", time_step = time_step)
  
  
  # Update lag and first occurrence columns
  m_TR[next_row, "a1c_lag"] <- m_TR[time_step, "a1c"]
  m_TR[next_row, "a1c_first"] <- m_TR[1, "a1c"]
  m_TR[next_row, "bmi_lag"] <- m_TR[time_step, "bmi"]
  m_TR[next_row, "bmi_lt18_5"] <- as.integer(m_TR[next_row, "bmi"] < 18.5)
  m_TR[next_row, "bmi_gte25"] <- as.integer(m_TR[next_row, "bmi"] >= 25)
  m_TR[next_row, "bmi_first"] <- m_TR[1, "bmi"]
  

  m_TR[next_row, "hdl_lag"] <- m_TR[time_step, "hdl_real"]
  m_TR[next_row, "hdl_first"] <- m_TR[1, "hdl_real"]
  
  m_TR[next_row, "heart_rate_lag"] <- m_TR[time_step, "heart_rate_real"]
  m_TR[next_row, "heart_rate_first"] <- m_TR[1, "heart_rate_real"]
  # check if this is functioning as a spline
  m_TR[next_row, "ldl_gt35"] <- as.integer(m_TR[next_row, "ldl_real"] > 35) /10
  m_TR[next_row, "ldl_lag"] <- m_TR[time_step, "ldl_real"]
  m_TR[next_row, "ldl_first"] <- m_TR[1, "ldl_real"]
  m_TR[next_row, "sbp_lag"] <- m_TR[time_step, "sbp_real"]
  m_TR[next_row, "sbp_first"] <- m_TR[1, "sbp_real"]
  m_TR[next_row, "wbc_lag"] <- m_TR[time_step, "wbc"]
  m_TR[next_row, "wbc_first"] <- m_TR[1, "wbc"]
  m_TR[next_row, "heamo_first"] <- m_TR[1, "heamo"]
  
  # Update additional values
  m_TR[next_row, "egfr"] <- m_TR[1, "egfr"]
  m_TR[next_row, "egfr_real"] <- m_TR[1, "egfr_real"]
  m_TR[next_row, "egfr_lt60"] <- m_TR[1, "egfr_lt60"]
  m_TR[next_row, "egfr_gte60"] <- m_TR[1, "egfr_gte60"]
  m_TR[next_row, "albumin_mm"] <- m_TR[1, "albumin_mm"]
  
  # Return updated matrix
  return(m_TR)
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
#' The function updates the provided `m_TR` matrix with the event occurrence 
#' at the specified time step.
#'
#' @param m_TR A matrix containing patient characteristics over time.
#' @param a_coef_C A 3D array of coefficients used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_TR` to use for calculations.
#' 
#' @return The updated `m_TR` matrix with the event occurrence stored.
#' @export
weibull_event <- function(m_TR, a_coef_C, health_outcome, health_event, time_step) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (m_TR[time_step, 1:62] %*%  a_coef_C[1:62, health_outcome, 1] + 
                        as.vector(a_coef_C["lambda", health_outcome, 1]) )
  
  # Compute cumulative hazard at the current time step
  cum_hazard_t <- exp(patient_factors) * (m_TR[time_step, "diab_dur"]^(a_coef_C["rho", health_outcome, 1]) )
  
  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- exp(patient_factors) * ((m_TR[time_step, "diab_dur"] + 1)^(a_coef_C["rho", health_outcome, 1]) )
  
  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(1)
  
  # Return the updated matrix
  return(event)
}




#' Calculate Transition Probability Based on a Logistic Regression and Update Patient State 
#'
#' This function calculates patient-specific factors, cumulative hazards, 
#' and the transition probability for a given health outcome (e.g., "ihd"). 
#' The function updates the provided `m_TR` matrix with the event occurrence 
#' at the specified time step.
#'
#' @param m_TR A matrix containing patient characteristics over time.
#' @param a_coef_C A 3D array of coefficients used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_TR` to use for calculations.
#' 
#' @return The updated `m_TR` matrix with the event occurrence stored.
#' @export
logistic_event <- function(m_TR, a_coef_C, health_outcome, health_event, time_step) {
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (m_TR[time_step, 1:62] %*%  a_coef_C[1:62, health_outcome, 1] + 
                        as.vector(a_coef_C["lambda", health_outcome, 1]) )
  
  # Calculate transition probability
  trans_prob=1-(exp(-patient_factors)/(1+exp(-patient_factors)))^1
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(1)
  

  # Return the value
  return(event)
}

#
# u <- logistic_event(m_TR, a_coef_C, health_outcome = "ulcer", health_event = "ulcer_event", time_step = 1)
#m_TR <- weibull_event(m_TR, a_coef_C, health_outcome = "ihd", health_event = "ihd_event", time_step = 1)



# Step 6: Initialize event and history variables ####
#' @title Update Health Events Over Time Steps
#' @description This function updates health events in a patient data matrix (`m_TR`) by applying Weibull 
#' and logistic event functions in a randomized order across multiple time steps.
#'
#' @param m_TR A matrix containing patient-level data, including health event history.
#' @param a_coef_C A coefficient matrix used in Weibull and logistic event calculations.
#' @param time_step An integer indicating the current time step to update events.
#'
#' @return Updated `m_TR` matrix with event and history values updated for the given time step.
#'
#' @examples
#' \dontrun{
#' m_TR <- update_health_events(m_TR, a_coef_C, time_step = 1)
#' }
#' 
#' @export
update_health_events <- function(m_TR, a_coef_C, time_step) {
  # Ensure m_TR remains a matrix
  if (!is.matrix(m_TR)) {
    stop("m_TR must be a matrix.")
  }

  # Atrial Fib and PVD do not update
  m_TR[, "atria_fib"] <- m_TR[1, "atria_fib"]
  m_TR[, "pvd_event"] <- m_TR[1, "pvd_event"]
  
  # Initialize event variables and update history
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")
  
  
  for (event in events) {
    m_TR[time_step, paste0(event, "_event")] <- 0
    m_TR[time_step, paste0(event, "_hist")] <- max(
      m_TR[max(1, time_step-1), paste0(event, "_hist")],
      m_TR[max(1, time_step-1), paste0(event, "_event")]
    )
  }
  m_TR[time_step, "amp_event2"] <- 0
  # Randomize event order
  randomized_events <- sample(events)
  
  for (event in randomized_events) {
    if (event == "amp") {
      amp1_no_ulcer <- weibull_event(m_TR, a_coef_C, health_outcome = "amp1_no_ulcer", health_event = "amp_event", time_step = time_step)
      amp1_yes_ulcer <- weibull_event(m_TR, a_coef_C, health_outcome = "amp1_yes_ulcer", health_event = "amp_event", time_step = time_step)
      
      m_TR[time_step, "amp_event"] <- (amp1_no_ulcer * (m_TR[time_step, "ulcer_hist"] == 0)) + 
        (amp1_yes_ulcer * (m_TR[time_step, "ulcer_hist"] == 1)) 
      #ensure that this is a new event
      m_TR[time_step, "amp_event"] <- m_TR[time_step, "amp_event"] * (m_TR[time_step, "amp_hist"] == 0)
        
      amp2 <- weibull_event(m_TR, a_coef_C, health_outcome = "amp2", health_event = "amp_event2", time_step = time_step)
      m_TR[time_step, "amp_event2"] <- 0 
      m_TR[time_step, "amp_event2"] <- amp2 * (m_TR[time_step, "amp_hist"] == 1)
      
      
    } else if (event == "mi") {
      mi1_male <- weibull_event(m_TR, a_coef_C, health_outcome = "mi1_male", health_event = "mi_event", time_step = time_step)
      mi1_female <- weibull_event(m_TR, a_coef_C, health_outcome = "mi1_female", health_event = "mi_event", time_step = time_step)
      
      m_TR[time_step, "mi_event"] <- (mi1_male * (m_TR[time_step, "female"] == 0)) + 
        (mi1_female * (m_TR[time_step, "female"] == 1))
      m_TR[time_step, "mi_event"] <- m_TR[time_step, "mi_event"] * (m_TR[time_step, "mi_hist"] == 0)
      
      mi2 <- weibull_event(m_TR, a_coef_C, health_outcome = "mi2", health_event = "mi_event", time_step = time_step)
      m_TR[time_step, "mi_event"] <- (m_TR[time_step, "mi_hist"] == 0) * m_TR[time_step, "mi_event"] + 
        (m_TR[time_step, "mi_hist"] == 1) * mi2
      
    } else if (event == "stroke") {
      stroke1 <- weibull_event(m_TR, a_coef_C, health_outcome = "stroke_1", health_event = "stroke_event", time_step = time_step)
      stroke2 <- weibull_event(m_TR, a_coef_C, health_outcome = "stroke_2", health_event = "stroke_event", time_step = time_step)
      
      m_TR[time_step, "stroke_event"] <- (stroke1 * (m_TR[time_step, "stroke_hist"] == 0)) + 
        (stroke2 * (m_TR[time_step, "stroke_hist"] == 1))
      m_TR[time_step, "stroke_event"] <- m_TR[time_step, "stroke_event"] * (m_TR[time_step, "stroke_hist"] == 0)
      
      
    } else if (event == "ulcer") {
      m_TR[time_step, "ulcer_event"] <- logistic_event(m_TR, a_coef_C, health_outcome = "ulcer", health_event = "ulcer_event", time_step = time_step)
      m_TR[time_step, "ulcer_event"] <- m_TR[time_step, "ulcer_event"] * (m_TR[time_step, "ulcer_hist"] == 0)
      
      
    } else {
      m_TR[time_step, paste0(event, "_event")] <- weibull_event(m_TR, a_coef_C, health_outcome = event, health_event = paste0(event, "_event"), time_step = time_step)
      
      m_TR[time_step,  paste0(event, "_event")] <- m_TR[time_step,  paste0(event, "_event")] * (m_TR[time_step,  paste0(event, "_hist")] == 0)
      
    }
  }
  
  return(m_TR)
}







# Step 7: Define a mortality function ####
# Combine relevant event functions affecting mortality
# Estimate survival probability over time
#' Calculate Transition Probability Based on a Gompertz Model and Update Patient State
#'
#'
#' This function calculates patient-specific factors, cumulative hazards, 
#' and the transition probability for mortality. 
#' The function updates the provided `m_TR` matrix with the event occurrence 
#' at the specified time step.
#'
#' @param m_TR A matrix containing patient characteristics over time.
#' @param a_coef_C A 3D array of coefficients used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_TR` to use for calculations.
#' 
#' @return The updated `m_TR` matrix with the event occurrence stored.
#' @export
gompertz_event <- function(m_TR, a_coef_C, health_outcome, health_event, time_step) {
  
  
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (m_TR[time_step, 1:62] %*%  a_coef_C[1:62, health_outcome, 1] + 
                        as.vector(a_coef_C["lambda", health_outcome, 1]) )
  
  # Compute cumulative hazard at the current time step
  cum_hazard_t <- (1/a_coef_C["rho", health_outcome, 1])* exp(patient_factors) * (exp(m_TR[time_step, "age"]*(a_coef_C["rho", health_outcome, 1])) -1 )
  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- (1/a_coef_C["rho", health_outcome, 1])* exp(patient_factors) * (exp((m_TR[time_step, "age"]+1)*(a_coef_C["rho", health_outcome, 1])) -1 )
  
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
#' @param m_TR A matrix containing patient-level data, including health event history.
#' @param a_coef_C A coefficient matrix used in Gompertz and logistic event calculations.
#' @param time_step An integer specifying the time step at which mortality should be calculated.
#'
#' @return The updated `m_TR` matrix with the mortality status recorded for the specified time step.
#'
#' @examples
#' \dontrun{
#' m_TR <- mortality(m_TR, a_coef_C, time_step = 5)
#' }
#'
#' @export
mortality <- function(m_TR, a_coef_C, time_step) {
  
  # Calculate new health event occurrence
  new_event <- max(m_TR[time_step, "amp_event"], m_TR[time_step, "amp_event2"],
                   m_TR[time_step, "blindness_event"], m_TR[time_step, "chf_event"],
                   m_TR[time_step, "esrd_event"], m_TR[time_step, "ihd_event"],
                   m_TR[time_step, "mi_event"], m_TR[time_step, "stroke_event"],
                   m_TR[time_step, "ulcer_event"])
  
  # Calculate any prior history of health events
  any_history <- max(m_TR[time_step, "amp_hist"], 
                     m_TR[time_step, "blindness_hist"], m_TR[time_step, "chf_hist"],
                     m_TR[time_step, "esrd_hist"], m_TR[time_step, "ihd_hist"],
                     m_TR[time_step, "mi_hist"], m_TR[time_step, "stroke_hist"],
                     m_TR[time_step, "ulcer_hist"])
  
  # Determine event-history combinations
  nhne <- new_event == 0 & any_history == 0  # No history, no event
  yhne <- new_event == 0 & any_history == 1  # Yes history, no event
  nhye <- new_event == 1 & any_history == 0  # No history, new event
  yhye <- new_event == 1 & any_history == 1  # Yes history, new event
  
  # Mortality calculations using Gompertz and logistic models
  death_nhne <- gompertz_event(m_TR, a_coef_C, health_outcome = "death_nhne", 
                               health_event = "death_nhne", time_step = time_step)
  
  death_yhne <- gompertz_event(m_TR, a_coef_C, health_outcome = "death_yhne", 
                               health_event = "death_yhne", time_step = time_step)
  
  death_nhye <- logistic_event(m_TR, a_coef_C, health_outcome = "death_1st_event", 
                               health_event = "death_nhye", time_step = time_step)
  
  death_yhye <- logistic_event(m_TR, a_coef_C, health_outcome = "death_yhye", 
                               health_event = "death_yhye", time_step = time_step)
  
  # Calculate new mortality status
  new_death <- nhne * death_nhne + yhne * death_yhne + nhye * death_nhye + yhye * death_yhye
  
  # Update the mortality status in the matrix for the given time step
  m_TR[time_step, "death"] <- new_death + m_TR[max(time_step - 1, 1), "death"]
  
  return(m_TR)
}


# Step 8: Simulate disease progression and mortality for an individual patient ####
# Initialize patient data
# Loop through time points to update risk factors and events
# Store results
# discount rate 


# Step 9: Simulate disease progression and mortality for 999 additional patients
# Loop over 1000 patients
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
  dimnames = list(c("1","2","3","4"),column_names)  
)

for (patient in 1:4)  {
#create a patient population 
m_TR <- initialize_patients(patient, ukpds_pop, m_TR)

# carry forward time invariant characteristics 
  m_TR[ ,"age_diag"]<-m_TR[1 ,"age_diag"]
  m_TR[ ,"black"]<-m_TR[1 ,"black"]
  m_TR[ ,"indian"  ]<-m_TR[1 ,"indian"]
  m_TR[ ,"female" ]<-m_TR[1,"female" ]
  m_TR[ ,"smoke"]<- m_TR[1,"smoke"]



for (time_step in 1:num_cycles) {

  m_TR[time_step,"death"]<-m_TR[max(time_step-1,1),"death"]
  m_TR[time_step,"lambda"]<-1
  m_TR[time_step,"rho"]<- 1
  
  m_TR[time_step,"age"]<-m_TR[max(1,time_step-1),"age"] +1
  m_TR[time_step,"diab_dur"]<-m_TR[max(1,time_step-1),"diab_dur"]+1    
  m_TR[time_step,"diab_dur_log"]<- (log(m_TR[time_step,"diab_dur"]))
  
  
  # ready to simulate 
  # event prediction at t
  m_TR <- update_health_events(m_TR, a_coef_C, time_step = time_step)
  # mortality prediction at t
  m_TR <- mortality(m_TR, a_coef_C, time_step = time_step)
  #predict the risk factors for the next cycle (t+1) 

      m_TR<- all_biomarkers(m_TR, a_coef_C, time_step = time_step, next_row = time_step+1) 
 
}
  m_TR_new <- m_TR[-nrow(m_TR), ]
  

 

  m_summary <- matrix(   
    data = NA, 
    nrow = length(cycles), 
    ncol = 4,   
    dimnames = list(cycles,column_names)  
  ) 
  m_summary <- m_summary[-nrow(m_summary), ]


  for (time_step in 1:num_cycles) {
    m_summary[time_step, "cost"]<- as.double(m_TR_new[time_step,"death"]==0)*c_baseline + 
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step,"blindness_event"] * c_blindness_e +
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step,"blindness_hist"] * c_blindness_c +
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step, "amp_event"] * c_amp_e + 
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step, "amp_event2"] * c_amp_e + 
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step, "amp_hist"] * c_amp_c + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "chf_event"] * c_chf_e + 
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step, "chf_hist"] * c_chf_c + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "esrd_event"] * c_esrd_e + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "esrd_hist"] * c_esrd_c + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "ihd_event"] * c_ihd_e + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "ihd_hist"] * c_ihd_c + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "mi_event"] * c_mi_e + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "mi_hist"] * c_mi_c + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "stroke_event"] * c_stroke_e + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "stroke_hist"] * c_stroke_c + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "ulcer_event"] * c_ulcer_e + 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "ulcer_hist"] * c_ulcer_c
    
    m_summary[time_step, "qalys"]<- as.double(m_TR_new[time_step,"death"]==0)*q_baseline + min(
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step,"blindness_event"] * q_blindness ,
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step,"blindness_hist"] * q_blindness ,
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step, "amp_event"] * q_amp , 
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step, "amp_event2"] * q_amp , 
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step, "amp_hist"] * q_amp ,
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "chf_event"] * q_chf , 
      as.double(m_TR_new[time_step,"death"]==0)*m_TR[time_step, "chf_hist"] * q_chf ,
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "esrd_event"] * q_esrd , 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "esrd_hist"] * q_esrd , 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "ihd_event"] * q_ihd , 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "ihd_hist"] * q_ihd , 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "mi_event"] * q_mi , 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "mi_hist"] * q_mi ,
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "stroke_event"] * q_stroke , 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "stroke_hist"] * q_stroke , 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "ulcer_event"] * q_ulcer , 
      as.double(m_TR_new[time_step,"death"]==0)* m_TR[time_step, "ulcer_hist"] * q_ulcer )
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
