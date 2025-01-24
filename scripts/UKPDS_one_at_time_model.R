
library(pacman)
pacman::p_load(haven,
               Rcpp,
               RcppArmadillo,
               readr,
               dplyr,
               profvis)


# Step 1: Import the matrix of coefficients
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

print(dim(a_coef_C)) # to verify the dimensions, 61 rows, 25 columns, 1 slice
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
patient<-1


# Initialize baseline values for the patient
# note the as.matrix ensures that the matrix isn't converted to a list
m_TR[1, "age"] <- as.matrix(ukpds_pop[patient, "age"])
m_TR[1, "age_diag"] <- as.matrix(ukpds_pop[patient, "age_diag"])
m_TR[1, "black"] <- as.matrix(ukpds_pop[patient, "black"])
m_TR[1, "indian"] <- as.matrix(ukpds_pop[patient, "indian"])
m_TR[1, "female"] <- as.matrix(ukpds_pop[patient, "female"])
m_TR[1, "diab_dur"] <- as.matrix(ukpds_pop[patient, "diab_dur"])
m_TR[1, "diab_dur_log"] <- as.matrix(ukpds_pop[patient, "diab_dur_log"])
m_TR[1, "smoke"] <- as.matrix(ukpds_pop[patient, "smoke"])
m_TR[1, "a1c"] <- as.matrix(ukpds_pop[patient, "a1c"])
m_TR[1, "a1c_lag"] <- as.matrix(ukpds_pop[patient, "a1c_lag"])
m_TR[1, "a1c_first"] <- as.matrix(ukpds_pop[patient, "a1c_first"])
m_TR[1, "bmi"] <- as.matrix(ukpds_pop[patient, "bmi"])
m_TR[1, "bmi_lt18_5"] <- as.matrix(ukpds_pop[patient, "bmi_lt18_5"])
m_TR[1, "bmi_gte25"] <- as.matrix(ukpds_pop[patient, "bmi_gte25"])
m_TR[1, "bmi_lag"] <- as.matrix(ukpds_pop[patient, "bmi_lag"])
m_TR[1, "bmi_first"] <- as.matrix(ukpds_pop[patient, "bmi_first"])
m_TR[1, "egfr"] <- as.matrix(ukpds_pop[patient, "egfr"])
m_TR[1, "egfr_lt60"] <- as.matrix(ukpds_pop[patient, "egfr_lt60"])
m_TR[1, "egfr_gte60"] <- as.matrix(ukpds_pop[patient, "egfr_gte60"])
m_TR[1, "hdl"] <- as.matrix(ukpds_pop[patient, "hdl"])
m_TR[1, "hdl_lag"] <- as.matrix(ukpds_pop[patient, "hdl_lag"])
m_TR[1, "hdl_first"] <- as.matrix(ukpds_pop[patient, "hdl_first"])
m_TR[1, "heart_rate"] <- as.matrix(ukpds_pop[patient, "heart_rate"])
m_TR[1, "heart_rate_lag"] <- as.matrix(ukpds_pop[patient, "heart_rate_lag"])
m_TR[1, "heart_rate_first"] <- as.matrix(ukpds_pop[patient, "heart_rate_first"])
m_TR[1, "ldl"] <- as.matrix(ukpds_pop[patient, "ldl"])
m_TR[1, "ldl_gt35"] <- as.matrix(ukpds_pop[patient, "ldl_gt35"])
m_TR[1, "ldl_lag"] <- as.matrix(ukpds_pop[patient, "ldl_lag"])
m_TR[1, "ldl_first"] <- as.matrix(ukpds_pop[patient, "ldl_first"])
m_TR[1, "albumin_mm"] <- as.matrix(ukpds_pop[patient, "albumin_mm"])
m_TR[1, "sbp"] <- as.matrix(ukpds_pop[patient, "sbp"])
m_TR[1, "sbp_lag"] <- as.matrix(ukpds_pop[patient, "sbp_lag"])
m_TR[1, "sbp_first"] <- as.matrix(ukpds_pop[patient, "sbp_first"])
m_TR[1, "wbc"] <- as.matrix(ukpds_pop[patient, "wbc"])
m_TR[1, "wbc_lag"] <- as.matrix(ukpds_pop[patient, "wbc_lag"])
m_TR[1, "wbc_first"] <- as.matrix(ukpds_pop[patient, "wbc_first"])
#m_TR[1, "heamo"] <- as.matrix(ukpds_pop[patient, "heamo"])
#m_TR[1, "heamo_first"] <- as.matrix(ukpds_pop[patient, "heamo_first"])
# find these and fill these in
m_TR[1, "heamo"] <- 1
m_TR[1, "heamo_first"] <- 1

m_TR[1, "amp_event"] <- as.matrix(ukpds_pop[patient, "amp_event"])
m_TR[1, "amp_event2"] <- as.matrix(ukpds_pop[patient, "amp_event2"])
m_TR[1, "amp_hist"] <- as.matrix(ukpds_pop[patient, "amp_hist"])
m_TR[1, "atria_fib"] <- as.matrix(ukpds_pop[patient, "atria_fib"])
m_TR[1, "blindness_event"] <- as.matrix(ukpds_pop[patient, "blindness_event"])
m_TR[1, "blindness_hist"] <- as.matrix(ukpds_pop[patient, "blindness_hist"])
m_TR[1, "chf_event"] <- as.matrix(ukpds_pop[patient, "chf_event"])
m_TR[1, "chf_hist"] <- as.matrix(ukpds_pop[patient, "chf_hist"])
m_TR[1, "esrd_event"] <- as.matrix(ukpds_pop[patient, "esrd_event"])
m_TR[1, "esrd_hist"] <- as.matrix(ukpds_pop[patient, "esrd_hist"])
m_TR[1, "ihd_event"] <- as.matrix(ukpds_pop[patient, "ihd_event"])
m_TR[1, "ihd_hist"] <- as.matrix(ukpds_pop[patient, "ihd_hist"])
m_TR[1, "mi_event"] <- as.matrix(ukpds_pop[patient, "mi_event"])
m_TR[1, "mi_hist"] <- as.matrix(ukpds_pop[patient, "mi_hist"])
m_TR[1, "pvd_event"] <- as.matrix(ukpds_pop[patient, "pvd_event"])
m_TR[1, "stroke_event"] <- as.matrix(ukpds_pop[patient, "stroke_event"])
m_TR[1, "stroke_hist"] <- as.matrix(ukpds_pop[patient, "stroke_hist"])
m_TR[1, "ulcer_event"] <- as.matrix(ukpds_pop[patient, "ulcer_event"])
m_TR[1, "ulcer_hist"] <- as.matrix(ukpds_pop[patient, "ulcer_hist"])

dimnames(m_TR)
# Step 3: Define functions for risk factor progression
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
  updated_biomarker <- (m_TR[time_step, 1:57] %*%  a_coef_C[1:57,  biomarker_eq, 1] + 
                          as.vector(a_coef_C["lambda",  biomarker_eq, 1]) )
  
  return(updated_biomarker)
}



# Step 4: Create a function to apply all risk factor models
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
  # List of biomarkers to update
  biomarkers <- c("a1c", "sbp", "ldl", "hdl", "bmi", 
                  "heart_rate", "wbc", "heamo")
  
  # Corresponding biomarker equation names in a_coef_C
  biomarker_eqs <- c("hba1c", "sbp", "ldl", "hdl", "bmi", 
                     "heart_rate", "wbc", "haem")
  
  # Loop through each biomarker and update the matrix
  for (i in seq_along(biomarkers)) {
    m_TR[next_row, biomarkers[i]] <- biomarker(
      m_TR, a_coef_C, biomarker_eq = biomarker_eqs[i], time_step = time_step
    )
  }
  
  # Update lag and first occurrence columns
  m_TR[next_row, "a1c_lag"] <- m_TR[time_step, "a1c"]
  m_TR[next_row, "a1c_first"] <- m_TR[1, "a1c"]
  m_TR[next_row, "bmi_lag"] <- m_TR[time_step, "bmi"]
  m_TR[next_row, "bmi_lt18_5"] <- as.integer(m_TR[next_row, "bmi"] < 18.5)
  m_TR[next_row, "bmi_gte25"] <- as.integer(m_TR[next_row, "bmi"] >= 25)
  
  m_TR[next_row, "bmi_first"] <- m_TR[1, "bmi"]
  m_TR[next_row, "hdl_lag"] <- m_TR[time_step, "hdl"]
  m_TR[next_row, "hdl_first"] <- m_TR[1, "hdl"]
  m_TR[next_row, "heart_rate_lag"] <- m_TR[time_step, "heart_rate"]
  m_TR[next_row, "heart_rate_first"] <- m_TR[1, "heart_rate"]
  
  m_TR[next_row, "ldl_gt35"] <- as.integer(m_TR[next_row, "ldl"] > 35)
  m_TR[next_row, "ldl_lag"] <- m_TR[time_step, "ldl"]
  m_TR[next_row, "ldl_first"] <- m_TR[1, "ldl"]
  m_TR[next_row, "sbp_lag"] <- m_TR[time_step, "sbp"]
  m_TR[next_row, "sbp_first"] <- m_TR[1, "sbp"]
  m_TR[next_row, "wbc_lag"] <- m_TR[time_step, "wbc"]
  m_TR[next_row, "wbc_first"] <- m_TR[1, "wbc"]
  m_TR[next_row, "heamo_first"] <- m_TR[1, "heamo"]
  
  # Update additional values
  m_TR[next_row, "egfr"] <- m_TR[1, "egfr"]
  m_TR[next_row, "egfr_lt60"] <- m_TR[1, "egfr_lt60"]
  m_TR[next_row, "egfr_gte60"] <- m_TR[1, "egfr_gte60"]
  m_TR[next_row, "albumin_mm"] <- m_TR[1, "albumin_mm"]
  
  # Return updated matrix
  return(m_TR)
}

#m_TR<- all_biomarkers(m_TR, a_coef_C, time_step=1, next_row=2) 


# Step 5: Define event functions (Weibull/Exponential and Logistic)
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
  patient_factors <- (m_TR[time_step, 1:57] %*%  a_coef_C[1:57, health_outcome, 1] + 
                        as.vector(a_coef_C["lambda", health_outcome, 1]) )
  
  # Compute cumulative hazard at the current time step
  cum_hazard_t <- exp(patient_factors) * (m_TR[time_step, "diab_dur"]^(a_coef_C["rho", health_outcome, 1]) )
  
  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- exp(patient_factors) * ((m_TR[time_step, "diab_dur"] + 1)^(a_coef_C["rho", health_outcome, 1]) )
  
  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(1)
  
  # Store the event occurrence in the next time step of m_TR
  m_TR[time_step +1 , health_event ] <- event
  # Return the updated matrix
  return(m_TR)
}


## UPDATE THIS - it should predict based on this period

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
  patient_factors <- (m_TR[time_step, 1:57] %*%  a_coef_C[1:57, health_outcome, 1] + 
                        as.vector(a_coef_C["lambda", health_outcome, 1]) )
  
  # Calculate transition probability
  trans_prob=1-(exp(-patient_factors)/(1+exp(-patient_factors)))^1
  
  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(1)
  
  # Store the event occurrence in the next time step of m_TR
  m_TR[time_step +1 , health_event ] <- event
  # Return the updated matrix
  return(m_TR)
}

m_TR <- logistic_event(m_TR, a_coef_C, health_outcome = "ulcer", health_event = "ulcer_event", time_step = 1)
m_TR <- weibull_event(m_TR, a_coef_C, health_outcome = "ihd", health_event = "ihd_event", time_step = 1)




# Step 6: Create a function to call all event functions in random order
# Randomly determine event order
# Apply event functions sequentially to patient data



# Step 7: Define a mortality function
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
  patient_factors <- (m_TR[time_step, 1:57] %*%  a_coef_C[1:57, health_outcome, 1] + 
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

death_nhne <- gompertz_event(m_TR, a_coef_C, health_outcome = "death_nhne", health_event = "death_nhne", time_step = 1)

# Step 8: Simulate disease progression and mortality for an individual patient
# Initialize patient data
# Loop through time points to update risk factors and events
# Store results
dimnames(m_TR)
m_TR[ ,"age_diag"]<-m_TR[1 ,"age_diag"]
m_TR[ ,"black"]<-m_TR[1 ,"black"]
m_TR[ ,"indian"  ]<-m_TR[1 ,"indian"]
m_TR[ ,"female" ]<-m_TR[1,"female" ]
m_TR[ ,"smoke"]<- m_TR[1,"smoke"]


m_TR[2,"age"]<-m_TR[1,"age"] +1
m_TR[2,"diab_dur"]<-m_TR[1,"diab_dur"]+1      
m_TR[2,"diab_dur_log"]<- log(m_TR[2,"diab_dur"])




# Step 9: Simulate disease progression and mortality for 999 additional patients
# Loop over 1000 patients
# Store and summarize population-level results

# Step 10: Summarize and visualize results
# Calculate summary statistics
# Generate plots to visualize disease progression trends
